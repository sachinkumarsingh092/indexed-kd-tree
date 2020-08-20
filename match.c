#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <error.h>
#include <stdbool.h>

#include "kdtree.c"

#include <gnuastro/qsort.h>
#include <gnuastro/table.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
#include <gnuastro/polygon.h>
#include <gnuastro/statistics.h>


/* Internally used macro to help in the processing */
# ifndef M_PI
# define M_PI 3.14159265358979323846
# endif

#define max(a,b) \
   ({ __typeof__ (a) _a = (a);  \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

# define DEG2RAD(d) (d)*M_PI/180
# define RAD2DEG(d) (d)*180/M_PI



/* The grid structure which has the necessary
   information of the grid on the catalogue. */
struct grid
{
	size_t size;
	size_t dim[2];
	double min[2], step_size[2];
};



struct matched_points
{
  double x, y;
  double ra, dec;
  struct matched_points *next;
};



/* Structure to keep parameters for multithreaded worker function. */
struct params
{
  /* Inputs. */
  gal_data_t *ra;
  gal_data_t *dec;
  gal_data_t *magnitude;
  gal_data_t *x, *y;
  gal_data_t *qmagnitude;

  /* Internal. */
  double max_star_dis_in_quad;
  size_t *brightest_star_id;

  /* The number of preallocated elements is the number of threads. */
  struct matched_points **matched;

  /* Output. */
  gal_data_t *cx, *cy;
  gal_data_t *dx, *dy;
  gal_data_t *rel_brightness;
  gal_data_t *a_ind, *b_ind, *c_ind, *d_ind;
  gal_data_t *left, *right;
};




/* For internal use only. */
struct quad_vertex
{
  size_t index;
  double distance;
};







/***************Visualization*****************/

/* Make a ds9 complatiable polygon region file for easy visualizations.

   Arguments
   =========
   char *filename      - The filename of the region file.
   size_t num_polygons - Total number of polygons.
   double *polygon     - The array containg the points of the polygon.
   size_t num_vertices - Number of vertices of each polygon.
   char *color         - The color of the polygons in DS9.

   Return
   =======
   The region file having polygons of `num_vertices` points.
*/
void
gal_polygon_to_ds9reg(char *filename, size_t num_polygons, double *polygon,
                      size_t num_vertices, char *color)
{
  FILE *fileptr;
  size_t i, j, n;
  size_t *ordinds=NULL;
  double *temp_sorted_arr=NULL;

  /* If no color is selected, then set default color to green. */
  if (!color) color="green";

  /* Open the file. */
  fileptr = fopen(filename, "w+");

  /* Write the format of polygon region file. */
  fprintf(fileptr, "# Region file format: DS9 version 4.1\n");
  fprintf(fileptr, "global color=%s dashlist=8 3 width=1 "
                    "font=\"helvetica 10 normal roman\" select=1 "
                    "highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 "
                    "include=1 source=1\n", color);
  fprintf(fileptr, "fk5\n");

  /* Write the vertices in the file. */
  j=0;
  for (i=0; i<num_polygons; ++i)
    {
      /* Allocate temporary arrays. */
      ordinds=malloc(num_vertices*2*(sizeof(*ordinds)));
      temp_sorted_arr=malloc(num_vertices*2*sizeof(*temp_sorted_arr));

      /* Make a temporary array for sorting. */
      for (n=0; n<2*num_vertices; ++n)
        temp_sorted_arr[n]=polygon[j++];

      /* Sort the vertices in couterclockwise. */
      gal_polygon_vertices_sort(temp_sorted_arr, num_vertices, ordinds);

      /* Start writing the polygons to file. */
      fprintf (fileptr, "polygon(");
      for (n=0; n<num_vertices; ++n)
        {
          if (n%num_vertices==0)
            fprintf(fileptr, "%lf,%lf", temp_sorted_arr[2*ordinds[n]  ],
                                        temp_sorted_arr[2*ordinds[n]+1]);
          else
            fprintf(fileptr, ",%lf,%lf", temp_sorted_arr[2*ordinds[n]  ],
                                         temp_sorted_arr[2*ordinds[n]+1]);
        }
      fprintf(fileptr, ")\n");

      /* free the temporary arrays. */
      free(temp_sorted_arr);
      free(ordinds);
    }

  /* Close the file. */
  fclose(fileptr);
}








/********************Grid formation*************************/

/* Make a grid of boxes on the catalogue dividing it in ranges
   of RA and DEC. */
void
grid_make(gal_data_t *ra, gal_data_t *dec, size_t ra_numbins,
          size_t dec_numbins, struct grid* out_grid)
{
	gal_data_t *ra_min=NULL, *ra_max=NULL;
	gal_data_t *dec_min=NULL, *dec_max=NULL;
	double *max_ra_val=NULL, *min_ra_val=NULL;
	double *max_dec_val=NULL, *min_dec_val=NULL;

  /* Assign the dimensions of the grid. */
	out_grid->dim[0]=ra_numbins;
	out_grid->dim[1]=dec_numbins;

  /* The size of the grid(=dim1xdim2). */
	out_grid->size=out_grid->dim[0]*out_grid->dim[1];

	/* Check range for RA. */
	ra_max=gal_statistics_maximum (ra);
  ra_min=gal_statistics_minimum (ra);

  max_ra_val=ra_max->array;
  min_ra_val=ra_min->array;

  /* Assign the minimum RA value and the size of steps=(max-min/no of grid)
     to reach the last grid fro the forst step. */
	out_grid->min[0]=min_ra_val[0];
	out_grid->step_size[0]=(max_ra_val[0]-min_ra_val[0])/ra_numbins;

	/* Check range for DEC. */
	dec_max=gal_statistics_maximum (dec);
  dec_min=gal_statistics_minimum (dec);

  max_dec_val=dec_max->array;
  min_dec_val=dec_min->array;

  /* Assign the minimum DEC value and the size of steps=(max-min/no of grid)
     to reach the last grid fro the forst step. */
	out_grid->min[1]=min_dec_val[0];
	out_grid->step_size[1]=(max_dec_val[0]-min_dec_val[0])/dec_numbins;

  /* Free allocated data. */
	gal_data_free(ra_max);
	gal_data_free(ra_min);

	gal_data_free(dec_max);
	gal_data_free(dec_min);
}







/* Return the index of the grid box. */
size_t
grid_findindex(double ra, double dec, struct grid *in_grid)
{
	size_t x=(ra-in_grid->min[0])/in_grid->step_size[0];
	size_t y=(dec-in_grid->min[1])/in_grid->step_size[1];

  /* If the value is maximum in ra and dec, the division returns
     the next integer instead of previous one. This only occurs
     for edge values. Hence only for the edge values x and y are
     checked and shifted to the last index rather than next to last. */
  x = ( x==in_grid->dim[0]
         ? in_grid->dim[0]-1
         : x );
  y = ( y==in_grid->dim[1]
         ? in_grid->dim[1]-1
         : y );


	return y*in_grid->dim[0]+x;
}






/****************finding brightest stars****************************/


/* Returns a 2x2 array with first index having the index of the
   box in the grid and the second index having the IDs of the stars. */
size_t *
find_brightest_stars(gal_data_t *ra, gal_data_t *dec, gal_data_t *magnitude,
                    size_t num_quads, struct grid *in_grid, size_t num_stars_per_gpixel)
{
  size_t i, index;
  size_t *sorted_id=NULL;
  double *ra_arr=ra->array;
	double *dec_arr=dec->array;
  float *magnitude_arr=magnitude->array;
  size_t *brightest_star_id=gal_pointer_allocate (GAL_TYPE_SIZE_T, num_quads, 0,
                                                  __func__, "brightest_star_id");
  size_t *bsi_counter=gal_pointer_allocate (GAL_TYPE_SIZE_T, in_grid->size, 1,
                                            __func__, "bsi_counter");

  /* Initialise output array. */
  for(i=0;i<num_quads;++i)
    brightest_star_id[i]=GAL_BLANK_SIZE_T;


  /* Allocate object id array and initialize. */
  sorted_id=gal_pointer_allocate (GAL_TYPE_SIZE_T, magnitude->size, 0,
                                  __func__, "sorted_id");

  for(i=0;i<magnitude->size; ++i)
    sorted_id[i]=i;


  /* Pointer to an array magnitude_arr to be used as a reference. */
  gal_qsort_index_single=magnitude_arr;


 /* Use brightness column as a reference to sort stars ID array
      in incresing order(brightness is inverse of number's
      magnitudes).  */
  qsort (sorted_id, magnitude->size, sizeof(size_t),
        gal_qsort_index_single_float32_d);


  /* Iterate through all the stars in the grid and save IDs
     of the brightest stars. The array of magnitudes of brightness
     is already sorted so for each box in grid only the brightest
     star indexes will be stored. */
  for (i=0;i<magnitude->size;++i)
    {
      /* Index of the grid box corresponding to particular RA and DEC
         values in the sorted order. */
      index = grid_findindex (ra_arr[sorted_id[i]],
                              dec_arr[sorted_id[i]], in_grid);


      /* If there are less number of star ids for the box than required,
         find more stars in that particular box. */
      if (bsi_counter[index] < num_stars_per_gpixel)
        {
          // printf("i = %zu index = %zu bsi_couter = %zu ind =  %zu\n",i , index, bsi_counter[index],
          //                                                   index*num_stars_per_gpixel+bsi_counter[index]);
          brightest_star_id[index*num_stars_per_gpixel+bsi_counter[index]]=sorted_id[i];
          bsi_counter[index]++;
        }
    }

  /* Clean and return. */
  free(sorted_id);
  free(bsi_counter);
  return brightest_star_id;
}








/**************making hashes************************/

/* Return the angle in range [0, 360) degrees between
   points A, O and B where O is the vertex of angle AOB. */
double
find_angle_aob(double ax, double ay, double ox, double oy, double bx, double by)
{
  double dir_o_to_a=0, dir_o_to_b=0, angle_aob=0;

  dir_o_to_a = atan2(ay - oy, ax - ox);
  dir_o_to_b = atan2(by - oy, bx - ox);

  angle_aob = dir_o_to_a - dir_o_to_b;

  /* atan2 returns anngle between [-pi, +pi] radians.
     Convert it [0, 360) degrees. */
  angle_aob = RAD2DEG(angle_aob);

  return angle_aob>=0
         ? angle_aob
         : angle_aob+360;
}





/* Make the final hash codes using angle between line AB and AC and AD.
   Finally return a array of size 4 cotaining the code - {Cx, Cy, Dx, Dy}.
   Also, after this the quads will be sorted by brightness, so this
   further removes redundencies. */
void
hash_geometric_fix_cd(size_t qindex, struct quad_vertex* sorted_vertices, struct params *p
                      /*TODO: double *cx, double *cy ...*/)
{
  size_t i;
  double scale=0;
  double angle_car=0, angle_dar=0;
  double mag_ab=0, mag_ac=0, mag_ad=0;
  double *cx=p->cx->array, *cy=p->cy->array;
  double *dx=p->dx->array, *dy=p->dy->array;
  uint32_t *c_ind=p->c_ind->array;
  uint32_t *d_ind=p->d_ind->array;
  double ra[4], dec[4], *ra_arr=p->ra->array, *dec_arr=p->dec->array;
  double ref_ra, ref_dec;

  /* Fill the ra and dec arrays. Here indexes of a, b, c, d
     are 0, 1, 2, 3 respectively. */
  for(i=0;i<4;++i)
    {
      ra[i]=ra_arr[sorted_vertices[i].index];
      dec[i]=dec_arr[sorted_vertices[i].index];
    }

  /* Firstly we find a point on the x-axis of the local coordinate
     system and use it as a reference(point R) so that we can find angles
     between CAR and DAR. We assume AB to be angle bisector of the local
     coordinate system.

     Point R will be the component of AB in x-axis.
     Therefore point R = (|AB|cos45, 0). */
  mag_ab=sqrt((ra[0]-ra[1])*(ra[0]-ra[1])
               +(dec[0]-dec[1])*(dec[0]-dec[1]));

  /* Make a scaled-down value for the new coordinate system. We want that
                |AB|=1, and hence scale = 1/|AB|
     where |AB|=sqrt(|AB.ra|^2+|AB.dec|^2).
     This value can be used to scale-down the other distaces mainly,
     AC and AD. */

  scale=1/mag_ab;
  // printf("scale = %g\n", scale);

  ref_ra =scale*mag_ab*cos(45);
  ref_dec=0;

  /* Find angles CAR and DAR. */
  angle_car=find_angle_aob(ra[2], dec[2], ra[0], dec[0], ref_ra, ref_dec);

  angle_dar=find_angle_aob(ra[3], dec[3], ra[0], dec[0], ref_ra, ref_dec);

  /* For a check:
  printf("angle_car = %g, angle_dar = %g\n", angle_car, angle_dar);
  */


  /* Find |AC| and |AD|. Then use distance to find cos and sin thetas to
     give ra and dec coordinates. */
  mag_ac=sqrt((ra[2]-ra[0])*(ra[2]-ra[0])
               +(dec[2]-dec[0])*(dec[2]-dec[0]));

  mag_ad=sqrt((ra[3]-ra[0])*(ra[3]-ra[0])
               +(dec[3]-dec[0])*(dec[3]-dec[0]));


  /* Make the final hash-code. */
  cx[qindex]=mag_ac*scale*cos(angle_car); /* Cx */
  cy[qindex]=mag_ac*scale*sin(angle_car); /* Cy */
  dx[qindex]=mag_ad*scale*cos(angle_dar); /* Dx */
  dy[qindex]=mag_ad*scale*sin(angle_dar); /* Dy */

  /* If the positions of star dont follow the given conditions,
   swap them. */
  if(cx[qindex] > dx[qindex] || cx[qindex]+dx[qindex]>1)
    {
      /* Swap c_ind[qindex] and d_ind[qindex]. */
      c_ind[qindex]=c_ind[qindex]+d_ind[qindex];
      d_ind[qindex]=c_ind[qindex]-d_ind[qindex];
      c_ind[qindex]=c_ind[qindex]-d_ind[qindex];
    }
}






/* Find the relative brightness of the stars in the quad
   and make a unique number to represtent the relative
   brightness of the quad to make it further unique
   while detection. */
void
hash_brightness(size_t qindex, struct quad_vertex* sorted_vertices, struct params *p)
{
  size_t i, j;
  float tmp;
  float magnitude[4];
  size_t rel_magnitude_index[4];  /* a, b, c, d values. */
  float *mag_arr=p->magnitude->array;
  /* Initialise to the nearest power of 2 such that all previous multiple of 4 bits are 0.*/
  uint16_t a_bit=1, b_bit=16, c_bit=256, d_bit=4096;
  uint16_t *rel_brightness=p->rel_brightness->array;

  /* Make an array for the values of the magnitude. */
  for(i=0;i<4;++i)
    magnitude[i]=mag_arr[sorted_vertices[i].index];

  /* Sort these 4 values in ascending order using bubble sort. */
  for(i=0; i<4; ++i)
    for(j=0; j<4; ++j)
      if(magnitude[i] < magnitude[j])
        {
          /* Swap these values. */
          tmp = magnitude[i];
				  magnitude[i] = magnitude[j];
				  magnitude[j] = tmp;
        }

  /* Assign values indexes for relative magnitudes. */
  for(i=0; i<4; ++i)
    for(j=0; j<4; ++j)
      if(mag_arr[sorted_vertices[j].index]==magnitude[i])
        {
          rel_magnitude_index[i]=j;
          break;
        }

  /* For eg, star A is the a-th star and is either of {0, 1, 2, 3}
     0 being the lowest value and 3 being the highest value. */

  /* Shift the bits wrt the relative brightness in the quad. */
  a_bit <<= rel_magnitude_index[0];
  b_bit <<= rel_magnitude_index[1];
  c_bit <<= rel_magnitude_index[2];
  d_bit <<= rel_magnitude_index[3];

  /* After we have 4 16-bits numbers representing the relative
     brigtness of each star, we do a bitiwise-or to join them
     together to give a unique 16-bit number to the quad.
     For eg, if a=3, b=2, c=1, d=0 then the exected output in
     binary system is:
     1 0 0 0   0 1 0 0   0 0 1 0   0 0 0 1
     which is 8737 in decimal system. */
  rel_brightness[qindex] = (a_bit | b_bit | c_bit | d_bit);

  /* For a check:
    printf("rel_brightness = %u\n", rel_brightness[qindex]);
  */
}




/***************making quads*************************/


/* Sort by distance as refrerence in ascending order. */
static int
sort_by_distance(const void *a, const void *b)
{
  struct quad_vertex *p1 = (struct quad_vertex *)a;
  struct quad_vertex *p2 = (struct quad_vertex *)b;

  return ( p1->distance==p2->distance
           ? 0
           : (p1->distance<p2->distance ? -1 : 1) );
}






/* Find indexes of stars A, B, C, D in the given quads.

   After finding the stars w.r.t a given star, we need to
   find the pair that is the most distant and assign them
   as star A and B. To further remove the redumdencies for
   C, D and to firmly ensure the position of A and B, we ensure
   that the following conditions hold true:
                     Cx<=Dx         (1)
                     Cx+Dx<=1       (2)
   The first condition assigns the correct star for C and D
   while the second condition assigns correct stars for  A and B.
*/
void
hash_build_write(size_t qindex, struct quad_vertex* sorted_vertices, struct params *p
                /*TODO: double *cx, double *cy ...*/)
{
  size_t i, j;
  double distance;
  int perm_set[4][4]={0}, c_assigned=0;
  double current_max_dis=DBL_MIN;
  struct quad_vertex temp_vertices[4]={0};
  uint32_t *a_ind=p->a_ind->array;
  uint32_t *b_ind=p->b_ind->array;
  uint32_t *c_ind=p->c_ind->array;
  uint32_t *d_ind=p->d_ind->array;
  double ra[4], dec[4], *ra_arr=p->ra->array, *dec_arr=p->dec->array;

  /* Fill the ra and dec arrays. */
  for(i=0;i<4;++i)
    {
      ra[i]=ra_arr[sorted_vertices[i].index];
      dec[i]=dec_arr[sorted_vertices[i].index];
    }

  /* For every star in the sorted_verices. */
  for(i=0;i<4;++i)
    for(j=0;j<4;++j)
      {
        /* Use a set to remove repeated pairs like {AB, BA} etc.
           Further remove cases where the same stars are used like AA etc.
           */
        if(!perm_set[i][j] && i!=j)
          {
            /* If this combination was unseen, increase its value to 1.*/
            perm_set[i][j]++;
            perm_set[j][i]=perm_set[i][j];

            /* Find the distance between the two stars. */
            distance=( (ra[i]-ra[j])*(ra[i]-ra[j])
                        +(dec[i]-dec[j])*(dec[i]-dec[j]) );

            /* If the calculated distance is greater than the
               current maximum, make the current distance equal to
               current maximum and save the indexes of the current
               stars as the indexes of A and B. */
            if(current_max_dis <= distance)
              {
                /* Make the current distance equal to current maximum */
                current_max_dis=distance;

                /* Assign the indexes of A and B.
                   Note that due to sorted nature of the input
                   `struct quad_vertex` and due to the above conditions,
                   i < j and hence the index of A will always be greater
                   than the index of B. */
                a_ind[qindex]=sorted_vertices[i].index;
                b_ind[qindex]=sorted_vertices[j].index;
              }
          }
      }


  /* Now we have the indexes of star A and B. The remaining vertices are
      randomly assigned as C and D. We will correctly determine the
      star C and D after this. */
  for(i=0;i<4;++i)
    /* If the remaining index is none of those assigned to A or B,
        assign them to C and D. */
    if( sorted_vertices[i].index != a_ind[qindex] &&
        sorted_vertices[i].index != b_ind[qindex] )
      {
        if(!c_assigned)
          {
            c_ind[qindex]=sorted_vertices[i].index;
            c_assigned=1;
          }
        else
          d_ind[qindex]=sorted_vertices[i].index;
      }

  /* We will use a temporary `struct quad_vertex` array to store the
     indexes of the four stars and then make hashes with those stars
     and firmly determine the true positions of C and D. */
  temp_vertices[0].index=a_ind[qindex];
  temp_vertices[1].index=b_ind[qindex];
  temp_vertices[2].index=c_ind[qindex];
  temp_vertices[3].index=d_ind[qindex];

  /* For a check
  printf("A = (%g, %g) B = (%g, %g) C = (%g, %g) D = (%g, %g)\n",
         ra_arr[a_ind[qindex]], dec_arr[a_ind[qindex]], ra_arr[b_ind[qindex]], dec_arr[b_ind[qindex]],
         ra_arr[c_ind[qindex]], dec_arr[c_ind[qindex]], ra_arr[d_ind[qindex]], dec_arr[d_ind[qindex]] );
  */

  /* Make the hash codes with this configuration of stars. */
  hash_geometric_fix_cd(qindex, temp_vertices, p /*TODO: double *cx, double *cy ...*/);

  /* Find the relative brightness of the quad. */
  hash_brightness(qindex, temp_vertices, p);
}








/* Make the quads from the top 5 brightest stars in each box in the
   grid. This is the worker function used by `gal_thread_spin_off` to
   be used on multiple threads to make the quads independently. */
void *
make_quads_worker(void *in_prm)
{
  /* Low-level definitions to be done first. */
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct params *p=(struct params *)tprm->params;

  /* Subsequent definitions. */
  size_t nstars;
  double ref_ra, ref_dec;
  size_t i, j, sid, qindex;
  struct quad_vertex *qvertices=NULL;
  struct quad_vertex good_vertices[4]={0};
  gal_list_sizet_t *good_stars=NULL, *tstar;

  /* For easy reading. */
  double *ra=p->ra->array;
  double *dec=p->dec->array;
  uint32_t *a_ind=p->a_ind->array;
  uint32_t *b_ind=p->b_ind->array;
  uint32_t *c_ind=p->c_ind->array;
  uint32_t *d_ind=p->d_ind->array;
  size_t *bsi=p->brightest_star_id;
  double max_dist=p->max_star_dis_in_quad;
  double *cx=p->cx->array, *cy=p->cy->array;
  double *dx=p->dx->array, *dy=p->dy->array;
  uint16_t *rel_brightness=p->rel_brightness->array;

  /* For polygon region file. Remove after testing. */
  int testing_ds9=1;
  size_t poly_counter=0;
  size_t total_quads=125;
  double *ra_arr=p->ra->array;
  double *dec_arr=p->dec->array;
  double *polygon=malloc(total_quads*4*2*sizeof(*polygon));

  /* Go over all the quads that were assigned to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* Extract this quad's index. */
      qindex = tprm->indexs[i];
      ref_ra = ra[bsi[qindex]];
      ref_dec = dec[bsi[qindex]];

      /* Go over all the stars and keep the ones around this
         quad's first star (same index as brightest_star_id). */
      good_stars=NULL;
      for(j=0;j<p->magnitude->size;++j)
        if (     ra[j]  <= ref_ra+max_dist
              && ra[j]  >= ref_ra-max_dist
              && dec[j] <= ref_dec+max_dist
              && dec[j] >= ref_dec-max_dist )
            gal_list_sizet_add(&good_stars, j);

      /* Count the number of good stars and if they are less
         than 4, set blank values in the columns. */
      nstars=gal_list_sizet_number(good_stars);
      if(nstars < 4)
        {
          cx[qindex]=cy[qindex]=dx[qindex]=dy[qindex]=NAN;
          rel_brightness[qindex]=GAL_BLANK_UINT16;
          a_ind[qindex]=b_ind[qindex]=c_ind[qindex]=d_ind[qindex]=GAL_BLANK_UINT32;
          continue;
        }

      /* Allocate array of qvertex. */
      errno=0;
      qvertices=malloc(nstars*sizeof(*qvertices));
      if(!qvertices)
        error(EXIT_FAILURE, errno, "%s: failed to allocate %zu "
              " bytes for 'qvertices'.", __func__,
              nstars*sizeof(*qvertices));

      /* Loop over the list extract inforamation for sorting. */
      j=0;
      for(tstar=good_stars; tstar!=NULL; tstar=tstar->next)
        {
          /* Assign the value of tstars as Star Id. */
          sid=tstar->v;

          /* Store the relevant details for the possible quad member star.
             We find distances from current star to later sort the stars
             on basis of diatance and form qvertices based on the distance.
             As only the magnitude of distance is required and not the
             absolute value, we use distance square(distance^2 = x*x+y*y)
             rather than using sqrt() to find actual distance. */
          qvertices[j].index=sid;
          qvertices[j].distance=(  (ra[sid]-ref_ra)*(ra[sid]-ref_ra)
                                  +(dec[sid]-ref_dec)*(dec[sid]-ref_dec) );

          /* Increment j. */
          ++j;
        }

      /* Clean the temporary space for the list. */
      gal_list_sizet_free(good_stars);

      /* Sort on the basis of distance and select stars at relative
         positions of n, 0.5 quartile and 0.75 quartile from the current star. */
      qsort(qvertices, nstars, sizeof(struct quad_vertex), sort_by_distance);

      /* 1 -> 23, 51
         2 -> 30, 55
         3 -> 25, 80
         4 -> 35, 90

         So find distance between all the stars in the quad wrt each other
         and select the 2 stars that are most distant as A or B.
         In the above ex the most distant star are 1 and 4. After we define star
         A to be the one which is more nearer to star C and D(or the other 2 stars)
         i.e if dis(1, 2) < dis(1, 3)

      */


      /* Make the temporary quad and assign relative positions
         of sorted stars to it. For temporaay C and D, we use the stars
         at 0.5 and 0.75 quantiles respectively. */
      switch(nstars)
        {
          case 4:
            good_vertices[0] = qvertices[0];         /* A */
            good_vertices[1] = qvertices[1];         /* C */
            good_vertices[2] = qvertices[2];         /* D */
            good_vertices[3] = qvertices[3];         /* B */
          case 5:
            good_vertices[0] = qvertices[0];         /* A */
            good_vertices[1] = qvertices[2];         /* C */
            good_vertices[2] = qvertices[3];         /* D */
            good_vertices[3] = qvertices[4];         /* B */
          default:
            good_vertices[0] = qvertices[0];           /* A */
            good_vertices[1] = qvertices[nstars/2];    /* C */
            good_vertices[2] = qvertices[nstars*3/4];  /* D */
            good_vertices[3] = qvertices[nstars-1];    /* B */
        }

      /* For a test:
      printf("good_vertices = (%g, %g) (%g, %g) (%g, %g) (%g, %g)\n",
             good_vertices[0].ra, good_vertices[0].dec,
             good_vertices[1].ra, good_vertices[1].dec,
             good_vertices[2].ra, good_vertices[2].dec,
             good_vertices[3].ra, good_vertices[3].dec  );
      */

      // /* Write a ds9 region file(on 1 thread only). Remove this section
      //    after testing.
      if (testing_ds9)
        {
          polygon[8*poly_counter+0]=ra_arr [good_vertices[0].index];
          polygon[8*poly_counter+1]=dec_arr[good_vertices[0].index];
          polygon[8*poly_counter+2]=ra_arr [good_vertices[1].index];
          polygon[8*poly_counter+3]=dec_arr[good_vertices[1].index];
          polygon[8*poly_counter+4]=ra_arr [good_vertices[2].index];
          polygon[8*poly_counter+5]=dec_arr[good_vertices[2].index];
          polygon[8*poly_counter+6]=ra_arr [good_vertices[3].index];
          polygon[8*poly_counter+7]=dec_arr[good_vertices[3].index];
          poly_counter++;
        }
      // */

      /* Find indexes of a,b,c,d and calculate and write the hashes. */
      hash_build_write(qindex, good_vertices, p /*TODO: double *cx, double *cy ...*/);

      /* For a check:
      printf("cx = %g, cy = %g, dx = %g, dy = %g\n", cx[qindex], cy[qindex], dx[qindex], dy[qindex]);
      */
    //  exit(0);

      /* Clean up. */
      free(qvertices);
    }
  /* Make region file. Not thread safe. Remove after testing. */
  if(testing_ds9)
    {
      gal_polygon_to_ds9reg("polygon.reg", total_quads, polygon, 4, NULL);
      free(polygon);
    }

  /* Wait for all the other threads to finish, then return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





void
quad_allocate_output(struct params *p, size_t num_quads)
{
  int quitemmap=1;
	size_t minmapsize=-1;

  /* Allocate all the output columns. */
  p->cx=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &num_quads, NULL, 0,
	                    minmapsize, quitemmap, "Cx", "position",
											"relative position in matching coordinates");
  p->cy=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &num_quads, NULL, 0,
	                    minmapsize, quitemmap, "Cy", "position",
											"relative position in matching coordinates");
  p->dx=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &num_quads, NULL, 0,
	                    minmapsize, quitemmap, "Dx", "position",
											"relative position in matching coordinates");
  p->dy=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &num_quads, NULL, 0,
	                    minmapsize, quitemmap, "Dy", "position",
											"relative position in matching coordinates");

  p->rel_brightness=gal_data_alloc(NULL, GAL_TYPE_UINT16, 1, &num_quads, NULL, 0,
	                                minmapsize, quitemmap, "rel-brightness", "none",
													         "relative brightness stored as bit-flags");

  p->a_ind=gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, &num_quads, NULL, 0,
	                        minmapsize, quitemmap, "A-index", "counter",
													"index of star A in the quad");
  p->b_ind=gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, &num_quads, NULL, 0,
	                        minmapsize, quitemmap, "B-index", "counter",
													"index of star B in the quad");
  p->c_ind=gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, &num_quads, NULL, 0,
	                        minmapsize, quitemmap, "C-index", "counter",
													"index of star C in the quad");
  p->d_ind=gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, &num_quads, NULL, 0,
	                        minmapsize, quitemmap, "D-index", "counter",
													"index of star D in the quad");

  /* Define them as a list for the final output. */
  p->cx->next=p->cy;
  p->cy->next=p->dx;
  p->dx->next=p->dy;
  p->dy->next=p->rel_brightness;
  p->rel_brightness->next=p->a_ind;
  p->a_ind->next=p->b_ind;
  p->b_ind->next=p->c_ind;
  p->c_ind->next=p->d_ind;
}




static size_t
high_level_quads(char *reference_name, char *kdtree_name)
{
  size_t kdtree_root;
  size_t num_threads=1;
  size_t x_numbin=5, y_numbin=5, num_stars_per_gpixel=5;

  /* Internal. */
  struct params p={0};
	struct grid in_grid={0};
	size_t num_quads=x_numbin*y_numbin*num_stars_per_gpixel;

  /* Choose columns to read. */
  gal_list_str_t *cols=NULL;
  gal_list_str_add(&cols, "ra", 0);
  gal_list_str_add(&cols, "dec", 0);
  gal_list_str_add(&cols, "phot_g_mean_mag", 0);
  gal_list_str_reverse(&cols);

  /* Read the columns. */
  gal_data_t *ref=gal_table_read (reference_name, "1",
                                  NULL, cols, GAL_TABLE_SEARCH_NAME,
                                  0, -1, 0, NULL);

  /* Seperate columns. */
  p.ra= gal_data_copy_to_new_type (ref, GAL_TYPE_FLOAT64);
  p.dec=gal_data_copy_to_new_type (ref->next, GAL_TYPE_FLOAT64);
  p.magnitude=gal_data_copy_to_new_type (ref->next->next, GAL_TYPE_FLOAT32);

  /* make a box-grid */
	grid_make(p.ra, p.dec, x_numbin, y_numbin, &in_grid);

	/* Find the top 5 star ids using structure. */
  p.brightest_star_id=find_brightest_stars(p.ra, p.dec, p.magnitude, num_quads, &in_grid,
                                          num_stars_per_gpixel);

  /* Allocate output columns. */
  quad_allocate_output(&p, num_quads);

  /* Finally config p and spin-off the threads. */
  p.max_star_dis_in_quad=max(in_grid.step_size[0],
                             in_grid.step_size[1]);
  gal_threads_spin_off(make_quads_worker, &p, num_quads, num_threads);

  /* Construct a tree and fix the column pointers. */
  p.dy->next=NULL;
  p.left=gal_kdtree_create(p.cx, &kdtree_root);
  p.dy->next=p.rel_brightness;
  p.d_ind->next=p.left;
  p.right=p.left->next;

  /* Write the final table with all the quad inforamations
     and the respective kd-tree. */
  gal_table_write(p.cx, NULL, GAL_TABLE_FORMAT_BFITS,
                  kdtree_name, "quad-kdtree", 0);

  /* Clean up and return. */
	gal_list_data_free (ref);
  return kdtree_root;
}



static void
match_prepare(struct params *p,  char *query_name, char *kdtree_name)
{
  gal_data_t *query;
  gal_data_t *ref_quad_kdtree;

  /* Read the query columns. */
  gal_list_str_t *query_cols=NULL;
  gal_list_str_add(&query_cols, "X", 0);
  gal_list_str_add(&query_cols, "Y", 0);
  gal_list_str_add(&query_cols, "Magnitude", 0);
  gal_list_str_reverse(&query_cols);
  query=gal_table_read (query_name, "1", NULL, query_cols,
                        GAL_TABLE_SEARCH_NAME, 0, -1, 0, NULL);
  
  /* Prepare the inputs and make sure only 3 columns are returned. */
  if(gal_list_data_number(query) != 3)
    error(EXIT_FAILURE, 0, "%s: %s should atleast have 3 columns "
          "(x, y, magnitude) but %zu columns matched these names",
          __func__, query_name, gal_list_data_number(query));

  p->y=query->next;
  p->qmagnitude=p->y->next;
  p->x=gal_data_copy_to_new_type_free(query, GAL_TYPE_FLOAT64);
  p->y=gal_data_copy_to_new_type_free(p->y, GAL_TYPE_FLOAT64);
  p->qmagnitude=gal_data_copy_to_new_type_free(p->qmagnitude, GAL_TYPE_FLOAT32);

  /* Read the input kd-tree. */
  ref_quad_kdtree=gal_table_read(kdtree_name, "1", NULL, 
                         NULL, 0, 0, -1, 0, NULL);

  if(gal_list_data_number(ref_quad_kdtree) != 11)
    error(EXIT_FAILURE, 0, "%s: %s should be 11 columns but it is %zu columns",
          __func__, kdtree_name, gal_list_data_number(ref_quad_kdtree));

  p->cx=ref_quad_kdtree;
  if(p->cx->type != GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "%s: %s 1st column should be float64 but it is %s",
          __func__, kdtree_name, gal_type_name(p->cx->type, 1));

  p->cy=p->cx->next;
  if(p->cy->type != GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "%s: %s 2nd column should be float64 but it is %s",
          __func__, kdtree_name, gal_type_name(p->cy->type, 1));
  
  p->dx=p->cy->next;
  if(p->dx->type != GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "%s: %s 3rd column should be float64 but it is %s",
          __func__, kdtree_name, gal_type_name(p->dx->type, 1));

  p->dy=p->dx->next;
  if(p->dy->type != GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "%s: %s 4th column should be float64 but it is %s",
          __func__, kdtree_name, gal_type_name(p->dy->type, 1));

  p->rel_brightness=p->dy->next;
  if(p->rel_brightness->type != GAL_TYPE_UINT16)
    error(EXIT_FAILURE, 0, "%s: %s 5th column should be uint16 but it is %s",
          __func__, kdtree_name, gal_type_name(p->rel_brightness->type, 1));
  
  p->a_ind=p->rel_brightness->next;
  if(p->a_ind->type != GAL_TYPE_UINT32)
    error(EXIT_FAILURE, 0, "%s: %s 6th column should be uint32 but it is %s",
          __func__, kdtree_name, gal_type_name(p->a_ind->type, 1));

  p->b_ind=p->a_ind->next;
  if(p->b_ind->type != GAL_TYPE_UINT32)
    error(EXIT_FAILURE, 0, "%s: %s 7th column should be uint32 but it is %s",
          __func__, kdtree_name, gal_type_name(p->b_ind->type, 1));

  p->c_ind=p->b_ind->next;
  if(p->c_ind->type != GAL_TYPE_UINT32)
    error(EXIT_FAILURE, 0, "%s: %s 8th column should be uint32 but it is %s",
          __func__, kdtree_name, gal_type_name(p->c_ind->type, 1));

  p->d_ind=p->c_ind->next;
  if(p->d_ind->type != GAL_TYPE_UINT32)
    error(EXIT_FAILURE, 0, "%s: %s 9th column should be uint32 but it is %s",
          __func__, kdtree_name, gal_type_name(p->d_ind->type, 1));
  
  p->left=p->d_ind->next;
  if(p->left->type != GAL_TYPE_UINT32)
    error(EXIT_FAILURE, 0, "%s: %s 10th column should be uint32 but it is %s",
          __func__, kdtree_name, gal_type_name(p->left->type, 1));
  
  p->right=p->left->next;
  if(p->right->type != GAL_TYPE_UINT32)
    error(EXIT_FAILURE, 0, "%s: %s 11th column should be uint32 but it is %s",
          __func__, kdtree_name, gal_type_name(p->right->type, 1));
}





static void
high_level_match(char *query_name, char *kdtree_name,
                 char *output_name, size_t kdtree_root)
{
  struct params p={0};
  size_t num_threads=1;
  struct grid in_grid={0};
  size_t x_numbin=5, y_numbin=5, num_stars_per_gpixel=5;
	size_t num_quads=x_numbin*y_numbin*num_stars_per_gpixel;

  /* Read and do sanity checks. */
  match_prepare(&p, query_name, kdtree_name);

  /* make a box-grid */
	grid_make(p.x, p.y, x_numbin, y_numbin, &in_grid);

	/* Find the top 5 star ids using structure. */
  p.brightest_star_id=find_brightest_stars(p.x, p.y, p.magnitude, num_quads, &in_grid,
                                          num_stars_per_gpixel);

  /* Allocate output columns. 
     X Y RA DEC
  */
  match_allocate_output(&p, num_quads);

  /* Finally config p and spin-off the threads. */
  p.max_star_dis_in_quad=max(in_grid.step_size[0],
                             in_grid.step_size[1]);
  gal_threads_spin_off(make_quads_worker, &p, num_quads, num_threads);
}





#if 0
static void
high_level_match(char *query_name, char *ref_quads_name,
                 char *kdtree_name, char *output_name)
{
  gal_data_t *cx, *cy, *dx, *dy;
  gal_list_str_t *query_cols=NULL;
  gal_list_str_add(&query_cols, "Cx", 0);
  gal_list_str_add(&query_cols, "Cy", 0);
  gal_list_str_add(&query_cols, "Dx", 0);
  gal_list_str_add(&query_cols, "Dy", 0);
  gal_list_str_reverse(&query_cols);

  /* Read the columns. */
  gal_data_t *query_table=gal_table_read (query_name, "1",
                                  NULL, query_cols, GAL_TABLE_SEARCH_NAME,
                                  0, -1, 0, NULL);

  /* Seperate columns. */
  cx=gal_data_copy_to_new_type (query_table, GAL_TYPE_FLOAT64);
  cy=gal_data_copy_to_new_type (query_table->next, GAL_TYPE_FLOAT64);
  dx=gal_data_copy_to_new_type (query_table->next->next, GAL_TYPE_FLOAT64);
  dy=gal_data_copy_to_new_type (query_table->next->next->next, GAL_TYPE_FLOAT64);

  double *cx_val=cx->array;
  double *cy_val=cy->array;
  double *dx_val=dx->array;
  double *dy_val=dy->array;

  size_t nearest_index;

  /* Find the nearest neighbour of the point. */
  for(size_t i=0;i<query_table->size; ++i)
  {
    double point[4]={cx_val[i], cy_val[i], dx_val[i], dy_val[i]};
    nearest_index=gal_kdtree_nearest_neighbour(coords, output, root, point);

    /* For a check */
    double *a=coords->array;
    double *b=coords->next->array;
    double *c=coords->next->next->array;
    double *d=coords->next->next->next->array;
    printf("i=%zu:(%g, %g, %g, %g) nearest_index=%zu->(%g, %g, %g, %g)\n\n",
            i, point[0], point[1], point[2], point[3],
            nearest_index, a[nearest_index], b[nearest_index], c[nearest_index], d[nearest_index]);
  }
  
}
#endif



/***********FILE STRUCTURE**************/
/* 
  RAW INPUTS
  ==========
  [i0] reference.fits -> ra, dec, magnitude 
  [i1] query.fits     -> x, y, magnitude

  REFERENCE QUAD-HASHES
  =====================
  [b0]<-[i0] ref-quads.fits (HDU=1) -> cx, cy, dx, dy
  [b1]<-[i0] ref-quads.fits (HDU=2) -> rel_brightness

  KD-TREE
  =======
  [b2]<-[b0] kdtree.fits  -> left, right


  OUTPUT
  ======
  [o0]<-[b0],[b1],[b2] matched-out.fits -> x, y, ra, dec
*/

int main()
{
  size_t kdtree_root;
  /* High-level file names. */
  char *query_name="./input/query.fits";
  char *reference_name="./input/reference.fits";
  char *output_name = "./build/matched-out.fits";

  char *kdtree_name="./build/kdtree.fits";

  /* Quad-hashes creation */
  kdtree_root=high_level_quads(reference_name, kdtree_name);

  /* Find quads on the query image and match them. */
  high_level_match(query_name, kdtree_name, output_name, kdtree_root);

  return EXIT_SUCCESS;
}