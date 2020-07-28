#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <error.h>

#include <gnuastro/qsort.h>
#include <gnuastro/table.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
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



/* For internal use only. */
struct quad_vertex
{
  size_t index;
  double ra, dec;
  double distance;
  float  brightness;
};





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
make_hash_codes(struct quad_vertex* sorted_vertices, double *cx, double *cy, 
                double *dx, double *dy, uint16_t *rel_brightness)
{
  double scale=0;
  double angle_cab=0, angle_dab=0;
  double mag_ab=0, mag_ac=0, mag_ad=0;


  /* First calculate the angles and add and subtract 45 degrees to make it
     in frame with the current coordinate system with AB as its angle
     bisector. */
  angle_cab=find_angle_aob(sorted_vertices[1].ra, sorted_vertices[1].dec,
                          sorted_vertices[0].ra, sorted_vertices[0].dec,
                          sorted_vertices[3].ra, sorted_vertices[3].dec);


  angle_dab=find_angle_aob(sorted_vertices[2].ra, sorted_vertices[2].dec,
                          sorted_vertices[0].ra, sorted_vertices[0].dec,
                          sorted_vertices[3].ra, sorted_vertices[3].dec);

  /* For a check:
  printf("angle_cab = %g, angle_dab = %g\n", angle_cab, angle_dab);
  */

  /* Make a scaled-down value for the new coordinate system. We want that
                |AB|=1, and hence scale = 1/|AB|
     where |AB|=sqrt(|AB.ra|^2+|AB.dec|^2).
     This value can be used to scale-down the other distaces mainly,
     AC and AD. */
  mag_ab=sqrt((sorted_vertices[0].ra-sorted_vertices[3].ra)
              *(sorted_vertices[0].ra-sorted_vertices[3].ra)
              +(sorted_vertices[0].dec-sorted_vertices[3].dec)
              *(sorted_vertices[0].dec-sorted_vertices[3].dec));

  scale=1/mag_ab;


  /* Find |AC| and |AD|. Then use distance to find cos and sin thetas to
     give ra and dec coordinates. */
  mag_ac=sqrt((sorted_vertices[1].ra-sorted_vertices[3].ra)
              *(sorted_vertices[1].ra-sorted_vertices[3].ra)
              +(sorted_vertices[1].dec-sorted_vertices[3].dec)
              *(sorted_vertices[1].dec-sorted_vertices[3].dec));

  mag_ad=sqrt((sorted_vertices[2].ra-sorted_vertices[3].ra)
              *(sorted_vertices[2].ra-sorted_vertices[3].ra)
              +(sorted_vertices[2].dec-sorted_vertices[3].dec)
              *(sorted_vertices[2].dec-sorted_vertices[3].dec));


  /* Make the final hash-code. */
  *cx=mag_ac*scale*cos(45+angle_cab); /* Cx */
  *cy=mag_ac*scale*sin(45+angle_cab); /* Cy */
  *dx=mag_ad*scale*cos(45-angle_dab); /* Dx */
  *dy=mag_ad*scale*sin(45-angle_dab); /* Dy */

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



/* Sort by brghtness as refrerence in ascending order. */
static int
sort_by_brightness(const void *a, const void *b)
{
  struct quad_vertex *p1 = (struct quad_vertex *)a;
  struct quad_vertex *p2 = (struct quad_vertex *)b;

  return ( p1->brightness==p2->brightness
           ? 0
           : (p1->brightness<p2->brightness ? -1 : 1) );
}




/* Structure to keep parameters for multithreaded worker function. */
struct params
{
  /* Inputs. */
  gal_data_t *ra;
  gal_data_t *dec;
  gal_data_t *magnitude;

  /* Internal. */
  double max_star_dis_in_quad;
  size_t *brightest_star_id;

  /* Output. */
  gal_data_t *cx, *cy;
  gal_data_t *dx, *dy;
  gal_data_t *rel_brightness;
  gal_data_t *a_ind, *b_ind, *c_ind, *d_ind;
  gal_data_t *left, *right;
};





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

  /* For easy reading*/
  double *ra=p->ra->array;
  double *dec=p->dec->array;
  size_t *a_ind=p->a_ind->array;
  size_t *b_ind=p->b_ind->array;
  size_t *c_ind=p->c_ind->array;
  size_t *d_ind=p->d_ind->array;
  size_t *bsi=p->brightest_star_id;
  float  *magnitude=p->magnitude->array;  
  double max_dist=p->max_star_dis_in_quad;
  double *cx=p->cx->array, *cy=p->cy->array;
  double *dx=p->dx->array, *dy=p->dy->array;
  uint16_t *rel_brightness=p->rel_brightness->array;

  /* Go over all the quads that were assigned to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* Extract this quad's index. */
      qindex = tprm->indexs[i];
      ref_ra = ra[bsi[qindex]];
      ref_dec = dec[bsi[qindex]];
      // printf("qindex = %zu\n", qindex);

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
          a_ind[qindex]=b_ind[qindex]=c_ind[qindex]=d_ind[qindex]=GAL_BLANK_SIZE_T;
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
          /* We find distances from current star to later sort the stars
             on basis of diatance and form qvertices based on the distance.
             As only the magnitude of distance is required and not the
             absolute value, we use distance square(distance^2 = x*x+y*y)
             rather than using sqrt() to find actual distance. */
          sid=tstar->v;

          /* Store the relevant details for the possible quad member star. */
          qvertices[j].index=sid;
          qvertices[j].ra=ra[sid];
          qvertices[j].dec=dec[sid];
          qvertices[j].brightness=magnitude[sid];
          qvertices[j].distance=(  (ra[sid]-ref_ra)*(ra[sid]-ref_ra)
                                  +(dec[sid]-ref_dec)*(dec[sid]-ref_dec) );

          /* Increment j*/
          ++j;
        }

      /* Clean the temporary space for the list. */
      gal_list_sizet_free(good_stars);

      /* Sort on the basis of distance and select stars at relative
         positions of n, n-2 and n-4 from the current star. */
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
      /* Make the final quad and assign relative positions
         of sorted stars to it. */
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
            good_vertices[0] = qvertices[0];         /* A */
            good_vertices[1] = qvertices[nstars-5];  /* C */
            good_vertices[2] = qvertices[nstars-3];  /* D */
            good_vertices[3] = qvertices[nstars-1];  /* B */
        }
      /* Write a ds9 region file(on 1 thread only). */

      /* Check if C and D are within the circle of radius 0.5. 
         Use middle point of AB as centre and check the distance.*/

      /* Sort according to brightness. */
      qsort(good_vertices, 4, sizeof(struct quad_vertex), sort_by_brightness);

      /* Run make hash function and fill the outputs columns for this quad. 
       
         void 
         make_hash_codes(struct quad_vertex* sorted_vertices, double *cx, double *cy, 
                         double *dx, double *dy, uint16_t *rel_brightness)
        */
      make_hash_codes(good_vertices, &cx[qindex], &cy[qindex], &dx[qindex], &dy[qindex],
                      &rel_brightness[qindex]);

      /*
      printf("cx = %g, cy = %g, dx = %g, dy = %g\n",  cx[qindex], cy[qindex], dx[qindex], dy[qindex]);
      */

      a_ind[qindex]=good_vertices[0].index;
      b_ind[qindex]=good_vertices[3].index;
      c_ind[qindex]=good_vertices[1].index;
      d_ind[qindex]=good_vertices[2].index;

      /* Clean up. */
      free(qvertices);
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

  p->a_ind=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &num_quads, NULL, 0,
	                        minmapsize, quitemmap, "A-index", "counter",
													"index of star A in the quad");
  p->b_ind=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &num_quads, NULL, 0,
	                        minmapsize, quitemmap, "B-index", "counter",
													"index of star B in the quad");
  p->c_ind=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &num_quads, NULL, 0,
	                        minmapsize, quitemmap, "C-index", "counter",
													"index of star C in the quad");
  p->d_ind=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &num_quads, NULL, 0,
	                        minmapsize, quitemmap, "D-index", "counter",
													"index of star D in the quad");

  p->left=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &num_quads, NULL, 0,
	                       minmapsize, quitemmap, "KD-left-index", "counter",
												 "index of left subtree");  
  p->right=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &num_quads, NULL, 0,
	                       minmapsize, quitemmap, "KD-right-index", "counter",
												 "index of right subtree");  

  /* Define them as a list for the final output. */
  p->cx->next=p->cy;
  p->cy->next=p->dx;
  p->dx->next=p->dy;
  p->dy->next=p->rel_brightness;
  p->rel_brightness->next=p->a_ind;
  p->a_ind->next=p->b_ind;
  p->b_ind->next=p->c_ind;
  p->c_ind->next=p->d_ind;
  p->d_ind->next=p->left;
  p->left->next=p->right;
}





int main()
{
  /* Input arguments. */
	size_t x_numbin=5, y_numbin=5, num_stars_per_gpixel=5;
  size_t num_threads=4;
  char *kd_tree_outname="kd-tree-output.fits";


  /* Internal. */
  struct params p;
	// size_t i,j, index;
	struct grid in_grid={0};
  // size_t num_gpix=x_numbin*y_numbin;
	size_t num_quads=x_numbin*y_numbin*num_stars_per_gpixel;

	/* Choose columns to read. */
  gal_list_str_t *cols=NULL;
  gal_list_str_add(&cols, "ra", 0);
  gal_list_str_add(&cols, "dec", 0);
  gal_list_str_add(&cols, "phot_g_mean_mag", 0);
  gal_list_str_reverse(&cols);

  /* Read the columns. */
  gal_data_t *ref=gal_table_read ("gaia-in-img.fits", "1",
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

  /* Write the quad calculation into a file. */
  p.dy->next=NULL;
  gal_table_write(p.cx, NULL, GAL_TABLE_FORMAT_BFITS,
                  kd_tree_outname, "quad-info.fits", 0);


  /* Calculate kd-tree. 
     struct kd_node
     {
       size_t coord_index; //index of this node in the input coordinate table
       struct kd_node *left, *right;
     }

     void
     kdtree_make(gal_data_t *coordinates, size_t ndim)
     {
       // Read 'quad-info.fits' into memory.(like ref)

       // Create a kd_node for every row of ref table.

       // Continue with old-implementation(instead of 4 values, use index to read from the input table.)
     }

     cx cy dx  dy           index
     1   1  1  1              0
     1   1  1  1              1
     1   1  1  1              2   
  */


  /* Write all information to a file. */
  // p.dy->next=p.rel_brightness;
  // gal_table_write(p.cx, NULL, GAL_TABLE_FORMAT_BFITS,
  //                 kd_tree_outname, "kd-tree.fits", 0);

  /* Clean up. */
	gal_list_data_free (ref);
}