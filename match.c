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










/***********************************************************/
/********          Basic macros/structures          ********/
/***********************************************************/

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



/* To sort quads. */
struct quad_vertex
{
  size_t index;
  double distance;
};



/* To keep reference of matched points. */
struct matched_points
{
  size_t qindex;
  size_t rindex;
  struct matched_points *next;
};



/* Structure to keep parameters for multithreaded worker function. */
struct params
{
  /* Inputs. */
  gal_data_t *x, *y;
  gal_data_t *ra, *dec;
  gal_data_t *r_mag, *q_mag;
  float max_mag_diff;
  char *ds9regprefix;

  /* Internal info for building quads. */
  size_t *brightest_star_id;
  double max_star_dis_in_quad;

  /* Reference table quad information. */
  size_t kdtree_root;
  gal_data_t *cx, *cy;
  gal_data_t *dx, *dy;
  gal_data_t *left, *right;
  gal_data_t *rel_brightness;
  gal_data_t *a_ind, *b_ind, *c_ind, *d_ind;

  /* These will either point to the respective query or reference
     catalog depending on the caller. They are defined for avoiding
     complexity (they are set once by the top function, and the proper
     arrays will be used in all steps afterwards). */
  gal_data_t *mag;
  gal_data_t *c1, *c2;

  /* Finally matched info. */
  struct matched_points **matched;
};




















/***********************************************************/
/**************          DS9 regions         ***************/
/***********************************************************/

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





static FILE *
match_quad_as_ds9_reg_init(char *regname, char *regprefix, size_t id)
{
  FILE *regfile;

  /* Set the filename (with thread number). */
  sprintf(regname, "%s-%zu.reg", regprefix, id);

  /* Open the file. */
  errno=0;
  regfile = fopen(regname, "w");
  if(regfile==NULL) error(EXIT_FAILURE, errno, "%s", regname);

  /* Write basic information. */
  fprintf(regfile, "# Region file format: DS9 version 4.1\n");
  fprintf(regfile, "fk5\n");

  /* Return the file pointer. */
  return regfile;
}





/* Add the given polygon into the given file. */
static void
match_quad_as_ds9_reg(FILE *regfile, double *c1, double *c2, size_t qindex,
		      struct quad_vertex *good_vertices)
{
  size_t ordinds[8];
  double ds9polygon[8]={
    c1[good_vertices[0].index], c2[good_vertices[0].index],
    c1[good_vertices[1].index], c2[good_vertices[1].index],
    c1[good_vertices[2].index], c2[good_vertices[2].index],
    c1[good_vertices[3].index], c2[good_vertices[3].index]};
  gal_polygon_vertices_sort(ds9polygon, 4, ordinds);
  fprintf(regfile, "polygon(%g,%g,%g,%g,%g,%g,%g,%g) # text={%zu}\n",
	  ds9polygon[ordinds[0]*2], ds9polygon[ordinds[0]*2+1],
	  ds9polygon[ordinds[1]*2], ds9polygon[ordinds[1]*2+1],
	  ds9polygon[ordinds[2]*2], ds9polygon[ordinds[2]*2+1],
	  ds9polygon[ordinds[3]*2], ds9polygon[ordinds[3]*2+1],
	  qindex);
}



















/***********************************************************/
/**************            Grid              ***************/
/***********************************************************/

/* Make a grid of boxes on the given list of coordinates. */
static void
grid_make(gal_data_t *c1, gal_data_t *c2, size_t c1_numbins,
          size_t c2_numbins, struct grid *grid)
{
  double c1min, c1max, c2min, c2max;
  gal_data_t *c1_min=NULL, *c1_max=NULL, *c2_min=NULL, *c2_max=NULL;

  /* Assign the dimensions of the grid. */
  grid->dim[0]=c1_numbins;
  grid->dim[1]=c2_numbins;

  /* The size of the grid(=dim1xdim2). */
  grid->size=grid->dim[0]*grid->dim[1];

  /* Extract the range of values. */
  c1_min=gal_statistics_minimum(c1);
  c1_max=gal_statistics_maximum(c1);
  c2_min=gal_statistics_minimum(c2);
  c2_max=gal_statistics_maximum(c2);
  c1min=((double *)(c1_min->array))[0];
  c1max=((double *)(c1_max->array))[0];
  c2min=((double *)(c2_min->array))[0];
  c2max=((double *)(c2_max->array))[0];

  /* Assign the minimum RA value and the size of steps=(max-min/no of grid)
     to reach the last grid from the first step. */
  grid->min[0]=c1min;
  grid->min[1]=c2min;
  grid->step_size[0]=(c1max-c1min)/c1_numbins;
  grid->step_size[1]=(c2max-c2min)/c2_numbins;

  /* Clean up. */
  gal_data_free(c1_min);
  gal_data_free(c1_max);
  gal_data_free(c2_max);
  gal_data_free(c2_min);
}







/* Return the index of the grid box. */
static size_t
grid_findindex(double X, double Y, struct grid *grid)
{
  size_t x=(X-grid->min[0])/grid->step_size[0];
  size_t y=(Y-grid->min[1])/grid->step_size[1];

  /* If the value is maximum in ra and dec, the division returns the
     next integer instead of previous one. This only occurs for edge
     values. Hence only for the edge values x and y are checked and
     shifted to the last index rather than next to last. */
  x = x==grid->dim[0] ? grid->dim[0]-1 : x;
  y = y==grid->dim[1] ? grid->dim[1]-1 : y;

  /* Return the index of the grid-element. */
  return y*grid->dim[0]+x;
}




















/***********************************************************/
/*********        Selecting brightest stars        *********/
/***********************************************************/

/* Select the brightest stars in each grid element. These will be used
   as a starting point for constructing quads. */
static size_t *
find_brightest_stars(gal_data_t *x_data, gal_data_t *y_data,
                     gal_data_t *mag_data, size_t num_quads,
                     struct grid *in_grid, size_t num_in_gpixel)
{
  size_t i, g;
  float *mag=mag_data->array;
  size_t *brightest_star_id, *bsi_counter;
  double *x=x_data->array, *y=y_data->array;
  size_t *sorted_id, npoints=mag_data->size;

  /* Allocate the necessary arrays. */
  sorted_id=gal_pointer_allocate (GAL_TYPE_SIZE_T, npoints, 0,
                                  __func__, "sorted_id");
  bsi_counter=gal_pointer_allocate(GAL_TYPE_SIZE_T, in_grid->size, 1,
                                   __func__, "bsi_counter");
  brightest_star_id=gal_pointer_allocate(GAL_TYPE_SIZE_T, num_quads, 0,
                                         __func__, "brightest_star_id");

  /* Initialise output and sorted_id arrays. */
  for(i=0;i<npoints; ++i) sorted_id[i]=i;
  for(i=0;i<num_quads;++i) brightest_star_id[i]=GAL_BLANK_SIZE_T;

  /* Set magnitude column as a reference to sort stars ID array in
     incresing order (note that the "magnitude" in astronomy is the
     inverse of brightness: as magnitude increases, the object becomes
     fainter).  */
  gal_qsort_index_single=mag;
  qsort(sorted_id, npoints, sizeof(size_t),
        gal_qsort_index_single_float32_i);

  /* Parse through the 'sorted_id' array (where the brightest star is
     the first), and put each star in the respective grid element's
     set of stars, until the capacity is reached. */
  for (i=0;i<npoints;++i)
    {
      /* Index of the grid box corresponding to particular RA and DEC
         values in the sorted order. */
      g = grid_findindex(x[sorted_id[i]], y[sorted_id[i]], in_grid);

      /* If there are less number of star ids for the box than
         required, find more stars in that particular box. */
      if (bsi_counter[g] < num_in_gpixel)
        {
          brightest_star_id[g*num_in_gpixel+bsi_counter[g]]=sorted_id[i];
          bsi_counter[g]++;
        }
    }

  /* Clean and return. */
  free(sorted_id);
  free(bsi_counter);
  return brightest_star_id;
}




















/***********************************************************/
/*********              Quad hashes                *********/
/***********************************************************/
/* Return the angle in range [0, 360) degrees between
   points A, O and B where O is the vertex of angle AOB. */
static double
find_angle_aob(double *a, double *o, double *b)
{
  double dir_o_to_a=0, dir_o_to_b=0, angle_aob=0;

  dir_o_to_a = atan2(a[1] - o[1], a[0] - o[0]);
  dir_o_to_b = atan2(b[1] - o[1], b[0] - o[0]);

  angle_aob = dir_o_to_a - dir_o_to_b;

  /* atan2 returns anngle between [-pi, +pi] radians.
     Convert it [0, 360) degrees. */
  angle_aob = RAD2DEG(angle_aob);

  return angle_aob>=0 ? angle_aob : (angle_aob+360);
}





/* An initial separation of the quad vertices was done before such
   that before calling this function, we know the two vertices that
   are most distant from each other, and the two that are closer. The
   job of this function is to uniquely identify A and B from the two
   that are farther, and uniquely identify C and D from the two that
   are closer. */
static void
hash_geometric_finalize(struct params *p, size_t *abcd, double *geohash)
{
  size_t tmpind;
  double tmpdbl, angle_c0d, angle_c1d;
  double *c1_arr=p->c1->array, *c2_arr=p->c2->array;

  /* Set the four points for easy reading. recall that before (in
     'hash_build_write'), we separated A (index 0 in 'abcd') and B
     (index 1) from C (index 2) and D (index 3). Just note that we
     have still not distinguished between A & B or C & D. */
  double a[2]={c1_arr[ abcd[0] ], c2_arr[ abcd[0] ]};
  double b[2]={c1_arr[ abcd[1] ], c2_arr[ abcd[1] ]};
  double c[2]={c1_arr[ abcd[2] ], c2_arr[ abcd[2] ]};
  double d[2]={c1_arr[ abcd[3] ], c2_arr[ abcd[3] ]};

  /* For a check:
  printf("\n============================================================="
	 "\nInitial (A&B separated from C&D, but not yet from each other)"
	 "\n=============================================================\n");
  printf("A: %-8zu -> (%g, %g)\n", abcd[0], a[0], a[1]);
  printf("B: %-8zu -> (%g, %g)\n", abcd[1], b[0], b[1]);
  printf("C: %-8zu -> (%g, %g)\n", abcd[2], c[0], c[1]);
  printf("D: %-8zu -> (%g, %g)\n", abcd[3], d[0], d[1]);
  */

  /* First we need to uniquely identify A from the pair A or B. We'll
     define A to be the one where the angle from C to D through that
     point is less. So if in the the candidate A&B, CBD angle is less
     than CAD, we should flip the two.. */
  angle_c0d=find_angle_aob(c, a, d);
  angle_c1d=find_angle_aob(c, b, d);
  if(angle_c0d > angle_c1d)
    {
      tmpind=abcd[0];    abcd[0]=abcd[1];     abcd[1]=tmpind;
      tmpdbl=a[0];       a[0]=b[0];           b[0]=tmpdbl;
      tmpdbl=a[1];       a[1]=b[1];           b[1]=tmpdbl;
    }

  /* For a check:
  printf("CAD angle: %.3f (deg)\nCBD angle: %.3f (deg)\n", angle_c0d, angle_c1d);
  if(angle_c0d > angle_c1d)
    printf("+++++++ A & B SWAPPED!\n\n");
  printf("\n=========================================================="
	 "\n A is now uniquely identified (based on smaller angle CXD)"
	 "\n==========================================================\n");
  printf("A: %-8zu -> (%g, %g)\n", abcd[0], a[0], a[1]);
  printf("B: %-8zu -> (%g, %g)\n", abcd[1], b[0], b[1]);
  printf("C: %-8zu -> (%g, %g)\n", abcd[2], c[0], c[1]);
  printf("D: %-8zu -> (%g, %g)\n\n", abcd[3], d[0], d[1]);
  printf("\n==================="
	 "\n Scaling parameters"
	 "\n===================\n");
  printf("Ax-Bx: %g\n", b[0]-a[0]);
  printf("Ay-By: %g\n", b[1]-a[1]);
  */

  /* Now that we uniquely know which point is A and which point is B,
     we need to scale C and D to be in a coordinate system where A is
     on (0,0) and B is on (1,1). Infact the Cx, Cy, Dx, Dy values are
     the ultimate hashes that we want. But we aren't finished yet! We
     still haven't uniquely identified C and D.*/
  geohash[0]=(c[0]-a[0])/(b[0]-a[0]);
  geohash[1]=(c[1]-a[1])/(b[1]-a[1]);
  geohash[2]=(d[0]-a[0])/(b[0]-a[0]);
  geohash[3]=(d[1]-a[1])/(b[1]-a[1]);

  /* We will define C as the point that has a smaller zero-th
     dimension length. So if Cx>Dx, then we'll flip the hashes of C
     and D (and their index in 'abcd'). */
  if(geohash[0]>geohash[2])
    {
      tmpind=abcd[2];     abcd[2]=abcd[3];        abcd[3]=tmpind;
      tmpdbl=geohash[0];  geohash[0]=geohash[2];  geohash[2]=tmpdbl;
      tmpdbl=geohash[1];  geohash[1]=geohash[3];  geohash[3]=tmpdbl;
    }

  /* For a check, note that this is happening after the actual
     flipping on the final hashes, so we need to use the opposite
     values.
  if(geohash[2]>geohash[0])
    {
      printf("Scaled-Cx=%g\n", geohash[2]);
      printf("Scaled-Dx=%g\n", geohash[0]);
      tmpdbl=c[0];   c[0]=d[0];   d[0]=tmpdbl;
      tmpdbl=c[1];   c[1]=d[1];   d[1]=tmpdbl;
      printf("+++++++ C & D SWAPPED!\n\n");
    }
  printf("\n=========================================================="
	 "\n FINAL A, B, C & D, "
	 "\n The coordinates of C and D when A is (0,0) and B is (1,1)"
	 "\n==========================================================\n");
  printf("A: %-8zu -> (%g, %g)\n", abcd[0], a[0], a[1]);
  printf("B: %-8zu -> (%g, %g)\n", abcd[1], b[0], b[1]);
  printf("C: %-8zu -> (%g, %g)\n", abcd[2], c[0], c[1]);
  printf("D: %-8zu -> (%g, %g)\n", abcd[3], d[0], d[1]);
  printf("Cx,Cy: %g, %g\n", geohash[0], geohash[1]);
  printf("Dx,Dy: %g, %g\n", geohash[2], geohash[3]);
  exit(0);
  */
}





/* Find the relative brightness of the stars in the quad and make a
   unique number to represtent the relative brightness of the quad to
   make it further unique while detection. */
static uint16_t
hash_brightness(struct params *p, size_t *abcd)
{
  size_t i, j;
  float tmpflt;
  float magnitude[4];
  uint16_t rel_brightness;
  size_t rel_mag_index[4];
  float *mag_arr=p->mag->array;

  /* Initialise to the nearest power of 2 such that all previous
     multiple of 4 bits are 0.*/
  uint16_t a_bit=1, b_bit=16, c_bit=256, d_bit=4096;

  /* Make an array for the values of the magnitude and their index. */
  for(i=0;i<4;++i)
    {
      rel_mag_index[i]=i;
      magnitude[i]=mag_arr[ abcd[i] ];
    }

  /* For a check:
  printf("\n========================"
	 "\n Magnitudes of A,B,C & D."
	 "\n========================\n");
  printf("A: %g\n", magnitude[0]);
  printf("B: %g\n", magnitude[1]);
  printf("C: %g\n", magnitude[2]);
  printf("D: %g\n", magnitude[3]);
  */

  /* Sort these 4 values in ascending order using bubble sort. Its
     only four points, so its not worth calling more complex
     systems.*/
  for(i=0; i<4; ++i)
    for(j=0; j<4; ++j)
      if(magnitude[i] < magnitude[j])
        {
          /* Swap the values. */
          tmpflt = magnitude[i];
          magnitude[i] = magnitude[j];
          magnitude[j] = tmpflt;
        }

  /* Assign values indexes for relative magnitudes. */
  for(i=0; i<4; ++i)
    for(j=0; j<4; ++j)
      if(mag_arr[ abcd[j] ]==magnitude[i])
        {
          rel_mag_index[i]=j;
          break;
        }

  /* For eg, star A is the a-th star and is either of {0, 1, 2, 3} 0
     being the lowest value and 3 being the highest value. Shift the
     bits wrt the relative brightness in the quad. */
  a_bit <<= rel_mag_index[0];
  b_bit <<= rel_mag_index[1];
  c_bit <<= rel_mag_index[2];
  d_bit <<= rel_mag_index[3];

  /* After we have 4 16-bits numbers representing the relative
     brigtness of each star, we do a bitiwise-or to join them
     together to give a unique 16-bit number to the quad.
     For eg, if a=3, b=2, c=1, d=0 then the exected output in
     binary system is:
     1 0 0 0   0 1 0 0   0 0 1 0   0 0 0 1
     which is 8737 in decimal system. */
  rel_brightness = (a_bit | b_bit | c_bit | d_bit);

  /* For a check:
  printf("\n============================================================"
	 "\n Magnitudes of A,B,C & D"
	 "\n (be careful with the bit-sequence in little-endian systems)"
	 "\n============================================================\n");
  printf("A: %zu: %s\n", rel_mag_index[0], gal_type_bit_string(&a_bit, 2));
  printf("B: %zu: %s\n", rel_mag_index[1], gal_type_bit_string(&b_bit, 2));
  printf("C: %zu: %s\n", rel_mag_index[2], gal_type_bit_string(&c_bit, 2));
  printf("D: %zu: %s\n", rel_mag_index[3], gal_type_bit_string(&d_bit, 2));
  printf("Final: %u: %s\n", rel_brightness,
	 gal_type_bit_string(&rel_brightness, 2));
  */

  /* Return the 16-bit integer. */
  return rel_brightness;
}





/* Find indexes of stars A, B, C, D (geometrically defined) in the
   given quads.

               ^
               |      C------- B
               |     /        /
               |    /        /
               |   A--------D
               |-------------------->

   After finding the vertices of the polygon, we need to find the pair
   that is the most distant. In a separate function we will then
   uniquely identify which one is A and which one is B and of the
   other two points, which one is uniquely C and D. */
static uint16_t
hash_build_write(struct params *p, size_t qindex,
                 struct quad_vertex* sorted_vertices, size_t *abcd,
		 double *geohash)
{
  size_t i, j;
  uint16_t rel_brightness;
  double c1[4], c2[4], distance;
  double current_max_dis=DBL_MIN;
  int perm_set[4][4]={0}, c_assigned=0;
  double *c1_arr=p->c1->array, *c2_arr=p->c2->array;

  /* Fill the polygon positions. */
  for(i=0;i<4;++i)
    {
      c1[i]=c1_arr[sorted_vertices[i].index];
      c2[i]=c2_arr[sorted_vertices[i].index];
    }

  /* Find the stars that are most distant from each other (A & B). */
  for(i=0;i<4;++i)
    for(j=0;j<4;++j)
      {
        /* Use a set to remove repeated pairs like {AB, BA} etc.
           Further remove cases where the same stars are used like AA etc.*/
        if(!perm_set[i][j] && i!=j)
          {
            /* If this combination was unseen, increase its value to 1.*/
            perm_set[i][j]++;
            perm_set[j][i]=perm_set[i][j];

            /* Find the distance between the two stars. */
            distance=( ( c1[i]-c1[j])*(c1[i]-c1[j])
                       +(c2[i]-c2[j])*(c2[i]-c2[j]) );

            /* If the calculated distance is greater than the
               current maximum, make the current distance equal to
               current maximum and save the indexes of the current
               stars as the indexes of A and B. */
            if(current_max_dis <= distance)
              {
                /* Make the current distance equal to current maximum */
                current_max_dis=distance;

                /* Assign the indexes of A and B. */
                abcd[0]=sorted_vertices[i].index;
                abcd[1]=sorted_vertices[j].index;
              }
          }
      }


  /* Now that we have the indexes of star A and B. The remaining
     vertices are simply initialized as C and D in the order that they
     appear. But THIS ISN'T THE FINAL ASSIGNMENT, we will finalize
     them after this.*/
  for(i=0;i<4;++i)
    if( sorted_vertices[i].index != abcd[0]
        && sorted_vertices[i].index != abcd[1] )
      {
        if(!c_assigned)
          {
            abcd[2]=sorted_vertices[i].index;
            c_assigned=1;
          }
        else
          abcd[3]=sorted_vertices[i].index;
      }

  /* Make the hash codes with this configuration of stars. */
  hash_geometric_finalize(p, abcd, geohash);

  /* Find the quad's relative brightness value. */
  rel_brightness=hash_brightness(p, abcd);

  /* Return the relative brightness. */
  return rel_brightness;
}





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




/* Build the quad structure of the desired star. */
static int
quad_from_point(struct params *p, size_t qindex,
		struct quad_vertex *good_vertices)
{
  float ref_mag;
  double ref_c1, ref_c2;
  struct quad_vertex *qvertices=NULL;
  gal_list_sizet_t *good_stars=NULL, *tstar;
  size_t i, sid, ngood, npoints = p->c1->size;

  /* For easy reading. */
  size_t *bsi=p->brightest_star_id;
  double max_dist=p->max_star_dis_in_quad;

  /* To properly deal with the reference or query catalogs. */
  double *c1 = p->c1->array;
  double *c2 = p->c2->array;
  float *mag = p->mag->array;

  /* Get the properties of the bright star used to construct the
     quad. */
  ref_c1  = c1[  bsi[qindex] ];
  ref_c2  = c2[  bsi[qindex] ];
  ref_mag = mag[ bsi[qindex] ];

  /* Go over all the stars and keep a list of candidates to constuct
     the quad. The candidates have to be within a certain spatial
     range, and have to be fainter than the "brightest_star", but not
     too fainter, only by 'max_mag_diff' (usually 2 or 3). */
  good_stars=NULL;
  for(i=0;i<npoints;++i)
    if ( c1[i]  <= ref_c1+max_dist
	 && c1[i]  >= ref_c1-max_dist
	 && c2[i]  <= ref_c2+max_dist
	 && c2[i]  >= ref_c2-max_dist
	 && mag[i] >  ref_mag
	 && mag[i] <= ref_mag+p->max_mag_diff )
      gal_list_sizet_add(&good_stars, i);

  /* Count the number of good stars and if they are less than 4,
     set blank values in the respective output columns and go to
     the next quad. */
  ngood=gal_list_sizet_number(good_stars);
  if(ngood < 4) return 0;

  /* Allocate array of qvertex to sort the vertices. */
  errno=0;
  qvertices=malloc(ngood*sizeof(*qvertices));
  if(!qvertices)
    error(EXIT_FAILURE, errno, "%s: failed to allocate %zu "
	  " bytes for 'qvertices'.", __func__,
	  ngood*sizeof(*qvertices));

  /* Loop over the list and extract information for sorting. */
  i=0;
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
      qvertices[i].index=sid;
      qvertices[i].distance=( ( c1[sid]-ref_c1)*(c1[sid]-ref_c1)
			      +(c2[sid]-ref_c2)*(c2[sid]-ref_c2) );

      /* Increment i. */
      ++i;
    }

  /* Clean the temporary space for the list. */
  gal_list_sizet_free(good_stars);

  /* Sort on the basis of distance from the brightest. */
  qsort(qvertices, ngood, sizeof(struct quad_vertex), sort_by_distance);

  /* Make the temporary quad and assign relative positions
     of sorted stars to it. For temporaay C and D, we use the stars
     at 0.5 and 0.75 quantiles respectively. */
  switch(ngood)
    {
    case 4:
      good_vertices[0] = qvertices[0];
      good_vertices[1] = qvertices[1];
      good_vertices[2] = qvertices[2];
      good_vertices[3] = qvertices[3];
      break;
    case 5:
      good_vertices[0] = qvertices[0];
      good_vertices[1] = qvertices[2];
      good_vertices[2] = qvertices[3];
      good_vertices[3] = qvertices[4];
      break;
    default:
      good_vertices[0] = qvertices[0];
      good_vertices[1] = qvertices[ngood/2];
      good_vertices[2] = qvertices[ngood*3/4];
      good_vertices[3] = qvertices[ngood-1];
      break;
    }

  /* Clean up and return a success value (1). */
  free(qvertices);
  return 1;
}





/* A quad has been found over the query catalog, we now want to match
   it with the reference quads and keep the indexs. */
void
match_quad_to_ref(struct params *p, struct matched_points *matched,
		  size_t *abcd, double *geohash, uint16_t rel_b)
{
  size_t qindex;
  double least_dist;

  /* For easy reading. */
  double *x=p->x->array;
  double *y=p->y->array;
  double *ra=p->ra->array;
  double *cx=p->cx->array;
  double *cy=p->cy->array;
  double *dx=p->dx->array;
  double *dy=p->dy->array;
  double *dec=p->dec->array;
  uint32_t *a_ind=p->a_ind->array;
  uint32_t *b_ind=p->b_ind->array;
  uint32_t *c_ind=p->c_ind->array;
  uint32_t *d_ind=p->d_ind->array;
  uint16_t *rel_brightness=p->rel_brightness->array;

  /* Find the matching quad from the reference catalog. */
  qindex=gal_kdtree_nearest_neighbour(p->cx, p->left, p->kdtree_root,
				       geohash, &least_dist);

  /* For a check: */
  printf("refquad: %zu (dist: %g)\n", qindex, least_dist);
  printf("%-10s%-30s%s\n", "Check", "Reference", "Query");
  printf("------------------------------------------------\n");
  printf("%-10s%-30u%u\n", "rel_b", rel_brightness[qindex], rel_b);
  printf("%-10s%-10.3f%-20.3f%-10.3f%.3f\n", "Cx,Cy", cx[qindex], cy[qindex],
	 geohash[0], geohash[1]);
  printf("%-10s%-10.3f%-20.3f%-10.3f%.3f\n", "Dx,Dy", dx[qindex], dy[qindex],
	 geohash[1], geohash[2]);
  printf("\n");
  printf("%-10s%-10g%-20g%-10g%g\n", "A", ra[a_ind[qindex]], dec[a_ind[qindex]],
	 x[abcd[0]], y[abcd[0]]);
  printf("%-10s%-10g%-20g%-10g%g\n", "B", ra[b_ind[qindex]], dec[b_ind[qindex]],
	 x[abcd[1]], y[abcd[1]]);
  printf("%-10s%-10g%-20g%-10g%g\n", "C", ra[c_ind[qindex]], dec[c_ind[qindex]],
	 x[abcd[2]], y[abcd[2]]);
  printf("%-10s%-10g%-20g%-10g%g\n", "D", ra[d_ind[qindex]], dec[d_ind[qindex]],
	 x[abcd[3]], y[abcd[3]]);


  printf("\n...%s...\n", __func__);
  exit(0);
}








/* Make the quads from the top 5 brightest stars in each box in the
   grid. This is the worker function used by `gal_thread_spin_off` to
   be used on multiple threads to make the quads independently. */
static void *
make_quads_worker(void *in_prm)
{
  /* Low-level definitions to be done first. */
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct params *p=(struct params *)tprm->params;

  /* Subsequent definitions. */
  uint16_t rel_b;
  size_t abcd[4];
  size_t i, qindex;
  double geohash[4];
  char regname[100];
  FILE *regfile=NULL;
  int make_ds9_reg=1, quadfound;
  struct quad_vertex good_vertices[4]={0};

  /* For easy reading. */
  double *c1=p->c1->array, *c2=p->c2->array;
  double *cx=p->cx->array, *cy=p->cy->array;
  double *dx=p->dx->array, *dy=p->dy->array;
  uint16_t *rel_brightness=p->rel_brightness->array;
  uint32_t *a_ind=p->a_ind->array, *b_ind=p->b_ind->array;
  uint32_t *c_ind=p->c_ind->array, *d_ind=p->d_ind->array;
  struct matched_points *matched= p->matched ? p->matched[tprm->id] : NULL;

  /* Visualize the quads as a ds9 region file on an image with an
     existing WCS. If you want the check, simply set 'makereg' to a
     value of 1. Note that this is thread-safe: each thread will make
     its own region file, you can then load them together into DS9. */
  if(make_ds9_reg)
    regfile=match_quad_as_ds9_reg_init(regname, p->ds9regprefix,
				       tprm->id);

  /* Go over all the quads that were assigned to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* Extract this quad's index in the final table. */
      qindex=tprm->indexs[i];

      /**************FOR TESTS*******************/
      if(p->c1==p->ra) qindex=1;
      /******************************************/

      /* Based on this bright star, build a quad. If there aren't
	 enough vertices to build a quad, just ignore this entry. */
      quadfound=quad_from_point(p, qindex, good_vertices);
      if(quadfound==0)
	{
	  if(p->c1==p->ra)	/* In the reference catalog, we need to */
	    {			/* set the elements to blank. */
	      a_ind[qindex]=GAL_BLANK_UINT32;
	      rel_brightness[qindex]=GAL_BLANK_UINT16;
	      cx[qindex]=cy[qindex]=dx[qindex]=dy[qindex]=NAN;
	      b_ind[qindex]=c_ind[qindex]=d_ind[qindex]=GAL_BLANK_UINT32;
	    }
	  continue;
	}

      /* If necessary (regfile isn't NULL), add a line for this
	 polygon into the ds9 region file. */
      if(regfile)
	match_quad_as_ds9_reg(regfile, c1, c2, qindex, good_vertices);

      /* Identify which vertices are A, B, C and D (based on special
	 geometric definitions described above) and calculate the
	 hashes.  */
      rel_b=hash_build_write(p, qindex, good_vertices, abcd, geohash);

      /* For a check:
      printf("\n======\nQuad hashes:\n");
      printf("cx,cy,dx,dy: %g, %g, %g, %g\n", geohash[0], geohash[1],
             geohash[2], geohash[3]);
      printf("a_ind,b_ind,c_ind,d_ind: %u, %u, %u, %u\n", abcd[0], abcd[1],
	     abcd[2], abcd[3]);
      printf("rel_brightness: %u\n", rel_b);
      exit(0);
      */

      /* If we are building reference catalog quads, then write the
	 values into the respective column. */
      if(p->c1==p->ra)
	{
	  cx[qindex]=geohash[0];    cy[qindex]=geohash[1];
	  dx[qindex]=geohash[2];    dy[qindex]=geohash[3];
	  a_ind[qindex]=abcd[0];    b_ind[qindex]=abcd[1];
	  c_ind[qindex]=abcd[2];    d_ind[qindex]=abcd[3];
	  rel_brightness[qindex]=rel_b;
	}
      else
	match_quad_to_ref(p, matched, abcd, geohash, rel_b);
    }

  /* Close the region file (if it was created). */
  if(regfile)
    if(fclose(regfile)==EOF)
      error(EXIT_FAILURE, errno, "%s", regname);

  /* Wait for all the other threads to finish, then return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





static void
read_ref_and_query(struct params *p, int r1q0, char *filename)
{
  gal_list_str_t *cols=NULL;
  gal_data_t *data, *c1, *c2, *mag;
  char *c1name  = r1q0 ? "ra"  : "X";
  char *c2name  = r1q0 ? "dec" : "Y";
  char *magname = r1q0 ? "phot_g_mean_mag" : "MAGNITUDE";

  /* Choose columns to read. */
  gal_list_str_add(&cols, c1name, 0);
  gal_list_str_add(&cols, c2name, 0);
  gal_list_str_add(&cols, magname, 0);
  gal_list_str_reverse(&cols);

  /* Read the columns and make sure that only three are found. */
  data=gal_table_read(filename, "1", NULL, cols,
		      GAL_TABLE_SEARCH_NAME, 0, -1, 0, NULL);
  if(gal_list_data_number(data) != 3)
    error(EXIT_FAILURE, 0, "%s: %s should atleast have 3 columns "
          "(x, y, magnitude) but %zu columns matched these names",
          __func__, filename, gal_list_data_number(data));

  /* Seperate columns. */
  c1=data;
  c2=c1->next;
  mag=c2->next;
  c1  = gal_data_copy_to_new_type_free(c1,  GAL_TYPE_FLOAT64);
  c2  = gal_data_copy_to_new_type_free(c2,  GAL_TYPE_FLOAT64);
  mag = gal_data_copy_to_new_type_free(mag, GAL_TYPE_FLOAT32);

  /* Write them into the parameters structure. */
  if(r1q0) { p->ra=c1; p->dec=c2; p->r_mag=mag; }
  else     { p->x=c1;  p->y=c2;   p->q_mag=mag; }
}





static void
prepare_reference(struct params *p, char *reference_name, size_t num_quads)
{
  int quitemmap=1;
  size_t minmapsize=-1;

  /* Read the three necessary references columns. */
  read_ref_and_query(p, 1, reference_name);

  /* For the core inputs, set the pointers to the reference arrays. */
  p->c1  = p->ra;
  p->c2  = p->dec;
  p->mag = p->r_mag;

  /* Allocate all the necessary columns. */
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

  p->rel_brightness=gal_data_alloc(NULL, GAL_TYPE_UINT16, 1, &num_quads,
				   NULL, 0, minmapsize, quitemmap,
				   "rel-brightness", "none",
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
highlevel_reference(char *reference_name, char *kdtree_name,
		    char *ds9regprefix)
{
  size_t kdtree_root;
  float max_mag_diff=2;
  size_t num_threads=1;
  size_t x_numbin=5, y_numbin=5, num_in_gpixel=5;

  /* Internal. */
  struct params p={0};
  struct grid in_grid={0};
  size_t num_quads=x_numbin*y_numbin*num_in_gpixel;

  /* Read the input, allocate output. */
  p.max_mag_diff=max_mag_diff;
  p.ds9regprefix=ds9regprefix;
  prepare_reference(&p, reference_name, num_quads);

  /* Build a box-grid, and set 'max_star_dis_in_quad'. */
  grid_make(p.ra, p.dec, x_numbin, y_numbin, &in_grid);
  p.max_star_dis_in_quad=max(in_grid.step_size[0],
                             in_grid.step_size[1]);

  /* Find the top 5 star ids in each grid element. */
  p.brightest_star_id=find_brightest_stars(p.ra, p.dec, p.r_mag,
					   num_quads, &in_grid,
                                           num_in_gpixel);

  /* Spin-off the threads to calculate quad hashes. */
  gal_threads_spin_off(make_quads_worker, &p, num_quads, num_threads);

  /* Construct a tree and fix the column pointers. Note that the
     kd-tree is ignorant to our higher-level columns. We only want to
     build the tree with the quad geometric hashes (Cx, Cy, Dx, Dy),
     so we'll set the 'dy->next' to NULL before calling the kdtree
     function. Afterwards, we'll set the 'next' pointers again to
     write them all into one final table. */
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
  gal_data_free(p.ra);    p.ra=NULL;
  gal_data_free(p.dec);   p.dec=NULL;
  gal_data_free(p.r_mag); p.r_mag=NULL;
  gal_list_data_free(p.cx);
  return kdtree_root;
}




















/***********************************************************/
/********      Match X,Y quads with ref quads       ********/
/***********************************************************/
static void
match_prepare(struct params *p, char *reference_name, char *kdtree_name,
	      char *query_name, size_t num_threads)
{
  gal_data_t *ref_quad_kdtree;

  /* Read the reference and query tables. */
  read_ref_and_query(p, 0, query_name);
  read_ref_and_query(p, 1, reference_name);

  /* Set the pointers to the query arrays. */
  p->c1  = p->x;
  p->c2  = p->y;
  p->mag = p->q_mag;

  /* Read the input kd-tree. */
  ref_quad_kdtree=gal_table_read(kdtree_name, "1", NULL,
                                 NULL, 0, 0, -1, 0, NULL);

  /* Make sure there are 11 columns. */
  if(gal_list_data_number(ref_quad_kdtree) != 11)
    error(EXIT_FAILURE, 0, "%s: %s should be 11 columns but it is %zu columns",
          __func__, kdtree_name, gal_list_data_number(ref_quad_kdtree));

  /* Cx */
  p->cx=ref_quad_kdtree;
  if(p->cx->type != GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "%s: %s 1st column should be float64 but it is %s",
          __func__, kdtree_name, gal_type_name(p->cx->type, 1));

  /* Cy */
  p->cy=p->cx->next;
  if(p->cy->type != GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "%s: %s 2nd column should be float64 but it is %s",
          __func__, kdtree_name, gal_type_name(p->cy->type, 1));

  /* Dx */
  p->dx=p->cy->next;
  if(p->dx->type != GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "%s: %s 3rd column should be float64 but it is %s",
          __func__, kdtree_name, gal_type_name(p->dx->type, 1));

  /* Dy */
  p->dy=p->dx->next;
  if(p->dy->type != GAL_TYPE_FLOAT64)
    error(EXIT_FAILURE, 0, "%s: %s 4th column should be float64 but it is %s",
          __func__, kdtree_name, gal_type_name(p->dy->type, 1));

  /* Relative brightness. */
  p->rel_brightness=p->dy->next;
  if(p->rel_brightness->type != GAL_TYPE_UINT16)
    error(EXIT_FAILURE, 0, "%s: %s 5th column should be uint16 but it is %s",
          __func__, kdtree_name, gal_type_name(p->rel_brightness->type, 1));

  /* Index of quad's A vertice */
  p->a_ind=p->rel_brightness->next;
  if(p->a_ind->type != GAL_TYPE_UINT32)
    error(EXIT_FAILURE, 0, "%s: %s 6th column should be uint32 but it is %s",
          __func__, kdtree_name, gal_type_name(p->a_ind->type, 1));

  /* Index of quad's B vertice */
  p->b_ind=p->a_ind->next;
  if(p->b_ind->type != GAL_TYPE_UINT32)
    error(EXIT_FAILURE, 0, "%s: %s 7th column should be uint32 but it is %s",
          __func__, kdtree_name, gal_type_name(p->b_ind->type, 1));

  /* Index of quad's C vertice */
  p->c_ind=p->b_ind->next;
  if(p->c_ind->type != GAL_TYPE_UINT32)
    error(EXIT_FAILURE, 0, "%s: %s 8th column should be uint32 but it is %s",
          __func__, kdtree_name, gal_type_name(p->c_ind->type, 1));

  /* Index of quad's D vertice */
  p->d_ind=p->c_ind->next;
  if(p->d_ind->type != GAL_TYPE_UINT32)
    error(EXIT_FAILURE, 0, "%s: %s 9th column should be uint32 but it is %s",
          __func__, kdtree_name, gal_type_name(p->d_ind->type, 1));

  /* kd-tree left subtrees. */
  p->left=p->d_ind->next;
  if(p->left->type != GAL_TYPE_UINT32)
    error(EXIT_FAILURE, 0, "%s: %s 10th column should be uint32 but it is %s",
          __func__, kdtree_name, gal_type_name(p->left->type, 1));

  /* kd-tree right subtrees. */
  p->right=p->left->next;
  if(p->right->type != GAL_TYPE_UINT32)
    error(EXIT_FAILURE, 0, "%s: %s 11th column should be uint32 but it is %s",
          __func__, kdtree_name, gal_type_name(p->right->type, 1));

  /* We don't need all of the columns to remain as lists any more so
     remove the extra 'next' pointers (that can cause problems when
     feeding to kd-tree functions for example). */
  p->dy->next=NULL;
  p->a_ind->next=p->b_ind->next=p->c_ind->next=p->d_ind->next=NULL;


  /* Allocate array keeping a separate list for each thread. */
  errno=0;
  p->matched = malloc(num_threads * sizeof *(p->matched) );
  if(p->matched==NULL)
    error(EXIT_FAILURE, errno, "%s: couldn't allocate %zu bytes for "
	  "'p->matched'", __func__, num_threads * sizeof *(p->matched));
}





static void
highlevel_query(char *reference_name, char *kdtree_name, char *query_name,
		char *output_name, size_t kdtree_root)
{
  struct params p={0};

  size_t num_threads=1;
  float max_mag_diff=2;
  struct grid in_grid={0};
  size_t x_numbin=5, y_numbin=5, num_in_gpixel=10;
  size_t num_quads=x_numbin*y_numbin*num_in_gpixel;

  /* Read and do sanity checks. */
  match_prepare(&p, reference_name, kdtree_name, query_name, num_threads);

  /* make a box-grid */
  grid_make(p.x, p.y, x_numbin, y_numbin, &in_grid);

  /* Find the top 5 star ids using structure. */
  p.brightest_star_id=find_brightest_stars(p.x, p.y, p.q_mag,
                                           num_quads, &in_grid,
                                           num_in_gpixel);

  /* Add last configuration and spin-off the threads. */
  p.kdtree_root=kdtree_root;
  p.max_mag_diff=max_mag_diff;
  p.max_star_dis_in_quad=max(in_grid.step_size[0],
                             in_grid.step_size[1]);
  gal_threads_spin_off(make_quads_worker, &p, num_quads, num_threads);
}




















/***********************************************************/
/********               User interface              ********/
/***********************************************************/
int
main()
{
  size_t kdtree_root;

  /* Reference catalog names. */
  char ds9regprefix[]="./build/quads";
  char *kdtree_name="./build/kdtree.fits";
  char *reference_name="./input/reference.fits";

  /* Query catalog name. */
  char *query_name="./input/query.fits";
  char *output_name = "./build/matched-out.fits";

  /* Quad-hashes creation */
  kdtree_root=highlevel_reference(reference_name, kdtree_name, ds9regprefix);
  printf("kdtree_root: %zu\n", kdtree_root);

  /* Find quads on the query image and match them. */
  highlevel_query(reference_name, kdtree_name, query_name,
		  output_name, kdtree_root);

  /* Finish the program successfully. */
  return EXIT_SUCCESS;
}
