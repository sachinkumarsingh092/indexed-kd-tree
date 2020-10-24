#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <error.h>
#include <stdbool.h>

#include <gnuastro/qsort.h>
#include <gnuastro/table.h>
#include <gnuastro/kdtree.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
#include <gnuastro/polygon.h>
#include <gnuastro/statistics.h>

#include "match.h"








/***********************************************************/
/********          Basic macros/structures          ********/
/***********************************************************/

/* Internally used macro to help in the processing */
# ifndef M_PI
# define M_PI 3.14159265358979323846
# endif

#define MATCH_FLT_ERROR 1e-6

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





/* Structure to keep parameters for multithreaded worker function. */
struct params
{
  /* Input parameters. */
  float max_mag_diff;
  char *ds9regprefix;

  /* Inputs datasets. */
  gal_data_t *x, *y;
  gal_data_t *ra, *dec;
  gal_data_t *r_mag, *q_mag;

  /* These will either point to the respective query or reference
     catalog depending on the caller. They are defined for avoiding
     complexity (they are set once by the top function, and the proper
     arrays will be used in all steps afterwards). */
  gal_data_t *mag;
  gal_data_t *c1, *c2;

  /* Internal info for building quads. */
  size_t *brightest_star_id;
  double max_star_dis_in_quad;

  /* Reference table quad information. */
  size_t kdtree_root;
  gal_data_t *left, *right;
  gal_data_t *bm, *cm, *dm;
  gal_data_t *cx, *cy, *dx, *dy;
  gal_data_t *a_ind, *b_ind, *c_ind, *d_ind;

  /* Query quad information. */
  gal_data_t *ref_qry;
  gal_data_t *a_ind_qry, *b_ind_qry, *c_ind_qry, *d_ind_qry;

  /* Same as above, but for keeping info generated in each thread. */
  gal_data_t **ref_qry_th;
  gal_data_t **a_ind_th, **b_ind_th;
  gal_data_t **c_ind_th, **d_ind_th;
  gal_data_t **bm_th, **cm_th, **dm_th;
  gal_data_t **cx_th, **cy_th, **dx_th, **dy_th;
};





struct params_on_thread
{
  /* Basic settings. */
  struct params *p;
  size_t thread_index;

  /* Pointers to thread-specific structures. */
  gal_data_t **th_ref_qry_d;
  gal_data_t **th_bm_d, **th_cm_d, **th_dm_d;
  gal_data_t **th_cx_d, **th_cy_d, **th_dx_d, **th_dy_d;
  gal_data_t **th_a_ind_d, **th_b_ind_d, **th_c_ind_d, **th_d_ind_d;

  /* Pointers to the arrays of the quads under study in the thread at
     every moment. */
  uint64_t *th_ref_qry;
  uint32_t *th_a_ind, *th_b_ind, *th_c_ind, *th_d_ind;
  double *th_bm, *th_cm, *th_dm, *th_cx, *th_cy, *th_dx, *th_dy;

  /* For the matching of query quads with reference quads. */
  size_t *refquads;
  gal_list_sizet_t **th_refquads_d;
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
match_quad_as_ds9_reg_init(char *regname, char *regprefix, size_t id,
			   int r1q0m2)
{
  FILE *regfile;

  /* Set the filename (with thread number). */
  sprintf(regname, "%s-%zu-%s.reg", regprefix, id,
	  r1q0m2?(r1q0m2==1?"ref":"match"):"query");

  /* Open the file. */
  errno=0;
  regfile = fopen(regname, "w");
  if(regfile==NULL) error(EXIT_FAILURE, errno, "%s", regname);

  /* Write basic information. */
  fprintf(regfile, "# Region file format: DS9 version 4.1\n");
  fprintf(regfile, "%s\n", r1q0m2?"fk5":"image");

  /* Return the file pointer. */
  return regfile;
}





/* Add the given polygon into the given file. */
static void
match_quad_as_ds9_reg(FILE *regfile, size_t refind, size_t quadind,
		      size_t *vertices, double *c1, double *c2, int r1q0)
{
  size_t ordinds[8];
  double polygon[8]={ c1[vertices[0]], c2[vertices[0]],
                      c1[vertices[1]], c2[vertices[1]],
                      c1[vertices[2]], c2[vertices[2]],
                      c1[vertices[3]], c2[vertices[3]] };

  /* Sort the vertices, so the polygon lines don't collide. */
  gal_polygon_vertices_sort(polygon, 4, ordinds);

  /* Print the ds9 definition. */
  fprintf(regfile, "polygon(%g,%g,%g,%g,%g,%g,%g,%g) # text={%s-%zu-%zu}\n",
	  polygon[ordinds[0]*2],  polygon[ordinds[0]*2+1],
	  polygon[ordinds[1]*2],  polygon[ordinds[1]*2+1],
	  polygon[ordinds[2]*2],  polygon[ordinds[2]*2+1],
	  polygon[ordinds[3]*2],  polygon[ordinds[3]*2+1],
	  r1q0?"r":"q", refind, quadind);
}





static void
match_quads_as_ds9_reg_matched(size_t refqid, double *ra, double *dec,
			       uint32_t *rfa, uint32_t *rfb, uint32_t *rfc,
			       uint32_t *rfd, size_t qryqid, double *x,
			       double *y, uint32_t *qa, uint32_t *qb,
			       uint32_t *qc, uint32_t *qd,
			       double least_dist, char *ds9regprefix)
{
  FILE *regfile;
  size_t ordinds[8];
  static int counter=0;
  char regname[1000], qname[1000];
  double ref[8]={ ra[rfa[refqid]], dec[rfa[refqid]],
                  ra[rfb[refqid]], dec[rfb[refqid]],
                  ra[rfc[refqid]], dec[rfc[refqid]],
                  ra[rfd[refqid]], dec[rfd[refqid]] };
  double query[8]={ x[qa[qryqid]], y[qa[qryqid]],
                    x[qb[qryqid]], y[qb[qryqid]],
                    x[qc[qryqid]], y[qc[qryqid]],
                    x[qd[qryqid]], y[qd[qryqid]] };

  /* Initialize the file and put the 'fk5' marker (for the reference
     catalog) and draw the reference quad. */
  regfile=match_quad_as_ds9_reg_init(regname, ds9regprefix, counter++, 2);
  gal_polygon_vertices_sort(ref, 4, ordinds);
  fprintf(regfile, "polygon(%g,%g,%g,%g,%g,%g,%g,%g) # text={%zu-%s}\n",
	  ref[ordinds[0]*2], ref[ordinds[0]*2+1],
	  ref[ordinds[1]*2], ref[ordinds[1]*2+1],
	  ref[ordinds[2]*2], ref[ordinds[2]*2+1],
	  ref[ordinds[3]*2], ref[ordinds[3]*2+1], refqid, "r");

  /* Set the query name (including the least distance). */
  sprintf(qname, "query-%g", least_dist);

  /* Add the query quad. */
  fprintf(regfile, "image\n");
  gal_polygon_vertices_sort(query, 4, ordinds);
  fprintf(regfile, "polygon(%g,%g,%g,%g,%g,%g,%g,%g) # text={%s}\n",
	  query[ordinds[0]*2], query[ordinds[0]*2+1],
	  query[ordinds[1]*2], query[ordinds[1]*2+1],
	  query[ordinds[2]*2], query[ordinds[2]*2+1],
	  query[ordinds[3]*2], query[ordinds[3]*2+1], qname);

  /* Close the file. */
  if(fclose(regfile)==EOF)
    error(EXIT_FAILURE, errno, "%s", regname);
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
                     gal_data_t *mag_data, size_t *num_quads,
                     struct grid *in_grid, size_t num_in_gpixel)
{
  size_t i, g;
  float *mag=mag_data->array;
  double *x=x_data->array, *y=y_data->array;
  size_t *brightest_star_id, *bsi_counter=NULL;
  size_t *sorted_id=NULL, npoints=mag_data->size;

  /************************************
   To temporarily disable sorting and using all stars. */
  *num_quads=npoints;
  /************************************/

  /* If the number of points is less than the (theoretical) number of
     quads, then just use all the stars (so we don't need to sort). */
  if(npoints > *num_quads)
    {
      sorted_id=gal_pointer_allocate (GAL_TYPE_SIZE_T, npoints, 0,
				      __func__, "sorted_id");
      bsi_counter=gal_pointer_allocate(GAL_TYPE_SIZE_T, in_grid->size,
				       1, __func__, "bsi_counter");
    }
  else *num_quads=npoints;

  /* Allocate the final output. */
  brightest_star_id=gal_pointer_allocate(GAL_TYPE_SIZE_T, *num_quads, 0,
                                         __func__, "brightest_star_id");

  /* If it was necessary to sort the points, do it. */
  if( sorted_id )
    {
      /* Initialise output and sorted_id arrays. */
      for(i=0;i<npoints; ++i) sorted_id[i]=i;
      for(i=0;i<*num_quads;++i) brightest_star_id[i]=GAL_BLANK_SIZE_T;

      /* Set magnitude column as a reference to sort stars ID array in
	 incresing order (note that the "magnitude" in astronomy is
	 the inverse of brightness: as magnitude increases, the object
	 becomes fainter).  */
      gal_qsort_index_single=mag;
      qsort(sorted_id, npoints, sizeof(size_t),
	    gal_qsort_index_single_float32_i);

      /* Parse through the 'sorted_id' array (where the brightest star
	 is the first), and put each star in the respective grid
	 element's set of stars, until the capacity is reached. */
      for (i=0;i<npoints;++i)
	{
	  /* Index of the grid box corresponding to particular RA and
             DEC values in the sorted order. */
	  g = grid_findindex(x[sorted_id[i]], y[sorted_id[i]], in_grid);

	  /* If there are less number of star ids for the box than
	     required, find more stars in that particular box. */
	  if (bsi_counter[g] < num_in_gpixel)
	    {
	      brightest_star_id[g*num_in_gpixel+bsi_counter[g]]=
		sorted_id[i];
	      bsi_counter[g]++;
	    }
	}

      /* Clean and return. */
      free(sorted_id);
      free(bsi_counter);
    }

  /* No sorting was necessary, just use all the stars. */
  else
    for(i=0;i<npoints; ++i) brightest_star_id[i]=i;

  /* Return the  */
  return brightest_star_id;
}




















/***********************************************************/
/*********              Quad hashes                *********/
/***********************************************************/

/* Allocate all the necessary columns for information on all the quads
   build from this bright star (which is shared in all of
   them). Within each thread, at every moment we will be filling the
   quad information on the quads from one bright star. So to simplify
   things, set the pointer arrays here. */
void
quads_allocate_for_thread(struct params_on_thread *p_th,
			  size_t refindinthread, size_t nquads)
{
  struct params *p=p_th->p;
  int quietmmap=p->ra->quietmmap;
  size_t minmapsize=p->ra->minmapsize;

  /* We need the indexs columns for both the references catalog and
     the query catalog. */
  p_th->th_a_ind_d[refindinthread]=gal_data_alloc(NULL, GAL_TYPE_UINT32,
		1, &nquads, NULL, 0, minmapsize, quietmmap, "A-index",
	        "index", "Index of polygon point A in input.");
  p_th->th_b_ind_d[refindinthread]=gal_data_alloc(NULL, GAL_TYPE_UINT32,
		1, &nquads, NULL, 0, minmapsize, quietmmap, "B-index",
		"index", "Index of polygon point B in input.");
  p_th->th_c_ind_d[refindinthread]=gal_data_alloc(NULL, GAL_TYPE_UINT32,
	        1, &nquads, NULL, 0, minmapsize, quietmmap, "C-index",
		"index", "Index of polygon point C in input.");
  p_th->th_d_ind_d[refindinthread]=gal_data_alloc(NULL, GAL_TYPE_UINT32,
		1, &nquads, NULL, 0, minmapsize, quietmmap, "D-index",
		"index", "Index of polygon point D in input.");
  p_th->th_a_ind = p_th->th_a_ind_d[ refindinthread ]->array;
  p_th->th_b_ind = p_th->th_b_ind_d[ refindinthread ]->array;
  p_th->th_c_ind = p_th->th_c_ind_d[ refindinthread ]->array;
  p_th->th_d_ind = p_th->th_d_ind_d[ refindinthread ]->array;


  /* Allocate the reference/query specific columns. */
  if(p->c1==p->ra)		/* Only for reference. */
    {
      p_th->th_bm_d[refindinthread]=gal_data_alloc(NULL, GAL_TYPE_FLOAT64,
		      1, &nquads, NULL, 0, minmapsize, quietmmap,
		      "B-mag-frac", "frac",
		      "Relative mag of B: (mag(B)-mag(A))/mag(A)");
      p_th->th_cm_d[refindinthread]=gal_data_alloc(NULL, GAL_TYPE_FLOAT64,
		      1, &nquads, NULL, 0, minmapsize, quietmmap,
		      "C-mag-frac", "frac",
		      "Relative mag of C: (mag(C)-mag(A))/mag(A)");
      p_th->th_dm_d[refindinthread]=gal_data_alloc(NULL, GAL_TYPE_FLOAT64,
		      1, &nquads, NULL, 0, minmapsize, quietmmap,
		      "D-mag-frac", "frac",
		      "Relative mag of D: (mag(D)-mag(A))/mag(A)");
      p_th->th_cx_d[refindinthread]=gal_data_alloc(NULL, GAL_TYPE_FLOAT64,
		      1, &nquads, NULL, 0, minmapsize, quietmmap,
		      "Cx", "frac", "C position in scaled dimension 1.");
      p_th->th_cy_d[refindinthread]=gal_data_alloc(NULL, GAL_TYPE_FLOAT64,
		      1, &nquads, NULL, 0, minmapsize, quietmmap,
		      "Cy", "frac", "C position in scaled dimension 2.");
      p_th->th_dx_d[refindinthread]=gal_data_alloc(NULL, GAL_TYPE_FLOAT64,
		      1, &nquads, NULL, 0, minmapsize, quietmmap,
		      "Dx", "frac", "D position in scaled dimension 1.");
      p_th->th_dy_d[refindinthread]=gal_data_alloc(NULL, GAL_TYPE_FLOAT64,
		      1, &nquads, NULL, 0, minmapsize, quietmmap,
		      "Dy", "frac", "D position in scaled dimension 2.");
      p_th->th_bm = p_th->th_bm_d [ refindinthread ]->array;
      p_th->th_cm = p_th->th_cm_d [ refindinthread ]->array;
      p_th->th_dm = p_th->th_dm_d [ refindinthread ]->array;
      p_th->th_cx = p_th->th_cx_d [ refindinthread ]->array;
      p_th->th_cy = p_th->th_cy_d [ refindinthread ]->array;
      p_th->th_dx = p_th->th_dx_d [ refindinthread ]->array;
      p_th->th_dy = p_th->th_dy_d [ refindinthread ]->array;
    }
  else				/* only for query. */
    {
      p_th->th_ref_qry_d[refindinthread]=gal_data_alloc(NULL, GAL_TYPE_UINT64,
		    1, &nquads, NULL, 0, minmapsize, quietmmap,
		    "ref_qry", "index", "Ref. quad matching this query quad");
      p_th->th_ref_qry = p_th->th_ref_qry_d[ refindinthread ]->array;
    }
}





/* Build the quad structure of the desired star. */
static void
quads_vertices_find(struct params_on_thread *p_th, size_t refind,
		    size_t refindinthread)
{
  struct params *p=p_th->p;
  gal_list_sizet_t *candidate_list=NULL;
  size_t i, nq, nquads, npoints, *candidates;

  /* For easy reading. */
  double      *c1 = p->c1->array;
  double      *c2 = p->c2->array;
  float      *mag = p->mag->array;
  double max_dist = p->max_star_dis_in_quad;

  /* Get the properties of the bright star used to construct this
     series of quads. */
  double ref_c1 = c1[  refind ];
  double ref_c2 = c2[  refind ];
  float ref_mag = mag[ refind ];

  if(refind==104) printf("ref: %f, %f\n", ref_c1, ref_c2);

  /* Go over all the stars of the input catalog and keep a list of
     candidates to constuct the quad. The candidates have to be within
     a certain spatial range around the bright star source of the
     quad, but also and have to be fainter, but not too fainter, only
     by 'max_mag_diff' (usually 2 or 3). Note that we don't want the
     other stars to have a magnitude too close to the reference/bright
     star because it may be due to */
  for(i=0;i<p->c1->size;++i)
    if ( c1[i]     <= ref_c1+max_dist
	 && c1[i]  >= ref_c1-max_dist
	 && c2[i]  <= ref_c2+max_dist
	 && c2[i]  >= ref_c2-max_dist
	 && mag[i] >  ref_mag
	 && mag[i] <= ref_mag + p->max_mag_diff )
      gal_list_sizet_add(&candidate_list, i);

  /* Having added all the other stars, add the actual reference star
     that must be present in all the polygons. We are putting it in
     the end, so it is the first popped element and becomes the
     zero-th element of the array later. */
  gal_list_sizet_add(&candidate_list, refind);

  /***********************************************
   **** Useful when debugging a certain quad. ****
   **** Just set it to one grid element, and  ****
   ****      one star per grid element        ****
   **********************************************
  gal_list_sizet_free(candidate_list);
  candidate_list=NULL;
  if(p->c1==p->ra)
    {
      gal_list_sizet_add(&candidate_list, XXX);
      gal_list_sizet_add(&candidate_list, XXX);
      gal_list_sizet_add(&candidate_list, XXX);
      gal_list_sizet_add(&candidate_list, XXX);
    }
  else
    {
      gal_list_sizet_add(&candidate_list, XXX);
      gal_list_sizet_add(&candidate_list, XXX);
      gal_list_sizet_add(&candidate_list, XXX);
      gal_list_sizet_add(&candidate_list, XXX);
    }
    *********************************************/

  /* Count the number of good stars and if they are less than 4, just
     return (this bright star is not useful for quads given the
     parameters). */
  npoints=gal_list_sizet_number(candidate_list);
  if(npoints < 4) return;

  /* Calculate the total number of quads possible from all the good
     points. Here is the logic: All the quads should share the bright
     star as one vertice (which is also the 0-th element in the
     'candidates' array below). Then we start filling the other 3
     vertices with the other indexs in this manner:

         0 1 2 3
         0 X 2 3        (where 3 < X < npoints: npoints-4 quads)
         0 1 X 3        (where 3 < X < npoints: npoints-4 quads)
         0 1 2 X        (where 3 < X < npoints: npoints-4 quads)

     So in total, we have '(npoints-4)*3 + 1' possible quads. */
  nquads = (npoints-4) * 3 + 1;

  /* Allocate the necessary datasets on each thread. Also convert the
     linked list to an array and free the list. Don't worry about the
     pointer to 'npoints', 'gal_list_size_t_to_array' will overwrite it
     with the same value. */
  candidates=gal_list_sizet_to_array(candidate_list, 0, &npoints);
  quads_allocate_for_thread(p_th, refindinthread, nquads);
  gal_list_sizet_free(candidate_list);

  /* Fill the indexs of all quad vertices into the "X_ind' arrays. In
     the end, these arrays should have geometric meaning and we will
     re-arrange the quad vertices based on geometry later. For now we
     just need to keep the indexs. As shown above, the first quad's
     vertices will be the first three in 'candidates' array. Then we
     can set the first vertice of all quads to 'candidate[0]' (the
     reference star) which should be shared by all quads. */
  nq=0;
  p_th->th_a_ind[ nq ]=candidates[ 0 ];
  p_th->th_b_ind[ nq ]=candidates[ 1 ];
  p_th->th_c_ind[ nq ]=candidates[ 2 ];
  p_th->th_d_ind[ nq ]=candidates[ 3 ];
  for(i=0; i<nquads; ++i) p_th->th_a_ind[i]=candidates[0];
  for(i=4; i<npoints; ++i)
    {
      p_th->th_b_ind[ ++nq ] = candidates[ i ];
      p_th->th_c_ind[   nq ] = candidates[ 2 ];
      p_th->th_d_ind[   nq ] = candidates[ 3 ];

      p_th->th_b_ind[ ++nq ] = candidates[ 1 ];
      p_th->th_c_ind[   nq ] = candidates[ i ];
      p_th->th_d_ind[   nq ] = candidates[ 3 ];

      p_th->th_b_ind[ ++nq ] = candidates[ 1 ];
      p_th->th_c_ind[   nq ] = candidates[ 2 ];
      p_th->th_d_ind[   nq ] = candidates[ i ];
    }

  /* We now have all the vertices of all the possible quads, so let's
     clean up and return. */
  free(candidates);
}





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





/* Given a set of vertice indexs and the full coordinate rows,
   re-arrange the vertices so point A is the first, point B is the
   second, point C is the third and point D is the fourth. */
static void
hash_geometric(double *c1_arr, double *c2_arr, size_t *vertices,
	       double *hash)
{
  size_t i, j, tmpind;
  int perm_set[4][4]={0};
  double tmpdbl, angle_a0b, angle_a1b, angle_c0d, angle_c1d;
  double a[2], b[2], c[2], d[2], distance, current_max_dis=0.0f;
  double c1[4]={ c1_arr[ vertices[0] ], c1_arr[ vertices[1] ],
                 c1_arr[ vertices[2] ], c1_arr[ vertices[3] ] };
  double c2[4]={ c2_arr[ vertices[0] ], c2_arr[ vertices[1] ],
                 c2_arr[ vertices[2] ], c2_arr[ vertices[3] ] };
  size_t abcd[4]={ GAL_BLANK_SIZE_T, GAL_BLANK_SIZE_T,
                   GAL_BLANK_SIZE_T, GAL_BLANK_SIZE_T };

  /* Find the stars that are most distant from each other (A & B). */
  for(i=0;i<4;++i)
    for(j=0;j<4;++j)
      {
	/* Remove repeated pairs like {AB, BA} etc. Also remove
	   cases where the same stars are used like AA.*/
	if(perm_set[i][j]==0 && i!=j)
	  {
	    /* This combination was unseen so flag it.*/
	    perm_set[j][i] = perm_set[i][j] = 1;

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
		abcd[0]=vertices[i];
		abcd[1]=vertices[j];
	      }
	  }
      }

  /* We have the indexes of A/B (not yet knowing which is which). For
     now, the remaining vertices are simply initialized as C and D in
     the order that they appear.*/
  for(i=0;i<4;++i)
    if( vertices[i] != abcd[0] && vertices[i] != abcd[1] )
      {
	if(abcd[2]==GAL_BLANK_SIZE_T) abcd[2]=vertices[i];
	else                          abcd[3]=vertices[i];
      }

  /* To help in readability. */
  a[0]=c1_arr[ abcd[0] ];      a[1]=c2_arr[ abcd[0] ];
  b[0]=c1_arr[ abcd[1] ];      b[1]=c2_arr[ abcd[1] ];
  c[0]=c1_arr[ abcd[2] ];      c[1]=c2_arr[ abcd[2] ];
  d[0]=c1_arr[ abcd[3] ];      d[1]=c2_arr[ abcd[3] ];

  /* For a check.
  printf("\n========================================================="
	 "\nInitial (A&B separated from C&D, not yet from each other)"
	 "\n=========================================================\n");
  printf("A: %-8zu -> (%g, %g)\n", abcd[0], a[0], a[1]);
  printf("B: %-8zu -> (%g, %g)\n", abcd[1], b[0], b[1]);
  printf("X: %-8zu -> (%g, %g)\n", abcd[2], c[0], c[1]);
  printf("X: %-8zu -> (%g, %g)\n", abcd[3], d[0], d[1]);
  */

  /* First we need to uniquely identify A from the pair A or B. We'll
     define A to be the one where the angle from C to D through that
     point is less. But if the two angles are almost identical, just
     ignore this quad. */
  angle_c0d=find_angle_aob(c, a, d);
  angle_c1d=find_angle_aob(c, b, d);
  if( fabs(angle_c0d - angle_c1d)<MATCH_FLT_ERROR )
    {
      /* For a check.
      printf("\nThe angle between C0D and C1D are equal!"
	     "\n .... DISCARDING QUAD ...\n");
      */
      hash[0]=hash[1]=hash[2]=hash[3]=NAN;
      vertices[0]=vertices[1]=vertices[2]=vertices[3]=GAL_BLANK_SIZE_T;
      return;
    }
  if(angle_c0d > angle_c1d)	        /* We need to flip A & B. */
    {
      tmpind=abcd[0];    abcd[0]=abcd[1];     abcd[1]=tmpind;
      tmpdbl=a[0];       a[0]=b[0];           b[0]=tmpdbl;
      tmpdbl=a[1];       a[1]=b[1];           b[1]=tmpdbl;
    }

  /* For a check:
  printf("CAD angle: %.3f (deg)\nCBD angle: %.3f (deg)\n",
	 angle_c0d, angle_c1d);
  if(angle_c0d > angle_c1d)
    printf("+++++++ A & B SWAPPED!\n\n");
  printf("\n=========================================================="
	 "\n A is now uniquely identified (based on smaller angle CXD)"
	 "\n==========================================================\n");
  printf("A: %-8zu -> (%g, %g)\n", abcd[0], a[0], a[1]);
  printf("B: %-8zu -> (%g, %g)\n", abcd[1], b[0], b[1]);
  printf("C: %-8zu -> (%g, %g)\n", abcd[2], c[0], c[1]);
  printf("D: %-8zu -> (%g, %g)\n\n", abcd[3], d[0], d[1]);
  */

  /* We will identify C from D using a similar approach as
     indentifying A from B: we will compare the angles ACB and ADB and
     define C to be the one that has the smaller angle. */
  angle_a0b=find_angle_aob(a, c, b);
  angle_a1b=find_angle_aob(a, d, b);
  if( fabs(angle_a0b - angle_a1b)<MATCH_FLT_ERROR )
    {
      /* To help in debugging.
      printf("\nThe angle between A0B and A1B are equal!\n"
	     " .... DISCARDING QUAD ...\n");
      */
      hash[0]=hash[1]=hash[2]=hash[3]=NAN;
      vertices[0]=vertices[1]=vertices[2]=vertices[3]=GAL_BLANK_SIZE_T;
      return;
    }
  if(angle_a0b > angle_a1b)	        /* We need to flip C & D. */
    {
      tmpind=abcd[2];    abcd[2]=abcd[3];     abcd[3]=tmpind;
      tmpdbl=c[0];       c[0]=d[0];           d[0]=tmpdbl;
      tmpdbl=c[1];       c[1]=d[1];           d[1]=tmpdbl;
    }

  /* For a check.
  printf("\n=========================================================="
	 "\n FINAL A, B, C & D, "
	 "\n==========================================================\n");
  printf("A: %-8zu -> (%g, %g)\n", abcd[0], a[0], a[1]);
  printf("B: %-8zu -> (%g, %g)\n", abcd[1], b[0], b[1]);
  printf("C: %-8zu -> (%g, %g)\n", abcd[2], c[0], c[1]);
  printf("D: %-8zu -> (%g, %g)\n", abcd[3], d[0], d[1]);
  */

  /* The four points are now found geometrically, so we can write the
     final proper order inside the vertices. */
  vertices[0]=abcd[0];
  vertices[1]=abcd[1];
  vertices[2]=abcd[2];
  vertices[3]=abcd[3];

  /* Now that we uniquely know which point is A and which point is B,
     we need to scale C and D to be in a coordinate system where A is
     on (0,0) and B is on (1,1). Infact the Cx, Cy, Dx, Dy values are
     the ultimate hashes that we want. */
  hash[0]=( c[0]-a[0] )/( b[0]-a[0] );
  hash[1]=( c[1]-a[1] )/( b[1]-a[1] );
  hash[2]=( d[0]-a[0] )/( b[0]-a[0] );
  hash[3]=( d[1]-a[1] )/( b[1]-a[1] );

  /* For a check.
  printf("Cx/Cy: %g, %g\n", hash[0], hash[1]);
  printf("Dx/Dy: %g, %g\n", hash[2], hash[3]);
  */
}





/* When a quad can't be used for any reason, we need to set the full
   row to blank. */
static void
results_set_row_to_blank(struct params_on_thread *p_th, size_t qind)
{
  /* Both on the query and reference catalog we need to set to the
     vertice indexs to blank. */
  p_th->th_a_ind[qind]=p_th->th_b_ind[qind]=GAL_BLANK_UINT32;
  p_th->th_c_ind[qind]=p_th->th_d_ind[qind]=GAL_BLANK_UINT32;

  /* Set the other values to blank. */
  if(p_th->p->c1 == p_th->p->ra)
    {
      p_th->th_dy[qind]=NAN;
      p_th->th_bm[qind]=p_th->th_cm[qind]=p_th->th_dm[qind]=NAN;
      p_th->th_cx[qind]=p_th->th_cy[qind]=p_th->th_dx[qind]=NAN;
    }
  else
    p_th->th_ref_qry[qind]=GAL_BLANK_UINT64;
}





/* A quad has been found over the query catalog, we now want to match
   it with the reference quads and keep the indexs. */
void
match_query_to_ref(struct params_on_thread *p_th, double *hash,
		   size_t qryqid)
{
  size_t refqid;
  double least_dist;
  struct params *p=p_th->p;

  /* For easy reading. */
  double   *x=p->x->array;
  double   *y=p->y->array;
  double   *ra=p->ra->array;
  double   *dec=p->dec->array;

  uint32_t *rfa=p->a_ind->array,    *qa=p_th->th_a_ind;
  uint32_t *rfb=p->b_ind->array,    *qb=p_th->th_b_ind;
  uint32_t *rfc=p->c_ind->array,    *qc=p_th->th_c_ind;
  uint32_t *rfd=p->d_ind->array,    *qd=p_th->th_d_ind;


  /* Find the matching quad from the reference catalog. */
  refqid=gal_kdtree_nearest_neighbour(p->cx, p->left, p->kdtree_root,
				      hash, &least_dist);

  /* For a check:
  {
    double *cx=p->cx->array;
    double *cy=p->cy->array;
    double *dx=p->dx->array;
    double *dy=p->dy->array;
    double *bm=p->bm->array;
    double *cm=p->cm->array;
    double *dm=p->dm->array;
    printf("\n------------------------------------------------\n");
    printf("ref-q-id: %zu (dist: %g)\n", refqid, least_dist);
    printf("%-10s%-30s%s\n", "Check", "Reference", "Query");
    printf("------------------------------------------------\n");
    printf("%-10s%-10.3f%-20.3f%-10.3f%.3f\n", "Cx,Cy", cx[refqid],
	   cy[refqid], hash[0], hash[1]);
    printf("%-10s%-10.3f%-20.3f%-10.3f%.3f\n", "Dx,Dy", dx[refqid],
	   dy[refqid], hash[2], hash[3]);
    printf("%-10s%-10.3f%-10.3f\n", "Bm", bm[refqid], hash[4]);
    printf("%-10s%-10.3f%-10.3f\n", "Cm", cm[refqid], hash[5]);
    printf("%-10s%-10.3f%-10.3f\n", "Dm", dm[refqid], hash[6]);
    printf("\n");
    printf("%-10s%-10g%-20g%-10g%g\n", "A", ra[rfa[refqid]],
	   dec[rfa[refqid]], x[qa[qryqid]], y[qa[qryqid]]);
    printf("%-10s%-10g%-20g%-10g%g\n", "B", ra[rfb[refqid]],
	   dec[rfb[refqid]], x[qb[qryqid]], y[qb[qryqid]]);
    printf("%-10s%-10g%-20g%-10g%g\n", "C", ra[rfc[refqid]],
	   dec[rfc[refqid]], x[qc[qryqid]], y[qc[qryqid]]);
    printf("%-10s%-10g%-20g%-10g%g\n", "D", ra[rfd[refqid]],
	   dec[rfd[refqid]], x[qd[qryqid]], y[qd[qryqid]]);
  }
  */

  /* If we don't have a match then set all this quad's indexs to NULL. */
  if(least_dist>0.01) results_set_row_to_blank(p_th, qryqid);
  else
    {
      /* For a check.
      printf("Matched points\n---------\n");
      printf("A: %u (query) --> %u (ref)\n", qa[qryqid], rfa[refqid]);
      printf("B: %u (query) --> %u (ref)\n", qb[qryqid], rfb[refqid]);
      printf("C: %u (query) --> %u (ref)\n", qc[qryqid], rfc[refqid]);
      printf("D: %u (query) --> %u (ref)\n", qd[qryqid], rfd[refqid]);
      */

      /* Put the ID of the matching reference quad in the proper place. */
      p_th->th_ref_qry[qryqid]=refqid;

      /* Make a DS9 region file of the matched quads if necessary. */
      if(1)
	match_quads_as_ds9_reg_matched(refqid, ra, dec, rfa, rfb, rfc, rfd,
				       qryqid, x, y, qa, qb, qc, qd,
				       least_dist, p->ds9regprefix);
    }
}





/* Calculate the hashs. */
static void
hashs_calculate(struct params_on_thread *p_th, size_t refind, size_t qind,
		size_t refindinthread, FILE *regfile)
{
  size_t vertices[4];
  double tmpdbl, hash[7];
  struct params *p=p_th->p;
  float *mag=p_th->p->mag->array;

  /* To make things simple, we'll put all the vertice indexs in one
     array. This is just based on the creation ordering. */
  vertices[0]=p_th->th_a_ind[qind];    vertices[1]=p_th->th_b_ind[qind];
  vertices[2]=p_th->th_c_ind[qind];    vertices[3]=p_th->th_d_ind[qind];


  /* Order the vertices based on geometry and calculate the geometric
     hash (first four numbers in 'hash'). If the quad's vertices can't
     be geometrically distinguished, this function will put NaNs in
     the respective hash elements. */
  hash_geometric(p->c1->array, p->c2->array, vertices, hash);
  if( isnan(hash[0]) ) { results_set_row_to_blank(p_th, qind); return; }


  /* Define the X and Y axis. So far, the A, B, C and D points have
     been uniquely defined. However, the X and Y axises haven't been
     uniquely defined yet. So, we'll define the X axis to be the one
     where C has a smaller value, or, where Cx<Cy. Therefore, if the
     current Cx>Cy, we should flip the first two (Cx <--> Cy) and the
     second two (Dx <--> Dy) hash values. But in case Cx is too close
     to Cy this can cause errors, so its safest to just discard the
     quad. */
  if( fabs( hash[0]-hash[1] ) < MATCH_FLT_ERROR )
    { results_set_row_to_blank(p_th, qind); return; }
  if( hash[0] > hash[1] )
    {
      tmpdbl=hash[0];     hash[0]=hash[1];       hash[1]=tmpdbl;
      tmpdbl=hash[2];     hash[2]=hash[3];       hash[3]=tmpdbl;
    }


  /* Make a ds9 region if requested. */
  if(regfile)
    match_quad_as_ds9_reg(regfile, refind, qind, vertices, p->c1->array,
			  p->c2->array, p->c1==p->ra);


  /* Calculate the relative brightness hash: We know that none of the
     magnitude differences are more than 'max_mag_diff', so we'll
     define the hash to be around unity (similar to the Cx and Cys!)
     by dividing the different with the magnitude of star A with
     'max_mag_diff'. Note that we want it to be "around" unity, not
     less than 1 (it is perfectly fine if its larger: A isn't the
     brightest star). */
  hash[4] = (mag[vertices[1]] - mag[vertices[0]]) / p->max_mag_diff;
  hash[5] = (mag[vertices[2]] - mag[vertices[0]]) / p->max_mag_diff;
  hash[6] = (mag[vertices[3]] - mag[vertices[0]]) / p->max_mag_diff;


  /* For a check:
  {
    double *c1=p->c1->array, *c2=p->c2->array;
    printf("\n%zu (%zu of %zu)\n", refind, qind+1,
	   p_th->th_a_ind_d[refindinthread]->size);
    printf(  "-----------------\n");
    printf("mags: %g, %g, %g, %g\n", mag[vertices[0]], mag[vertices[1]],
	   mag[vertices[2]], mag[vertices[3]]);
    printf("hashes: %-8.4g %-8.4g %-8.4g %-8.4g %-8.4g %-8.4g %-8.4g\n",
	   hash[0], hash[1], hash[2], hash[3], hash[4], hash[5],
	   hash[6]);
    printf("coords: (%g, %g), (%g, %g), (%g, %g), (%g, %g)\n",
	   c1[vertices[0]], c2[vertices[0]], c1[vertices[1]], c2[vertices[1]],
	   c1[vertices[2]], c2[vertices[2]], c1[vertices[3]], c2[vertices[3]] );
  }
  */


  /* If we are on the reference catalog, write the hash values
     into the table. */
  if(p->c1==p->ra)
    {
      p_th->th_cx[qind] = hash[0];      p_th->th_bm[qind] = hash[4];
      p_th->th_cy[qind] = hash[1];      p_th->th_cm[qind] = hash[5];
      p_th->th_dx[qind] = hash[2];      p_th->th_dm[qind] = hash[6];
      p_th->th_dy[qind] = hash[3];
    }
  else
    match_query_to_ref(p_th, hash, qind);
}





static gal_data_t *
gal_data_array_ptr_append_noblank_inplace(gal_data_t **dataarr, size_t size)
{
  void *optr;
  int quietmmap;
  gal_data_t *out=NULL;
  char *name, *unit, *comment;
  size_t i, nout=0, minmapsize;
  uint8_t type=GAL_TYPE_INVALID;

  /* Go through the components and use the first one that is not
     NULL. */
  for(i=0;i<size;++i)
    if(dataarr[i])
      {
	type=dataarr[i]->type;
	name=dataarr[i]->name;
	unit=dataarr[i]->unit;
	comment=dataarr[i]->comment;
	quietmmap=dataarr[i]->quietmmap;
	minmapsize=dataarr[i]->minmapsize;
	break;
      }

  /* Get the total number of points in each dataset (if something
     actually exists for it). */
  for(i=0;i<size;++i)
    if( dataarr[i] )
      {
	/* Make sure that they all have the same numerical data type
	   (all are the same type as the first). */
	if( dataarr[i]->type != type )
	  error(EXIT_FAILURE, 0, "%s: all inputs must have the same type, "
		"but the input atleast contains these two different types: "
		"%s and %s", __func__, gal_type_name(type, 1),
		gal_type_name(dataarr[i]->type, 1));

	/* Remove all blanks and add the number of non-blank elements to
	   'out_num'. */
	gal_blank_remove(dataarr[i]);
	nout += dataarr[i]->size;
      }

  /* Only continue if there is actually anything. */
  if(nout)
    {
      /* Allocate the output dataset in the same type. */
      out=gal_data_alloc(NULL, type, 1, &nout, NULL, 0,
			 minmapsize, quietmmap, name, unit, comment);

      /* Copy the array element of each dataset into the output array and
	 return it. */
      nout=0;
      for(i=0;i<size;++i)
	if(dataarr[i])
	  {
	    optr=gal_pointer_increment(out->array, nout, type);
	    memcpy(optr, dataarr[i]->array,
		   dataarr[i]->size*gal_type_sizeof(type));
	    nout+=dataarr[i]->size;
	  }

      /* For a check.
      {
	double *oarr=out->array;
	for(i=0;i<out->size;++i) printf("%zu: %f\n", i, oarr[i]);
	for(i=0;i<size;++i) printf("partial: %zu\n", dataarr[i]->size);
	printf("total: %zu\n", out->size);
      }
      exit(0);
      */
    }

  /* Return the output dataset. */
  return out;
}





/* Merge all the results obtained in this thread (as separate arrays
   of 'gal_data_t' for each star with varying number of good cases)
   into one gal_data_t, then clean up. */
static void
unite_thread_results_clean(struct params_on_thread *p_th, size_t starnum)
{
  struct params *p=p_th->p;
  size_t th_id=p_th->thread_index;

  p->a_ind_th[th_id]=
    gal_data_array_ptr_append_noblank_inplace(p_th->th_a_ind_d, starnum);
  p->b_ind_th[th_id]=
    gal_data_array_ptr_append_noblank_inplace(p_th->th_b_ind_d, starnum);
  p->c_ind_th[th_id]=
    gal_data_array_ptr_append_noblank_inplace(p_th->th_c_ind_d, starnum);
  p->d_ind_th[th_id]=
    gal_data_array_ptr_append_noblank_inplace(p_th->th_d_ind_d, starnum);
  if(p->c1==p->ra)
    {
      p->cx_th[th_id]=
	gal_data_array_ptr_append_noblank_inplace(p_th->th_cx_d, starnum);
      p->cy_th[th_id]=
	gal_data_array_ptr_append_noblank_inplace(p_th->th_cy_d, starnum);
      p->dx_th[th_id]=
	gal_data_array_ptr_append_noblank_inplace(p_th->th_dx_d, starnum);
      p->dy_th[th_id]=
	gal_data_array_ptr_append_noblank_inplace(p_th->th_dy_d, starnum);
      p->bm_th[th_id]=
	gal_data_array_ptr_append_noblank_inplace(p_th->th_bm_d, starnum);
      p->cm_th[th_id]=
	gal_data_array_ptr_append_noblank_inplace(p_th->th_cm_d, starnum);
      p->dm_th[th_id]=
	gal_data_array_ptr_append_noblank_inplace(p_th->th_dm_d, starnum);
    }
  else
    {
      p->ref_qry_th[th_id]=
	gal_data_array_ptr_append_noblank_inplace(p_th->th_ref_qry_d, starnum);
    }

  /* Now that everything is copied outside, clean up all the allocated
     space in this thread. */
  gal_data_array_ptr_free(p_th->th_a_ind_d, starnum, 1);
  gal_data_array_ptr_free(p_th->th_b_ind_d, starnum, 1);
  gal_data_array_ptr_free(p_th->th_c_ind_d, starnum, 1);
  gal_data_array_ptr_free(p_th->th_d_ind_d, starnum, 1);
  if(p->c1==p->ra)
    {
      gal_data_array_ptr_free(p_th->th_cx_d, starnum, 1);
      gal_data_array_ptr_free(p_th->th_cy_d, starnum, 1);
      gal_data_array_ptr_free(p_th->th_dx_d, starnum, 1);
      gal_data_array_ptr_free(p_th->th_dy_d, starnum, 1);
      gal_data_array_ptr_free(p_th->th_bm_d, starnum, 1);
      gal_data_array_ptr_free(p_th->th_cm_d, starnum, 1);
      gal_data_array_ptr_free(p_th->th_dm_d, starnum, 1);
    }
  else
    gal_data_array_ptr_free(p_th->th_ref_qry_d, starnum, 1);

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
  struct params_on_thread p_th={0};

  /* Subsequent definitions. */
  char regname[100];
  FILE *regfile=NULL;
  int make_ds9_reg=1;
  size_t i, j, refind, starnum=0;

  /* Initialize the thread-specific parameters. We need to allocate
     space to keep the information for the quads found in this
     thread. sS first count how many stars have been given to this
     thread, then allocate the necessary 'gal_data_t **' (to keep all
     the quads for each star). Note that for the query catalog, we
     don't need to store the hashes, only the vertice indexs and the
     indexs of the matching quads. */
  p_th.p=p;
  p_th.thread_index=tprm->id;
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i) ++starnum;
  p_th.th_a_ind_d = gal_data_array_ptr_calloc(starnum);
  p_th.th_b_ind_d = gal_data_array_ptr_calloc(starnum);
  p_th.th_c_ind_d = gal_data_array_ptr_calloc(starnum);
  p_th.th_d_ind_d = gal_data_array_ptr_calloc(starnum);
  if(p->c1==p->ra)
    {
      p_th.th_cx_d = gal_data_array_ptr_calloc(starnum);
      p_th.th_cy_d = gal_data_array_ptr_calloc(starnum);
      p_th.th_dx_d = gal_data_array_ptr_calloc(starnum);
      p_th.th_dy_d = gal_data_array_ptr_calloc(starnum);
      p_th.th_bm_d = gal_data_array_ptr_calloc(starnum);
      p_th.th_cm_d = gal_data_array_ptr_calloc(starnum);
      p_th.th_dm_d = gal_data_array_ptr_calloc(starnum);
    }
  else
    p_th.th_ref_qry_d = gal_data_array_ptr_calloc(starnum);

  /* Prepare for visualization of the quads as a ds9 region file on an
     image with an existing WCS. If you want the check, simply set
     'make_ds9_reg' to a value of 1. Note that this is thread-safe:
     each thread will make its own region file, you can then load them
     together into DS9 if you want, or see them separately. */
  if(make_ds9_reg)
    regfile=match_quad_as_ds9_reg_init(regname, p->ds9regprefix,
				       tprm->id, p->c1==p->ra);

  /* Go over all the bright stars that were assigned to this
     thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* Extract the input catalog index of this reference star. */
      refind = p->brightest_star_id[ tprm->indexs[i] ];

      /* Based on this bright star, find all the possible quad
	 vertices and allocate the necessary datasets.  */
      quads_vertices_find(&p_th, refind, i);

      /* If at least one quad could be found, then calculate the
	 hashes and match them if necessary. */
      if(p_th.th_a_ind_d[ i ])
	for(j=0; j<p_th.th_a_ind_d[ i ]->size; ++j)
	  hashs_calculate(&p_th, refind, j, i, regfile);
    }

  /* Unite all the results for each star into one array for this
     thread, then clean up all the no-more-needed spaces. */
  unite_thread_results_clean(&p_th, starnum);

  /* Close the region file for this thread (if it was created). */
  if(regfile)
    if(fclose(regfile)==EOF)
      error(EXIT_FAILURE, errno, "%s", regname);

  /* Wait for all the other threads to finish, then return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





static void
read_ref_and_query(struct params *p, int r1q0, char *filename, char *hdu)
{
  gal_list_str_t *cols=NULL;
  gal_data_t *data, *c1, *c2, *mag;
  char *c1name  = r1q0 ? "1" : "1";
  char *c2name  = r1q0 ? "2" : "2";
  char *magname = r1q0 ? "3" : "3";

  /* Choose columns to read. */
  gal_list_str_add(&cols, c1name, 0);
  gal_list_str_add(&cols, c2name, 0);
  gal_list_str_add(&cols, magname, 0);
  gal_list_str_reverse(&cols);

  /* Read the columns and make sure that only three are found. */
  data=gal_table_read(filename, hdu, NULL, cols,
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




/* Prepare the necessary inputs for each array. */
static void
prepare_reference(struct params *p, char *reference_name, size_t numthreads)
{
  char *refhdu="1";

  /* Read the three necessary references columns. */
  read_ref_and_query(p, 1, reference_name, refhdu);

  /* For the reference, set the three pointers to the reference arrays. */
  p->c1  = p->ra;
  p->c2  = p->dec;
  p->mag = p->r_mag;

  /* Allocate arrays of 'gal_data_t *'s to keep the results of each
     thread. */
  p->cx_th=gal_data_array_ptr_calloc(numthreads);
  p->cy_th=gal_data_array_ptr_calloc(numthreads);
  p->dx_th=gal_data_array_ptr_calloc(numthreads);
  p->dy_th=gal_data_array_ptr_calloc(numthreads);
  p->bm_th=gal_data_array_ptr_calloc(numthreads);
  p->cm_th=gal_data_array_ptr_calloc(numthreads);
  p->dm_th=gal_data_array_ptr_calloc(numthreads);
  p->a_ind_th=gal_data_array_ptr_calloc(numthreads);
  p->b_ind_th=gal_data_array_ptr_calloc(numthreads);
  p->c_ind_th=gal_data_array_ptr_calloc(numthreads);
  p->d_ind_th=gal_data_array_ptr_calloc(numthreads);
}





/* Merge the results on each thread into one. */
static void
unite_threads_all_clean(struct params *p, size_t nt)
{
  /* Unite the results depending on the reference/query phase, and
     clean the arrays that are unique to that phase. */
  if(p->c1==p->ra)
    {
      /* Merge all the arrays into a final array. */
      p->cx    = gal_data_array_ptr_append_noblank_inplace(p->cx_th,    nt);
      p->cy    = gal_data_array_ptr_append_noblank_inplace(p->cy_th,    nt);
      p->dx    = gal_data_array_ptr_append_noblank_inplace(p->dx_th,    nt);
      p->dy    = gal_data_array_ptr_append_noblank_inplace(p->dy_th,    nt);
      p->bm    = gal_data_array_ptr_append_noblank_inplace(p->bm_th,    nt);
      p->cm    = gal_data_array_ptr_append_noblank_inplace(p->cm_th,    nt);
      p->dm    = gal_data_array_ptr_append_noblank_inplace(p->dm_th,    nt);
      p->a_ind = gal_data_array_ptr_append_noblank_inplace(p->a_ind_th, nt);
      p->b_ind = gal_data_array_ptr_append_noblank_inplace(p->b_ind_th, nt);
      p->c_ind = gal_data_array_ptr_append_noblank_inplace(p->c_ind_th, nt);
      p->d_ind = gal_data_array_ptr_append_noblank_inplace(p->d_ind_th, nt);

      /* Clean up the arrays that were only for references.*/
      gal_data_array_ptr_free(p->cx_th, nt, 1);
      gal_data_array_ptr_free(p->cy_th, nt, 1);
      gal_data_array_ptr_free(p->dx_th, nt, 1);
      gal_data_array_ptr_free(p->dy_th, nt, 1);
      gal_data_array_ptr_free(p->bm_th, nt, 1);
      gal_data_array_ptr_free(p->cm_th, nt, 1);
      gal_data_array_ptr_free(p->dm_th, nt, 1);
    }
  else
    {
      /* Merge all the arrays into a final array. */
      p->ref_qry = gal_data_array_ptr_append_noblank_inplace(p->ref_qry_th, nt);
      p->a_ind_qry = gal_data_array_ptr_append_noblank_inplace(p->a_ind_th, nt);
      p->b_ind_qry = gal_data_array_ptr_append_noblank_inplace(p->b_ind_th, nt);
      p->c_ind_qry = gal_data_array_ptr_append_noblank_inplace(p->c_ind_th, nt);
      p->d_ind_qry = gal_data_array_ptr_append_noblank_inplace(p->d_ind_th, nt);

      /* Clean up the arrays that were only for references.*/
      gal_data_array_ptr_free(p->ref_qry_th, nt, 1);
    }

  /* Clean up common arrays. */
  gal_data_array_ptr_free(p->a_ind_th, nt, 1);
  gal_data_array_ptr_free(p->b_ind_th, nt, 1);
  gal_data_array_ptr_free(p->c_ind_th, nt, 1);
  gal_data_array_ptr_free(p->d_ind_th, nt, 1);
}





static void
highlevel_reference_write(struct params *p, char *in_filename,
			   char *out_filename)
{
  gal_data_t *outtable=NULL;
  gal_fits_list_key_t *keylist=NULL;

  /* If any quads could be found. */
  if(p->a_ind)
    {
      /* Define all columns as a list to be printed in order. */
      outtable       = p->a_ind;
      p->a_ind->next = p->b_ind;
      p->b_ind->next = p->c_ind;
      p->c_ind->next = p->d_ind;
      p->d_ind->next = p->cx;    /* Cx already has the hashes as 'next'.*/
      p->dm->next    = p->left;  /* 'right' is already 'next' to 'left'.*/
    }
  /* There weren't any quads: allocate empty datasets to have a table
     with metadata, but without any rows. */
  else
    {
      printf("\n%s: DO THIS PART: THERE AREN'T ANY DATASETS.\n",
	     __func__);
      exit(1);
    }

  /* Add the necessary keywords to write in the output. */
  gal_fits_key_list_title_add_end(&keylist, "Information on table", 0);
  gal_fits_key_list_add_end(&keylist, GAL_TYPE_SIZE_T,
			    GAL_MATCH_KDROOT_KEY_NAME, 0, &p->kdtree_root, 0,
			    "k-d tree root index (counting from 0)", 0,
			    "counter", 0);
  gal_fits_key_write_filename("INPUT", in_filename, &keylist, 0);

  /* Write the final table with all the quad information and the
     respective kd-tree. */
  gal_table_write(outtable, &keylist, NULL, GAL_TABLE_FORMAT_BFITS,
                  out_filename, "quad-kdtree", 0);

  /* To avoid confusion later, we'll remove all the extra 'next'
     pointers that were just added to write into a table. */
  p->dm->next = NULL;
  p->a_ind->next = p->b_ind->next = p->c_ind->next = p->d_ind->next = NULL;
}





static void
highlevel_reference(struct params *p, char *reference_name,
		    char *kdtree_name, size_t numthreads,
		    size_t x_numbin, size_t y_numbin,
		    size_t num_in_gpixel)
{
  /* Internal. */
  struct grid in_grid={0};
  size_t num_quads=x_numbin*y_numbin*num_in_gpixel;

  /* Read the input, allocate output. */
  prepare_reference(p, reference_name, numthreads);

  /* Build a box-grid, and set 'max_star_dis_in_quad'. */
  grid_make(p->ra, p->dec, x_numbin, y_numbin, &in_grid);
  p->max_star_dis_in_quad=max(in_grid.step_size[0],
			      in_grid.step_size[1]);

  /* Find the top 'num_in_gpixel' star ids in each grid element. */
  p->brightest_star_id=find_brightest_stars(p->ra, p->dec, p->r_mag,
					   &num_quads, &in_grid,
                                           num_in_gpixel);

  /* Spin-off the threads to calculate quad hashes. */
  gal_threads_spin_off(make_quads_worker, p, num_quads, numthreads);

  /* All the threads have processed the quads. Now, we can allocate
     the final colums and merge all the different results into one. */
  unite_threads_all_clean(p, numthreads);

  /* Construct a tree and fix the column pointers if there is
     any. Note that the kd-tree is ignorant to our higher-level
     columns. We only want to build the tree with the quad geometric
     hashes (Cx, Cy, Dx, Dy), so we'll set the 'dy->next' to NULL
     before calling the kdtree function.  */
  if(p->cx)
    {
      p->cx->next=p->cy;
      p->cy->next=p->dx;
      p->dx->next=p->dy;
      p->dy->next=p->bm;
      /***** To only have geo-hash *****
      p->dy->next=NULL;
      **********************************/
      p->bm->next=p->cm;
      p->cm->next=p->dm;
      p->dm->next=NULL;
      p->left=gal_kdtree_create(p->cx, &p->kdtree_root);
      /***** To only have geo-hash *****
      p->dy->next=p->bm;
      *********************************/
    }

  /* Write the k-d tree and all quad data into a table, and include
     the kdtree root index and input filename as keyword arguments. */
  if(kdtree_name)
    highlevel_reference_write(p, reference_name, kdtree_name);
}




















/***********************************************************/
/********      Match X,Y quads with ref quads       ********/
/***********************************************************/
static void
prepare_query(struct params *p, char *reference_name, char *kdtree_name,
	      char *query_name, size_t numthreads)
{
  char *refhdu="1", *qhdu="1";
  gal_data_t *ref_quad_kdtree;
  gal_data_t *keysll=gal_data_array_calloc(1);

  /* Read the query table and set the necessary pointers. */
  read_ref_and_query(p, 0, query_name, qhdu);
  p->c1  = p->x;
  p->c2  = p->y;
  p->mag = p->q_mag;

  /* If the reference catalog isn't loaded yet, then read it. */
  if(p->ra==NULL)
    read_ref_and_query(p, 1, reference_name, refhdu);

  /* If the k-d tree and reference quad information aren't ready yet,
     read them. */
  if(p->left==NULL)
    {
      /* Read the input kd-tree. */
      ref_quad_kdtree=gal_table_read(kdtree_name, "1", NULL,
				     NULL, 0, 0, -1, 0, NULL);

      /* Make sure there are 13 columns. */
      if(gal_list_data_number(ref_quad_kdtree) != 13)
	error(EXIT_FAILURE, 0, "%s: %s should be 13 columns but it is %zu "
	      "columns", __func__, kdtree_name,
	      gal_list_data_number(ref_quad_kdtree));

      /* Index of quad's A vertice */
      p->a_ind=ref_quad_kdtree;
      if(p->a_ind->type != GAL_TYPE_UINT32)
	error(EXIT_FAILURE, 0, "%s: %s 1st column should be uint32 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->a_ind->type, 1));

      /* Index of quad's B vertice */
      p->b_ind=p->a_ind->next;
      if(p->b_ind->type != GAL_TYPE_UINT32)
	error(EXIT_FAILURE, 0, "%s: %s 2nd column should be uint32 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->b_ind->type, 1));

      /* Index of quad's C vertice */
      p->c_ind=p->b_ind->next;
      if(p->c_ind->type != GAL_TYPE_UINT32)
	error(EXIT_FAILURE, 0, "%s: %s 3rd column should be uint32 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->c_ind->type, 1));

      /* Index of quad's D vertice */
      p->d_ind=p->c_ind->next;
      if(p->d_ind->type != GAL_TYPE_UINT32)
	error(EXIT_FAILURE, 0, "%s: %s 4th column should be uint32 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->d_ind->type, 1));

      /* Cx */
      p->cx=p->d_ind->next;
      if(p->cx->type != GAL_TYPE_FLOAT64)
	error(EXIT_FAILURE, 0, "%s: %s 5th column should be float64 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->cx->type, 1));

      /* Cy */
      p->cy=p->cx->next;
      if(p->cy->type != GAL_TYPE_FLOAT64)
	error(EXIT_FAILURE, 0, "%s: %s 6th column should be float64 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->cy->type, 1));

      /* Dx */
      p->dx=p->cy->next;
      if(p->dx->type != GAL_TYPE_FLOAT64)
	error(EXIT_FAILURE, 0, "%s: %s 7th column should be float64 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->dx->type, 1));

      /* Dy */
      p->dy=p->dx->next;
      if(p->dy->type != GAL_TYPE_FLOAT64)
	error(EXIT_FAILURE, 0, "%s: %s 8th column should be float64 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->dy->type, 1));

      /* B-mag-frac */
      p->bm=p->dy->next;
      if(p->bm->type != GAL_TYPE_FLOAT64)
	error(EXIT_FAILURE, 0, "%s: %s 9th column should be float64 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->dy->type, 1));

      /* C-mag-frac */
      p->cm=p->bm->next;
      if(p->cm->type != GAL_TYPE_FLOAT64)
	error(EXIT_FAILURE, 0, "%s: %s 10th column should be float64 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->dy->type, 1));

      /* D-mag-frac */
      p->dm=p->cm->next;
      if(p->dm->type != GAL_TYPE_FLOAT64)
	error(EXIT_FAILURE, 0, "%s: %s 11th column should be float64 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->dy->type, 1));

      /* kd-tree left subtrees. */
      p->left=p->dm->next;
      if(p->left->type != GAL_TYPE_UINT32)
	error(EXIT_FAILURE, 0, "%s: %s 12th column should be uint32 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->left->type, 1));

      /* kd-tree right subtrees. */
      p->right=p->left->next;
      if(p->right->type != GAL_TYPE_UINT32)
	error(EXIT_FAILURE, 0, "%s: %s 13th column should be uint32 but it is %s",
	      __func__, kdtree_name, gal_type_name(p->right->type, 1));

      /* We don't need all of the columns to remain as lists any more so
	 remove the extra 'next' pointers (that can cause problems when
	 feeding to kd-tree functions for example). */
      p->dm->next=NULL;
      p->a_ind->next=p->b_ind->next=p->c_ind->next=p->d_ind->next=NULL;

      /***** To only have geo-hash *****
	     p->dy->next=NULL;
      ***********************************/

      /* Read the k-d tree root node. */
      keysll->type=GAL_TYPE_SIZE_T;
      keysll->name=GAL_MATCH_KDROOT_KEY_NAME;
      gal_fits_key_read(kdtree_name, refhdu, keysll, 0, 0);
      if(keysll->status)
	error(EXIT_FAILURE, 0, "%s (hdu: %s): no '%s' keyword found",
	      kdtree_name, refhdu, GAL_MATCH_KDROOT_KEY_NAME);
      p->kdtree_root=((size_t *)(keysll->array))[0];
      keysll->name=NULL;
      gal_data_array_free(keysll, 1, 1);
    }

  /* Allocate space to keep the quad information in the query catalog
     generated on each thread (and later merge into one table). */
  p->a_ind_th=gal_data_array_ptr_calloc(numthreads);
  p->b_ind_th=gal_data_array_ptr_calloc(numthreads);
  p->c_ind_th=gal_data_array_ptr_calloc(numthreads);
  p->d_ind_th=gal_data_array_ptr_calloc(numthreads);
  p->ref_qry_th=gal_data_array_ptr_calloc(numthreads);
}





static void
highlevel_query(struct params *p, char *reference_name, char *kdtree_name,
		char *query_name, char *output_name, size_t numthreads,
		size_t x_numbin, size_t y_numbin, size_t num_in_gpixel)
{
  struct grid in_grid={0};
  size_t num_quads=x_numbin*y_numbin*num_in_gpixel;

  /* Read and do sanity checks. */
  prepare_query(p, reference_name, kdtree_name, query_name, numthreads);

  /* Make a box-grid. */
  grid_make(p->x, p->y, x_numbin, y_numbin, &in_grid);

  /* Find the necessary stars. */
  p->brightest_star_id=find_brightest_stars(p->x, p->y, p->q_mag,
                                           &num_quads, &in_grid,
                                           num_in_gpixel);

  /* Add last configuration and spin-off the threads. */
  p->max_star_dis_in_quad=max(in_grid.step_size[0],
			      in_grid.step_size[1]);
  gal_threads_spin_off(make_quads_worker, p, num_quads, numthreads);

  /* All the threads have processed their quads. Now, we can allocate
     the final colums and merge all the different results into one. */
  unite_threads_all_clean(p, numthreads);
}





static void
highlevel_wcs(struct params *p)
{
  /* Just to get things started: */
  size_t i;
  uint64_t *ref_qry_a=p->ref_qry->array;
  double *x_a=p->x->array, *y_a=p->y->array;
  uint32_t *a_ind_qry_a=p->a_ind_qry->array;
  uint32_t *b_ind_qry_a=p->b_ind_qry->array;
  uint32_t *c_ind_qry_a=p->c_ind_qry->array;
  uint32_t *d_ind_qry_a=p->d_ind_qry->array;

  uint32_t r_qid;
  uint32_t *a_ind_a=p->a_ind->array;
  double *ra_a=p->ra->array, *dec_a=p->dec->array;

  for(i=0;i<p->a_ind_qry->size;++i)
    {
      printf("%u(%g,%g), %u(%g, %g), %u(%g,%g), %u(%g,%g) matches with quad %lu in ref\n",
	     a_ind_qry_a[i], x_a[ a_ind_qry_a[i] ],  y_a[ a_ind_qry_a[i] ],
	     b_ind_qry_a[i], x_a[ b_ind_qry_a[i] ],  y_a[ b_ind_qry_a[i] ],
	     c_ind_qry_a[i], x_a[ c_ind_qry_a[i] ],  y_a[ c_ind_qry_a[i] ],
	     d_ind_qry_a[i], x_a[ d_ind_qry_a[i] ],  y_a[ d_ind_qry_a[i] ],
	     ref_qry_a[i]);

      r_qid=ref_qry_a[i];
      printf("starA: (%g, %g) --> (%g, %g)\n", x_a[ a_ind_qry_a[i] ],  y_a[ a_ind_qry_a[i] ]
                                             , ra_a[a_ind_a[r_qid]], dec_a[a_ind_a[r_qid]]);

      exit(0);
    }
}









/* We'll use homogenous matrix to find 'theta'(t) and 'scale'(s) value
 * with the coordinates given from the query image, the [X] matrix, and
 * the reference matrix with `ra` and `dec` values, matrix [R].
 * We assume counter-clock wise rotation when matching the images.
 *
 * We merge the rotation and scaling warpings by multiplying them.
 *
 *       [X]   =      [ROTATION]       *    [SCALING]   *   [R]
 *      -----    --------------------      -----------    -------
 *      | x |    | cos(t) -sin(t) 0 |      | s  0  0 |    | ra  |
 *      | y |  = | sin(t)  cos(t) 0 |      | 0  s  0 |    | dec |
 *      | 1 |    |   0       0    1 |      | 0  0  1 |    |  1  |
 *
 * Opening up the matrix multiplication, we finally get:
 *
 *	x = s*cos(t)*ra - s*sin(t)*dec				(1)
 *	y = s*sin(t)*ra + s*cos(t)*dec				(2)
 *
 * Multiplying eq. (1) with `ra` and eq. (2) with `dec`, we get:
 *
 * 	ra*x  = s*cos(t)*ra^2  - s*sin(t)*ra*dec			(3)
 * 	dec*y = s*cos(t)*dec^2 + s*sin(t)*ra*dec			(4)
 *
 * Adding (3) and (4), we finally get:
 *
 * 	s*cos(t) = (ra*x + dec*y)/(ra^2 + dec^2)			(5)
 *
 * Taking X = s*cos(t) and Y = s*sin(t), and putting X in (5):
 *
 * 	X = (ra*x + dec*y)/(ra^2 + dec^2)				(6)
 *
 * From eq. (1), we get Y:
 *
 * 	Y = (ra*X - x)/dec						(7)
 *
 * Dividing Y by X, we get theta(t):
 *
 * 	t = atan(Y/X)							(8)
 *
 * Using the value of theta(t), we get scale(s):
 *
 * 	s = X/cos(t)							(9)
 *
 * The sign of theta will reverse for any clockwise rotations.
 */
static void
calculate_rot_scale(struct params *p)
{
  size_t i, j;
  uint64_t *ref_qry_a=p->ref_qry->array;
  double *x_a=p->x->array, *y_a=p->y->array;
  uint32_t *a_ind_qry_a=p->a_ind_qry->array;
  uint32_t *b_ind_qry_a=p->b_ind_qry->array;
  uint32_t *c_ind_qry_a=p->c_ind_qry->array;
  uint32_t *d_ind_qry_a=p->d_ind_qry->array;

  uint32_t r_qid;
  uint32_t *a_ind_a=p->a_ind->array;
  uint32_t *b_ind_a=p->b_ind->array;
  uint32_t *c_ind_a=p->c_ind->array;
  uint32_t *d_ind_a=p->d_ind->array;
  double *ra_a=p->ra->array, *dec_a=p->dec->array;


  for(i=0;i<p->a_ind_qry->size;++i)
    {
      r_qid=ref_qry_a[i];
      double x[4]={ x_a[ a_ind_qry_a[i] ] - x_a[ a_ind_qry_a[i] ],
	            x_a[ b_ind_qry_a[i] ] - x_a[ a_ind_qry_a[i] ],
                    x_a[ c_ind_qry_a[i] ] - x_a[ a_ind_qry_a[i] ],
		    x_a[ d_ind_qry_a[i] ] - x_a[ a_ind_qry_a[i] ] };
      double y[4]={ y_a[ a_ind_qry_a[i] ] - y_a[ a_ind_qry_a[i] ],
	            y_a[ b_ind_qry_a[i] ] - y_a[ a_ind_qry_a[i] ],
                    y_a[ c_ind_qry_a[i] ] - y_a[ a_ind_qry_a[i] ],
		    y_a[ d_ind_qry_a[i] ] - y_a[ a_ind_qry_a[i] ] };

      double r[4]={ ra_a[a_ind_a[r_qid]] - ra_a[a_ind_a[r_qid]],
	            ra_a[b_ind_a[r_qid]] - ra_a[a_ind_a[r_qid]],
                    ra_a[c_ind_a[r_qid]] - ra_a[a_ind_a[r_qid]],
		    ra_a[d_ind_a[r_qid]] - ra_a[a_ind_a[r_qid]] };
      double d[4]={ dec_a[a_ind_a[r_qid]] - dec_a[a_ind_a[r_qid]],
	            dec_a[b_ind_a[r_qid]] - dec_a[a_ind_a[r_qid]],
                    dec_a[c_ind_a[r_qid]] - dec_a[a_ind_a[r_qid]],
		    dec_a[d_ind_a[r_qid]] - dec_a[a_ind_a[r_qid]] };
      for(j=1;j<4;++j)
        {
          double X=((x[j])*r[j]+(y[j])*d[j])/(r[j]*r[j]+d[j]*d[j]);
          double Y=(r[j]*X-x[j])/d[j];
          float theta=atan(Y/X);
          float scale=X/cos(theta);

          printf("\ttheta=%g, scale=%g\n", theta, scale);
        }
      printf("\t===============================\n");
    }
  exit(0);
}












/***********************************************************/
/********               User interface              ********/
/***********************************************************/
int
main()
{
  int buildref=1;
  size_t numthreads=1;
  float max_mag_diff=3;
  size_t x_numbin=1, y_numbin=1, num_in_gpixel=50;

  /* File names. */
  char *query_name="./input2/xy_small.txt";
  char *reference_name="./input2/wcs_small.txt";
  char *ds9regprefix="./build/reg/quads";
  char *kdtree_name="./build/kdtree.fits";
  char *output_name = "./build/matched-out.fits";

  /* Internal parameters. */
  struct params p={0};
  p.max_mag_diff=max_mag_diff;
  p.ds9regprefix=ds9regprefix;

  /* Process reference catalog. When 'kdtree_name==NULL', then the
     reference quad information won't be stored. */
  if(buildref)
    highlevel_reference(&p, reference_name, kdtree_name, numthreads,
			x_numbin, y_numbin, num_in_gpixel);

  /* Find quads on the query image and match them. */
  highlevel_query(&p, reference_name, kdtree_name, query_name,
		  output_name, numthreads, x_numbin, y_numbin,
		  num_in_gpixel);

  /* Using matched quads, find WCS. */
  // highlevel_wcs(&p);
    calculate_rot_scale(&p);


  /* Finish the program successfully. */
  return EXIT_SUCCESS;
}
