#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/qsort.h>
#include <gnuastro/table.h>
#include <gnuastro/threads.h>
#include <gnuastro/statistics.h>



/* The grid structure which has the necessary
   information of the grid on the catalogue. */
struct grid
{
	size_t size;
	size_t dim[2];
	double min[2], step_size[2];
};



/* For internal use only. */
struct brightestStarsId_t
{
  size_t size;
  size_t star_id[5];
};


/* For internal use only. */
struct quad_t
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

	return y*in_grid->dim[0]+x;
}






/****************finding brightest stars****************************/


/* Returns a 2x2 array with first index having the index of the
   box in the grid and the second index having the IDs of the stars. */
void
find_brightestStars(gal_data_t *ra, gal_data_t *dec, gal_data_t *magnitude,
                    struct grid *in_grid, size_t num_brightestStarsId,
                    size_t brightestStarId_arr[][5])
{
  size_t i, j, index;
  size_t *sortedId_arr=NULL;
  double *ra_arr=ra->array;
	double *dec_arr=dec->array;
  float *magnitude_arr=magnitude->array;
  struct brightestStarsId_t *brightestIDs=NULL;


  /* Allocate and initialize the structure. */
  brightestIDs=malloc(magnitude->size*sizeof(*brightestIDs));
  for(i=0;i<magnitude->size;++i)
  {
    for(j=0;j<num_brightestStarsId;j++) brightestIDs[i].star_id[j]=0;

    brightestIDs[i].size=0;
  }


  /* Allocate object id array and initialize. */
  sortedId_arr=malloc (magnitude->size*sizeof (*sortedId_arr));
  for(i=0;i<magnitude->size; ++i)
    sortedId_arr[i]=i;


  /* Pointer to an array magnitude_arr to be used as a reference. */
  gal_qsort_index_single=magnitude_arr;

 /* Use brightness column as a reference to sort stars ID array
      in incresing order(brightness is inverse of number's
      magnitudes).  */
  qsort (sortedId_arr, magnitude->size, sizeof(size_t),
        gal_qsort_index_single_float32_d);


  /* Iterate through all the stars in the grid and save IDs
     of the brightest stars. The array of magnitudes of brightness
     is already sorted so for each box in grid only the brightest
     star indexes will be stored. */
  for (i=0;i<magnitude->size;++i)
    {
      /* Index of the grid box corresponding to particular RA and DEC
         values in the sorted order. */
      index = grid_findindex (ra_arr[sortedId_arr[i]],
                              dec_arr[sortedId_arr[i]], in_grid);

      /* For a check:
      printf("magnitude = %g, magnitude_sorted = %g"
              " index = %zu sorted index = %zu\n",
              magnitude_arr[i], magnitude_arr[sortedId_arr[i]], i, sortedId_arr[i]);
      */

      /* If there are less number of star ids for the box than required,
         find more stars in that particular box. */
      if ( brightestIDs[index].size < num_brightestStarsId )
        {
          /* Store the ID of the brightest stars in the structure. */
          brightestIDs[index].star_id[brightestIDs[index].size]=sortedId_arr[i];

          /* Store the ID of the brightest stars in the 2d array that
             will be returned and finally used. First index having
             the index of the box in the grid and the second index
             having the IDs of the stars.*/
          brightestStarId_arr[index][brightestIDs[index].size]=sortedId_arr[i];

          /* For a check:
            printf("brightest_arr[%zu][%zu] = %zu\n",
                    index, brightestIDs[index].size, sortedId_arr[i]);
          */

          /* Increase the size, which signifies the number of IDs that are
             currently stored for the particular box. */
          brightestIDs[index].size++;
        }
    }

  /* Free the allocations. */
  free(sortedId_arr);
  free(brightestIDs);
}






/***************making quads*************************/


/* Sort by distance as refrerence in ascending order. */
static int
sortByDistance(const void *a, const void *b)
{
  struct quad_t *p1 = (struct quad_t *)a;
  struct quad_t *p2 = (struct quad_t *)b;

  return ( p1->distance==p2->distance
           ? 0
           : (p1->distance<p2->distance ? -1 : 1) );
}



/* Sort by brghtness as refrerence in ascending order. */
static int
sortByBrightness(const void *a, const void *b)
{
  struct quad_t *p1 = (struct quad_t *)a;
  struct quad_t *p2 = (struct quad_t *)b;

  return ( p1->brightness==p2->brightness
           ? 0
           : (p1->brightness<p2->brightness ? -1 : 1) );
}




/* Structure to keep parameters for multithreaded worker function. */
struct params
{
  gal_data_t *ra;
  gal_data_t *dec;
  gal_data_t *magnitude;
  double     maxStarDisInQuad;
  size_t     num_timesUsed;
  size_t     brightestStarId_arr[26][5];
};



/* Internally used macro to help in the processing */
#define max(a,b) \
   ({ __typeof__ (a) _a = (a);  \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })



void *
make_hashcodes_worker(void *in_prm)
{
  /* Low-level definitions to be done first. */
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct params *p=(struct params *)tprm->params;

  /* Subsequent definitions. */
  size_t i, j, k=0;
  double distance_sq;
  struct quad_t *quads=NULL;
  size_t *brightness_Id=NULL;
  size_t *timesUsed_count_arr=NULL;
  struct quad_t final_quad[4]={0};

  /* Array declerations. */
  double *ra_arr=p->ra->array;
  double *dec_arr=p->dec->array;
  float  *magnitude_arr=p->magnitude->array;

  /* Grid dimentions. */
  size_t row=sizeof(p->brightestStarId_arr)
              /sizeof(p->brightestStarId_arr[0]);
  size_t column=sizeof(p->brightestStarId_arr[0])
                /sizeof(p->brightestStarId_arr[0][0]);
  size_t num_totalStars=sizeof(p->brightestStarId_arr)
                        /sizeof(p->brightestStarId_arr[0][0]);

  /* Make an array that keeps the count of the number
     of times a star is used by increasing the value
     for the index of the star. */
  timesUsed_count_arr=malloc(p->magnitude->size*sizeof(timesUsed_count_arr));

  /* Make a 1d array from the 2d brightestStarId_arr
     for easy use of indexes. */
  brightness_Id=malloc(num_totalStars*sizeof(*brightness_Id));
  for(i=0;i<row;++i)
    for(j=0;j<column;++j)
    {
      brightness_Id[k++]=p->brightestStarId_arr[i][j];
      // printf("brightness ID = %zu, %zu\n", brightness_Id[k], p->brightestStarId_arr[i][j]);
    }


  /* Allocate memory for quads and initialise it. */
  quads=malloc(num_totalStars*sizeof(*quads));
  for (i=0;i<num_totalStars;++i)
    {
      quads[i].index=0;
      quads[i].ra=0;
      quads[i].dec=0;
      quads[i].distance=0;
      quads[i].brightness=0;
    }

  /* For every star,Niners do a complete search and select the stars
     for quads which are in a particular range from it. */
  for (i=0;i<num_totalStars;++i)
    {
      for (j=0,k=0;j<num_totalStars;++j)
        {
          /* Check the number of times the stars is used. */
          if( timesUsed_count_arr[brightness_Id[i]] < p->num_timesUsed )
            /* Search for stars in a particular range. We take stars
              which are at a maximum of `maxStarDisInQuad` distance
              from the current star. */
            if ( ra_arr[brightness_Id[j]] <= ra_arr[brightness_Id[i]]+p->maxStarDisInQuad &&
                ra_arr[brightness_Id[j]] >= ra_arr[brightness_Id[i]]-p->maxStarDisInQuad)
              if ( dec_arr[brightness_Id[j]] <= dec_arr[brightness_Id[i]]+p->maxStarDisInQuad &&
                  dec_arr[brightness_Id[j]] >= dec_arr[brightness_Id[i]]-p->maxStarDisInQuad)
                {
                  /* We find distances from current star to later sort the stars
                    on basis of diatance and form quads based on the distance.
                    As only the magnitude of distance is required and not the
                    absolute value, we use distance square(distance^2 = x*x+y*y)
                    rather than using sqrt() to find actual distance. */
                  distance_sq = ( ra_arr[brightness_Id[j]]-ra_arr[brightness_Id[i]])*
                                ( ra_arr[brightness_Id[j]]-ra_arr[brightness_Id[i]] )
                              +( dec_arr[brightness_Id[j]]-dec_arr[brightness_Id[i]])*
                                ( dec_arr[brightness_Id[j]]-dec_arr[brightness_Id[i]]);

                  /* Store the relevant details for the possible quad member star. */
                  quads[k].index=brightness_Id[j];
                  quads[k].ra=ra_arr[brightness_Id[j]];
                  quads[k].dec=dec_arr[brightness_Id[j]];
                  quads[k].brightness=magnitude_arr[brightness_Id[j]];
                  quads[k].distance=distance_sq;

                  /* Increase the index so that stars are stored sequentially
                    in the structure. */
                  ++k;
                }
        }

        /* Sort on the basis of distance and select stars at relative
            positions of n, n-2 and n-4 from the current star. */
        qsort(quads, k+1, sizeof(struct quad_t), sortByDistance);


        /* For n, n-2 and n-4 positions to be present, at least 6
            elements must be present inside the structure. */
        if ( k>=6 )
          {
            /* Make the final quad and assign relative positions
                of sorted stars to it. */
            final_quad[0] = quads[0];    /* A */
            final_quad[1] = quads[k-4];  /* C */
            final_quad[2] = quads[k-2];  /* D */
            final_quad[3] = quads[k];    /* B */

            /* Increase the count of number of times the stars
                are used. */
            timesUsed_count_arr[quads[0].index]++;
            timesUsed_count_arr[quads[k-4].index]++;
            timesUsed_count_arr[quads[k-2].index]++;
            timesUsed_count_arr[quads[k].index]++;
          }

        /* Reset the quads. */
        for(j=0;j<=k;++j)
          {
            quads[j].index=0;
            quads[j].ra=0;
            quads[j].dec=0;
            quads[j].distance=0;
            quads[j].brightness=0;
          }

    }

  /* Free allocated space. */
  free(quads);
  free(brightness_Id);
  free(timesUsed_count_arr);

  /* Wait for all the other threads to finish, then return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}







int main()
{
  struct params p;
	size_t i,j, index;
	int quitemmap=1;
	size_t minmapsize=-1;
	struct grid in_grid={0};
  size_t brightestStar_arr[26][5]={0};
	size_t x_numbin=5, y_numbin=5, num_stars_per_gpixel=5;

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
  gal_data_t *ra= gal_data_copy_to_new_type (ref, GAL_TYPE_FLOAT64);
  gal_data_t *dec=gal_data_copy_to_new_type (ref->next, GAL_TYPE_FLOAT64);
  gal_data_t *magnitude=gal_data_copy_to_new_type (ref->next->next, GAL_TYPE_FLOAT32);

  /* make a box-grid */
	grid_make(ra, dec, x_numbin, y_numbin, &in_grid);

	double *ra_arr=ra->array;
	double *dec_arr=dec->array;

  /* Find indexes. */
	for(i=0;i<ra->size;++i)
		index = grid_findindex(ra_arr[i], dec_arr[i], &in_grid);

	/* Find the top 5 star ids using structure. */
  find_brightestStars(ra, dec, magnitude, &in_grid, 5, brightestStar_arr);

	/* Total number of quads. */
	size_t num_quads=x_numbin*y_numbin*num_stars_per_gpixel;


  /* Configure p. */
  double maxStarDisInQuad=max(in_grid.step_size[0],
                              in_grid.step_size[1]);
  size_t num_timesUsed=10;

  p.ra=ra;
  p.dec=dec;
  p.magnitude=magnitude;
  p.num_timesUsed=num_timesUsed;
  p.maxStarDisInQuad=maxStarDisInQuad;
  for(i=0;i<26;++i)
    for(j=0;j<5;j++)
      p.brightestStarId_arr[i][j]=brightestStar_arr[i][j];


  /* Spin-off the threads and do the processing on each thread. */
  gal_threads_spin_off(make_hashcodes_worker, &p, 10, 1);

	gal_data_t *cx=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &num_quads, NULL, 0,
	                             minmapsize, quitemmap, "Cx", "position",
															 "relative position in matching coordinates");

	/* For all 9 columns. */
	/* Link the columns
	cx->next = cy;
	cy->next = dx;
	....
	*/

	/* Use multithreading for quad calculation. num_quads will be total no of jobs.
	   structure for option in worker.

		 gindex = qindex/in_grid->size

		 sindex=  qindex%in_grid->size

		 from gindex and sinde and map, you can extract the first star.
		 */

	gal_list_data_free (ref);
}