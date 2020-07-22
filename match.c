#include <stdio.h>

#include <gnuastro/qsort.h>
#include <gnuastro/table.h>
#include <gnuastro/statistics.h>



struct grid
{
	size_t size;
	size_t dim[2];
	double min[2], step_size[2];
};




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

	// printf("ingrid_min1 = %f, ingid_step_size = %f\n", in_grid->min[1], in_grid->step_size[1]);

	// printf("x = %zu y = %zu in_grid = %zu\n", x, y, in_grid->dim[0]);

	return y*in_grid->dim[0]+x;
}




int main()
{
	size_t i;
	int quitemmap=1;
	size_t minmapsize=-1;
	struct grid in_grid={0};
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

	for(i=0;i<ra->size;++i)
	{
		size_t index = grid_findindex(ra_arr[i], dec_arr[i], &in_grid);
		printf("index = %zu ra = %f, dec = %f\n", index, ra_arr[i], dec_arr[i]);
	}

	/* Find the top 5 star ids using structure. */

	/* Total number of quads. */
	size_t num_quads=x_numbin*y_numbin*num_stars_per_gpixel;

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