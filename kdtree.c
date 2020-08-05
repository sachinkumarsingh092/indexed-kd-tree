#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <error.h>

#include <gnuastro/table.h>
#include <gnuastro/data.h>


struct kdtree_node
{
	size_t index;
	struct kdtree_node *left, *right;
};



/* Find the median to seperate the hyperspace. Instead of randomly
   chossing the media-point, we use `quickselect alogorithm` to find
   median in average complexity of O(n). This also makes nodes partially
   sorted w.r.t each axis.
   See `https://en.wikipedia.org/wiki/Quickselect` for pseudocode and
   more information of the algorithm.

   return : median point(node) that splits the hyperspace.
*/
static struct kdtree_node*
kdtree_median_find(struct kdtree_node *left, struct kdtree_node *right,
									 double *coordinate)
{
  double pivotValue;
  struct kdtree_node *p=NULL, *storeIndex=NULL, *mid=NULL;

  /* False state, return null. */
  if (right <= left) return NULL;

  /* If the tree contains only one element, return that element. */
  if (right == left + 1) return left;

  /* The middle value(here used as pivot) between left and right. */
  mid=left + (right - left) / 2;

  /* Loop until the median(for the current axis) is returned. */
  while(1)
    {
      storeIndex = left;
      pivotValue = coordinate[mid->index];

      /* Select a pivotIndex between left and right. */
      swap(mid, right - 1, MAX_DIM);

      for (p = left; p < right; ++p)
        if (coordinate[p->index] < pivotValue)
          {
            if (p != storeIndex)
                swap(p, storeIndex, MAX_DIM);

            /* Increase the index of partition(pivot). */
            storeIndex++;
          }

      /* Move pivot to its final place. */
      swap(storeIndex, right - 1, MAX_DIM);

      /* If median is found, return it. */
      if (coordinate[storeIndex->index] == coordinate[mid->index])
        return mid;

      /* The pivot is now in its final sorted position. Assign left
         or right based on its current value w.r.t mid. */
      if (storeIndex > mid) right = storeIndex;
      else                  left  = storeIndex;
    }
}


static struct kdtree_node *
kdtree_fill_subtrees(struct kdtree_node *current_node, size_t remaining, 
								 size_t current_axis, gal_data_t **coords)
{
	struct kdtree_node *median=NULL;

	/* Base criteria for termination for recursion. */
	if(!remaining) return NULL;

	median=kdtree_median_find(current_node, current_node+remaining, coords[current_axis]->array);
	

}



gal_data_t *
gal_kdtree_create(gal_data_t *coords_raw)
{
	size_t i;
	struct kdtree_node *nodes;
	size_t *left_arr, *right_arr; 
	gal_data_t *left, *right, **coords, *tmp;
	size_t ndim=gal_list_data_number(coords_raw);

	/* Allocate the coordinate array. */
	errno=0;
	coords=malloc(ndim*sizeof(*coords));
	if(coords==NULL)
		error(EXIT_FAILURE, errno, "%s: couldn't allocate %zu bytes "
					"for 'coords'", __func__, ndim*sizeof(*coords));

	/* Allocate output. */
	left=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, coords_raw->dsize, NULL, 0,
	                    coords_raw->minmapsize, coords_raw->quietmmap, "left", "index",
											"index of left subtree in the kd-tree");
	right=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, coords_raw->dsize, NULL, 0,
	                    coords_raw->minmapsize, coords_raw->quietmmap, "right", "index",
											"index of right subtree in the kd-tree");

	/* Set right to be next of left. */
	left->next=right;

	/* Convert input to double type. */
	tmp=coords_raw;
	for(i=0; i<ndim; ++i)
		{
			if(tmp->type == GAL_TYPE_FLOAT64)
      	coords[i]=tmp;
			else
				coords[i]=gal_data_copy_to_new_type (tmp, GAL_TYPE_FLOAT64);

			/* Go to the next column list. */
			tmp=tmp->next;
		}

	/* Allocate and initialise the kd-tree nodes. */
	errno=0;
	nodes=malloc(coords_raw->size*sizeof(*nodes));
	if(nodes==NULL)
		error(EXIT_FAILURE, errno, "%s: couldn't allocate %zu bytes "
					"for 'nodes'", __func__, coords_raw->size*sizeof(*nodes));
	for(i=0; i<coords_raw->size; ++i)
		{
			nodes[i].index=i;
			nodes[i].left=nodes[i].right=NULL;
		}
	
	/* Fill the kd-tree*/
	kdtree_fill_subtrees(nodes, coords_raw->size, 0, coords);

	/* Write the left and right indexes into the final output. */
	for(i=0; i<coords_raw->size; ++i)
		{
			left_arr[i]=nodes[i].left->index;
			right_arr[i]=nodes[i].right->index;
		}

	/* Clean up. */
	tmp=coords_raw;
	for(i=0; i<ndim; ++i)
		{
			if(coords[i]!=tmp) gal_data_free(coords[i]);
			tmp=tmp->next;
		}
	free(coords);
	free(nodes);

	/* Return results. */
	return left;
}





int main()
{
  char *inputname="quad-out.fits";
  char *outputname="kdtree-out.fits";

	/* Read the input table. */
  gal_data_t *coords=gal_table_read (inputname, "1",
                                  NULL, NULL, GAL_TABLE_SEARCH_NAME,
                                  0, -1, 0, NULL);
	gal_kdtree_create(coords);

	/* Write output. */
  return EXIT_SUCCESS;
}