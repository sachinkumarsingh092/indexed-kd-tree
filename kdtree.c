#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <error.h>

#include <gnuastro/data.h>
#include <gnuastro/table.h>
#include <gnuastro/blank.h>


struct kdtree_node
{
	uint32_t index;
	struct kdtree_node *left, *right;
};




/* Swap 2 nodes of the tree. */
static void
kdtree_node_swap(struct kdtree_node *x, struct kdtree_node *y)
{
	uint32_t tmp_index=x->index;
	struct kdtree_node *tmp_left=x->left;
	struct kdtree_node *tmp_right=x->right;

	x->index=y->index;
	x->left=y->left;
	x->right=y->right;

	y->index=tmp_index;
	y->left=tmp_left;
	y->right=tmp_right;
}





/* Find the distance between 2 nodes of the tree.

   return : distance(squared) between 2 nodes of the tree.
*/
static double
kdtree_distance_find(struct kdtree_node *a, struct kdtree_node *b, gal_data_t **coords, size_t dim)
{
  size_t i;
  double t_dis, dis = 0;

  /* For all points in the node. */
  for(i=0; i<dim; ++i)
    {
      t_dis = (	(double *)coords[i]->array)[a->index]
								-((double *)coords[i]->array)[b->index];

      /* As only the relative magnitude of distance is required,
         it is efficient to return distance square rather than using
         `sqrt()` function to calculate actual distace value. */
      dis += t_dis * t_dis;
    }

  return dis;
}





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
      kdtree_node_swap(mid, right - 1);

      for (p = left; p < right; ++p)
        if (coordinate[p->index] < pivotValue)
          {
            if (p != storeIndex)
                kdtree_node_swap(p, storeIndex);

            /* Increase the index of partition(pivot). */
            storeIndex++;
          }

      /* Move pivot to its final place. */
      kdtree_node_swap(storeIndex, right - 1);

      /* If median is found, return it. */
      if (coordinate[storeIndex->index] == coordinate[mid->index])
        return mid;

      /* The pivot is now in its final sorted position. Assign left
         or right based on its current value w.r.t mid. */
      if (storeIndex > mid) right = storeIndex;
      else                  left  = storeIndex;
    }
}





/* Make a kd-tree from a given set of points.
   For tree construction, a median point is selected for each
   axis and the left and right branches are made by comparing
   points based on that axis.

   return : a balanced kd-tree. */
static struct kdtree_node *
kdtree_fill_subtrees(struct kdtree_node *current_node, size_t remaining, 
								     size_t current_axis, gal_data_t **coords, size_t dim)
{
	struct kdtree_node *median_node=NULL;

	/* Base criteria for termination of recursion. */
	if(!remaining) return NULL;

	/* Find the median node. */
	median_node=kdtree_median_find(current_node, current_node+remaining, coords[current_axis]->array);
	
	/* If median node is present, recursively make left and right
     subtree. */
  if(median_node)
    {
      current_axis = (current_axis + 1) % dim;
      median_node->left  = kdtree_fill_subtrees(current_node,
                                     						median_node-current_node,
                                     						current_axis, coords, dim);
      median_node->right = kdtree_fill_subtrees(median_node + 1,
                                     						current_node + remaining - (median_node + 1),
                                     						current_axis, coords, dim);
    }

  return median_node;

}






gal_data_t *
gal_kdtree_create(gal_data_t *coords_raw)
{
	size_t i;
	struct kdtree_node *nodes;
	uint32_t *left_arr, *right_arr; 
	gal_data_t *left, *right, **coords, *tmp;
	size_t ndim=gal_list_data_number(coords_raw);


	/* Allocate the coordinate array. */
	errno=0;
	coords=malloc(ndim*sizeof(**coords));
	if(coords==NULL)
		error(EXIT_FAILURE, errno, "%s: couldn't allocate %zu bytes "
					"for 'coords'", __func__, ndim*sizeof(**coords));


	/* Allocate output. */
	left=gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, coords_raw->dsize, NULL, 0,
	                    coords_raw->minmapsize, coords_raw->quietmmap, "left",
											"index", "index of left subtree in the kd-tree");
	right=gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, coords_raw->dsize, NULL, 0,
	                    coords_raw->minmapsize, coords_raw->quietmmap, "right",
											"index", "index of right subtree in the kd-tree");

	/* Set right to be next of left and array pointer. */
	left->next=right;
	left_arr=left->array;
	right_arr=right->array;

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
	kdtree_fill_subtrees(nodes, coords_raw->size, 0, coords, ndim);

	/* Write the left and right indexes into the final output. */
	for(i=0; i<coords_raw->size; ++i)
		{
			// printf("i=%zu, %s, %s\n", i, nodes[i].left?"Allocated":"Not", nodes[i].right?"Allocated":"Not");
			left_arr[nodes[i].index]=nodes[i].left ? nodes[i].left->index : GAL_BLANK_UINT32;
			right_arr[nodes[i].index]=nodes[i].right ? nodes[i].right->index : GAL_BLANK_UINT32;
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





/* Find the nearest neighbour of the `point`.
   See `https://en.wikipedia.org/wiki/K-d_tree#Nearest_neighbour_search`
   for more information.

  return:
  nearest : The nearest node to the search point.
  least_dist : The distance from the search point to the nearest point.
  */
uint32_t
gal_kdtree_nearest_neighbour(struct kdtree_node *current_node, struct kdtree_node *point,
												     size_t current_axis, size_t dim,
												     struct kdtree_node **nearest, double *least_dist,
												     double *coordinate, gal_data_t **coords)
{
  double d, dx, dx_square;

  /* If no tree/subtree present, don't search further. */
  if(!current_node) return GAL_BLANK_UINT32;

  /* The distance between search point to the current nearest.*/
  d = kdtree_distance_find(current_node, point, coords, dim);

  /* Distance between the splitting coordinate of the search
     point and current node*/
  dx = coordinate[current_node->index] - coordinate[point->index];
  dx_square = dx*dx;

  /* Check if the current node is nearer than the previous
     nearest node. */
  if(!*nearest || d < *least_dist)
    {
      *least_dist = d;
      *nearest = current_node;
    }

  /* If exact match found(least distance 0), return it. */
  if(!*least_dist) return (*nearest)->index;

  current_axis = (current_axis + 1) % dim;

  /* Recursively search in subtrees. */
  gal_kdtree_nearest_neighbour(dx > 0 ? current_node->left : current_node->right, point,
										           current_axis, dim, nearest, least_dist,
										           coordinate, coords);

  /* Since the hyperplanes are all axis-aligned, for checking
  if there is a node in other branch that is nearer to the
  search node,we do a simple comparison to see whether the
  distance between the splitting coordinate of the search
  point and current node is lesser(i.e on same side of hyperplane)
  than the distance (overall coordinates) from the search point to
  the current nearest. */
  if(dx_square >= *least_dist) return GAL_BLANK_UINT32;
  gal_kdtree_nearest_neighbour(dx > 0 ? current_node->right : current_node->left, point,
									             current_axis, dim, nearest, least_dist,
										           coordinate, coords);

	return (*nearest)->index;
}






int main()
{
  char *inputname="dummytable.txt";
  char *outputname="kdtree-out.fits";

	gal_data_t *output;
	/* Read the input table. */
  gal_data_t *coords=gal_table_read (inputname, "1",
                                  NULL, NULL, GAL_TABLE_SEARCH_NAME,
                                  0, -1, 0, NULL);

	/* Construct a tree. */
	output=gal_kdtree_create(coords);

	/* Write output. */
	gal_table_write(output, NULL, GAL_TABLE_FORMAT_BFITS,
                  outputname, "kdtree-info", 0);

  return EXIT_SUCCESS;
}