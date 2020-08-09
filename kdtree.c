#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <error.h>

#include <gnuastro/data.h>
#include <gnuastro/table.h>
#include <gnuastro/blank.h>
#include <gnuastro/pointer.h>
#include <gnuastro/permutation.h>




struct kdtree_params
{
  size_t ndim;
  size_t *input_row;
  gal_data_t **coords;
  uint32_t *left, *right;
};



/* Swap 2 nodes of the tree. */
static void
kdtree_node_swap(struct kdtree_params *p, size_t node1, size_t node2)
{
	uint32_t tmp_left=p->left[node1];
	uint32_t tmp_right=p->right[node1];
	size_t tmp_input_row=p->input_row[node1];

  /* No need to swap same node. */
  if(node1==node2) return;

  // printf("left = %u, right = %u\n", tmp_left, tmp_right);
  p->left[node1]=p->left[node2];
  p->right[node1]=p->right[node2];
  p->input_row[node1]=p->input_row[node2];

  p->left[node2]=tmp_left;
  p->right[node2]=tmp_right;
  p->input_row[node2]=tmp_input_row;

}





/* Find the median to seperate the hyperspace. Instead of randomly
   chossing the media-point, we use `quickselect alogorithm` to find
   median in average complexity of O(n). This also makes nodes partially
   sorted w.r.t each axis.
   See `https://en.wikipedia.org/wiki/Quickselect` for pseudocode and
   more information of the algorithm.

   return : median point(node) that splits the hyperspace.
*/
static uint32_t
kdtree_median_find(struct kdtree_params *p, size_t node_left, size_t node_right,
									 double *coordinate)
{
  double pivot_value;
  size_t i, node_store, node_mid;

  /* False state, return null. */
  if (node_right <= node_left) return GAL_BLANK_UINT32;

  /* If the tree contains only one element, return that element. */
  if (node_right == node_left + 1) return node_left;

  /* The middle value (here used as pivot) between node_left and node_right. */
  node_mid=node_left+(node_right-node_left)/2;

  printf("\n\n%s\n\tnl=%zu, nr=%zu, nm=%zu\n", __func__, node_left, node_right, node_mid); 

  for(i=0;i<p->coords[0]->size; ++i)
    printf("input row[%zu]=%zu  coordinates[%zu]=%f      (before)\n",
           i, p->input_row[i], i, coordinate[p->input_row[i]]);

  /* Loop until the median (for the current axis) is returned. */
  while(1)
    {
      node_store = node_left;
      pivot_value = coordinate[p->input_row[node_mid]];


      printf("\tpv=%f\n", pivot_value);

      // printf("\tBEFORE nl=%zu, nr=%zu, nm=%zu, p->input_row[nr]=%f, p->input_row[nm]=%f\n",
      //        node_left, node_right, node_mid, coordinate[p->input_row[node_right]], 
      //        coordinate[p->input_row[node_mid]]);   

      /* Move the pivot_value to the right/larger end. */
      kdtree_node_swap(p, node_mid, node_right);

      // printf("\tAFTER nl=%zu, nr=%zu, nm=%zu, p->input_row[nr]=%f, p->input_row[nm]=%f\n",
      //        node_left, node_right, node_mid, coordinate[p->input_row[node_right]], 
      //        coordinate[p->input_row[node_mid]]);

      /* Move all nodes smaller than pivot value to the left/smaller side of 
         the nodes. */
      for (i = node_left; i < node_right; ++i)
        if (coordinate[p->input_row[i]] < pivot_value)
          {
            /* Move ith-node to the left/smaller side. */
            kdtree_node_swap(p, i, node_store);

            /* Prepare the place of next smaller node. */
            // printf("\t\tIN COND i=%zu, node_store=%zu, p->input[ns]=%f, pivotV=%f\n",
            //         i, node_store, coordinate[p->input_row[node_store]], pivot_value);
            node_store++;
          }

      /* Move pivot to node_store. */
      kdtree_node_swap(p, node_store, node_right);

      // /* For a check
      printf("\tnode_store=%zu, node_mid=%zu, node_right=%zu\n\n", 
             node_store, node_mid, node_right);
      // */

      /* If median is found, break the loop and return the node_mid. */
      if (node_store == node_mid) break;

      /* Change the left or right node based on the position of node_store
         so we can continue the search/sort in the loop. */
      if (node_store > node_mid) node_right = node_store;
      else                       node_left  = node_store;
    }
  
  for(i=0;i<p->coords[0]->size; ++i)
    printf("input row[%zu]=%zu  coordinates[%zu]=%f      (after)\n",
           i, p->input_row[i], i, coordinate[p->input_row[i]]);
  printf("returned median = %zu\n", node_mid);

  return node_mid;
}





/* Make a kd-tree from a given set of points.
   For tree construction, a median point is selected for each
   axis and the left and right branches are made by comparing
   points based on that axis.

   return : a balanced kd-tree. */
static uint32_t
kdtree_fill_subtrees(struct kdtree_params *p, size_t remaining, 
                     size_t dim_current, size_t node_current)
{
  static size_t counter=0;
  printf("\n%s (START %zu)\nremaining=%zu\n", __func__, ++counter, remaining);

  double *carr=p->coords[dim_current]->array;
  for(size_t i=0;i<p->coords[0]->size; ++i)
    printf("%zu: %-10f, %-10u, %-10u\n",p->input_row[i], carr[p->input_row[i]], p->left[i], p->right[i]);

  printf("\n");
  /* node_median is a counter over the `input_row` array. 
     `input_row` array has the input_row(row number). */
	size_t node_median;

	/* Base criteria for termination of recursion. */
	if(!remaining) return GAL_BLANK_UINT32;

	/* Find the median node. */
	node_median=kdtree_median_find(p, node_current, node_current+remaining,
                                 p->coords[dim_current]->array);


  // /* For a check
	printf("node_current=%zu, node_median=%zu (input[c]=%f input[m]=%f)\n", 
         node_current, node_median,
         carr[p->input_row[node_current]],
         carr[p->input_row[node_median]]);
  // exit(0);
  // */

  /* If median index is present, set the left and right subtree. */
  if(node_median != GAL_BLANK_UINT32)
    {
      /* Increament the current dimension and make sure it doesn't
         exceed the total number of dimensions. */
      dim_current = (dim_current + 1) % p->ndim;

      /* Fill left and right subtrees by calling this functions again. */
      printf("\nSETTING LEFT subtree (in dim %zu, for med=%zu)\n", dim_current, node_median);
      p->left[node_median]  = kdtree_fill_subtrees(p, node_median-node_current,
                                     						    dim_current, node_current);
      printf("\nLEFT NODE[%zu]=%u\n",p->input_row[node_median], p->left[node_median]);

      // if(p->input_row[node_median]==4) exit(0);

      printf("\nSETTING RIGHT subtree (in dim %zu, for med=%zu)\n", dim_current, node_median);
      p->right[node_median] = kdtree_fill_subtrees(p, node_current + remaining - (node_median + 1),
                                     						    dim_current, node_median+1);
      printf("\nRIGHT NODE[%zu]=%u\n",p->input_row[node_median], p->right[node_median]);
    }

  printf("\n%s (END %zu)\n", __func__, counter);

  return p->input_row[node_median];
}






gal_data_t *
gal_kdtree_create(gal_data_t *coords_raw)
{
	size_t i;
  struct kdtree_params p;
	gal_data_t *left, *right, *tmp;

	/* Allocate output. */
	left=gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, coords_raw->dsize, NULL, 0,
	                    coords_raw->minmapsize, coords_raw->quietmmap, "left",
											"index", "index of left subtree in the kd-tree");
	right=gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, coords_raw->dsize, NULL, 0,
	                    coords_raw->minmapsize, coords_raw->quietmmap, "right",
											"index", "index of right subtree in the kd-tree");

	/* Fill the elements of the params structure. */
	left->next=right;
	p.left=left->array;
	p.right=right->array;
	p.ndim=gal_list_data_number(coords_raw);

	/* Allocate and initialise the kd-tree input_row. */
	p.input_row=gal_pointer_allocate(GAL_TYPE_SIZE_T, coords_raw->size, 0,
                                   __func__, "p.input_row");
	for(i=0; i<coords_raw->size; ++i)	p.input_row[i]=i;

	/* Allocate the coordinate array. */
	errno=0;
	p.coords=malloc(p.ndim*sizeof(**(p.coords)));
	if(p.coords==NULL)
		error(EXIT_FAILURE, errno, "%s: couldn't allocate %zu bytes "
					"for 'coords'", __func__, p.ndim*sizeof(**(p.coords)));

	/* Convert input to double type. */
	tmp=coords_raw;
	for(i=0; i<p.ndim; ++i)
		{
			if(tmp->type == GAL_TYPE_FLOAT64)
      	p.coords[i]=tmp;
			else
				p.coords[i]=gal_data_copy_to_new_type (tmp, GAL_TYPE_FLOAT64);

			/* Go to the next column list. */
			tmp=tmp->next;
		}
  
  for(i=0;i<coords_raw->size;++i)
    {
      p.left[i]=GAL_BLANK_UINT32;
      p.right[i]=GAL_BLANK_UINT32;
    }

	/* Fill the kd-tree*/
	kdtree_fill_subtrees(&p, coords_raw->size-1, 0, 0);

  /* Do a reverse permutation to sort the indexes back in the input order. */
  gal_permutation_apply_inverse (left, p.input_row);
  gal_permutation_apply_inverse (right, p.input_row);

	/* Clean up. */
	tmp=coords_raw;
	for(i=0; i<p.ndim; ++i)
		{
			if(p.coords[i]!=tmp) gal_data_free(p.coords[i]);
			tmp=tmp->next;
		}
	free(p.coords);
	free(p.input_row);

	/* Return results. */
	return left;
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







#if 0

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

#endif