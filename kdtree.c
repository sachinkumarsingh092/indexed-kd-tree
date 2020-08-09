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
static size_t
kdtree_median_find(struct kdtree_params *p, size_t node_left,
		   size_t node_right, double *coordinate)
{
  double pivot_value;
  size_t i, store_i, node_pivot, node_k;

  /* False state, this is a programming error. */
  if (node_right < node_left)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us to fix the problem! "
	  "For some reason, the node_right is smaller than node_left", __func__);

  /* If the two nodes are the same, just return the node. */
  if (node_right == node_left) return node_left;

  /* The middle value (here used as pivot) between node_left and node_right. */
  node_k=node_left+(node_right-node_left)/2;

  /* For a check. */
  printf("\n\n%s\n\tnl=%zu, nr=%zu, nm=%zu\n", __func__, node_left,
	 node_right, node_k);
  for(i=0;i<p->coords[0]->size; ++i)
    printf("input row[%zu]=%zu  coordinates[%zu]=%f      (before)\n",
           i, p->input_row[i], i, coordinate[p->input_row[i]]);

  /* Loop until the median (for the current axis) is returned (and in
     the mean-while, the underlying indexs are sorted). */
  while(1)
    {
      /* In every run, 'node_left' and 'node_right' change, so the
	 pivot index and pivot value change. */
      node_pivot = node_left+(node_right-node_left)/2;
      pivot_value = coordinate[p->input_row[node_pivot]];
      printf("\tpv=%f\n", pivot_value);

      // printf("\tBEFORE nl=%zu, nr=%zu, nm=%zu, p->input_row[nr]=%f, p->input_row[nm]=%f\n",
      //        node_left, node_right, node_k, coordinate[p->input_row[node_right]],
      //        coordinate[p->input_row[node_k]]);

      /* Move the pivot_value to the right/larger end. */
      kdtree_node_swap(p, node_pivot, node_right);

      // printf("\tAFTER nl=%zu, nr=%zu, nm=%zu, p->input_row[nr]=%f, p->input_row[nm]=%f\n",
      //        node_left, node_right, node_k, coordinate[p->input_row[node_right]],
      //        coordinate[p->input_row[node_k]]);

      /* Move all nodes smaller than pivot value to the left/smaller side of
         the nodes. */
      store_i=node_left;
      for (i = node_left; i < node_right; ++i)
        if (coordinate[p->input_row[i]] < pivot_value)
          {
            /* Move ith-node to the left/smaller side, and increment store_i */
            kdtree_node_swap(p, store_i, i);
	    store_i++;

            /* Prepare the place of next smaller node. */
            // printf("\t\tIN COND i=%zu, node_pivot=%zu, p->input[ns]=%f, pivotV=%f\n",
            //         i, node_pivot, coordinate[p->input_row[node_pivot]], pivot_value);
          }

      /* Move pivot, to be just after of all the nodes that are less
	 than it (because pivot was moved to 'node_right'). */
      kdtree_node_swap(p, node_right, store_i);

      /* Set node_pivot to be the store_i. */
      node_pivot=store_i;

      // /* For a check
      printf("\tnode_pivot=%zu, node_k=%zu, node_right=%zu\n\n",
             node_pivot, node_k, node_right);
      // */

      /* If median is found, break the loop and return the node_k. */
      if (node_k == node_pivot) break;

      /* Change the left or right node based on the position of node_pivot
         so we can continue the search/sort in the loop. */
      if (node_k < node_pivot) node_right = node_pivot - 1;
      else                     node_left  = node_pivot + 1;
    }

  /* For a check. */
  for(i=0;i<p->coords[0]->size; ++i)
    printf("input row[%zu]=%zu  coordinates[%zu]=%f      (after)\n",
           i, p->input_row[i], i, coordinate[p->input_row[i]]);
  printf("returned median = %zu\n", node_k);
  /**/

  /* Return the pivot node. */
  return node_pivot;
}





/* Make a kd-tree from a given set of points.
   For tree construction, a median point is selected for each
   axis and the left and right branches are made by comparing
   points based on that axis.

   return : a balanced kd-tree. */
static uint32_t
kdtree_fill_subtrees(struct kdtree_params *p, size_t node_left,
		     size_t node_right, size_t depth)
{
  /* Set the working axis. */
  size_t axis=depth % p->ndim;

  /********DEBUGGING*********/
  static size_t counter=0;
  printf("\n%s (START %zu)\nl-r=%zu\n", __func__, ++counter, node_right-node_left);

  double *carr=p->coords[axis]->array;
  for(size_t i=0;i<p->coords[0]->size; ++i)
    printf("%zu: %-10f, %-10u, %-10u %s\n",p->input_row[i], carr[p->input_row[i]],
	   p->left[i], p->right[i],
	   i>=node_left && i<=node_right ? "#": "");
  printf("\n");
  /*************************/

  /* node_median is a counter over the `input_row` array.
     `input_row` array has the input_row(row number). */
  size_t node_median;

  /* Recursion terminates when the left and right nodes are the
     same. */
  if(node_left==node_right) return GAL_BLANK_UINT32;

  /* Find the median node. */
  node_median=kdtree_median_find(p, node_left, node_right,
                                 p->coords[axis]->array);

  /* For a check */
  printf("node_left=%zu, node_median=%zu (input[c]=%f input[m]=%f)\n",
         node_left, node_median,
         carr[p->input_row[node_left]],
         carr[p->input_row[node_median]]);
  // exit(0);
  /**/

  /* If median index is present, set the left and right subtree. */
  if(node_median != GAL_BLANK_UINT32)
    {
      /* Fill left and right subtrees by calling this functions again. */
      printf("\nSETTING LEFT subtree (in axis %zu, for node_median=%zu)\n", axis, node_median);
      p->left[p->input_row[node_median]] = kdtree_fill_subtrees(p, node_left, node_median, depth+1);
      printf("\nLEFT NODE[%zu]=%u\n",p->input_row[node_median],
	     p->left[p->input_row[node_median]]);

      printf("\nSETTING RIGHT subtree (in axis %zu, for node_median=%zu)\n", axis, node_median);
      p->right[p->input_row[node_median]] = kdtree_fill_subtrees(p, node_median+1, node_right, depth+1);
      printf("\nRIGHT NODE[%zu]=%u\n",p->input_row[node_median],
	     p->right[p->input_row[node_median]]);
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

  /* Allocate output and initialize them. */
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
  for(i=0;i<coords_raw->size;++i)
    { p.left[i]=p.right[i]=GAL_BLANK_UINT32; }

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

  /* Fill the kd-tree*/
  kdtree_fill_subtrees(&p, 0, coords_raw->size-1, 0);

  /* For a check before permutation. */
  for(i=0;i<coords_raw->size;++i)
    printf("%-15zu%-15u%-15u\n", p.input_row[i], p.left[i], p.right[i]);

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
kdtree_distance_find(struct kdtree_node *a, struct kdtree_node *b,
		     gal_data_t **coords, size_t dim)
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
gal_kdtree_nearest_neighbour(struct kdtree_node *current_node,
			     struct kdtree_node *point,
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
