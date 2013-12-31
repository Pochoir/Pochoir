/*
 * ============================================================================
 *       Filename:  pochoir_modified_cuts.hpp
 *    Description:  Has routines 
 *					1. that implement the modified space/power of two
 *					time cuts.
 *					The code uses the time/space cut code framework in
 *					pochoir_walk_recursive.hpp
 *        Created:  10/02/2012
 *         Author:  Eka Palamadai, epn@mit.edu
 * ============================================================================
 */
#ifndef POCHOIR_MODIFIED_CUTS_HPP 
#define POCHOIR_MODIFIED_CUTS_HPP 

#include "pochoir_common.hpp"
#include "pochoir_walk.hpp"

#define push_queue(_dep, _level, _t0, _t1, _grid) \
do { \
    assert(queue_len_[_dep] < ALGOR_QUEUE_SIZE); \
    circular_queue_[_dep][queue_tail_[_dep]].level = _level; \
    circular_queue_[_dep][queue_tail_[_dep]].t0 = _t0; \
    circular_queue_[_dep][queue_tail_[_dep]].t1 = _t1; \
    circular_queue_[_dep][queue_tail_[_dep]].grid = _grid; \
    ++queue_len_[_dep]; \
    queue_tail_[_dep] = pmod((queue_tail_[_dep] + 1), ALGOR_QUEUE_SIZE); \
} while(0)

#define top_queue(_dep, _queue_elem) \
do { \
    assert(queue_len_[_dep] > 0); \
    _queue_elem = &(circular_queue_[_dep][queue_head_[_dep]]); \
} while(0)

#define pop_queue(_dep) \
do { \
    assert(queue_len_[_dep] > 0); \
    queue_head_[_dep] = pmod((queue_head_[_dep] + 1), ALGOR_QUEUE_SIZE); \
    --queue_len_[_dep]; \
} while(0)

#if defined (TRAP) || !defined (KERNEL_SELECTION)
//Coarsen base case wrt the bottom side
#define COARSEN_BASE_CASE_WRT_BOTTOM_SIDE
#endif

#ifdef COARSEN_BASE_CASE_WRT_BOTTOM_SIDE
//to do : merge the macros below and avoid redundance.

#define CHECK_WIDTH_LB (lb >= 2 * thres && lb > dx_recursive_[level])
#define CHECK_WIDTH_TB (tb >= 2 * thres && lb > dx_recursive_[level])
#define CHECK_WIDTH_LB_BOUNDARY (lb >= 2 * thres && lb > dx_recursive_boundary_[level])
#define CHECK_WIDTH_TB_BOUNDARY (tb >= 2 * thres && lb > dx_recursive_boundary_[level])

//#define CAN_CUT_I (cut_lb ? (lb >= 2 * thres && lb > dx_recursive_[level]) : (tb >= 2 * thres && lb > dx_recursive_[level]))

//#define CAN_CUT_B (cut_lb ? (l_touch_boundary ? (lb >= 2 * thres && lb > dx_recursive_boundary_[level]) : (lb >= 2 * thres && lb > dx_recursive_[level])) : (l_touch_boundary ? (tb >= 2 * thres && lb > dx_recursive_boundary_[level]) : (tb >= 2 * thres && lb > dx_recursive_[level])))

#define CAN_CUT_IN (cut_lb ? (lb >= 2 * thres && lb > dx_recursive_[i]) : (tb >= 2 * thres && lb > dx_recursive_[i]))
#define SIM_CAN_CUT_I (sim_can_cut || CAN_CUT_IN)

#define CAN_CUT_BO (cut_lb ? (l_touch_boundary ? (lb >= 2 * thres && lb > dx_recursive_boundary_[i]) : (lb >= 2 * thres && lb > dx_recursive_[i])) : (l_touch_boundary ? (tb >= 2 * thres && lb > dx_recursive_boundary_[i]) : (tb >= 2 * thres && lb > dx_recursive_[i])))
#define SIM_CAN_CUT_B (sim_can_cut || CAN_CUT_BO)

#else
#define CHECK_WIDTH_LB (lb >= 2 * thres && lb > dx_recursive_[level])
#define CHECK_WIDTH_TB (tb >= 2 * thres && tb > dx_recursive_[level])
#define CHECK_WIDTH_LB_BOUNDARY (lb >= 2 * thres && lb > dx_recursive_boundary_[level])
#define CHECK_WIDTH_TB_BOUNDARY (tb >= 2 * thres && tb > dx_recursive_boundary_[level])

//#define CAN_CUT_I (cut_lb ? (lb >= 2 * thres && lb > dx_recursive_[level]) : (tb >= 2 * thres && tb > dx_recursive_[level]))

//#define CAN_CUT_B (cut_lb ? (l_touch_boundary ? (lb >= 2 * thres && lb > dx_recursive_boundary_[level]) : (lb >= 2 * thres && lb > dx_recursive_[level])) : (l_touch_boundary ? (tb >= 2 * thres && tb > dx_recursive_boundary_[level]) : (tb >= 2 * thres && tb > dx_recursive_[level])))

#define CAN_CUT_IN (cut_lb ? (lb >= 2 * thres && lb > dx_recursive_[i]) : (tb >= 2 * thres && tb > dx_recursive_[i]))
#define SIM_CAN_CUT_I (sim_can_cut || CAN_CUT_IN)

#define CAN_CUT_BO (cut_lb ? (l_touch_boundary ? (lb >= 2 * thres && lb > dx_recursive_boundary_[i]) : (lb >= 2 * thres && lb > dx_recursive_[i])) : (l_touch_boundary ? (tb >= 2 * thres && tb > dx_recursive_boundary_[i]) : (tb >= 2 * thres && tb > dx_recursive_[i])))
#define SIM_CAN_CUT_B (sim_can_cut || CAN_CUT_BO)
#endif

#define CAN_CUT_I (cut_lb ? CHECK_WIDTH_LB : CHECK_WIDTH_TB)
#define CAN_CUT_B (cut_lb ? (l_touch_boundary ? CHECK_WIDTH_LB_BOUNDARY : CHECK_WIDTH_LB) : (l_touch_boundary ? CHECK_WIDTH_TB_BOUNDARY : CHECK_WIDTH_TB))

template <int N_RANK> template <typename F, typename BF>
inline void Algorithm<N_RANK>::space_cut_boundary(int t0, int t1, grid_info<N_RANK> const grid, F const & f, BF const & bf)
{
	//cout << "bhere" << endl ;
	int h = t1 - t0 ;
	//cout << "h " << h << endl ;
	int num_upright_zoids [N_RANK], num_inverted_zoids [N_RANK] ; 
	int offset [N_RANK], two_sigma_h [N_RANK] ;
	int l_slope [N_RANK] ;
	grid_info<N_RANK> grid2 = grid ;
	cout << "done" << endl ; 
	for (int i = N_RANK - 1 ; i >= 0 ; i--)
	{
		int lb = grid.x1[i] - grid.x0[i] ;
		int tb = grid.x1[i] + grid.dx1[i] * h - grid.x0[i] - grid.dx0[i] * h;
		const bool l_touch_boundary = touch_boundary(i, h, grid2);
//#ifndef NDEBUG
		cout << " x0 [" << i << "] " << grid.x0 [i] 
		<< " x1 [" << i << "] " << grid.x1 [i] 
		<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * h
		<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * h
		<< " h " << h << endl ; 
//#endif
		two_sigma_h [i] = 2 * slope_ [i] * h ;
		l_slope [i] = slope_ [i] ;
		//num_upright_zoids [i] = lb / (two_sigma_h [i]) ; 
		if (l_touch_boundary)
		{
			num_upright_zoids [i] = lb / (two_sigma_h [i]) ; 
		}
		else
		{
			num_upright_zoids [i] = max (lb / (two_sigma_h [i]), 
									(lb - dx_recursive_ [i]) / two_sigma_h [i]);
		}
		assert (num_upright_zoids [i] >= 0) ;
		if (lb <= tb)
		{
			//inverted trapezoid
			num_inverted_zoids [i] = num_upright_zoids [i] + 1 ;
			offset [i] = 0 ;
		}
		else
		{
			//upright trapezoid
			num_inverted_zoids [i] = max (num_upright_zoids [i] - 1, 0) ;
			offset [i]  = two_sigma_h [i] ;
		}
		//cout << " num_upright_zoids [" << i << "] " << num_upright_zoids [i] << endl ;
		//cout << " num_inverted_zoids [" << i << "] " << num_inverted_zoids [i] << endl ;
		assert (num_inverted_zoids [i] >= 0) ;

		if (lb == tb)
		{
			//projection trapezoid is a rectangle.
			l_slope [i] = 0 ;
		}
	}
	//coarsen the base case.
	int popcount = 0 ;
	grid_info<N_RANK> subgrid ;
	for (int j = 0 ; j < num_orientations ; j++)
	{
		unsigned int bits = space_cut_sequence [N_RANK - 1] [j] ;
		//cout << "bits " << bits << endl ;
		//find the total number of zoids to process.
		int total_num_zoids = 1 ; 
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
			int inverted = bits & 1 << i ;
			if (inverted)
			{
				total_num_zoids *= num_inverted_zoids [i] ;
			}
			else
			{
				total_num_zoids *= num_upright_zoids [i] ;
			}
		}
		//cout << "total_num_zoids " << total_num_zoids << endl ;
		if (total_num_zoids && __builtin_popcount(bits) != popcount)
		{
			//cout << "sync" << endl ;
			//sync when the popcount of the sequence changes.
			cilk_sync ;
			popcount = __builtin_popcount(bits) ;
		}
		//process each of the zoids
		for (int k = 0 ; k < total_num_zoids ; k++)
		{
			int volume = 1 ;
			//fix the zoid co-ordinates in each dimension
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				int index = -1 ;
				int inverted = bits & 1 << i ; 
				if (inverted)
				{
					subgrid.dx0 [i] = -l_slope [i] ;
					subgrid.dx1 [i] = l_slope [i] ;
					index = (k / volume) % num_inverted_zoids [i] ;
					int mid = num_inverted_zoids [i] / 2 ;
					//cout << "index " << index << endl ;
					if (index < mid)
					{
						//cout << "il " << endl ;
						//zoid lies to the left of midpoint
						subgrid.x0 [i] = grid.x0[i] + index * two_sigma_h [i] +
										offset [i] ;
						subgrid.x1 [i] = subgrid.x0 [i] ;
					}
					else if (index == mid && num_inverted_zoids [i] & 1)
					{
						//cout << "im " << endl ;
						//zoid lies in the middle
						subgrid.x0 [i] = grid.x0[i] + index * two_sigma_h [i] +
										offset [i] ;
						subgrid.x1 [i] = grid.x1[i] - index * two_sigma_h [i] -
										offset [i] ;
					}
					else
					{
						//cout << "ir " << endl ;
						//zoid lies to the right of midpoint
						subgrid.x0 [i] = grid.x1[i] - 
							(num_inverted_zoids [i] - index - 1) * 
							two_sigma_h [i] - offset [i] ;
						subgrid.x1 [i] = subgrid.x0 [i] ;
					}
					volume *= num_inverted_zoids [i] ;
				}
				else
				{
					subgrid.dx0 [i] = l_slope [i] ;
					subgrid.dx1 [i] = -l_slope [i] ;
					index = (k / volume) % num_upright_zoids [i] ;
					int mid = num_upright_zoids [i] / 2 ;
					//cout << "index " << index << endl ;
					if (index < mid)
					{
						//cout << "ul" << endl ;
						//zoid lies to the left of midpoint
						subgrid.x0 [i] = grid.x0[i] + index * two_sigma_h [i] ;
						subgrid.x1 [i] = subgrid.x0 [i] + two_sigma_h [i] ;
					}
					else if (index == mid && num_upright_zoids [i] & 1) 
					{
						//zoid lies in the middle
						//cout << "um" << endl ;
						subgrid.x0 [i] = grid.x0[i] + index * two_sigma_h [i] ;
						subgrid.x1 [i] = grid.x1[i] - index * two_sigma_h [i] ;
					}
					else
					{
						//cout << "ur" << endl ;
						//zoid lies to the right of midpoint
						subgrid.x0 [i] = grid.x1[i] - 
						(num_upright_zoids [i] - index) * two_sigma_h [i] ;
						subgrid.x1 [i] = subgrid.x0 [i] + two_sigma_h [i] ;
					}
					volume *= num_upright_zoids [i] ;
				}
#ifndef NDEBUG
		/*cout << "spawning " ;
		cout << " x0 [" << i << "] " << subgrid.x0 [i] 
		<< " x1 [" << i << "] " << subgrid.x1 [i] 
		<< " x2 [" << i << "] " << subgrid.x0[i] + subgrid.dx0[i] * h
		<< " x3 [" << i << "] " << subgrid.x1[i] + subgrid.dx1[i] * h
		<< " h " << h << endl ;*/
#endif
			}
			cilk_spawn space_time_cut_boundary(t0, t1, subgrid, f, bf) ;
		}
	}
}

template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::space_cut_interior(int t0, int t1, grid_info<N_RANK> const grid, F const & f)
{
	//cout << "here" << endl ;
	int h = t1 - t0 ;
	//cout << "h " << h << endl ;
	int num_upright_zoids [N_RANK], num_inverted_zoids [N_RANK] ; 
	int offset [N_RANK], two_sigma_h [N_RANK] ;
	//cout << "done" << endl ;
	for (int i = N_RANK - 1 ; i >= 0 ; i--)
	{
		int lb = grid.x1[i] - grid.x0[i] ;
		int tb = grid.x1[i] + grid.dx1[i] * h - grid.x0[i] - grid.dx0[i] * h;
#ifndef NDEBUG
		cout << " x0 [" << i << "] " << grid.x0 [i] 
		<< " x1 [" << i << "] " << grid.x1 [i] 
		<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * h
		<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * h
		<< " h " << h << endl ; 
#endif
		two_sigma_h [i] = 2 * slope_ [i] * h ;
		//num_upright_zoids [i] = lb / (two_sigma_h [i]) ; 
		num_upright_zoids [i] = max (lb / (two_sigma_h [i]), 
									(lb - dx_recursive_ [i]) / two_sigma_h [i]);
		//cout << " num_upright_zoids [" << i << "] " << num_upright_zoids [i] << endl ;
		assert (num_upright_zoids [i] >= 0) ;
		if (lb <= tb)
		{
			//inverted trapezoid
			num_inverted_zoids [i] = num_upright_zoids [i] + 1 ;
			offset [i] = 0 ;
		}
		else
		{
			//upright trapezoid
			num_inverted_zoids [i] = num_upright_zoids [i] - 1 ;
			offset [i]  = two_sigma_h [i] ;
		}
		assert (num_inverted_zoids [i] >= 0) ;
		//cout << " num_inverted_zoids [" << i << "] " << num_inverted_zoids [i] << endl ;
	}
	//coarsen the base case.
	grid_info<N_RANK> subgrid ;
	int popcount = 0 ;
	for (int j = 0 ; j < num_orientations ; j++)
	{
		unsigned int bits = space_cut_sequence [N_RANK - 1] [j] ;
		//cout << "bits " << bits << endl ;
		//find the total number of zoids to process.
		int total_num_zoids = 1 ; 
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
			int inverted = bits & 1 << i ;
			if (inverted)
			{
				total_num_zoids *= num_inverted_zoids [i] ;
			}
			else
			{
				total_num_zoids *= num_upright_zoids [i] ;
			}
		}
		//cout << "total_num_zoids " << total_num_zoids << endl ;
		if (total_num_zoids && __builtin_popcount(bits) != popcount)
		{
			//cout << "sync" << endl ;
			//sync when the popcount of the sequence changes.
			cilk_sync ;
			popcount = __builtin_popcount(bits) ;
		}
		//process each of the zoids
		for (int k = 0 ; k < total_num_zoids ; k++)
		{
			int volume = 1 ;
			//fix the zoid co-ordinates in each dimension
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				int index = -1 ;
				int inverted = bits & 1 << i ; 
				if (inverted)
				{
					subgrid.dx0 [i] = -slope_ [i] ;
					subgrid.dx1 [i] = slope_ [i] ;
					index = (k / volume) % num_inverted_zoids [i] ;
					int mid = num_inverted_zoids [i] / 2 ;
					//cout << "index " << index << endl ;
					if (index < mid)
					{
						//cout << "il " << endl ;
						//zoid lies to the left of midpoint
						subgrid.x0 [i] = grid.x0[i] + index * two_sigma_h [i] +
										offset [i] ;
						subgrid.x1 [i] = subgrid.x0 [i] ;
					}
					else if (index == mid && num_inverted_zoids [i] & 1)
					{
						//cout << "im " << endl ;
						//zoid lies in the middle
						subgrid.x0 [i] = grid.x0[i] + index * two_sigma_h [i] +
										offset [i] ;
						subgrid.x1 [i] = grid.x1[i] - index * two_sigma_h [i] -
										offset [i] ;
					}
					else
					{
						//cout << "ir " << endl ;
						//zoid lies to the right of midpoint
						subgrid.x0 [i] = grid.x1[i] - 
							(num_inverted_zoids [i] - index - 1) * 
							two_sigma_h [i] - offset [i] ;
						subgrid.x1 [i] = subgrid.x0 [i] ;
					}
					volume *= num_inverted_zoids [i] ;
				}
				else
				{
					subgrid.dx0 [i] = slope_ [i] ;
					subgrid.dx1 [i] = -slope_ [i] ;
					index = (k / volume) % num_upright_zoids [i] ;
					int mid = num_upright_zoids [i] / 2 ;
					//cout << "index " << index << endl ;
					if (index < mid)
					{
						//cout << "ul" << endl ;
						//zoid lies to the left of midpoint
						subgrid.x0 [i] = grid.x0[i] + index * two_sigma_h [i] ;
						subgrid.x1 [i] = subgrid.x0 [i] + two_sigma_h [i] ;
					}
					else if (index == mid && num_upright_zoids [i] & 1) 
					{
						//zoid lies in the middle
						//cout << "um" << endl ;
						subgrid.x0 [i] = grid.x0[i] + index * two_sigma_h [i] ;
						subgrid.x1 [i] = grid.x1[i] - index * two_sigma_h [i] ;
					}
					else
					{
						//cout << "ur" << endl ;
						//zoid lies to the right of midpoint
						subgrid.x0 [i] = grid.x1[i] - 
						(num_upright_zoids [i] - index) * two_sigma_h [i] ;
						subgrid.x1 [i] = subgrid.x0 [i] + two_sigma_h [i] ;
					}
					volume *= num_upright_zoids [i] ;
				}
#ifndef NDEBUG
		/*cout << "spawning " ;
		cout << " x0 [" << i << "] " << subgrid.x0 [i] 
		<< " x1 [" << i << "] " << subgrid.x1 [i] 
		<< " x2 [" << i << "] " << subgrid.x0[i] + subgrid.dx0[i] * h
		<< " x3 [" << i << "] " << subgrid.x1[i] + subgrid.dx1[i] * h
		<< " h " << h << endl ;*/
#endif
			}
			cilk_spawn space_time_cut_interior(t0, t1, subgrid, f) ;
		}
	}
}

/*
template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::space_cut_interior(int t0, int t1, grid_info<N_RANK> const grid, F const & f)
{
	//cout << "here" << endl ;
	int h = t1 - t0 ;
	int num_upright_zoids [N_RANK], num_inverted_zoids [N_RANK] ; 
	int num_zoids [N_RANK] ;
	int offset [N_RANK], sigma_h [N_RANK] ;
	for (int i = N_RANK - 1 ; i >= 0 ; i--)
	{
		int lb = grid.x1[i] - grid.x0[i] ;
		int tb = grid.x1[i] + grid.dx1[i] * h - grid.x0[i] - grid.dx0[i] * h;
		sigma_h [i] = slope_ [0] * h ;
		num_upright_zoids [i] = lb / (2 * sigma_h [i]) ; 
		if (lb < tb)
		{
			//inverted trapezoid
			num_inverted_zoids [i] = num_upright_zoids [i] + 1 ;
			offset [i] = 0 ;
		}
		else
		{
			//upright trapezoid
			num_inverted_zoids [i] = num_upright_zoids [i] - 1 ;
			offset [i]  = 2 * sigma_h [i] ;
		}
	}
	//seems to be working correctly.
	//coarsen the base case.
	//continue from here
	grid_info<N_RANK> subgrid [num_orientations] ;
	for (int j = 0 ; j < num_orientations ; j++)
	{
		int bits = space_cut_sequence [N_RANK - 1] [j] ;
		//find the total number of zoids to process.
		int total_num_zoids = 1 ; 
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
			num_zoids [i] = 0 ;
			int inverted = bits & 1 << i ;
			if (inverted)
			{
				total_num_zoids *= num_inverted_zoids [i] ;
			}
			else
			{
				total_num_zoids *= num_upright_zoids [i] ;
			}
		}
		int counter = 0 ;
		//process each of the zoids
		//this may not work. need to process all cartesian product of zoids
		//consider enumerating zoids and walking over each of them.
		while (counter < total_num_zoids)
		{
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
			int num_zoids_to_spawn = 1 ;
			int inverted = bits & 1 << i ; 
			int num_zoids_to_consider = 0 ;
			int x0_left = 0, x1_left = 0, x0_right = 0, x1_right = 0 ;
			if (inverted)
			{
				num_zoids_to_consider = num_inverted_zoids [i] ;
				left_end = 
			}
			else
			{
				num_zoids_to_consider = num_upright_zoids [i] ;
				x0_left = grid.x0[i] + num_zoids [i] * sigma_h [i] ;
				x1_left = x0_left + 2 * sigma_h [i] ;

				x1_right = grid.x1[i] - num_zoids [i] * sigma_h [i] ;
				x0_right = x1_right - 2 * sigma_h [i] ;
			
				slope_left = slope_ [i] ;
				slope_right = -slope_ [i] ;
			}
				//process the upright zoids.
				if (num_zoids_to_consider - num_zoids [i] > 1)
				{ 
					for (int k = 0 ; k < num_orientations ; k++)
					{
						subgrid [k].dx0 [i] = slope_left ;
						subgrid [k].dx1 [i] = slope_right ;
						int left = k & 1 << N_RANK - 1 ;
						if (left)
						{
							subgrid [k].x0 [i] = x0_left ;
							subgrid [k].x1 [i] = x1_left ;
						}
						else
						{
							subgrid [k].x0 [i] = x0_right ;
							subgrid [k].x1 [i] = x1_right ;
						}
					}
					num_zoids [i] += 2 ;
					num_zoids_to_spawn *= 2 ;
				}
				else if (num_zoids_to_consider - num_zoids[i] == 1)
				{
					//process the middle upright zoid
					for (int k = 0 ; k < num_orientations / 2 ; k++)
					{
						subgrid [k].dx0 [i] = slope_left ;
						subgrid [k].dx1 [i] = slope_right ;
						subgrid [k].x0 [i] = x0_left ;
						subgrid [k].x1 [i] = x1_right ;
					}
					num_zoids [i] += 1 ;
				}
			//process the upright zoids first	
			for ( ; num_upright_zoids [i] > 1 ; num_upright_zoids [i] -= 2)
			{
				//process the triangle at left
				subgrid.x0 [0] = left_end ;
				left_end = subgrid.x0 [1] = left_end + 2 * sigma_h [i] ;

				cilk_spawn space_time_cut_interior(t0, t1, subgrid, f) ;

				//process the triangle at right
				subgrid.x0 [0] = right_end - 2 * sigma_h [i] ;
				subgrid.x0 [1] = right_end ;
				right_end = subgrid.x0 [0] ;

				cilk_spawn space_time_cut_interior(t0, t1, subgrid, f) ;
			}

			//process the middle trapezoid
			subgrid.x0 [0] = left_end ;
			subgrid.x0 [1] = right_end ;
			if (num_upright_zoids [i] == 1) 
			{
				cilk_spawn space_time_cut_interior(t0, t1, subgrid, f) ;
			}
			cilk_sync ;
			
			left_end = grid.x0[0] + offset [i] ;
			right_end = grid.x1[0] - offset [i] ;
			subgrid.dx0 [0] = -slope_ [0] ;
			subgrid.dx1 [0] = slope_ [0] ;
			//process the inverted zoids next
			for ( ; num_inverted_zoids [i] > 1; num_inverted_zoids [i] -= 2)
			{
				//process the triangle at left
				subgrid.x0 [0] = subgrid.x0 [1] = left_end ;
				left_end += 2 * sigma_h [i] ;

				cilk_spawn space_time_cut_interior(t0, t1, subgrid, f) ;

				//process the triangle at right
				subgrid.x0 [0] = subgrid.x0 [1] = right_end  ;
				right_end -= 2 * sigma_h [i] ;

				cilk_spawn space_time_cut_interior(t0, t1, subgrid, f) ;
			}

			if (num_inverted_zoids [i] == 1)
			{
				subgrid.x0 [0] = left_end ;
				subgrid.x0 [1] = right_end ;
				cilk_spawn space_time_cut_interior(t0, t1, subgrid, f) ;
			}
			cilk_sync ;
		}
		}
	}
}
*/

/*
template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::space_cut_interior(int t0, int t1, grid_info<N_RANK> const grid, F const & f)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }

    // set up the initial grid 
    push_queue(0, N_RANK-1, t0, t1, grid);
    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
        const int curr_dep_pointer = (curr_dep & 0x1);
        while (queue_len_[curr_dep_pointer] > 0) {
            top_queue(curr_dep_pointer, l_father);
            if (l_father->level < 0) {
                // spawn all the grids in circular_queue_[curr_dep][] 
#if USE_CILK_FOR 
                // use cilk_for to spawn all the sub-grid 
// #pragma cilk_grainsize = 1
                cilk_for (int j = 0; j < queue_len_[curr_dep_pointer]; ++j) {
                    int i = pmod((queue_head_[curr_dep_pointer]+j), ALGOR_QUEUE_SIZE);
                    queue_info * l_son = &(circular_queue_[curr_dep_pointer][i]);
                    // assert all the sub-grid has done N_RANK spatial cuts 
                    assert(l_son->level == -1);
                    space_time_cut_interior(l_son->t0, 
						l_son->t1, l_son->grid, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0)
				{
                    space_time_cut_interior(l_father->t0,
 							l_father->t1, l_father->grid, f);
				}
                else
				{
                    cilk_spawn space_time_cut_interior(
						l_father->t0, l_father->t1, l_father->grid, f) ;
				}
#endif
            } else {
                // performing a space cut on dimension 'level' 
                pop_queue(curr_dep_pointer);
                const grid_info<N_RANK> l_father_grid = l_father->grid;
                const int t0 = l_father->t0, t1 = l_father->t1;
                const int lt = (t1 - t0);
                const int level = l_father->level;
                const int thres = slope_[level] * lt;
                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
                const bool cut_lb = (lb < tb);
				const bool can_cut = CAN_CUT_I ;
				if (! can_cut) {
                //if (num_zoids [level] == 0) {
                    // if we can't cut into this dimension, just directly push 
                    // it into the circular queue 
                    //
                    push_queue(curr_dep_pointer, level-1, t0, t1, 
							l_father_grid) ;
                } else  {
                    // can_cut! 
                    if (cut_lb) {
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        l_son_grid.x0[level] = l_start ;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
                        l_son_grid.x1[level] = l_start ;
                        l_son_grid.dx1[level] = slope_[level] ;
                        push_queue(next_dep_pointer, level-1, t0, t1, 
									l_son_grid) ;

                        l_son_grid.x0[level] = l_end ;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end ;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
                        push_queue(next_dep_pointer, level-1, t0, t1, 
									l_son_grid) ;
						//const int cut_more = ((lb - (thres << 2)) >= 0) ;
						tb -= (thres << 2) ;
						//const bool cut_more = (tb >= thres << 1)  && 
						//					  (lb > dx_recursive_[level]) ;
						const bool cut_more = CHECK_WIDTH_TB ;
						if (cut_more)
						{
							const int offset = (thres << 1) ;
							l_son_grid.x0[level] = l_start ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_start + offset ;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
							
							l_son_grid.x0[level] = l_end - offset ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_end ;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							l_son_grid.x0[level] = l_start + offset ;
	                        l_son_grid.dx0[level] = -slope_[level] ;
    	                    l_son_grid.x1[level] = l_end - offset ;
        	                l_son_grid.dx1[level] = slope_[level] ;
            	            push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						}
						else
						{
							l_son_grid.x0[level] = l_start ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_end;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						}
                    } // end if (cut_lb) 
                    else {
                        // cut_tb 
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);
						const int offset = (thres << 1) ;

                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start + offset;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, 
									l_son_grid) ;

                        l_son_grid.x0[level] = l_end - offset ;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, 
									l_son_grid);

                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
						//const int cut_more = ((tb - (thres << 2)) >= 0) ;
						lb -= (thres << 2) ;
						//const bool cut_more = (lb >= thres << 1)  && 
						//					  (lb > dx_recursive_[level]) ;
						const bool cut_more = CHECK_WIDTH_LB ;
						if (cut_more)
						{
							l_son_grid.x0[level] = l_start + offset ;
                            l_son_grid.dx0[level] = -slope_[level];
							l_son_grid.x1[level] = l_start + offset;
                        	l_son_grid.dx1[level] = slope_[level];
                        	push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							l_son_grid.x0[level] = l_end - offset ;
                        	l_son_grid.dx0[level] = -slope_[level];
							l_son_grid.x1[level] = l_end - offset ;
                            l_son_grid.dx1[level] = slope_[level];
                        	push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							l_son_grid.x0[level] = l_start + offset ;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end - offset ;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						}
						else
						{
                        	l_son_grid.x0[level] = l_start + offset ;
                        	l_son_grid.dx0[level] = -slope_[level];
                        	l_son_grid.x1[level] = l_end - offset ;
                        	l_son_grid.dx1[level] = slope_[level];
                        	push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						}
                    } // end else (cut_tb) 
                } // end if (can_cut) 
            } // end if (performing a space cut) 
        } // end while (queue_len_[curr_dep] > 0) 
#if !USE_CILK_FOR
        cilk_sync;
#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}
*/

// modified space cuts. 
/*
template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::space_cut_interior(int t0, int t1, grid_info<N_RANK> const grid, F const & f)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }

    // set up the initial grid 
    push_queue(0, N_RANK-1, t0, t1, grid);
    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
        const int curr_dep_pointer = (curr_dep & 0x1);
        while (queue_len_[curr_dep_pointer] > 0) {
            top_queue(curr_dep_pointer, l_father);
            if (l_father->level < 0) {
                // spawn all the grids in circular_queue_[curr_dep][] 
#if USE_CILK_FOR 
                // use cilk_for to spawn all the sub-grid 
// #pragma cilk_grainsize = 1
                cilk_for (int j = 0; j < queue_len_[curr_dep_pointer]; ++j) {
                    int i = pmod((queue_head_[curr_dep_pointer]+j), ALGOR_QUEUE_SIZE);
                    queue_info * l_son = &(circular_queue_[curr_dep_pointer][i]);
                    // assert all the sub-grid has done N_RANK spatial cuts 
                    assert(l_son->level == -1);
                    space_time_cut_interior(l_son->t0, 
						l_son->t1, l_son->grid, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0)
				{
                    space_time_cut_interior(l_father->t0,
 							l_father->t1, l_father->grid, f);
				}
                else
				{
                    cilk_spawn space_time_cut_interior(
						l_father->t0, l_father->t1, l_father->grid, f) ;
				}
#endif
            } else {
                // performing a space cut on dimension 'level' 
                pop_queue(curr_dep_pointer);
                const grid_info<N_RANK> l_father_grid = l_father->grid;
                const int t0 = l_father->t0, t1 = l_father->t1;
                const int lt = (t1 - t0);
                const int level = l_father->level;
                const int thres = slope_[level] * lt;
                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
                const bool cut_lb = (lb < tb);
				const bool can_cut = CAN_CUT_I ;
				if (! can_cut) {
                //if (num_zoids [level] == 0) {
                    // if we can't cut into this dimension, just directly push 
                    // it into the circular queue 
                    //
                    push_queue(curr_dep_pointer, level-1, t0, t1, 
							l_father_grid) ;
                } else  {
                    // can_cut! 
                    if (cut_lb) {
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        l_son_grid.x0[level] = l_start ;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
                        l_son_grid.x1[level] = l_start ;
                        l_son_grid.dx1[level] = slope_[level] ;
                        push_queue(next_dep_pointer, level-1, t0, t1, 
									l_son_grid) ;

                        l_son_grid.x0[level] = l_end ;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end ;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
                        push_queue(next_dep_pointer, level-1, t0, t1, 
									l_son_grid) ;
						//const int cut_more = ((lb - (thres << 2)) >= 0) ;
						tb -= (thres << 2) ;
						//const bool cut_more = (tb >= thres << 1)  && 
						//					  (lb > dx_recursive_[level]) ;
						const bool cut_more = CHECK_WIDTH_TB ;
						if (cut_more)
						{
							const int offset = (thres << 1) ;
							l_son_grid.x0[level] = l_start ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_start + offset ;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
							
							l_son_grid.x0[level] = l_end - offset ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_end ;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							l_son_grid.x0[level] = l_start + offset ;
	                        l_son_grid.dx0[level] = -slope_[level] ;
    	                    l_son_grid.x1[level] = l_end - offset ;
        	                l_son_grid.dx1[level] = slope_[level] ;
            	            push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						}
						else
						{
							l_son_grid.x0[level] = l_start ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_end;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						}
                    } // end if (cut_lb) 
                    else {
                        // cut_tb 
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);
						const int offset = (thres << 1) ;

                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start + offset;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, 
									l_son_grid) ;

                        l_son_grid.x0[level] = l_end - offset ;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, 
									l_son_grid);

                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
						//const int cut_more = ((tb - (thres << 2)) >= 0) ;
						lb -= (thres << 2) ;
						//const bool cut_more = (lb >= thres << 1)  && 
						//					  (lb > dx_recursive_[level]) ;
						const bool cut_more = CHECK_WIDTH_LB ;
						if (cut_more)
						{
							l_son_grid.x0[level] = l_start + offset ;
                            l_son_grid.dx0[level] = -slope_[level];
							l_son_grid.x1[level] = l_start + offset;
                        	l_son_grid.dx1[level] = slope_[level];
                        	push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							l_son_grid.x0[level] = l_end - offset ;
                        	l_son_grid.dx0[level] = -slope_[level];
							l_son_grid.x1[level] = l_end - offset ;
                            l_son_grid.dx1[level] = slope_[level];
                        	push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							l_son_grid.x0[level] = l_start + offset ;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end - offset ;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						}
						else
						{
                        	l_son_grid.x0[level] = l_start + offset ;
                        	l_son_grid.dx0[level] = -slope_[level];
                        	l_son_grid.x1[level] = l_end - offset ;
                        	l_son_grid.dx1[level] = slope_[level];
                        	push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						}
                    } // end else (cut_tb) 
                } // end if (can_cut) 
            } // end if (performing a space cut) 
        } // end while (queue_len_[curr_dep] > 0) 
#if !USE_CILK_FOR
        cilk_sync;
#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}
*/

/* Boundary space cut. Uses modified space cut.
 */

template <int N_RANK> template <typename F, typename BF>
inline void Algorithm<N_RANK>::space_cut_boundary_initial(int t0, int t1, grid_info<N_RANK> const grid, F const & f, BF const & bf)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }
	
    // set up the initial grid 
    push_queue(0, N_RANK-1, t0, t1, grid);
    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
        const int curr_dep_pointer = (curr_dep & 0x1);
        while (queue_len_[curr_dep_pointer] > 0) {
            top_queue(curr_dep_pointer, l_father);
            if (l_father->level < 0) {
                // spawn all the grids in circular_queue_[curr_dep][] 
#if USE_CILK_FOR 
                // use cilk_for to spawn all the sub-grid 
// #pragma cilk_grainsize = 1
				
                cilk_for (int j = 0; j < queue_len_[curr_dep_pointer]; ++j) {
                    int i = pmod((queue_head_[curr_dep_pointer]+j), ALGOR_QUEUE_SIZE);
                    queue_info * l_son = &(circular_queue_[curr_dep_pointer][i]);
                    // assert all the sub-grid has done N_RANK spatial cuts 
                    //assert(l_son->level == -1);
                    space_time_cut_boundary(l_son->t0, 
						l_son->t1, l_son->grid, f, bf);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0) {
                    space_time_cut_boundary(l_father->t0,
						l_father->t1, l_father->grid, f, bf);
                } else {
                    cilk_spawn space_time_cut_boundary(
						l_father->t0, l_father->t1, l_father->grid, f, bf);
                }
#endif
            } else {
                // performing a space cut on dimension 'level' 
                pop_queue(curr_dep_pointer);
                grid_info<N_RANK> l_father_grid = l_father->grid;
                const int t0 = l_father->t0, t1 = l_father->t1;
                const int lt = (t1 - t0);
                const int level = l_father->level;
                const int thres = slope_[level] * lt;
                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
                const bool cut_lb = (lb < tb);
                const bool l_touch_boundary = touch_boundary(level, lt, l_father_grid);
                //const bool l_touch_boundary = true ;
				const bool can_cut = CAN_CUT_B ;
                if (! can_cut) {
                //if (num_zoids [level] == 0) {
                    // if we can't cut into this dimension, just directly push
                    // it into the circular queue
                    //
                    push_queue(curr_dep_pointer, level-1, t0, t1, 
							l_father_grid) ;
                } else  {
                    /* can_cut */
                    if (cut_lb) {
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        l_son_grid.x0[level] = l_start ;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
                        l_son_grid.x1[level] = l_start ;
                        l_son_grid.dx1[level] = slope_[level] ;
                        push_queue(next_dep_pointer, level-1, t0, t1, 
											l_son_grid) ;

                        l_son_grid.x0[level] = l_end ;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end ;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
                        push_queue(next_dep_pointer, level-1, t0, t1, 
											l_son_grid) ;

						//const int cut_more = ((lb - (thres << 2)) >= 0) ;
						tb -= (thres << 2) ;
						//const bool cut_more = (tb >= thres << 1)  && 
						//				  (lb > dx_recursive_[level]) ;
						const bool cut_more = CHECK_WIDTH_TB ;
						if (cut_more)
						{
							const int offset = (thres << 1) ;
							l_son_grid.x0[level] = l_start ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_start + offset;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
									t1, l_son_grid) ;

							
							l_son_grid.x0[level] = l_end - offset ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_end ;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
									t1, l_son_grid) ;

							l_son_grid.x0[level] = l_start + offset;
	                        l_son_grid.dx0[level] = -slope_[level];
    	                    l_son_grid.x1[level] = l_end - offset ;
        	                l_son_grid.dx1[level] = slope_[level] ;
            	            push_queue(next_dep_pointer, level-1, t0, 
									t1, l_son_grid);
						}
						else
						{
							l_son_grid.x0[level] = l_start ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_end;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
									t1, l_son_grid);
						}
                    } /* end if (cut_lb) */
                    else { /* cut_tb */
                        if (lb == phys_length_[level] && l_father_grid.dx0[level] == 0 && l_father_grid.dx1[level] == 0) { /* initial cut on the dimension */
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
							//const int cut_more = ((lb - (thres << 2)) >= 0) ;
							tb -= (thres << 1) ;
							//const bool cut_more = (tb >= thres << 1)  && 
							//		  (lb > dx_recursive_boundary_[level]) ;
							const bool cut_more = CHECK_WIDTH_TB_BOUNDARY ;
							if (cut_more)
							{
								const int offset = (thres << 1) ;
                            	l_son_grid.x0[level] = l_start ;
                            	l_son_grid.dx0[level] = slope_[level];
                            	l_son_grid.x1[level] = l_start + offset ;
                            	l_son_grid.dx1[level] = -slope_[level];
                            	push_queue(curr_dep_pointer, level-1, 
									t0, t1, l_son_grid) ;

                            	l_son_grid.x0[level] = l_end - offset ;
                            	l_son_grid.dx0[level] = slope_[level];
                            	l_son_grid.x1[level] = l_end ;
                            	l_son_grid.dx1[level] = -slope_[level];
                            	push_queue(curr_dep_pointer, level-1, 
									t0, t1, l_son_grid) ;

                            	l_son_grid.x0[level] = l_start + offset ;
                            	l_son_grid.dx0[level] = -slope_[level];
                            	l_son_grid.x1[level] = l_end - offset ;
                            	l_son_grid.dx1[level] = slope_[level];
                            	push_queue(next_dep_pointer, level-1, 
									t0, t1, l_son_grid);
							}
							else
							{
                            	l_son_grid.x0[level] = l_start ;
                            	l_son_grid.dx0[level] = slope_[level];
                            	l_son_grid.x1[level] = l_end ;
                            	l_son_grid.dx1[level] = -slope_[level];
                            	push_queue(curr_dep_pointer, level-1, 
									t0, t1, l_son_grid);
							}
                            l_son_grid.x0[level] = l_end ;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = l_end ;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						} else  {
							/* cut_tb */
							grid_info<N_RANK> l_son_grid = l_father_grid;
							const int l_start = (l_father_grid.x0[level]);
							const int l_end = (l_father_grid.x1[level]);
							const int offset = (thres << 1) ;

							l_son_grid.x0[level] = l_start;
							l_son_grid.dx0[level] = l_father_grid.dx0[level];
							l_son_grid.x1[level] = l_start + offset;
							l_son_grid.dx1[level] = -slope_[level];
							push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							l_son_grid.x0[level] = l_end - offset ;
							l_son_grid.dx0[level] = slope_[level];
							l_son_grid.x1[level] = l_end;
							l_son_grid.dx1[level] = l_father_grid.dx1[level];
							push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							const int next_dep_pointer = (curr_dep + 1) & 0x1;
							//const int cut_more = ((tb - (thres << 2)) >= 0) ;
							lb -= (thres << 2) ;
							//const bool cut_more = (lb >= thres << 1)  && 
							//			  (lb > dx_recursive_[level]) ;
							const bool cut_more = CHECK_WIDTH_LB ;
							if (cut_more)
							{
								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = -slope_[level];
								l_son_grid.x1[level] = l_start + offset;
								l_son_grid.dx1[level] = slope_[level];
								push_queue(next_dep_pointer, level-1, 
									t0, t1, l_son_grid) ;

								l_son_grid.x0[level] = l_end - offset ;
								l_son_grid.dx0[level] = -slope_[level];
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = slope_[level];
								push_queue(next_dep_pointer, level-1, 
									t0, t1, l_son_grid) ;

								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = slope_[level];
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = -slope_[level];
								push_queue(curr_dep_pointer, level-1, 
									t0, t1, l_son_grid);
							}
							else
							{
								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = -slope_[level];
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = slope_[level];
								push_queue(next_dep_pointer, level-1, 
									t0, t1, l_son_grid);
							}
						}                   
                    } /* end if (cut_tb) */
                }// end if (can_cut) 
            } // end if (performing a space cut) 
        } // end while (queue_len_[curr_dep] > 0) 
#if !USE_CILK_FOR
        cilk_sync;
#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}

/* This is the version for interior region cut! */
template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::space_time_cut_interior(int t0, int t1, grid_info<N_RANK> const grid, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;
 
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);

        bool cut_lb = (lb < tb);
        thres = (slope_[i] * lt);
        sim_can_cut = SIM_CAN_CUT_I ;
#ifndef NDEBUG
		/*cout << " x0 [" << i << "] " << grid.x0 [i] 
		<< " x1 [" << i << "] " << grid.x1 [i] 
		<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
		<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
		<< " lt " << lt << endl ;*/
#endif
    }
    if (sim_can_cut) 
	{
        /* cut into space */
        space_cut_interior(t0, t1, grid, f) ;
    } 
	else if (lt > dt_recursive_) 
	{
        /* cut into time */
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
        space_time_cut_interior(t0, t0+halflt, 
									l_son_grid, f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        space_time_cut_interior(t0+halflt, t1, 
								l_son_grid, f);
    }
	else 
	{
        f(t0, t1, grid);
	}
}

/* This is the version for boundary region cut! */
template <int N_RANK> template <typename F, typename BF>
inline void Algorithm<N_RANK>::space_time_cut_boundary(int t0, int t1, grid_info<N_RANK> const grid, F const & f, BF const & bf)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;
	bool initial_cut = true ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);

        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);

		if (grid.dx0[i] != 0 || grid.dx1[i] != 0)
		{
			initial_cut = false ;
		}

        thres = slope_[i] * lt ;
		bool cut_lb = (lb < tb);
		sim_can_cut = SIM_CAN_CUT_B ;
        call_boundary |= l_touch_boundary;
#ifndef NDEBUG
		/*cout << " x0 [" << i << "] " << grid.x0 [i] 
		<< " x1 [" << i << "] " << grid.x1 [i] 
		<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
		<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
		<< " lt " << lt << endl ;*/
#endif
    }
    if (call_boundary)
	{
        l_dt_stop = dt_recursive_boundary_;
	}
    else
	{
        l_dt_stop = dt_recursive_;
	}
    if (sim_can_cut) 
	{
        //cut into space 
        if (call_boundary) 
		{
			//cout << "space cut b" << endl ;
			if (initial_cut)
			{
            	space_cut_boundary_initial(t0, t1, l_father_grid, f, bf) ;
			}
			else
			{
            	space_cut_boundary(t0, t1, l_father_grid, f, bf) ;
			}
        }
		else
		{
			//cout << "space cut i" << endl ;
            space_cut_interior(t0, t1, l_father_grid, 
										f) ;
		}
    } 
	else if (lt > l_dt_stop)  //time cut
	{
		//cout << "time cut " << endl ;
        // cut into time 
		
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
        if (call_boundary) {
            space_time_cut_boundary(t0, t0+halflt, 
								l_son_grid, f, bf);
        } else {
            space_time_cut_interior(t0, t0+halflt, 
								l_son_grid, f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            space_time_cut_boundary(t0+halflt, t1, 
								l_son_grid, f, bf);
        } else {
            space_time_cut_interior(t0+halflt, t1, 
								l_son_grid, f);
        }
    } 
	else
	{
        if (call_boundary) {
            base_case_kernel_boundary(t0, t1, l_father_grid, bf);
        } else {
            f(t0, t1, l_father_grid);
        }
	}
}

/* This is the version for boundary region cut! */
template <int N_RANK> template <typename F, typename BF>
inline void Algorithm<N_RANK>::abnormal_region_space_time_cut_boundary(int t0, int t1, grid_info<N_RANK> const grid, F const & f, BF const & bf)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);

        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);

        thres = slope_[i] * lt ;
		bool cut_lb = (lb < tb);
		sim_can_cut = SIM_CAN_CUT_B ;
        call_boundary |= l_touch_boundary;
#ifndef NDEBUG
		/*cout << " x0 [" << i << "] " << grid.x0 [i] 
		<< " x1 [" << i << "] " << grid.x1 [i] 
		<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
		<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
		<< " lt " << lt << endl ;*/
#endif
    }
    if (call_boundary)
	{
        l_dt_stop = dt_recursive_boundary_;
	}
    else
	{
        l_dt_stop = dt_recursive_;
	}
    if (sim_can_cut) 
	{
		//cout << "space cut " << endl ;
        //cut into space 
        if (call_boundary) 
		{
            abnormal_region_space_cut_boundary(t0, t1, l_father_grid, 
										f, bf) ;
        }
		else
		{
            abnormal_region_space_cut_interior(t0, t1, l_father_grid, 
										f) ;
		}
    } 
	else if (lt > l_dt_stop)  //time cut
	{
		//cout << "time cut " << endl ;
        // cut into time 
		
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
        if (call_boundary) {
            abnormal_region_space_time_cut_boundary(t0, t0+halflt, 
								l_son_grid, f, bf);
        } else {
            abnormal_region_space_time_cut_interior(t0, t0+halflt, 
								l_son_grid, f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            abnormal_region_space_time_cut_boundary(t0+halflt, t1, 
								l_son_grid, f, bf);
        } else {
            abnormal_region_space_time_cut_interior(t0+halflt, t1, 
								l_son_grid, f);
        }
    } 
	else
	{
        if (call_boundary) {
            base_case_kernel_boundary(t0, t1, l_father_grid, bf);
        } else {
            f(t0, t1, l_father_grid);
        }
	}
}

/* This is the version for interior region cut! */
template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::abnormal_region_space_time_cut_interior(int t0, int t1, grid_info<N_RANK> const grid, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;
 
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);

        bool cut_lb = (lb < tb);
        thres = (slope_[i] * lt);
        sim_can_cut = SIM_CAN_CUT_I ;
#ifndef NDEBUG
		/*cout << " x0 [" << i << "] " << grid.x0 [i] 
		<< " x1 [" << i << "] " << grid.x1 [i] 
		<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
		<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
		<< " lt " << lt << endl ;*/
#endif
    }
    if (sim_can_cut) 
	{
        /* cut into space */
        abnormal_region_space_cut_interior(t0, t1, grid, f) ;
    } 
	else if (lt > dt_recursive_) 
	{
        /* cut into time */
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
        abnormal_region_space_time_cut_interior(t0, t0+halflt, 
									l_son_grid, f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        abnormal_region_space_time_cut_interior(t0+halflt, t1, 
								l_son_grid, f);
    }
	else 
	{
        f(t0, t1, grid);
	}
}


template <int N_RANK> template <typename F, typename BF>
inline void Algorithm<N_RANK>::abnormal_region_space_cut_boundary(int t0, int t1, grid_info<N_RANK> const grid, F const & f, BF const & bf)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }
	
    // set up the initial grid 
    push_queue(0, N_RANK-1, t0, t1, grid);
    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
        const int curr_dep_pointer = (curr_dep & 0x1);
        while (queue_len_[curr_dep_pointer] > 0) {
            top_queue(curr_dep_pointer, l_father);
            if (l_father->level < 0) {
                // spawn all the grids in circular_queue_[curr_dep][] 
#if USE_CILK_FOR 
                // use cilk_for to spawn all the sub-grid 
// #pragma cilk_grainsize = 1
				
                cilk_for (int j = 0; j < queue_len_[curr_dep_pointer]; ++j) {
                    int i = pmod((queue_head_[curr_dep_pointer]+j), ALGOR_QUEUE_SIZE);
                    queue_info * l_son = &(circular_queue_[curr_dep_pointer][i]);
                    // assert all the sub-grid has done N_RANK spatial cuts 
                    //assert(l_son->level == -1);
                    abnormal_region_space_time_cut_boundary(l_son->t0, 
						l_son->t1, l_son->grid, f, bf);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0) {
                    abnormal_region_space_time_cut_boundary(l_father->t0,
						l_father->t1, l_father->grid, f, bf);
                } else {
                    //cilk_spawn abnormal_region_space_time_cut_boundary(
                    abnormal_region_space_time_cut_boundary(
						l_father->t0, l_father->t1, l_father->grid, f, bf);
                }
#endif
            } else {
                // performing a space cut on dimension 'level' 
                pop_queue(curr_dep_pointer);
                grid_info<N_RANK> l_father_grid = l_father->grid;
                const int t0 = l_father->t0, t1 = l_father->t1;
                const int lt = (t1 - t0);
                const int level = l_father->level;
                const int thres = slope_[level] * lt;
                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
                const bool cut_lb = (lb < tb);
                const bool l_touch_boundary = touch_boundary(level, lt, l_father_grid);
                //const bool l_touch_boundary = true ;
				const bool can_cut = CAN_CUT_B ;
                if (! can_cut) {
                //if (num_zoids [level] == 0) {
                    // if we can't cut into this dimension, just directly push
                    // it into the circular queue
                    //
                    push_queue(curr_dep_pointer, level-1, t0, t1, 
							l_father_grid) ;
                } else  {
                    /* can_cut */
                    if (cut_lb) {
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        l_son_grid.x0[level] = l_start ;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
                        l_son_grid.x1[level] = l_start ;
                        l_son_grid.dx1[level] = slope_[level] ;
                        push_queue(next_dep_pointer, level-1, t0, t1, 
											l_son_grid) ;

                        l_son_grid.x0[level] = l_end ;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end ;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
                        push_queue(next_dep_pointer, level-1, t0, t1, 
											l_son_grid) ;

						//const int cut_more = ((lb - (thres << 2)) >= 0) ;
						tb -= (thres << 2) ;
						//const bool cut_more = (tb >= thres << 1)  && 
						//				  (lb > dx_recursive_[level]) ;
						const bool cut_more = CHECK_WIDTH_TB ;
						if (cut_more)
						{
							const int offset = (thres << 1) ;
							l_son_grid.x0[level] = l_start ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_start + offset;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
									t1, l_son_grid) ;

							
							l_son_grid.x0[level] = l_end - offset ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_end ;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
									t1, l_son_grid) ;

							l_son_grid.x0[level] = l_start + offset;
	                        l_son_grid.dx0[level] = -slope_[level];
    	                    l_son_grid.x1[level] = l_end - offset ;
        	                l_son_grid.dx1[level] = slope_[level] ;
            	            push_queue(next_dep_pointer, level-1, t0, 
									t1, l_son_grid);
						}
						else
						{
							l_son_grid.x0[level] = l_start ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_end;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
									t1, l_son_grid);
						}
                    } /* end if (cut_lb) */
                    else { /* cut_tb */
                        if (lb == phys_length_[level] && l_father_grid.dx0[level] == 0 && l_father_grid.dx1[level] == 0) { /* initial cut on the dimension */
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
							//const int cut_more = ((lb - (thres << 2)) >= 0) ;
							tb -= (thres << 1) ;
							//const bool cut_more = (tb >= thres << 1)  && 
							//		  (lb > dx_recursive_boundary_[level]) ;
							const bool cut_more = CHECK_WIDTH_TB_BOUNDARY ;
							if (cut_more)
							{
								const int offset = (thres << 1) ;
                            	l_son_grid.x0[level] = l_start ;
                            	l_son_grid.dx0[level] = slope_[level];
                            	l_son_grid.x1[level] = l_start + offset ;
                            	l_son_grid.dx1[level] = -slope_[level];
                            	push_queue(curr_dep_pointer, level-1, 
									t0, t1, l_son_grid) ;

                            	l_son_grid.x0[level] = l_end - offset ;
                            	l_son_grid.dx0[level] = slope_[level];
                            	l_son_grid.x1[level] = l_end ;
                            	l_son_grid.dx1[level] = -slope_[level];
                            	push_queue(curr_dep_pointer, level-1, 
									t0, t1, l_son_grid) ;

                            	l_son_grid.x0[level] = l_start + offset ;
                            	l_son_grid.dx0[level] = -slope_[level];
                            	l_son_grid.x1[level] = l_end - offset ;
                            	l_son_grid.dx1[level] = slope_[level];
                            	push_queue(next_dep_pointer, level-1, 
									t0, t1, l_son_grid);
							}
							else
							{
                            	l_son_grid.x0[level] = l_start ;
                            	l_son_grid.dx0[level] = slope_[level];
                            	l_son_grid.x1[level] = l_end ;
                            	l_son_grid.dx1[level] = -slope_[level];
                            	push_queue(curr_dep_pointer, level-1, 
									t0, t1, l_son_grid);
							}
                            l_son_grid.x0[level] = l_end ;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = l_end ;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						} else  {
							/* cut_tb */
							grid_info<N_RANK> l_son_grid = l_father_grid;
							const int l_start = (l_father_grid.x0[level]);
							const int l_end = (l_father_grid.x1[level]);
							const int offset = (thres << 1) ;

							l_son_grid.x0[level] = l_start;
							l_son_grid.dx0[level] = l_father_grid.dx0[level];
							l_son_grid.x1[level] = l_start + offset;
							l_son_grid.dx1[level] = -slope_[level];
							push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							l_son_grid.x0[level] = l_end - offset ;
							l_son_grid.dx0[level] = slope_[level];
							l_son_grid.x1[level] = l_end;
							l_son_grid.dx1[level] = l_father_grid.dx1[level];
							push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							const int next_dep_pointer = (curr_dep + 1) & 0x1;
							//const int cut_more = ((tb - (thres << 2)) >= 0) ;
							lb -= (thres << 2) ;
							//const bool cut_more = (lb >= thres << 1)  && 
							//			  (lb > dx_recursive_[level]) ;
							const bool cut_more = CHECK_WIDTH_LB ;
							if (cut_more)
							{
								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = -slope_[level];
								l_son_grid.x1[level] = l_start + offset;
								l_son_grid.dx1[level] = slope_[level];
								push_queue(next_dep_pointer, level-1, 
									t0, t1, l_son_grid) ;

								l_son_grid.x0[level] = l_end - offset ;
								l_son_grid.dx0[level] = -slope_[level];
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = slope_[level];
								push_queue(next_dep_pointer, level-1, 
									t0, t1, l_son_grid) ;

								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = slope_[level];
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = -slope_[level];
								push_queue(curr_dep_pointer, level-1, 
									t0, t1, l_son_grid);
							}
							else
							{
								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = -slope_[level];
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = slope_[level];
								push_queue(next_dep_pointer, level-1, 
									t0, t1, l_son_grid);
							}
						}                   
                    } /* end if (cut_tb) */
                }// end if (can_cut) 
            } // end if (performing a space cut) 
        } // end while (queue_len_[curr_dep] > 0) 
#if !USE_CILK_FOR
        //cilk_sync;
#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}


template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::abnormal_region_space_cut_interior(int t0, int t1, grid_info<N_RANK> const grid, F const & f)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }

    // set up the initial grid 
    push_queue(0, N_RANK-1, t0, t1, grid);
    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
        const int curr_dep_pointer = (curr_dep & 0x1);
        while (queue_len_[curr_dep_pointer] > 0) {
            top_queue(curr_dep_pointer, l_father);
            if (l_father->level < 0) {
                // spawn all the grids in circular_queue_[curr_dep][] 
#if USE_CILK_FOR 
                // use cilk_for to spawn all the sub-grid 
// #pragma cilk_grainsize = 1
                cilk_for (int j = 0; j < queue_len_[curr_dep_pointer]; ++j) {
                    int i = pmod((queue_head_[curr_dep_pointer]+j), ALGOR_QUEUE_SIZE);
                    queue_info * l_son = &(circular_queue_[curr_dep_pointer][i]);
                    // assert all the sub-grid has done N_RANK spatial cuts 
                    assert(l_son->level == -1);
                    abnormal_region_space_time_cut_interior(l_son->t0, 
						l_son->t1, l_son->grid, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0)
				{
                    abnormal_region_space_time_cut_interior(l_father->t0,
 							l_father->t1, l_father->grid, f);
				}
                else
				{
                    //cilk_spawn abnormal_region_space_time_cut_interior(
                    abnormal_region_space_time_cut_interior(
						l_father->t0, l_father->t1, l_father->grid, f) ;
				}
#endif
            } else {
                // performing a space cut on dimension 'level' 
                pop_queue(curr_dep_pointer);
                const grid_info<N_RANK> l_father_grid = l_father->grid;
                const int t0 = l_father->t0, t1 = l_father->t1;
                const int lt = (t1 - t0);
                const int level = l_father->level;
                const int thres = slope_[level] * lt;
                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
                const bool cut_lb = (lb < tb);
				const bool can_cut = CAN_CUT_I ;
				if (! can_cut) {
                //if (num_zoids [level] == 0) {
                    // if we can't cut into this dimension, just directly push 
                    // it into the circular queue 
                    //
                    push_queue(curr_dep_pointer, level-1, t0, t1, 
							l_father_grid) ;
                } else  {
                    /* can_cut! */
                    if (cut_lb) {
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        l_son_grid.x0[level] = l_start ;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
                        l_son_grid.x1[level] = l_start ;
                        l_son_grid.dx1[level] = slope_[level] ;
                        push_queue(next_dep_pointer, level-1, t0, t1, 
									l_son_grid) ;

                        l_son_grid.x0[level] = l_end ;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end ;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
                        push_queue(next_dep_pointer, level-1, t0, t1, 
									l_son_grid) ;
						//const int cut_more = ((lb - (thres << 2)) >= 0) ;
						tb -= (thres << 2) ;
						//const bool cut_more = (tb >= thres << 1)  && 
						//					  (lb > dx_recursive_[level]) ;
						const bool cut_more = CHECK_WIDTH_TB ;
						if (cut_more)
						{
							const int offset = (thres << 1) ;
							l_son_grid.x0[level] = l_start ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_start + offset ;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
							
							l_son_grid.x0[level] = l_end - offset ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_end ;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							l_son_grid.x0[level] = l_start + offset ;
	                        l_son_grid.dx0[level] = -slope_[level] ;
    	                    l_son_grid.x1[level] = l_end - offset ;
        	                l_son_grid.dx1[level] = slope_[level] ;
            	            push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						}
						else
						{
							l_son_grid.x0[level] = l_start ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_end;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						}
                    } /* end if (cut_lb) */
                    else {
                        /* cut_tb */
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);
						const int offset = (thres << 1) ;

                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start + offset;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, 
									l_son_grid) ;

                        l_son_grid.x0[level] = l_end - offset ;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, 
									l_son_grid);

                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
						//const int cut_more = ((tb - (thres << 2)) >= 0) ;
						lb -= (thres << 2) ;
						//const bool cut_more = (lb >= thres << 1)  && 
						//					  (lb > dx_recursive_[level]) ;
						const bool cut_more = CHECK_WIDTH_LB ;
						if (cut_more)
						{
							l_son_grid.x0[level] = l_start + offset ;
                            l_son_grid.dx0[level] = -slope_[level];
							l_son_grid.x1[level] = l_start + offset;
                        	l_son_grid.dx1[level] = slope_[level];
                        	push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							l_son_grid.x0[level] = l_end - offset ;
                        	l_son_grid.dx0[level] = -slope_[level];
							l_son_grid.x1[level] = l_end - offset ;
                            l_son_grid.dx1[level] = slope_[level];
                        	push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;

							l_son_grid.x0[level] = l_start + offset ;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end - offset ;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						}
						else
						{
                        	l_son_grid.x0[level] = l_start + offset ;
                        	l_son_grid.dx0[level] = -slope_[level];
                        	l_son_grid.x1[level] = l_end - offset ;
                        	l_son_grid.dx1[level] = slope_[level];
                        	push_queue(next_dep_pointer, level-1, t0, 
								t1, l_son_grid) ;
						}
                    } /* end else (cut_tb) */
                } // end if (can_cut) 
            } // end if (performing a space cut) 
        } // end while (queue_len_[curr_dep] > 0) 
#if !USE_CILK_FOR
        //cilk_sync;
#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}

template <int N_RANK> template <typename F, typename BF>
inline void Algorithm<N_RANK>::power_of_two_time_cut(int t0, int t1, 
				grid_info<N_RANK> const grid, F const & f, BF const & bf) 
{
	int T = t1 - t0 ;
	int W = 0 ;  //max_width among all dimensions
	int slope ;
	for (int i = 0 ; i < N_RANK ; i++)
	{
		if (phys_length_ [i] > W)
		{
			W = phys_length_ [i] ;
			slope = slope_ [i] ;
		}		
	}
	//cout << " ALGOR_QUEUE_SIZE " << ALGOR_QUEUE_SIZE << endl ;
	//find index of most significant bit that is set
	int Wn = W / (slope << 1) ;
	int index_msb = (sizeof(int) << 3) - __builtin_clz(Wn) - 1 ;
	//h1 = 2^floor(lg(Wn)). The zoid with height h1 undergoes a space cut.
	int h1 = 1 << index_msb ;
	for (int i = 0 ; i < T / h1 ; i++)
	{	
		cout << "t0 " << t0 << " t1 " << t1 << 
			" h1 " << h1 << " t0 + h1 " <<
			t0 + h1 << endl ;
		space_time_cut_boundary(t0, t0 + h1, grid, f, bf) ;
		t0 += h1 ;
	}

	int h2 = t1 - t0 ;
	//time cuts happen only if height > dt_recursive_
	//while (h2 > dt_recursive_)
	while (h2 >= 1 )
	{
		//find index of most significant bit that is set
		index_msb = (sizeof(int) << 3) - __builtin_clz(h2) - 1 ;
		int h = 1 << index_msb ;
		cout << "t0 " << t0 << " t1 " << t1 << 
			" h " << h << " t0 + h " <<
			t0 + h << endl ;
		space_time_cut_boundary(t0, t0 + h, grid, f, bf) ;
		t0 += h ;
		h2 = t1 - t0 ;
	}
	/*while (h2 > 1)
	{
		//find index of most significant bit that is set
		index_msb = (sizeof(int) << 3) - __builtin_clz(h2) - 1 ;
		int h = 1 << index_msb ;
		cout << "t0 " << t0 << " t1 " << t1 << 
			" h " << h << " t0 + h " <<
			t0 + h << endl ;
		bool abnormal = false ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			int n = dx_recursive_ [i] / ((slope_ [i] * h) << 1) ;
			num_triangles [i] = max(1, n) ;
		}
		abnormal_region_space_time_cut_boundary(t0, t0 + h, grid, f, bf) ;
		//space_time_cut_boundary(t0, t0 + h, grid, f, bf) ;
		t0 += h ;
		h2 = t1 - t0 ;
	}
	if (h2 == 1)
	{
		cout << "h = 1 t0 " << t0 << " t1 " << t1 << 
			 " t0 + h " << t0 + h2 << endl ;
		//base_case_kernel_boundary(t0, t0 + h2, grid, bf);
		shorter_duo_sim_obase_bicut_p(t0, t0 + h2, grid, f, bf) ;
	}*/
}
#endif
