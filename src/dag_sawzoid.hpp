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
#ifndef DAG_SAWZOID_HPP 
#define DAG_SAWZOID_HPP 

#include "dag_header.hpp"

#define dx_recursive_boundary_  (m_algo.dx_recursive_boundary_)
#define dx_recursive_ (m_algo.dx_recursive_)
#define dt_recursive_boundary_ (m_algo.dt_recursive_boundary_)
#define dt_recursive_ (m_algo.dt_recursive_)
#define slope_ m_algo.slope_
#define touch_boundary m_algo.touch_boundary
#define phys_length_ m_algo.phys_length_
#define base_case_kernel_boundary m_algo.base_case_kernel_boundary
#define uub_boundary m_algo.uub_boundary
#define ulb_boundary m_algo.ulb_boundary
#define lub_boundary m_algo.lub_boundary

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::symbolic_sawzoid_space_cut_interior(int t0,
int t1, grid_info<N_RANK> const & grid, unsigned long parent_index, 
F const & f, int * num_zoids, double & redundant_time, double & projected_time)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }

	int child_index = 0 ;
	unsigned int num_children_per_level = 0 ;
    // set up the initial grid 
    push_queue(0, N_RANK-1, t0, t1, grid);
    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
		num_children_per_level = child_index ;
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
                    symbolic_sawzoid_space_time_cut_interior(l_son->t0, 
						l_son->t1, l_son->grid, parent_index, child_index, 
						redundant_time, projected_time, f);
					child_index++ ;
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
				//dependency++ ;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				symbolic_sawzoid_space_time_cut_interior(l_father->t0, 
					l_father->t1, l_father->grid, parent_index, child_index, 
					redundant_time, projected_time, f);
				child_index++ ;
				//dependency++ ;
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
                //if (num_zoids [level] == 0) {
                if (! can_cut) {
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
						if (num_zoids [level] == 5)
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
						if (num_zoids [level] == 5)
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
//        cilk_sync;
#endif
		num_children_per_level = child_index - num_children_per_level ;
		//if (dependency > 0)
		//{
		//assert (dependency) ;
		assert (num_children_per_level <= ALGOR_QUEUE_SIZE) ;
		assert (num_children_per_level <= m_zoids [parent_index].num_children) ;
		//m_zoids [parent_index].dependency <<= zoid_type::NUM_BITS_DEPENDENCY;
		//m_zoids [parent_index].dependency |= num_children_per_level ;
		m_zoids [parent_index].dependency |= num_children_per_level << 
									curr_dep * zoid_type::NUM_BITS_DEPENDENCY ;
		assert (m_zoids [parent_index].dependency) ;
		//}
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
	assert (child_index == m_zoids [parent_index].num_children) ;
}

/* Boundary space cut. Uses sawzoid space cut.
 */
template <int N_RANK> template <typename F, typename BF> 
inline void auto_tune<N_RANK>::symbolic_sawzoid_space_cut_boundary(int t0,
		int t1, grid_info<N_RANK> const & grid, unsigned long parent_index, 
		F const & f, BF const & bf, int * num_zoids, double & redundant_time, 
		double & projected_time)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }
	
	int child_index = 0 ;
	unsigned int num_children_per_level = 0 ;
    // set up the initial grid 
    push_queue(0, N_RANK-1, t0, t1, grid);
    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
		num_children_per_level = child_index ;
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
                    symbolic_sawzoid_space_time_cut_boundary(l_son->t0, 
						l_son->t1, l_son->grid, parent_index, child_index, 
						redundant_time, projected_time, f, bf);
					child_index++ ; //this can be a race.
					//dependency++ ;
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				symbolic_sawzoid_space_time_cut_boundary(l_father->t0, 
					l_father->t1, l_father->grid, parent_index, child_index,
					redundant_time, projected_time, f, bf);
				child_index++ ; 
				//dependency++ ;
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

						if (num_zoids [level] == 5)
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
							if (num_zoids [level] == 4)
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
							if (num_zoids [level] == 5)
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
//#if !USE_CILK_FOR
//        cilk_sync;
//#endif
		num_children_per_level = child_index - num_children_per_level ;
		//if (dependency > 0)
		//{
		//assert (dependency) ;
		assert (num_children_per_level <= ALGOR_QUEUE_SIZE) ;
		assert (num_children_per_level <= m_zoids [parent_index].num_children) ;
		//m_zoids [parent_index].dependency <<= zoid_type::NUM_BITS_DEPENDENCY ;
		//m_zoids [parent_index].dependency |= num_children_per_level ;
		m_zoids [parent_index].dependency |= num_children_per_level << 
									curr_dep * zoid_type::NUM_BITS_DEPENDENCY ;
		assert (m_zoids [parent_index].dependency) ;
		//}
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
	assert (child_index == m_zoids [parent_index].num_children) ;
}

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::symbolic_sawzoid_space_time_cut_interior(
		int t0, int t1, grid_info<N_RANK> const & grid, 
		unsigned long parent_index, int child_index, 
		double & redundant_time, double & projected_time, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;

	struct timespec start, end;
	struct timespec start1, end1 ;
	clock_gettime(CLOCK_MONOTONIC, &start);

	int centroid = 0, width = 1 ; 
	int offset = 0, total_num_subzoids = 1 ;
	int num_subzoids [N_RANK] ;
	unsigned long key = 0 ;
	//unsigned short decision = 0 ;
	unsigned char decision = 0 ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);

		assert (grid.x0[i] < phys_length_ [i]);
		assert (grid.x1[i] < phys_length_ [i]);

		centroid = pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) * width + 
					centroid ;
		width = phys_length_ [i] ;

		unsigned long dim_key = (unsigned long) lb << num_bits_width | tb ;
		key = key << offset | dim_key ;
		offset += num_bits_dim ;
		/*cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
			<< " lt " << lt << endl ;*/
        thres = (slope_[i] * lt);
		int short_side ;
		bool space_cut = false ;
		if (lb < tb)
		{
			short_side = lb ;
			//set if projection trapezoid is inverted
			//decision |= 1 << i + 1 + N_RANK ;
		}
		else
		{
			short_side = tb ;
		}
		num_subzoids [i] = 0 ;
		if (short_side >= (thres << 1) && lb > dx_recursive_[i])
		{
			space_cut = true ;
			//set if a space cut can happen in dimension i
			//decision |= 1 << i + 1 ;
			num_subzoids [i] = 3 ;
			if (short_side - (thres << 2) >= 0)
			{
				num_subzoids [i]= 5 ;
				//set if space cut yields 5 pieces
				//decision |= 1 << i + 1 + 2 * N_RANK ;
			}
			total_num_subzoids *= num_subzoids [i] ;
		}
        sim_can_cut |= space_cut ;
    }
	unsigned long index ;
	bool projection_exists = check_and_create_projection (key, lt, 
											centroid, index, grid) ;
	zoid_type & z = m_zoids [index];
	zoid_type & parent = m_zoids [parent_index] ;
	//add the zoid as a child of the parent
	parent.add_child(&z, child_index, index) ;
	assert (index < m_zoids.size()) ;
	
	if (projection_exists)
	{ 
		//Add the projected time of the zoid to that of the parent.
		clock_gettime(CLOCK_MONOTONIC, &end) ;
		redundant_time += tdiff2(&end, &start) ;
		projected_time += z.time  ;
		//a zoid with the projection already exists. return
		return ;
	}
	//print_dag() ;

	bool divide_and_conquer = false ;
	double necessary_time = 0, projected_time1 = 0 ;
    if (sim_can_cut) 
	{
		//assert (decision) ;
		divide_and_conquer = true ; 
		double elapsed_time = 0, rtime = 0, ptime = 0 ;
		z.resize_children(total_num_subzoids) ;
        /* cut into space */
		clock_gettime(CLOCK_MONOTONIC, &start1);
		symbolic_sawzoid_space_cut_interior(t0, t1, grid, index, f, 
							num_subzoids, rtime, ptime) ;
		clock_gettime(CLOCK_MONOTONIC, &end1);
		elapsed_time = tdiff2(&end1, &start1) - rtime ;
		assert (elapsed_time >= 0.) ;
		assert (ptime >= 0.) ;

		necessary_time = elapsed_time ;
		projected_time1 = ptime ;
		decision = 2 ;
		m_zoids [index].decision |= decision ;
#ifndef NDEBUG
		m_zoids [index].stime = elapsed_time + ptime ;
#endif
    }
	else if (lt > dt_recursive_)
	{
		divide_and_conquer = true ;
		z.resize_children(2) ;
		double elapsed_time = 0, rtime = 0, ptime = 0 ;
        /* cut into time */
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
		clock_gettime(CLOCK_MONOTONIC, &start1);
        symbolic_sawzoid_space_time_cut_interior(t0, t0+halflt, l_son_grid, 
				index, 0, rtime, ptime, f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        symbolic_sawzoid_space_time_cut_interior(t0+halflt, t1, l_son_grid, 
				index, 1, rtime, ptime, f);

		clock_gettime(CLOCK_MONOTONIC, &end1) ;
		elapsed_time = tdiff2(&end1, &start1) - rtime ;
		assert (elapsed_time >= 0.) ;
		assert (ptime >= 0.) ;

		//if space cut didn't happen or 
		//   space cut happened and time for time cut < time for space cut
		//then reset the necessary and projected times.
		if ((sim_can_cut && 
			elapsed_time + ptime < necessary_time + projected_time1)
			|| !sim_can_cut)
		{
			necessary_time = elapsed_time ;
			projected_time1 = ptime ;
			decision = 1 ;
			m_zoids [index].decision |= decision ;
		}
#ifndef NDEBUG
		m_zoids [index].ttime = elapsed_time + ptime ;
#endif
    }
	else
	{
    	//base case
		assert (m_zoids [index].decision == 0) ;
		z.resize_children(1) ;
		//add the grid as a child of the leaf
		m_zoids [index].add_child((zoid_type *) 0, 0, m_leaves.size()) ;
		m_leaves.push_back(grid) ;
	}
	clock_gettime(CLOCK_MONOTONIC, &end);
	double total_time = tdiff2(&end, &start) ;
	redundant_time += total_time - necessary_time ;
}

template <int N_RANK> template <typename F, typename BF>
inline void auto_tune<N_RANK>::symbolic_sawzoid_space_time_cut_boundary(
int t0, int t1,	grid_info<N_RANK> const & grid, 
unsigned long parent_index, int child_index, 
double & redundant_time, double & projected_time, F const & f, BF const & bf)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;

	int centroid = 0, width = 1 ; 
	int offset = 0, total_num_subzoids = 1 ;
	int num_subzoids [N_RANK] ;
	unsigned long key = 0 ;
	//unsigned short decision = 0 ;
	unsigned char decision = 0 ;

	struct timespec start, end;
	struct timespec start1, end1 ;
	clock_gettime(CLOCK_MONOTONIC, &start);

    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
		if (! l_touch_boundary)
		{
			assert (l_father_grid.x0 [i] < phys_length_ [i]) ;
			assert (l_father_grid.x1 [i] < phys_length_ [i]) ;
		}
        lb = (grid.x1[i] - grid.x0[i]);
		centroid = pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) * width + 
					centroid ;
		width = phys_length_ [i] ;
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
		/*cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
			<< " lt " << lt << endl ;*/
        thres = slope_[i] * lt ;
		bool initial_space_cut = false ;
        if (lb == phys_length_[i] && grid.dx0[i] == 0 && grid.dx1[i] == 0) 
		{ 
			//set if initial cut on the dimension 
			//decision |= 1 << i + 1 + 3 * N_RANK ;
			initial_space_cut = true ;
		}
		int short_side ;
		bool space_cut = false ;
		if (lb < tb)
		{
			short_side = lb ;
			//set if projection trapezoid is inverted
			//decision |= 1 << i + 1 + N_RANK ;
		}
		else
		{
			short_side = tb ;
		}
		int limit ;
		if (l_touch_boundary)
		{
			limit = dx_recursive_boundary_[i] ;
		}
		else
		{
			limit = dx_recursive_[i] ;
		}
		num_subzoids [i] = 0 ;
		if (short_side >= (thres << 1) && lb > limit)
		{
			space_cut = true ;
			//set if a space cut can be done in dimension i
			//decision |= 1 << i + 1 ;
			num_subzoids [i] = initial_space_cut ? 2 : 3 ;
			if (short_side - (thres << 2) >= 0)
			{
				//set if space cut yields 5 pieces
				//decision |= 1 << i + 1 + 2 * N_RANK ;
				num_subzoids [i] = initial_space_cut ? 4 : 5 ;
			}
			total_num_subzoids *= num_subzoids [i] ;
		}
        sim_can_cut |= space_cut ;
		unsigned long dim_key = (unsigned long) lb << num_bits_width | tb ;
		key = key << offset | dim_key ;
		offset += num_bits_dim ;
        call_boundary |= l_touch_boundary;
    }
	unsigned long index ;
	bool projection_exists = check_and_create_projection (key, lt, 
									centroid, index, l_father_grid) ;
	zoid_type & z = m_zoids [index] ;
	zoid_type & parent = m_zoids [parent_index] ;
	//add the zoid as a child of the parent
	parent.add_child(&z, child_index, index) ;
	assert (index < m_zoids.size()) ;
	if (projection_exists)
	{ 
		//Add the time of the zoid to the parent.
		clock_gettime(CLOCK_MONOTONIC, &end) ;
		redundant_time += tdiff2(&end, &start) ;
		projected_time += z.time  ;
		//a zoid with the projection already exists. return
		return ;
	}
	//print_dag() ;

    if (call_boundary)
	{
		//z.decision |= (unsigned short) 1 << 
		z.decision |= (unsigned char) 1 << 
					  zoid<N_RANK>::NUM_BITS_DECISION - 1 ;
        l_dt_stop = dt_recursive_boundary_;
	}
    else
	{
        l_dt_stop = dt_recursive_;
	}
	double projected_time1 = 0, necessary_time = 0 ;
	bool divide_and_conquer = false ;
    if (sim_can_cut) 
	{
		//assert (decision) ;
		divide_and_conquer = true ;
		z.resize_children(total_num_subzoids) ;
		double elapsed_time = 0, rtime = 0, ptime = 0 ;
		//cout << "space cut " << endl ;
        //cut into space 
		clock_gettime(CLOCK_MONOTONIC, &start1) ;
    	/*for (int i = N_RANK-1; i >= 0; --i) {
        	touch_boundary(i, lt, l_father_grid) ;
    	}*/
        if (call_boundary) 
		{
            symbolic_sawzoid_space_cut_boundary(t0, t1, l_father_grid, index, 
					 f, bf, num_subzoids, rtime, ptime);
        }
		else
		{
            symbolic_sawzoid_space_cut_interior(t0, t1, l_father_grid, index, 
						f, num_subzoids, rtime, ptime);
		}
		clock_gettime(CLOCK_MONOTONIC, &end1) ;
		elapsed_time = tdiff2(&end1, &start1) - rtime ;
		assert (elapsed_time >= 0.) ;
		assert (ptime >= 0.) ;

		necessary_time = elapsed_time ;
		projected_time1 = ptime ;
		decision = 2 ;
		m_zoids [index].decision |= decision ;
#ifndef NDEBUG
		m_zoids [index].stime = elapsed_time + ptime ;
#endif
    } 
	else if (lt > l_dt_stop)  //time cut
	{
		divide_and_conquer = true ;
		double elapsed_time = 0, rtime = 0, ptime = 0 ;
		z.resize_children(2) ;
        // cut into time 
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
		clock_gettime(CLOCK_MONOTONIC, &start1) ;
    	/*for (int i = N_RANK-1; i >= 0; --i) {
        	touch_boundary(i, lt, l_father_grid) ;
    	}*/
        if (call_boundary) {
            symbolic_sawzoid_space_time_cut_boundary(t0, t0+halflt, l_son_grid, 							index, 0, rtime, ptime , f, bf);
        } else {
            symbolic_sawzoid_space_time_cut_interior(t0, t0+halflt, l_son_grid,
								index, 0, rtime, ptime, f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            symbolic_sawzoid_space_time_cut_boundary(t0+halflt, t1, l_son_grid,
							index, 1, rtime, ptime, f, bf);
        } else {
            symbolic_sawzoid_space_time_cut_interior(t0+halflt, t1, l_son_grid,
							index, 1, rtime, ptime, f);
        }
		clock_gettime(CLOCK_MONOTONIC, &end1) ;
		elapsed_time = tdiff2(&end1, &start1) - rtime ;
		assert (elapsed_time >= 0.) ;
		assert (ptime >= 0.) ;

		//if space cut didn't happen or 
		//   space cut happened and time for time cut < time for space cut
		//then reset the necessary and projected times.
		if ((sim_can_cut && 
			elapsed_time + ptime < necessary_time + projected_time1)
			|| !sim_can_cut)
		{
			necessary_time = elapsed_time ;
			projected_time1 = ptime ;
			decision = 1 ;
			m_zoids [index].decision |= decision ;
		}
#ifndef NDEBUG
		m_zoids [index].ttime = elapsed_time + ptime ;
#endif
    } 
	else
	{
		//base case
#ifndef NDEBUG
		if (call_boundary)
		{
			assert (m_zoids [index].decision == (unsigned char) 1 <<
                      zoid<N_RANK>::NUM_BITS_DECISION - 1) ;
		}
		else
		{
			assert (m_zoids [index].decision == 0) ;
		}
#endif
		z.resize_children(1) ;
		//add the grid as a child of the leaf
		m_zoids [index].add_child((zoid_type *) 0, 0, m_leaves.size()) ;
		m_leaves.push_back(l_father_grid) ;
		
	}
	clock_gettime(CLOCK_MONOTONIC, &end) ;
	double total_time = tdiff2(&end, &start) ;
	redundant_time += total_time - necessary_time ;
}

// sawzoid space cuts. 
//template <int N_RANK> template <typename F>
//inline void auto_tune<N_RANK>::sawzoid_space_cut_interior(
//				int t0, int t1, grid_info<N_RANK> const & grid, 
//				simple_zoid_type * projection_zoid, F const & f)
//{
//    queue_info *l_father;
//    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
//    int queue_head_[2], queue_tail_[2], queue_len_[2];
//	const int lt = t1 - t0 ;
//
//    for (int i = 0; i < 2; ++i) {
//        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
//    }
//
//	int child_index = 0 ;
//	assert (projection_zoid) ;
//    // set up the initial grid 
//    push_queue(0, N_RANK-1, t0, t1, grid);
//    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
//        const int curr_dep_pointer = (curr_dep & 0x1);
//        while (queue_len_[curr_dep_pointer] > 0) {
//            top_queue(curr_dep_pointer, l_father);
//            if (l_father->level < 0) {
//                // spawn all the grids in circular_queue_[curr_dep][] 
//#if USE_CILK_FOR 
//                // use cilk_for to spawn all the sub-grid 
//// #pragma cilk_grainsize = 1
//                cilk_for (int j = 0; j < queue_len_[curr_dep_pointer]; ++j) {
//                    int i = pmod((queue_head_[curr_dep_pointer]+j), ALGOR_QUEUE_SIZE);
//                    queue_info * l_son = &(circular_queue_[curr_dep_pointer][i]);
//                    // assert all the sub-grid has done N_RANK spatial cuts 
//                    assert(l_son->level == -1);
//					//assert(child_index < projection_zoid->num_children) ;
//					unsigned long index = projection_zoid->children [child_index] ;
//					simple_zoid_type * child = &(m_simple_zoids [index]) ;
//					//zoid_type * child = &(m_zoids [index]) ;
//					child_index++ ; //looks like a race
//                    sawzoid_space_time_cut_interior(l_son->t0, 
//						l_son->t1, l_son->grid, child, f);
//                } // end cilk_for 
//                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
//                queue_len_[curr_dep_pointer] = 0;
//#else
//                // use cilk_spawn to spawn all the sub-grid 
//                pop_queue(curr_dep_pointer);
//				//assert(child_index < projection_zoid->num_children) ;
//				unsigned long index = projection_zoid->children [child_index] ;
//				simple_zoid_type * child = &(m_simple_zoids [index]) ;
//				//zoid_type * child = &(m_zoids [index]) ;
//				child_index++ ;
//                if (queue_len_[curr_dep_pointer] == 0)
//				{
//                    sawzoid_space_time_cut_interior(l_father->t0,
// 							l_father->t1, l_father->grid, child, f);
//				}
//                else
//				{
//                    cilk_spawn sawzoid_space_time_cut_interior(
//						l_father->t0, l_father->t1, l_father->grid, child, f) ;
//				}
//#endif
//            } else {
//                // performing a space cut on dimension 'level' 
//                pop_queue(curr_dep_pointer);
//                const grid_info<N_RANK> l_father_grid = l_father->grid;
//                //const int t0 = l_father->t0, t1 = l_father->t1;
//                //const int lt = (t1 - t0);
//                const int level = l_father->level;
//                const int thres = slope_[level] * lt;
//                //const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
//                //const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
//                //const bool cut_lb = (lb < tb);
//				//const bool can_cut = CAN_CUT_I ;
//				//if (! can_cut) {
//				if ((projection_zoid->decision & 1 << level + 1) == 0) {
//                    // if we can't cut into this dimension, just directly push 
//                    // it into the circular queue 
//                    //
//                    push_queue(curr_dep_pointer, level-1, t0, t1, 
//							l_father_grid) ;
//                } else  {
//					assert ((projection_zoid->decision & 1 << level + 1) != 0) ;
//                    /* can_cut! */
//                    //if (cut_lb) {
//					if (projection_zoid->decision & 1 << level + 1 + N_RANK) {
//                        grid_info<N_RANK> l_son_grid = l_father_grid;
//                        const int l_start = (l_father_grid.x0[level]);
//                        const int l_end = (l_father_grid.x1[level]);
//
//                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
//                        l_son_grid.x0[level] = l_start ;
//                        l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
//                        l_son_grid.x1[level] = l_start ;
//                        l_son_grid.dx1[level] = slope_[level] ;
//                        push_queue(next_dep_pointer, level-1, t0, t1, 
//									l_son_grid) ;
//
//                        l_son_grid.x0[level] = l_end ;
//                        l_son_grid.dx0[level] = -slope_[level];
//                        l_son_grid.x1[level] = l_end ;
//                        l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
//                        push_queue(next_dep_pointer, level-1, t0, t1, 
//									l_son_grid) ;
//						//const int cut_more = ((lb - (thres << 2)) >= 0) ;
//						//if (cut_more)
//						if (projection_zoid->decision & 
//									1 << level + 1 + 2 * N_RANK)
//						{
//							const int offset = (thres << 1) ;
//							l_son_grid.x0[level] = l_start ;
//	                        l_son_grid.dx0[level] = slope_[level];
//    	                    l_son_grid.x1[level] = l_start + offset ;
//        	                l_son_grid.dx1[level] = -slope_[level] ;
//            	            push_queue(curr_dep_pointer, level-1, t0, 
//								t1, l_son_grid) ;
//							
//							l_son_grid.x0[level] = l_end - offset ;
//	                        l_son_grid.dx0[level] = slope_[level];
//    	                    l_son_grid.x1[level] = l_end ;
//        	                l_son_grid.dx1[level] = -slope_[level] ;
//            	            push_queue(curr_dep_pointer, level-1, t0, 
//								t1, l_son_grid) ;
//
//							l_son_grid.x0[level] = l_start + offset ;
//	                        l_son_grid.dx0[level] = -slope_[level] ;
//    	                    l_son_grid.x1[level] = l_end - offset ;
//        	                l_son_grid.dx1[level] = slope_[level] ;
//            	            push_queue(next_dep_pointer, level-1, t0, 
//								t1, l_son_grid) ;
//						}
//						else
//						{
//							l_son_grid.x0[level] = l_start ;
//	                        l_son_grid.dx0[level] = slope_[level];
//    	                    l_son_grid.x1[level] = l_end;
//        	                l_son_grid.dx1[level] = -slope_[level] ;
//            	            push_queue(curr_dep_pointer, level-1, t0, 
//								t1, l_son_grid) ;
//						}
//                    } /* end if (cut_lb) */
//                    else {
//                        /* cut_tb */
//                        grid_info<N_RANK> l_son_grid = l_father_grid;
//                        const int l_start = (l_father_grid.x0[level]);
//                        const int l_end = (l_father_grid.x1[level]);
//						const int offset = (thres << 1) ;
//
//                        l_son_grid.x0[level] = l_start;
//                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
//                        l_son_grid.x1[level] = l_start + offset;
//                        l_son_grid.dx1[level] = -slope_[level];
//                        push_queue(curr_dep_pointer, level-1, t0, t1, 
//									l_son_grid) ;
//
//                        l_son_grid.x0[level] = l_end - offset ;
//                        l_son_grid.dx0[level] = slope_[level];
//                        l_son_grid.x1[level] = l_end;
//                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
//                        push_queue(curr_dep_pointer, level-1, t0, t1, 
//									l_son_grid);
//
//                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
//						//const int cut_more = ((tb - (thres << 2)) >= 0) ;
//						//if (cut_more)
//						if (projection_zoid->decision & 
//									1 << level + 1 + 2 * N_RANK)
//						{
//							l_son_grid.x0[level] = l_start + offset ;
//                            l_son_grid.dx0[level] = -slope_[level];
//							l_son_grid.x1[level] = l_start + offset;
//                        	l_son_grid.dx1[level] = slope_[level];
//                        	push_queue(next_dep_pointer, level-1, t0, 
//								t1, l_son_grid) ;
//
//							l_son_grid.x0[level] = l_end - offset ;
//                        	l_son_grid.dx0[level] = -slope_[level];
//							l_son_grid.x1[level] = l_end - offset ;
//                            l_son_grid.dx1[level] = slope_[level];
//                        	push_queue(next_dep_pointer, level-1, t0, 
//								t1, l_son_grid) ;
//
//							l_son_grid.x0[level] = l_start + offset ;
//                            l_son_grid.dx0[level] = slope_[level];
//                            l_son_grid.x1[level] = l_end - offset ;
//                            l_son_grid.dx1[level] = -slope_[level];
//                            push_queue(curr_dep_pointer, level-1, t0, 
//								t1, l_son_grid) ;
//						}
//						else
//						{
//                        	l_son_grid.x0[level] = l_start + offset ;
//                        	l_son_grid.dx0[level] = -slope_[level];
//                        	l_son_grid.x1[level] = l_end - offset ;
//                        	l_son_grid.dx1[level] = slope_[level];
//                        	push_queue(next_dep_pointer, level-1, t0, 
//								t1, l_son_grid) ;
//						}
//                    } /* end else (cut_tb) */
//                } // end if (can_cut) 
//            } // end if (performing a space cut) 
//        } // end while (queue_len_[curr_dep] > 0) 
//#if !USE_CILK_FOR
//        cilk_sync;
//#endif
//        assert(queue_len_[curr_dep_pointer] == 0);
//    } // end for (curr_dep < N_RANK+1) 
//}
//
//template <int N_RANK> template <typename F, typename BF>
//inline void auto_tune<N_RANK>::sawzoid_space_cut_boundary(
//	int t0, int t1,	grid_info<N_RANK> const & grid, 
//	simple_zoid_type * projection_zoid, F const & f, BF const & bf)
//{
//    queue_info *l_father;
//    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
//    int queue_head_[2], queue_tail_[2], queue_len_[2];
//    const int lt = t1 - t0 ;
//
//    for (int i = 0; i < 2; ++i) {
//        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
//    }
//	
//	int child_index = 0 ;
//	assert (projection_zoid) ;
//    // set up the initial grid 
//    push_queue(0, N_RANK-1, t0, t1, grid);
//    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
//        const int curr_dep_pointer = (curr_dep & 0x1);
//        while (queue_len_[curr_dep_pointer] > 0) {
//            top_queue(curr_dep_pointer, l_father);
//            if (l_father->level < 0) {
//                // spawn all the grids in circular_queue_[curr_dep][] 
//#if USE_CILK_FOR 
//                // use cilk_for to spawn all the sub-grid 
//// #pragma cilk_grainsize = 1
//				
//                cilk_for (int j = 0; j < queue_len_[curr_dep_pointer]; ++j) {
//                    int i = pmod((queue_head_[curr_dep_pointer]+j), ALGOR_QUEUE_SIZE);
//                    queue_info * l_son = &(circular_queue_[curr_dep_pointer][i]);
//                    // assert all the sub-grid has done N_RANK spatial cuts 
//                    //assert(l_son->level == -1);
//					//assert(child_index < projection_zoid->num_children) ;
//					unsigned long index = projection_zoid->children [child_index] ;
//					simple_zoid_type * child = &(m_simple_zoids [index]) ;
//					//zoid_type * child = &(m_zoids [index]) ;
//					child_index++ ;
//                    sawzoid_space_time_cut_boundary(l_son->t0, 
//						l_son->t1, l_son->grid, child, f, bf);
//                } // end cilk_for 
//                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
//                queue_len_[curr_dep_pointer] = 0;
//#else
//                // use cilk_spawn to spawn all the sub-grid 
//                pop_queue(curr_dep_pointer);
//				//assert(child_index < projection_zoid->num_children) ;
//				unsigned long index = projection_zoid->children [child_index] ;
//				simple_zoid_type * child = &(m_simple_zoids [index]) ;
//				//zoid_type * child = &(m_zoids [index]) ;
//				child_index++ ;
//                if (queue_len_[curr_dep_pointer] == 0) {
//                    sawzoid_space_time_cut_boundary(l_father->t0,
//						l_father->t1, l_father->grid, child, f, bf);
//                } else {
//                    cilk_spawn sawzoid_space_time_cut_boundary(
//						l_father->t0, l_father->t1, l_father->grid, child, f, bf);
//                }
//#endif
//            } else {
//                // performing a space cut on dimension 'level' 
//                pop_queue(curr_dep_pointer);
//                grid_info<N_RANK> l_father_grid = l_father->grid;
//                //const int t0 = l_father->t0, t1 = l_father->t1;
//                //const int lt = (t1 - t0);
//                const int level = l_father->level;
//                const int thres = slope_[level] * lt;
//                //const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
//                //const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
//                //const bool cut_lb = (lb < tb);
////#define touch_boundary m_algo.touch_boundary
//                //const bool l_touch_boundary = touch_boundary(level, lt, l_father_grid);
////#undef touch_boundary
//				//const bool can_cut = CAN_CUT_B ;
//                //if (! can_cut) {
//				if ((projection_zoid->decision & 1 << level + 1) == 0) {
//                //if (num_zoids [level] == 0) {
//                    // if we can't cut into this dimension, just directly push
//                    // it into the circular queue
//                    //
//                    push_queue(curr_dep_pointer, level-1, t0, t1, 
//							l_father_grid) ;
//                } else  {
//					assert ((projection_zoid->decision & 1 << level + 1) != 0) ;
//                    /* can_cut */
//                    //if (cut_lb) {
//					if (projection_zoid->decision & 1 << level + 1 + N_RANK) {
//                        grid_info<N_RANK> l_son_grid = l_father_grid;
//                        const int l_start = (l_father_grid.x0[level]);
//                        const int l_end = (l_father_grid.x1[level]);
//
//                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
//                        l_son_grid.x0[level] = l_start ;
//                        l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
//                        l_son_grid.x1[level] = l_start ;
//                        l_son_grid.dx1[level] = slope_[level] ;
//                        push_queue(next_dep_pointer, level-1, t0, t1, 
//											l_son_grid) ;
//
//                        l_son_grid.x0[level] = l_end ;
//                        l_son_grid.dx0[level] = -slope_[level];
//                        l_son_grid.x1[level] = l_end ;
//                        l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
//                        push_queue(next_dep_pointer, level-1, t0, t1, 
//											l_son_grid) ;
//
//						//const int cut_more = ((lb - (thres << 2)) >= 0) ;
//						//if (cut_more)
//						if (projection_zoid->decision & 
//									1 << level + 1 + 2 * N_RANK)
//						{
//							const int offset = (thres << 1) ;
//							l_son_grid.x0[level] = l_start ;
//	                        l_son_grid.dx0[level] = slope_[level];
//    	                    l_son_grid.x1[level] = l_start + offset;
//        	                l_son_grid.dx1[level] = -slope_[level] ;
//            	            push_queue(curr_dep_pointer, level-1, t0, 
//									t1, l_son_grid) ;
//
//							
//							l_son_grid.x0[level] = l_end - offset ;
//	                        l_son_grid.dx0[level] = slope_[level];
//    	                    l_son_grid.x1[level] = l_end ;
//        	                l_son_grid.dx1[level] = -slope_[level] ;
//            	            push_queue(curr_dep_pointer, level-1, t0, 
//									t1, l_son_grid) ;
//
//							l_son_grid.x0[level] = l_start + offset;
//	                        l_son_grid.dx0[level] = -slope_[level];
//    	                    l_son_grid.x1[level] = l_end - offset ;
//        	                l_son_grid.dx1[level] = slope_[level] ;
//            	            push_queue(next_dep_pointer, level-1, t0, 
//									t1, l_son_grid);
//						}
//						else
//						{
//							l_son_grid.x0[level] = l_start ;
//	                        l_son_grid.dx0[level] = slope_[level];
//    	                    l_son_grid.x1[level] = l_end;
//        	                l_son_grid.dx1[level] = -slope_[level] ;
//            	            push_queue(curr_dep_pointer, level-1, t0, 
//									t1, l_son_grid);
//						}
//                    } /* end if (cut_lb) */
//                    else { /* cut_tb */
//                        //if (lb == phys_length_[level] && l_father_grid.dx0[level] == 0 && l_father_grid.dx1[level] == 0) { /* initial cut on the dimension */
//						if (projection_zoid->decision & 
//									1 << level + 1 + 3 * N_RANK) {
//                            grid_info<N_RANK> l_son_grid = l_father_grid;
//                            const int l_start = (l_father_grid.x0[level]);
//                            const int l_end = (l_father_grid.x1[level]);
//                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
//							//const int cut_more = ((lb - (thres << 2)) >= 0) ;
//							//if (cut_more)
//							if (projection_zoid->decision & 
//									1 << level + 1 + 2 * N_RANK)
//							{
//								const int offset = (thres << 1) ;
//                            	l_son_grid.x0[level] = l_start ;
//                            	l_son_grid.dx0[level] = slope_[level];
//                            	l_son_grid.x1[level] = l_start + offset ;
//                            	l_son_grid.dx1[level] = -slope_[level];
//                            	push_queue(curr_dep_pointer, level-1, 
//									t0, t1, l_son_grid) ;
//
//                            	l_son_grid.x0[level] = l_end - offset ;
//                            	l_son_grid.dx0[level] = slope_[level];
//                            	l_son_grid.x1[level] = l_end ;
//                            	l_son_grid.dx1[level] = -slope_[level];
//                            	push_queue(curr_dep_pointer, level-1, 
//									t0, t1, l_son_grid) ;
//
//                            	l_son_grid.x0[level] = l_start + offset ;
//                            	l_son_grid.dx0[level] = -slope_[level];
//                            	l_son_grid.x1[level] = l_end - offset ;
//                            	l_son_grid.dx1[level] = slope_[level];
//                            	push_queue(next_dep_pointer, level-1, 
//									t0, t1, l_son_grid);
//							}
//							else
//							{
//                            	l_son_grid.x0[level] = l_start ;
//                            	l_son_grid.dx0[level] = slope_[level];
//                            	l_son_grid.x1[level] = l_end ;
//                            	l_son_grid.dx1[level] = -slope_[level];
//                            	push_queue(curr_dep_pointer, level-1, 
//									t0, t1, l_son_grid);
//							}
//                            l_son_grid.x0[level] = l_end ;
//                            l_son_grid.dx0[level] = -slope_[level];
//                            l_son_grid.x1[level] = l_end ;
//                            l_son_grid.dx1[level] = slope_[level];
//                            push_queue(next_dep_pointer, level-1, t0, 
//								t1, l_son_grid) ;
//						} else  {
//							/* cut_tb */
//							grid_info<N_RANK> l_son_grid = l_father_grid;
//							const int l_start = (l_father_grid.x0[level]);
//							const int l_end = (l_father_grid.x1[level]);
//							const int offset = (thres << 1) ;
//
//							l_son_grid.x0[level] = l_start;
//							l_son_grid.dx0[level] = l_father_grid.dx0[level];
//							l_son_grid.x1[level] = l_start + offset;
//							l_son_grid.dx1[level] = -slope_[level];
//							push_queue(curr_dep_pointer, level-1, t0, 
//								t1, l_son_grid) ;
//
//							l_son_grid.x0[level] = l_end - offset ;
//							l_son_grid.dx0[level] = slope_[level];
//							l_son_grid.x1[level] = l_end;
//							l_son_grid.dx1[level] = l_father_grid.dx1[level];
//							push_queue(curr_dep_pointer, level-1, t0, 
//								t1, l_son_grid) ;
//
//							const int next_dep_pointer = (curr_dep + 1) & 0x1;
//							//const int cut_more = ((tb - (thres << 2)) >= 0) ;
//							//if (cut_more)
//							if (projection_zoid->decision & 
//									1 << level + 1 + 2 * N_RANK)
//							{
//								l_son_grid.x0[level] = l_start + offset ;
//								l_son_grid.dx0[level] = -slope_[level];
//								l_son_grid.x1[level] = l_start + offset;
//								l_son_grid.dx1[level] = slope_[level];
//								push_queue(next_dep_pointer, level-1, 
//									t0, t1, l_son_grid) ;
//
//								l_son_grid.x0[level] = l_end - offset ;
//								l_son_grid.dx0[level] = -slope_[level];
//								l_son_grid.x1[level] = l_end - offset ;
//								l_son_grid.dx1[level] = slope_[level];
//								push_queue(next_dep_pointer, level-1, 
//									t0, t1, l_son_grid) ;
//
//								l_son_grid.x0[level] = l_start + offset ;
//								l_son_grid.dx0[level] = slope_[level];
//								l_son_grid.x1[level] = l_end - offset ;
//								l_son_grid.dx1[level] = -slope_[level];
//								push_queue(curr_dep_pointer, level-1, 
//									t0, t1, l_son_grid);
//							}
//							else
//							{
//								l_son_grid.x0[level] = l_start + offset ;
//								l_son_grid.dx0[level] = -slope_[level];
//								l_son_grid.x1[level] = l_end - offset ;
//								l_son_grid.dx1[level] = slope_[level];
//								push_queue(next_dep_pointer, level-1, 
//									t0, t1, l_son_grid);
//							}
//						}                   
//                    } /* end if (cut_tb) */
//                }// end if (can_cut) 
//            } // end if (performing a space cut) 
//        } // end while (queue_len_[curr_dep] > 0) 
//#if !USE_CILK_FOR
//        cilk_sync;
//#endif
//        assert(queue_len_[curr_dep_pointer] == 0);
//    } // end for (curr_dep < N_RANK+1) 
//}
//
//template <int N_RANK> template <typename F>
//inline void auto_tune<N_RANK>::
//sawzoid_space_time_cut_interior(int t0,
//			int t1, grid_info<N_RANK> const & grid, 
//			simple_zoid_type * pzoid, F const & f)
//{
//    const int lt = t1 - t0;
//    bool sim_can_cut = false;
//    grid_info<N_RANK> l_son_grid;
//
//	assert (pzoid) ;
//	
//    for (int i = N_RANK-1; i >= 0; --i) {
//        int lb, thres, tb;
//        lb = (grid.x1[i] - grid.x0[i]);
//        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
//
//#ifndef NDEBUG
//		simple_zoid_type * z = pzoid ;
//		grid_info <N_RANK> & grid2 = z->grid ;
//		int x0 = pmod(grid.x0 [i], phys_length_ [i]) ;
//		int x1 = pmod(grid.x1 [i], phys_length_ [i]) ;
//
//		//int x0_ = pmod(grid2.x0 [i], phys_length_ [i]) ;
//		//int x1_ = pmod(grid2.x1 [i], phys_length_ [i]) ;
//		int x0_ = grid2.x0 [i] ;
//		int x1_ = grid2.x1 [i] ;
//
//		int x2 = pmod(grid.x0[i] + grid.dx0[i] * lt, phys_length_ [i]) ;
//		int x3 = pmod(grid.x1[i] + grid.dx1[i] * lt, phys_length_ [i]) ;
//
//		//int x2_ = pmod(grid2.x0[i] + grid2.dx0[i] * lt, phys_length_ [i]) ;
//		//int x3_ = pmod(grid2.x1[i] + grid2.dx1[i] * lt, phys_length_ [i]) ;
//		int x2_ = grid2.x0[i] + grid2.dx0[i] * lt ;
//		int x3_ = grid2.x1[i] + grid2.dx1[i] * lt ;
//
//		if (x0 != x0_ || x1 != x1_ || x2 != x2_ || x3 != x3_)
//		{
//			cout << "zoid and proj zoid differ " << endl ;
//			//cout << "\nzoid " << z->id << " height " << z->height <<
//				//" num children " << z->children.size() <<
//			cout <<	" num children " << z->num_children << endl ;
//			//	" num_parents " << z->parents.size() << " geneity " ;
//			//print_bits(&(z->geneity), sizeof(word_type) * 8);
//			int h = lt ;
//			for (int i = N_RANK - 1 ; i >= 0 ; i--)
//			{
//				cout << " x0 [" << i << "] " << grid.x0 [i] 
//				 << " x1 [" << i << "] " << grid.x1 [i] 
//				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
//				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
//				<< " lt " << lt << endl ;
//				cout << " x0 [" << i << "] " << grid2.x0 [i]
//				 << " x1 [" << i << "] " << grid2.x1 [i]
//				<< " x2 [" << i << "] " << grid2.x0[i] + grid2.dx0[i] * h
//				<< " x3 [" << i << "] " << grid2.x1[i] + grid2.dx1[i] * h
//				<< " h " << h << endl ;
//			}
//			assert(0) ;
//		}
//#endif
//        bool cut_lb = (lb < tb);
//        thres = (slope_[i] * lt);
//        sim_can_cut = SIM_CAN_CUT_I ;
//    }
//    //if (pzoid->decision & 2) 
//    if (pzoid->decision & m_space_cut_mask) 
//	{
//        /* cut into space */
//        sawzoid_space_cut_interior(t0, t1, grid, pzoid, f) ;
//    } 
//	else if (pzoid->decision & 1)
//	{
//        /* cut into time */
//		assert (pzoid->children [0]) ;
//		assert (pzoid->children [1]) ;
//        assert(lt > dt_recursive_);
//        int halflt = lt / 2;
//        l_son_grid = grid;
//		unsigned long index = pzoid->children [0] ;
//        sawzoid_space_time_cut_interior(t0, t0+halflt, 
//									l_son_grid, &(m_simple_zoids [index]), f);
//
//        for (int i = 0; i < N_RANK; ++i) {
//            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
//            l_son_grid.dx0[i] = grid.dx0[i];
//            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
//            l_son_grid.dx1[i] = grid.dx1[i];
//        }
//		index = pzoid->children [1] ;
//        sawzoid_space_time_cut_interior(t0+halflt, t1, 
//								l_son_grid, &(m_simple_zoids [index]), f);
//    }
//	else
//	{
//		assert (pzoid->decision == 0) ;
//		//loop
//		f(t0, t1, grid);
//		return ;
//	}
//}

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::sawzoid_space_cut_interior(
	int t0, int t1,	
	simple_zoid_type * pzoid, F const & f)
{
	assert (pzoid) ;
	int n = (int) pzoid->num_children ;
	unsigned int dependency = pzoid->dependency ;
	int child_index = 0 ;
#ifndef NDEBUG
	int level = 0 ;
#endif
	while (n > 0)
	{
		//cout << "dependency " << dependency << endl ;
		int num_children_level = dependency & zoid_type::DEPENDENCY_MASK;
		//assert (num_children_level > 0) ;
		assert (num_children_level <= pzoid->num_children) ;
		dependency >>= zoid_type::NUM_BITS_DEPENDENCY ;
		//iterate over children in each dependency level and spawn them
		for (int i = 0 ; i < num_children_level ; i++)
		{
			unsigned long index = pzoid->children [child_index] ;
			child_index++ ;
			cilk_spawn sawzoid_space_time_cut_interior(t0, t1,
						&(m_simple_zoids [index]), f);
		}
		if (num_children_level > 0)
		{
			cilk_sync ;
		}
		n -= num_children_level ;
#ifndef NDEBUG
		level++ ;
#endif
	}
	assert (level <= N_RANK + 1) ;
	/*for (int i = 0 ; i < N_RANK + 1 ; i++)
	{
		//cout << "dependency " << dependency << endl ;
		int num_children_level = dependency & zoid_type::DEPENDENCY_MASK;
		//assert (num_children_level > 0) ;
		assert (num_children_level <= pzoid->num_children) ;
		dependency >>= zoid_type::NUM_BITS_DEPENDENCY ;
		//iterate over children in each dependency level and spawn them
		for (int i = 0 ; i < num_children_level ; i++)
		{
			unsigned long index = pzoid->children [child_index] ;
			child_index++ ;
			cilk_spawn sawzoid_space_time_cut_interior(t0, t1,
						&(m_simple_zoids [index]), f);
		}
		if (num_children_level > 0)
		{
			cilk_sync ;
		}
	}*/
	assert (child_index == pzoid->num_children) ;
}

template <int N_RANK> template <typename F, typename BF>
inline void auto_tune<N_RANK>::sawzoid_space_cut_boundary(
	int t0, int t1,	
	simple_zoid_type * pzoid, F const & f, BF const & bf)
{
	assert (pzoid) ;
	int n = (int) pzoid->num_children ;
	unsigned int dependency = pzoid->dependency ;
	int child_index = 0 ;
#ifndef NDEBUG
	int level = 0 ;
#endif
	while (n > 0)
	{
		//cout << "dependency " << dependency << endl ;
		int num_children_level = dependency & zoid_type::DEPENDENCY_MASK;
		//assert (num_children_level > 0) ;
		assert (num_children_level <= pzoid->num_children) ;
		dependency >>= zoid_type::NUM_BITS_DEPENDENCY ;
		//iterate over children in each dependency level and spawn them
		for (int i = 0 ; i < num_children_level ; i++)
		{
			unsigned long index = pzoid->children [child_index] ;
			child_index++ ;
			cilk_spawn sawzoid_space_time_cut_boundary(t0, t1,
						&(m_simple_zoids [index]), f, bf);
		}
		if (num_children_level > 0)
		{
			cilk_sync ;
		}
		n -= num_children_level ;
#ifndef NDEBUG
		level++ ;
#endif
	}
	assert (level <= N_RANK + 1) ;

	/*for (int i = 0 ; i < N_RANK + 1 ; i++)
	{
		//cout << "dependency " << dependency << endl ;
		int num_children_level = dependency & zoid_type::DEPENDENCY_MASK;
		//assert (num_children_level > 0) ;
		assert (num_children_level <= pzoid->num_children) ;
		dependency >>= zoid_type::NUM_BITS_DEPENDENCY ;
		//iterate over children in each dependency level and spawn them
		for (int i = 0 ; i < num_children_level ; i++)
		{
			unsigned long index = pzoid->children [child_index] ;
			child_index++ ;
			cilk_spawn sawzoid_space_time_cut_boundary(t0, t1,
						&(m_simple_zoids [index]), f, bf);
		}
		if (num_children_level > 0)
		{
			cilk_sync ;
		}
	}*/
	assert (child_index == pzoid->num_children) ;
}

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::
sawzoid_space_time_cut_interior(int t0, int t1, 
			simple_zoid_type * pzoid, F const & f)
{
    const int lt = t1 - t0;

	assert (pzoid) ;
	
    if (pzoid->decision & m_space_cut_mask) 
	{
        /* cut into space */
        sawzoid_space_cut_interior(t0, t1, pzoid, f) ;
    } 
	else if (pzoid->decision & 1)
	{
        /* cut into time */
		assert (pzoid->children [0]) ;
		assert (pzoid->children [1]) ;
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
		unsigned long index = pzoid->children [0] ;
        sawzoid_space_time_cut_interior(t0, t0+halflt, 
									&(m_simple_zoids [index]), f);

		index = pzoid->children [1] ;
        sawzoid_space_time_cut_interior(t0+halflt, t1, 
								&(m_simple_zoids [index]), f);
    }
	else
	{
		assert (pzoid->decision == 0) ;
		//loop
		//f(t0, t1, pzoid->grid);
		f(t0, t1, m_leaves [pzoid->children [0]]);
		return ;
	}
}

template <int N_RANK> template <typename F, typename BF>
inline void auto_tune<N_RANK>::
sawzoid_space_time_cut_boundary(int t0, int t1,	
				simple_zoid_type * pzoid, F const & f, BF const & bf)
{
    const int lt = t1 - t0;
	assert (pzoid) ;

	unsigned short call_boundary = pzoid->decision >> 
					zoid<N_RANK>::NUM_BITS_DECISION - 1 ;
	//grid_info<N_RANK> l_father_grid = pzoid->grid ;

	if (pzoid->decision & m_space_cut_mask)
	{
		//cout << "space cut " << endl ;
		//cut into space 
		if (call_boundary) 
		{
			sawzoid_space_cut_boundary(t0, t1, pzoid, f, bf) ;
		}
		else
		{
			sawzoid_space_cut_interior(t0, t1, pzoid, f) ;
		}
	} 
	else if (pzoid->decision & 1)
	{
		//cout << "time cut " << endl ;
		// cut into time 
		assert (pzoid->children [0]) ;
		assert (pzoid->children [1]) ;
		
		int halflt = lt / 2;
		unsigned long index = pzoid->children [0] ;
		if (call_boundary) {
			sawzoid_space_time_cut_boundary(t0, t0+halflt, 
								&(m_simple_zoids [index]), f, bf);
		} else {
			sawzoid_space_time_cut_interior(t0, t0+halflt, 
								&(m_simple_zoids [index]), f);
		}
		index = pzoid->children [1] ;
		if (call_boundary) {
			sawzoid_space_time_cut_boundary(t0+halflt, t1, 
								&(m_simple_zoids [index]), f, bf);
		} else {
			sawzoid_space_time_cut_interior(t0+halflt, t1, 
								&(m_simple_zoids [index]), f);
		}
	}
	else
	{
		//loop
		if (call_boundary) {
			assert (pzoid->decision == 
					1 << zoid<N_RANK>::NUM_BITS_DECISION - 1) ;
            //base_case_kernel_boundary(t0, t1, pzoid->grid, bf);
            base_case_kernel_boundary(t0, t1, m_leaves [pzoid->children [0]], 
									  bf);
        } else { 
			assert (pzoid->decision == 0) ;
            //f(t0, t1, pzoid->grid);
			f(t0, t1, m_leaves [pzoid->children [0]]);
        }
	}
}


//template <int N_RANK> template <typename F, typename BF>
//inline void auto_tune<N_RANK>::
//sawzoid_space_time_cut_boundary(int t0,
//				int t1,	grid_info<N_RANK> const & grid, 
//				simple_zoid_type * pzoid, F const & f, BF const & bf)
//{
//    const int lt = t1 - t0;
//    bool sim_can_cut = false, call_boundary = false;
//    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
//
//	assert (pzoid) ;
//    for (int i = N_RANK-1; i >= 0; --i) {
//        int lb, thres, tb;
//        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
//        lb = (grid.x1[i] - grid.x0[i]);
//        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
//#ifndef NDEBUG
//		simple_zoid_type * z = pzoid ;
//		grid_info <N_RANK> & grid2 = z->grid ;
//		int x0 = pmod(grid.x0 [i], phys_length_ [i]) ;
//		int x1 = pmod(grid.x1 [i], phys_length_ [i]) ;
//
//		int x0_ = pmod(grid2.x0 [i], phys_length_ [i]) ;
//		int x1_ = pmod(grid2.x1 [i], phys_length_ [i]) ;
//
//		int x2 = pmod(grid.x0[i] + grid.dx0[i] * lt, phys_length_ [i]) ;
//		int x3 = pmod(grid.x1[i] + grid.dx1[i] * lt, phys_length_ [i]) ;
//
//		int x2_ = pmod(grid2.x0[i] + grid2.dx0[i] * lt, phys_length_ [i]) ;
//		int x3_ = pmod(grid2.x1[i] + grid2.dx1[i] * lt, phys_length_ [i]) ;
//
//		if (x0 != x0_ || x1 != x1_ || x2 != x2_ || x3 != x3_)
//		{
//			cout << "zoid and proj zoid differ " << endl ;
//			//cout << "\nzoid " << z->id << " height " << z->height <<
//				//" num children " << z->children.size() <<
//			cout <<	" num children " << z->num_children << endl ;
//				//" num_parents " << z->parents.size() ;
//				//" num_parents " << z->parents.size() << " geneity " ;
//			//print_bits(&(z->geneity), sizeof(word_type) * 8);
//			int h = lt ;
//			for (int i = N_RANK - 1 ; i >= 0 ; i--)
//			{
//				cout << " x0 [" << i << "] " << grid.x0 [i] 
//				 << " x1 [" << i << "] " << grid.x1 [i] 
//				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
//				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
//				<< " lt " << lt << endl ;
//				cout << " x0 [" << i << "] " << grid2.x0 [i]
//				 << " x1 [" << i << "] " << grid2.x1 [i]
//				<< " x2 [" << i << "] " << grid2.x0[i] + grid2.dx0[i] * h
//				<< " x3 [" << i << "] " << grid2.x1[i] + grid2.dx1[i] * h
//				<< " h " << h << endl ;
//			}
//			assert(0) ;
//		}
//#endif
//        thres = slope_[i] * lt ;
//        bool cut_lb = (lb < tb);
//
//        sim_can_cut = SIM_CAN_CUT_B ;
//        call_boundary |= l_touch_boundary;
//    }
//    /*for (int i = N_RANK-1; i >= 0; --i) {
//        touch_boundary(i, lt, l_father_grid) ;
//    }*/
//	//unsigned short call_boundary = pzoid->decision >> 
//	//				zoid<N_RANK>::NUM_BITS_DECISION - 1 ;
//
//	//if (pzoid->decision & 2)
//	if (pzoid->decision & m_space_cut_mask)
//	{
//		//cout << "space cut " << endl ;
//		//cut into space 
//		if (call_boundary) 
//		{
//			sawzoid_space_cut_boundary(t0, t1, l_father_grid, 
//										pzoid, f, bf) ;
//		}
//		else
//		{
//			sawzoid_space_cut_interior(t0, t1, l_father_grid, 
//										pzoid, f) ;
//		}
//	} 
//	else if (pzoid->decision & 1)
//	{
//		//cout << "time cut " << endl ;
//		// cut into time 
//		assert (pzoid->children [0]) ;
//		assert (pzoid->children [1]) ;
//		
//		int halflt = lt / 2;
//		l_son_grid = l_father_grid;
//		unsigned long index = pzoid->children [0] ;
//		if (call_boundary) {
//			sawzoid_space_time_cut_boundary(t0, t0+halflt, 
//								l_son_grid, &(m_simple_zoids [index]), f, bf);
//		} else {
//			sawzoid_space_time_cut_interior(t0, t0+halflt, 
//								l_son_grid, &(m_simple_zoids [index]), f);
//		}
//
//		for (int i = 0; i < N_RANK; ++i) {
//			l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
//			l_son_grid.dx0[i] = l_father_grid.dx0[i];
//			l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
//			l_son_grid.dx1[i] = l_father_grid.dx1[i];
//		}
//		index = pzoid->children [1] ;
//		if (call_boundary) {
//			sawzoid_space_time_cut_boundary(t0+halflt, t1, 
//								l_son_grid, &(m_simple_zoids [index]), f, bf);
//		} else {
//			sawzoid_space_time_cut_interior(t0+halflt, t1, 
//								l_son_grid, &(m_simple_zoids [index]), f);
//		}
//	}
//	else
//	{
//		//loop
//		if (call_boundary) {
//			assert (pzoid->decision == 
//					1 << zoid<N_RANK>::NUM_BITS_DECISION - 1) ;
//            base_case_kernel_boundary(t0, t1, l_father_grid, bf);
//        } else { 
//			assert (pzoid->decision == 0) ;
//            f(t0, t1, l_father_grid);
//        }
//	}
//}


#undef dx_recursive_boundary_  
#undef dx_recursive_ 
#undef dt_recursive_boundary_ 
#undef dt_recursive_ 
#undef slope_ 
#undef touch_boundary 
#undef base_case_kernel_boundary

#endif
