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
#ifndef KERNEL_SELECTION_SAWZOID_HPP
#define KERNEL_SELECTION_SAWZOID_HPP 

#include "kernel_selection_sawzoid_header.hpp"

#define dx_recursive_boundary_  (m_algo.dx_recursive_boundary_)
#define dx_recursive_ (m_algo.dx_recursive_)
#define dt_recursive_boundary_ (m_algo.dt_recursive_boundary_)
#define dt_recursive_ (m_algo.dt_recursive_)
#define slope_ m_algo.slope_
#define touch_boundary m_algo.touch_boundary
#define phys_length_ m_algo.phys_length_
#define base_case_kernel_boundary m_algo.base_case_kernel_boundary

template <int N_RANK> template <typename F>
inline void kernel_selection_sawzoid<N_RANK>::symbolic_modified_space_cut_interior(int t0,
			int t1, grid_info<N_RANK> const & grid, F const & f)
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
                    symbolic_modified_space_time_cut_interior(l_son->t0, 
						l_son->t1, l_son->grid, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				symbolic_modified_space_time_cut_interior(l_father->t0, 
					l_father->t1, l_father->grid, f);
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
						const bool cut_more = CHECK_WIDTH_TB ;
						//const bool cut_more = (tb >= thres << 1)  && 
						//					  (lb > dx_recursive_[level]) ;
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
//        cilk_sync;
#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}

/* Boundary space cut. Uses modified space cut.
 */
template <int N_RANK> template <typename F>
inline void kernel_selection_sawzoid<N_RANK>::symbolic_modified_space_cut_boundary(int t0,
		int t1, grid_info<N_RANK> const & grid, F const & f)
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
                    symbolic_modified_space_time_cut_boundary(l_son->t0, 
						l_son->t1, l_son->grid, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				symbolic_modified_space_time_cut_boundary(l_father->t0, 
					l_father->t1, l_father->grid, f);
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
                //if (num_zoids [level] == 0) {
                if (! can_cut) {
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
						//			(lb > dx_recursive_[level]) ;
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
//#if !USE_CILK_FOR
//        cilk_sync;
//#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}

template <int N_RANK> template <typename F>
inline void kernel_selection_sawzoid<N_RANK>::symbolic_modified_space_time_cut_interior(
		int t0, int t1, grid_info<N_RANK> const & grid, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;
	int p_lb [N_RANK], p_tb [N_RANK] ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        lb = p_lb [i] = (grid.x1[i] - grid.x0[i]);
        tb = p_tb [i] = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);

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
		symbolic_modified_space_cut_interior(t0, t1, grid, f); 
    } 
	else if (lt > dt_recursive_) 
	{
        /* cut into time */
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
        symbolic_modified_space_time_cut_interior(t0, t0+halflt, l_son_grid, f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        symbolic_modified_space_time_cut_interior(t0+halflt, t1, l_son_grid, f);
    } 
	else 
	{
		// base case
		//determine the geneity of the leaf
		word_type geneity = 0 ;
		compute_geneity(lt, grid, geneity, f) ;
		
		//int key = (p_lb [0] >= p_tb) ; 
		int kernel = -1 ;
		if (geneity == 0)
		{
			kernel = 0 ;
		}
		else if (__builtin_popcount(geneity) == 1)
		{
			//zoid is homogeneous
			kernel = __builtin_ffs(geneity) ;
			//cout << "zoid is homogeneous" << endl ;
		}
		else
		{
			kernel = m_clone_array->clones.size() - 1 ;
		}
		
		assert (kernel != -1) ;
		assert (kernel < m_clone_array->clones.size()) ;
		
		//print_bits(&(z->geneity), sizeof(word_type) * 8) ;
		int index = 0, width = 1; 
		//if (lt == 1)
		if (lt < dt_recursive_)
		{
			int offset = 0 ;
			unsigned long key = 0 ; 
	    	for (int i = N_RANK-1; i >= 0; --i) 
			{
				//index = (grid.x0[i] + (p_lb [i] >> 1)) * 
				//			width + index ;
				index = pmod(grid.x0[i] + (p_lb [i] >> 1), phys_length_ [i]) * 
							width + index ;
				//index = pmod(grid.x0[i], phys_length_ [i]) * width + 
				//			index ;
				width = phys_length_ [i] ;
				//width = m_extended_length [i] ;
				unsigned long dim_key = (unsigned long) p_lb [i] << 
										num_bits_width | p_tb [i] ;
				key = key << offset | dim_key ;
				//key = key << offset | p_lb [i] ;
				offset += num_bits_dim ;
			}
			leaf_kernel_table & table = 
							m_arr_leaf_kernel_table [index] ;
			std::pair<lk_table_iterator, lk_table_iterator> p = 
											table.equal_range (key) ;
			if (p.first == p.second)
			{
				//key doesn't exist.
#ifndef NDEBUG
				zoid_type z ;
				z.kernel = kernel ;
				z.info = grid ;
				z.height = lt ;
				z.t0 = t0 , z.t1 = t1 ;
				//insert (key,zoid) pair
				table.insert(std::pair<unsigned long, zoid_type>(key, z)) ;
#else
				//insert (key,kernel) pair
				table.insert(std::pair<unsigned long, char>(key, (char) kernel)) ;
#endif
			}
		}
		else
		{
			int key ;
	    	for (int i = N_RANK-1; i >= 0; --i) 
			{
				//index = (grid.x0[i] + (p_lb [i] >> 1)) * 
				//			width + index ;
				index = pmod(grid.x0[i] + (p_lb [i] >> 1), phys_length_ [i]) * 
							width + index ;
				//index = pmod(grid.x0[i], phys_length_ [i]) * width + 
				//			index ;
				width = phys_length_ [i] ;
				//width = m_extended_length [i] ;
				//key = key << 1 | p_lb [i] >= p_tb [i] ;
			}
			key = p_lb [0] >= p_tb [0] ;
			//assert (key >= 0 && key < 1 << N_RANK) ;
#ifndef NDEBUG
			//zoid_type &  z = m_kernel_map [(index << N_RANK) + key] ;
			zoid_type &  z = m_kernel_map [(index << 1) + key] ;
			if (z.kernel != (char) 0  && 
				kernel != 0 && z.kernel != (char) kernel)
			{
				cout << "z.kernel  "  << (int) z.kernel << " kernel " << 
						kernel << endl ; 
				cout << " t0 " << t0 << " t1 " << t1 << endl ;
				grid_info <N_RANK> g2 = z.info ;
				int h = z.height ;
				for (int i = N_RANK - 1 ; i >= 0 ; i--)
				{
					cout << " x0 [" << i << "] " << grid.x0 [i] 
					<< " x1 [" << i << "] " << grid.x1 [i] 
					<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
					<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
					<< " lt " << lt << endl ;
				}
				cout << " z.t0 " << z.t0 << " z.t1 " << z.t1 << endl ;
				for (int i = N_RANK - 1 ; i >= 0 ; i--)
				{
					cout << " x0 [" << i << "] " << g2.x0 [i]
					<< " x1 [" << i << "] " << g2.x1 [i]
					<< " x2 [" << i << "] " << g2.x0[i] + g2.dx0[i] * h
					<< " x3 [" << i << "] " << g2.x1[i] + g2.dx1[i] * h
					<< " h " << h << endl ;
				}
				assert (z.kernel == (char) kernel) ;
			}
			else
			{
				z.info = grid ;
				z.height = lt ;
				z.kernel = (char) kernel ;
				z.t0 = t0 , z.t1 = t1 ;
			}
#else
			//m_kernel_map [(index << N_RANK) + key] = (char) kernel ;
			m_kernel_map [(index << 1) + key] = (char) kernel ;
			//char * c = & (m_kernel_map [index << 1]) ;
			//c [key] = (char) kernel ;
#endif
		}
	}
}

template <int N_RANK> template <typename F>
inline void kernel_selection_sawzoid<N_RANK>::symbolic_modified_space_time_cut_boundary(
		int t0, int t1,	grid_info<N_RANK> const & grid, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;

	int p_lb [N_RANK], p_tb [N_RANK] ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);

        lb = p_lb [i] = (grid.x1[i] - grid.x0[i]);
        tb = p_tb [i] = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
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
            symbolic_modified_space_cut_boundary(t0, t1, l_father_grid, f);
        }
		else
		{
            symbolic_modified_space_cut_interior(t0, t1, l_father_grid, f);
		}
    } 
	else if (lt > l_dt_stop)  //time cut
	{
		//cout << "time cut " << endl ;
        // cut into time 
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
        if (call_boundary) {
            symbolic_modified_space_time_cut_boundary(t0, t0+halflt, l_son_grid, 													f);
        } else {
            symbolic_modified_space_time_cut_interior(t0, t0+halflt, l_son_grid,
													f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            symbolic_modified_space_time_cut_boundary(t0+halflt, t1, l_son_grid,
													f);
        } else {
            symbolic_modified_space_time_cut_interior(t0+halflt, t1, l_son_grid,
													f);
        }
    } 
	else
	{
		// base case
		//determine the geneity of the leaf
		word_type geneity = 0 ;
		compute_geneity(lt, l_father_grid, geneity, f) ;
		//print_bits(&(z->geneity), sizeof(word_type) * 8) ;
		
		int kernel = -1 ;
		if (geneity == 0)
		{
			kernel = 0 ;
		}
		else if (__builtin_popcount(geneity) == 1)
		{
			//zoid is homogeneous
			kernel = __builtin_ffs(geneity) ;
			//cout << "zoid is homogeneous" << endl ;
		}
		else
		{
			//kernel_map [index][centroid] = m_clone_array->clones.size() - 1 ;
			kernel = m_clone_array->clones.size() - 1 ;
		}
		
		assert (kernel != -1) ;
		assert (kernel < m_clone_array->clones.size()) ;
		
		int index = 0, width = 1, key = 0 ; 
	    //if (lt == 1)
		if (lt < dt_recursive_)
		{
			int offset = 0 ;
			unsigned long key = 0 ; 
	    	for (int i = N_RANK-1; i >= 0; --i) 
			{
				//index = (grid.x0[i] + (p_lb [i] >> 1)) * 
				//			width + index ;
				index = pmod(grid.x0[i] + (p_lb [i] >> 1), phys_length_ [i]) * 
							width + index ;
				//index = pmod(grid.x0[i], phys_length_ [i]) * width + 
				//			index ;
				width = phys_length_ [i] ;
				//width = m_extended_length [i] ;
				unsigned long dim_key = (unsigned long) p_lb [i] << 
										num_bits_width | p_tb [i] ;
				key = key << offset | dim_key ;
				//key = key << offset | p_lb [i] ;
				offset += num_bits_dim ;
			}
			leaf_kernel_table & table = 
							m_arr_leaf_kernel_table [index] ;
			std::pair<lk_table_iterator, lk_table_iterator> p = 
											table.equal_range (key) ;
			if (p.first == p.second)
			{
				//key doesn't exist.
#ifndef NDEBUG
				zoid_type z ;
				z.kernel = kernel ;
				z.info = grid ;
				z.height = lt ;
				z.t0 = t0 , z.t1 = t1 ;
				//insert (key,zoid) pair
				table.insert(std::pair<unsigned long, zoid_type>(key, z)) ;
#else
				//insert (key,kernel) pair
				table.insert(std::pair<unsigned long, char>(key, (char) kernel)) ;
#endif
			}
		}
		else
		{
			int key = 0 ; 
	    	for (int i = N_RANK-1; i >= 0; --i) 
			{
				//index = (grid.x0[i] + (p_lb [i] >> 1)) * 
							width + index ;
				index = pmod(grid.x0[i] + (p_lb [i] >> 1), phys_length_ [i]) * 
							width + index ;
				//index = pmod(grid.x0[i], phys_length_ [i]) * width + 
				//			index ;
				width = phys_length_ [i] ;
				//width = m_extended_length [i] ;
				//key = key << 1 | p_lb [i] >= p_tb [i] ;
			}
			key = p_lb [0] >= p_tb [0] ;
			//assert (key >= 0 && key < 1 << N_RANK) ;
#ifndef NDEBUG
			//zoid_type &  z = m_kernel_map [(index << N_RANK) + key] ;
			zoid_type &  z = m_kernel_map [(index << 1) + key] ;
			if (z.kernel != (char) 0 && 
				kernel != 0 && z.kernel != (char) kernel)
			{
				cout << "z.kernel  "  << (int) z.kernel << " kernel " << 
						kernel << endl ; 
				cout << "index " << index << " key  " << key << endl ;
				cout << " t0 " << t0 << " t1 " << t1 << endl ;
				grid_info <N_RANK> g2 = z.info ;
				int h = z.height ;
				index = 0, width = 1, key = 0 ; 
				for (int i = N_RANK - 1 ; i >= 0 ; i--)
				{
					cout << " x0 [" << i << "] " << grid.x0 [i] 
					<< " x1 [" << i << "] " << grid.x1 [i] 
					<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
					<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
					<< " lt " << lt << endl ;
				}
				cout << " z.t0 " << z.t0 << " z.t1 " << z.t1 << endl ;
				for (int i = N_RANK - 1 ; i >= 0 ; i--)
				{
					cout << " x0 [" << i << "] " << g2.x0 [i]
					<< " x1 [" << i << "] " << g2.x1 [i]
					<< " x2 [" << i << "] " << g2.x0[i] + g2.dx0[i] * h
					<< " x3 [" << i << "] " << g2.x1[i] + g2.dx1[i] * h
					<< " h " << h << endl ;
					int lb = (g2.x1[i] - g2.x0[i]);
					int tb = (g2.x1[i] + g2.dx1[i] * h - g2.x0[i] - g2.dx0[i] * h);
					//index = (g2.x0[i] + (lb >> 1)) * 
					//		width + index ;
					index = pmod(g2.x0[i] + (lb >> 1), phys_length_ [i]) * 
							width + index ;
					width = phys_length_ [i] ;
					//width = m_extended_length [i] ;
					key = key << 1 | lb >= tb ;
				}
				cout << "index " << index << " key  " << key << endl ;
				assert (z.kernel == (char) kernel) ;
			}
			else
			{
				z.info = grid ;
				z.height = lt ;
				z.kernel = (char) kernel ;
				z.t0 = t0 , z.t1 = t1 ;
			}
#else
			//m_kernel_map [(index << N_RANK) + key] = (char) kernel ;
			m_kernel_map [(index << 1) + key] = (char) kernel ;
			//char * c = & (m_kernel_map [index << 1]) ;
			//c [key] = (char) kernel ;
#endif
		}
	}
}

// modified space cuts. 
template <int N_RANK> template <typename F>
inline void kernel_selection_sawzoid<N_RANK>::heterogeneous_modified_space_cut_interior(
				int t0, int t1, grid_info<N_RANK> const & grid, F const & f)
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
                    heterogeneous_modified_space_time_cut_interior(l_son->t0, 
						l_son->t1, l_son->grid, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0)
				{
                    heterogeneous_modified_space_time_cut_interior(l_father->t0,
 							l_father->t1, l_father->grid, f);
				}
                else
				{
                    cilk_spawn heterogeneous_modified_space_time_cut_interior(
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
        cilk_sync;
#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}


template <int N_RANK> template <typename F, typename BF>
inline void kernel_selection_sawzoid<N_RANK>::heterogeneous_modified_space_cut_boundary(
	int t0, int t1,	grid_info<N_RANK> const & grid, F const & f, BF const & bf)
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
                    heterogeneous_modified_space_time_cut_boundary(l_son->t0, 
						l_son->t1, l_son->grid, f, bf);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0) {
                    heterogeneous_modified_space_time_cut_boundary(l_father->t0,
						l_father->t1, l_father->grid, f, bf);
                } else {
                    cilk_spawn heterogeneous_modified_space_time_cut_boundary(
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

template <int N_RANK> template <typename F>
inline void kernel_selection_sawzoid<N_RANK>::
heterogeneous_modified_space_time_cut_interior(int t0,
				int t1, grid_info<N_RANK> const & grid, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;
	int p_lb [N_RANK], p_tb [N_RANK] ;
 
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        lb = p_lb [i] = (grid.x1[i] - grid.x0[i]);
        tb = p_tb [i] = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);

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
        heterogeneous_modified_space_cut_interior(t0, t1, grid, f) ;
    } 
	else if (lt > dt_recursive_) 
	{
        /* cut into time */
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
        heterogeneous_modified_space_time_cut_interior(t0, t0+halflt, 
									l_son_grid, f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        heterogeneous_modified_space_time_cut_interior(t0+halflt, t1, 
								l_son_grid, f);
    }
	else 
	{
		// base case
		int kernel = -1 ;
		//if (lt == 1)
		if (lt < dt_recursive_)
		{
			int index = 0, width = 1 ;
			int offset = 0 ;
			unsigned long key = 0 ; 
	    	for (int i = N_RANK-1; i >= 0; --i) 
			{
				//index = (grid.x0[i] + (p_lb [i] >> 1)) * 
				//			width + index ;
				index = pmod(grid.x0[i] + (p_lb [i] >> 1), phys_length_ [i]) * 
							width + index ;
				//index = pmod(grid.x0[i], phys_length_ [i]) * width + 
				//			index ;
				width = phys_length_ [i] ;
				//width = m_extended_length [i] ;
				unsigned long dim_key = (unsigned long) p_lb [i] << 
										num_bits_width | p_tb [i] ;
				key = key << offset | dim_key ;
				//key = key << offset | p_lb [i] ;
				offset += num_bits_dim ;
			}
			leaf_kernel_table & table = 
							m_arr_leaf_kernel_table [index] ;
			std::pair<lk_table_iterator, lk_table_iterator> p = 
											table.equal_range (key) ;
			assert (p.first != p.second) ;
			assert (p.first->first == key) ;
#ifndef NDEBUG
			kernel = p.first->second.kernel ;
#else
			kernel = p.first->second ;
#endif
		}
		else
		{
			int index = 0, width = 1, key = 0 ; 
			for (int i = N_RANK-1; i >= 0; --i) 
			{
				//index = (grid.x0[i] + (p_lb [i] >> 1)) * 
				//			width + index ;
				index = pmod(grid.x0[i] + (p_lb [i] >> 1), phys_length_ [i]) * 
							width + index ;
				//index = pmod(grid.x0[i], phys_length_ [i]) * width + 
				//			index ;
				width = phys_length_ [i] ;
				//width = m_extended_length [i] ;
				//key = key << 1 | p_lb [i] >= p_tb [i] ;
			}
			key = p_lb [0] >= p_tb [0] ;
			//assert (key >= 0 && key < 1 << N_RANK) ;
#ifndef NDEBUG
			//kernel = m_kernel_map [(index << N_RANK) + key].kernel ; 
			kernel = m_kernel_map [(index << 1) + key].kernel ; 
#else
			//kernel = m_kernel_map [(index << N_RANK) + key] ;
			kernel = m_kernel_map [(index << 1) + key] ;
#endif
		}
		//char * c = & (m_kernel_map [index << 1]) ;
		//kernel = (int) c [key] ;
		
		assert (kernel != -1) ;
		assert (kernel < m_clone_array->clones.size()) ;
		(*m_clone_array) [kernel] (t0, t1, grid);
	}
}

template <int N_RANK> template <typename F, typename BF>
inline void kernel_selection_sawzoid<N_RANK>::
heterogeneous_modified_space_time_cut_boundary(int t0,
					int t1,	grid_info<N_RANK> const & grid, 
					F const & f, BF const & bf)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;
	int p_lb [N_RANK], p_tb [N_RANK] ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);

        lb = p_lb [i] = (grid.x1[i] - grid.x0[i]);
        tb = p_tb [i] = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);

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
            heterogeneous_modified_space_cut_boundary(t0, t1, l_father_grid, 
										f, bf) ;
        }
		else
		{
            heterogeneous_modified_space_cut_interior(t0, t1, l_father_grid, 
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
            heterogeneous_modified_space_time_cut_boundary(t0, t0+halflt, 
								l_son_grid, f, bf);
        } else {
            heterogeneous_modified_space_time_cut_interior(t0, t0+halflt, 
								l_son_grid, f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            heterogeneous_modified_space_time_cut_boundary(t0+halflt, t1, 
								l_son_grid, f, bf);
        } else {
            heterogeneous_modified_space_time_cut_interior(t0+halflt, t1, 
								l_son_grid, f);
        }
    } 
	else
	{
		int kernel = -1 ;
		//if (lt == 1)
		if (lt < dt_recursive_)
		{
			int index = 0, width = 1 ;
			int offset = 0 ;
			unsigned long key = 0 ; 
	    	for (int i = N_RANK-1; i >= 0; --i) 
			{
				//index = (grid.x0[i] + (p_lb [i] >> 1)) * 
				//			width + index ;
				index = pmod(grid.x0[i] + (p_lb [i] >> 1), phys_length_ [i]) * 
							width + index ;
				//index = pmod(grid.x0[i], phys_length_ [i]) * width + 
				//			index ;
				width = phys_length_ [i] ;
				//width = m_extended_length [i] ;
				unsigned long dim_key = (unsigned long) p_lb [i] << 
										num_bits_width | p_tb [i] ;
				key = key << offset | dim_key ;
				//key = key << offset | p_lb [i] ;
				offset += num_bits_dim ;
			}
			leaf_kernel_table & table = 
							m_arr_leaf_kernel_table [index] ;
			std::pair<lk_table_iterator, lk_table_iterator> p = 
											table.equal_range (key) ;
			assert (p.first != p.second) ;
			assert (p.first->first == key) ;
#ifndef NDEBUG
			kernel = p.first->second.kernel ;
#else
			kernel = p.first->second ;
#endif
		}
		else
		{
			// base case
			int index = 0, width = 1, key = 0 ; 
			for (int i = N_RANK-1; i >= 0; --i) 
			{
				//index = (grid.x0[i] + (p_lb [i] >> 1)) * 
				//			width + index ;
				index = pmod(grid.x0[i] + (p_lb [i] >> 1), phys_length_ [i]) * 
							width + index ;
				//index = pmod(grid.x0[i], phys_length_ [i]) * width + 
				//			index ;
				width = phys_length_ [i] ;
				//width = m_extended_length [i] ;
				//key = key << 1 | p_lb [i] >= p_tb [i] ;
			}
			key = p_lb [0] >= p_tb [0] ;
			//assert (key >= 0 && key < 1 << N_RANK) ;
#ifndef NDEBUG
			//kernel = m_kernel_map [(index << N_RANK) + key].kernel ; 
			kernel = m_kernel_map [(index << 1) + key].kernel ; 
#else
			//kernel = m_kernel_map [(index << N_RANK) + key] ;
			kernel = m_kernel_map [(index << 1) + key] ;
#endif
		}
		//char * c = & (m_kernel_map [index << 1]) ;
		//kernel = (int) c [key] ;
		
		assert (kernel != -1) ;
		assert (kernel < m_clone_array->clones.size()) ;
		if (call_boundary) {
			base_case_kernel_boundary(t0, t1, l_father_grid, 
									(*m_clone_array) [kernel]);
		} else { 
			(*m_clone_array) [kernel] (t0, t1, l_father_grid);
		}
	}
}



#undef dx_recursive_boundary_  
#undef dx_recursive_ 
#undef dt_recursive_boundary_ 
#undef dt_recursive_ 
#undef slope_ 
#undef touch_boundary 
#undef base_case_kernel_boundary

#endif
