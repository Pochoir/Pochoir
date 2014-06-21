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
#ifndef AUTO_TUNING_ARBITRARY_CUTS_TRAP_HPP 
#define AUTO_TUNING_ARBITRARY_CUTS_TRAP_HPP 

//#include "auto_tuning_arbitrary_cuts_header.hpp"
#include "auto_tuning_homogeneous_header.hpp"
//#include "auto_tuning_bottomup_header.hpp"

#define dx_recursive_boundary_  (m_algo.dx_recursive_boundary_)
#define dx_recursive_ (m_algo.dx_recursive_)
#define dt_recursive_boundary_ (m_algo.dt_recursive_boundary_)
#define dt_recursive_ (m_algo.dt_recursive_)
#define slope_ m_algo.slope_
#define touch_boundary m_algo.touch_boundary
#define phys_length_ m_algo.phys_length_
#define base_case_kernel_boundary m_algo.base_case_kernel_boundary
#define base_case_kernel_boundary_rectangle m_algo.base_case_kernel_boundary_rectangle
#define uub_boundary m_algo.uub_boundary
#define ulb_boundary m_algo.ulb_boundary
#define lub_boundary m_algo.lub_boundary

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::symbolic_trap_space_cut_interior(int t0,
	int t1, grid_info<N_RANK> const & grid, unsigned long parent_index, 
	F const & f, int * num_zoids, time_type & linkage_time, 
	time_type & child_time, time_type & max_loop_time, int bits)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }

	int child_index = 0 ;
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
                    symbolic_trap_space_time_cut_interior(l_son->t0, 
						l_son->t1, l_son->grid, parent_index, child_index, 
						linkage_time, f);
					child_index++ ;
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                pop_queue(curr_dep_pointer);
				time_type time = 0 ;
				symbolic_trap_space_time_cut_interior(l_father->t0, 
					l_father->t1, l_father->grid, parent_index, child_index, 
					linkage_time, child_time, f, time);
#ifdef SUBSUMPTION_TIME 
				max_loop_time = max(time, max_loop_time) ;
#endif
				child_index++ ;
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
				int cut_dim = bits & 1 << level ;
                if (num_zoids [level] == 0 || ! cut_dim) {
                    // if we can't cut into this dimension, just directly push 
                    // it into the circular queue 
                    push_queue(curr_dep_pointer, level-1, t0, t1, 
							l_father_grid) ;
                } else {
                    // can_cut! 
                    if (cut_lb) {
                        const int mid = (lb/2);
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        // push the middle triangular minizoid (gray) into 
                        // circular queue of (curr_dep) 
                        //
                        l_son_grid.x0[level] = l_start + mid - thres;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_start + mid + thres;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        // cilk_sync 
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        // push the left big trapezoid (black)
                        // into circular queue of (curr_dep + 1)
                        //
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start + mid - thres;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

                        // push the right big trapezoid (black)
                        // into circular queue of (curr_dep + 1)
                        //
                        l_son_grid.x0[level] = l_start + mid + thres;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

                    } // end if (cut_lb) 
                    else {
                        // cut_tb 
                        const int mid = (tb/2);
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);
                        const int ul_start = (l_father_grid.x0[level] + l_father_grid.dx0[level] * lt);

                        // push left black sub-grid into circular queue of (curr_dep) 
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = ul_start + mid;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        // push right black sub-grid into circular queue of (curr_dep) 
                        l_son_grid.x0[level] = ul_start + mid;;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        // cilk_sync 
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        // push the middle gray triangular minizoid into 
                        // circular queue of (curr_dep + 1)
                        //
                        l_son_grid.x0[level] = ul_start + mid;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = ul_start + mid;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                    } // end else (cut_tb) 
                }// end if (can_cut) 
            } // end if (performing a space cut) 
        } // end while (queue_len_[curr_dep] > 0) 
#if !USE_CILK_FOR
//        cilk_sync;
#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}

/* Boundary space cut. Uses trap space cut.
 */
template <int N_RANK> template <typename F, typename BF> 
inline void auto_tune<N_RANK>::symbolic_trap_space_cut_boundary(int t0,
		int t1, grid_info<N_RANK> const & grid, unsigned long parent_index, 
		F const & f, BF const & bf, int * num_zoids, time_type & linkage_time,
		time_type & child_time, time_type & max_loop_time, int bits)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }
	
	int child_index = 0 ;
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
                    symbolic_trap_space_time_cut_boundary(l_son->t0, 
						l_son->t1, l_son->grid, parent_index, child_index, 
						linkage_time, f, bf);
					child_index++ ; //this can be a race.
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                pop_queue(curr_dep_pointer);
				time_type time = 0 ;
				symbolic_trap_space_time_cut_boundary(l_father->t0, 
					l_father->t1, l_father->grid, parent_index, child_index,
					linkage_time, child_time, f, bf, time) ;
#ifdef SUBSUMPTION_TIME 
				max_loop_time = max(time, max_loop_time) ;
#endif
				child_index++ ; 
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
				int cut_dim = bits & 1 << level ;
				//cout << "level " << level << "cut_dim " << cut_dim << endl ;
                if (num_zoids [level] == 0 || ! cut_dim) {
                    // if we can't cut into this dimension, just directly push
                    // it into the circular queue
                    push_queue(curr_dep_pointer, level-1, t0, t1, 
							l_father_grid) ;
                } else  {
                    // can_cut 
                    if (cut_lb) {
                        // if cutting lb, there's no initial cut! 
                        assert(lb != phys_length_[level] || l_father_grid.dx0[level] != 0 || l_father_grid.dx1[level] != 0);
                        const int mid = lb/2;
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        // push the middle gray minizoid
                        // into circular queue of (curr_dep) 
                        //
                        l_son_grid.x0[level] = l_start + mid - thres;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_start + mid + thres;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        // cilk_sync 
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        // push one sub-grid into circular queue of (curr_dep + 1)
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start + mid - thres;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

                        // push one sub-grid into circular queue of (curr_dep + 1)
                        l_son_grid.x0[level] = l_start + mid + thres;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                    } // end if (cut_lb) 
                    else { // cut_tb 
                        if (lb == phys_length_[level] && l_father_grid.dx0[level] == 0 && l_father_grid.dx1[level] == 0) { /* initial cut on the dimension */
						   /* initial cut on the dimension */
                            assert(l_father_grid.dx0[level] == 0);
                            assert(l_father_grid.dx1[level] == 0);
                            //const int mid = tb/2;
							//const int dx = slope_ [level] * lt ;
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);

                            l_son_grid.x0[level] = l_start ;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end ;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
							
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            // push triangle into circular queue of (curr_dep + 1)
                            l_son_grid.x0[level] = l_end ;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = l_end ;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
							/*
                            //draw a triangle with a vertex at midpoint of 
							//top base.
                            l_son_grid.x0[level] = mid - dx ;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = mid + dx ;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            // cilk_sync
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            // push trapezoid into circular queue of (curr_dep + 1)
                            l_son_grid.x0[level] = mid + dx ;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = l_end + mid - dx ;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
							*/
                        } else { /* NOT the initial cut! */
                            const int mid = tb/2;
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);
                            const int ul_start = (l_father_grid.x0[level] + l_father_grid.dx0[level] * lt);
                            // push one sub-grid into circular queue of (curr_dep) 
                            l_son_grid.x0[level] = l_start;
                            l_son_grid.dx0[level] = l_father_grid.dx0[level];
                            l_son_grid.x1[level] = ul_start + mid;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            // push one sub-grid into circular queue of (curr_dep) 
                            l_son_grid.x0[level] = ul_start + mid;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end;
                            l_son_grid.dx1[level] = l_father_grid.dx1[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            // cilk_sync 
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            // push one sub-grid into circular queue of (curr_dep + 1)
                            l_son_grid.x0[level] = ul_start + mid;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = ul_start + mid;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                        }                    
                    } // end if (cut_tb) 
                } // end if (can_cut) 
            } // end if (performing a space cut) 
        } // end while (queue_len_[curr_dep] > 0) 
//#if !USE_CILK_FOR
//        cilk_sync;
//#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::symbolic_trap_space_time_cut_interior(
		int t0, int t1, grid_info<N_RANK> const & grid, 
		unsigned long parent_index, int child_index, time_type & linkage_time,
		time_type & child_time, F const & f, time_type & max_loop_time)
{
	time_type t = 0 ;
	stopwatch * ptr = &m_stopwatch ;
	//stopwatch_stop_and_start(ptr, t) ;
	stopwatch_stop(ptr) ;
	stopwatch_get_elapsed_time(ptr, t) ;
	linkage_time += t ;
	//linkage_time += stopwatch_stop_and_start(&m_stopwatch) ;
    const int lt = t1 - t0;

    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;

	int centroid = 0, width = 1 ; 
	int total_num_subzoids = 1 ;
	int num_subzoids [N_RANK] ;
	
	unsigned long key = 0 ;
	decision_type decision = 0 ;
    for (int i = N_RANK-1; i >= 0; --i) {
        unsigned long lb, tb;
        int thres ;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
		centroid = pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) * width + 
					centroid ;
		/*	cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
			<< " lt " << lt << endl ;
		*/
		assert (lb >= 0) ;
		assert (tb >= 0) ;
		assert (centroid >= 0) ;
		width *= phys_length_ [i] ;

		key <<= num_bits_width ;
		key |= lb ;
		key <<= num_bits_width ;
		key |= tb ;
        thres = slope_[i] * lt ;
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
		if (short_side >= 2 * thres)
		{
			space_cut = true ;
			//set if a space cut can happen in dimension i
			decision |= 1 << (i + 1) ;
			num_subzoids [i] = 3 ;
			total_num_subzoids *= num_subzoids [i] ;
		}
        sim_can_cut |= space_cut ;
    }
	unsigned long index ;
#ifdef TIME_INVARIANCE_INTERIOR
	bool projection_exists = check_and_create_time_invariant_replica (key,
					lt, centroid, centroid, index, grid, (1 << N_RANK) - 1);
#else
	bool projection_exists = check_and_create_space_time_invariant_replica (
								key, lt, centroid, index, grid) ;
#endif
	zoid_type & z = m_zoids [index];
	
	zoid_type & parent = m_zoids [parent_index] ;
	//add the zoid as a child of the parent
	parent.add_child(&z, child_index, index) ;
	assert (index < m_zoids.size()) ;
	
	if (projection_exists)
	{ 
		//Add the projected time of the zoid to that of the parent.
#ifdef SUBSUMPTION_TIME
		max_loop_time = z.max_loop_time ;
#endif
		//a zoid with the projection already exists. return
		child_time += z.time ;
		//start measuring linkage time
		//stopwatch_stop_and_start(ptr, t) ;
		stopwatch_start(ptr) ;
		return ;
	}
	bool time_cut = false ;
	bool divide_and_conquer = false ;
	time_type time_cut_elapsed_time = 0, space_cut_elapsed_time = LONG_MAX ; 
#ifdef FIXED_TIME_CUT
	//cut in time only when space cut is not possible
	if (lt > dt_recursive_ && ! sim_can_cut) 
#else
	if (lt > dt_recursive_)
#endif
	{
		divide_and_conquer = true ;
		time_cut = true ;
		m_zoids [index].set_capacity(max (2, total_num_subzoids)) ;
		m_zoids [index].resize_children(2) ;
	
		time_type time1 = 0, time2 = 0, ltime = 0, ctime = 0 ;
        /* cut into time */
		//stopwatch_stop_and_start(ptr, t) ;
		stopwatch_start(ptr) ;
        int halflt = lt / 2;
        l_son_grid = grid;
        symbolic_trap_space_time_cut_interior(t0, t0+halflt, l_son_grid, 
				index, 0, ltime, ctime, f, time1);
        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        symbolic_trap_space_time_cut_interior(t0+halflt, t1, l_son_grid, 
				index, 1, ltime, ctime, f, time2);
		//calculate the linkage time.
		//stopwatch_stop_and_start(ptr, t) ;
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t) ;
		time_cut_elapsed_time = t + ltime ;
		//time_cut_elapsed_time = stopwatch_stop_and_start(&m_stopwatch) + 
		//												ltime ;
		time1 = max(time1, time2) ;
		max_loop_time = max(time1, max_loop_time) ;
		assert (ltime >= 0) ;
		assert (time_cut_elapsed_time >= 0) ;
		assert (m_zoids [index].num_children == 2) ;
		time_cut_elapsed_time += ctime ;
		/*for (int i = 0 ; i < m_zoids [index].num_children ; i++)
		{
			unsigned long child_index = m_zoids [index].children [i] ;
			//add the time of children.
			time_cut_elapsed_time += m_zoids [child_index].time ;
		}*/
		assert (time_cut_elapsed_time >= 0) ;
		
#ifndef NDEBUG
		m_zoids [index].ttime = time_cut_elapsed_time ;
#endif
    }
	zoid_type bak ;
    if (sim_can_cut) 
	{
		assert (decision) ;
		divide_and_conquer = true ; 
		if (time_cut)
		{
			//back up the time cut children data
			bak = m_zoids [index] ;
			assert(bak.num_children == 2) ;
		}
		m_zoids [index].set_capacity(total_num_subzoids) ;
	}
	
	//set num_children correctly before recursing down to child.
	if (sim_can_cut) 
	{
		zoid_type bak2 ;
		int num_cases = 1 << N_RANK ;
		int best_case = 0, num_children_best_case = 1 ;
#ifdef FIXED_SPACE_CUT
		//do a hyper space cut
		for (int i = (decision >> 1) ; i == (decision >> 1) ; i--)
#else
		for (int i = num_cases - 1 ; i > 0 ; i--)
#endif
		{
			int invalid_case = 0 ;
#ifdef FIXED_SPACE_CUT
			int num_children = total_num_subzoids ;
#else
			int num_children = 1 ;
			for (int j = 0 ; j < N_RANK && ! invalid_case ; j++)
			{
				int like_to_cut_dim_j = i & 1 << j ;
				int can_cut_dim_j = decision & 1 << (j + 1) ;
				//if we like to cut dim j but cannot cut dim j,
				//the case is invalid.
				invalid_case = like_to_cut_dim_j && ! can_cut_dim_j ;
				//find # of subzoids in dim j
				if (like_to_cut_dim_j && can_cut_dim_j)
				{
					num_children *= num_subzoids [j] ;
				}
			}
#endif
			if (invalid_case)
			{
				//invalid case.
				continue ;
			}
			assert (num_children <= total_num_subzoids) ;
			m_zoids [index].resize_children(num_children) ;
			//stopwatch_stop_and_start(ptr, t) ;
			stopwatch_start(ptr) ;
			/* cut into space */
			time_type time = 0, ltime = 0, ctime = 0, elapsed_time = 0  ;
			symbolic_trap_space_cut_interior(t0, t1, grid, index, f, 
							num_subzoids, ltime, ctime, time, i) ;
			//calculate the linkage time
			//stopwatch_stop_and_start(ptr, t) ;
			stopwatch_stop(ptr) ;
			stopwatch_get_elapsed_time(ptr, t) ;
			elapsed_time = t + ltime ;
			//elapsed_time = stopwatch_stop_and_start(&m_stopwatch) + 
			//										ltime ;
			max_loop_time = max(time, max_loop_time) ;
			assert (ltime >= 0) ;
			assert (elapsed_time >= 0) ;
			assert (m_zoids [index].num_children == num_children) ;
			assert (m_zoids [index].num_children <= total_num_subzoids) ;
			elapsed_time += ctime ;
			/*for (int i = 0 ; i < m_zoids [index].num_children ; i++)
			{
				unsigned long child_index = m_zoids [index].children [i] ;
				//add the time of children.
				elapsed_time += m_zoids [child_index].time ;
			}*/
			assert (elapsed_time >= 0) ;
			if (elapsed_time < space_cut_elapsed_time)
			{
				space_cut_elapsed_time = elapsed_time ;
				best_case = i ;
				//back up the zoid with its children.
				bak2 = m_zoids [index] ;
				num_children_best_case = num_children ;
			}
		}
		assert (space_cut_elapsed_time >= 0) ;

		//set the decision with the best case found
		for (int j = 0 ; j < N_RANK ; j++)
		{
			int bit = best_case & 1 << j ;
			decision = (decision & ~(1 << (j + 1))) | ((bit != 0) << (j + 1)) ;
#ifndef NDEBUG
			int decision_bit = decision & 1 << (j + 1) ;
			assert (decision_bit == bit << 1) ;
#endif
		}
		//restore the back up.
		m_zoids [index] = bak2 ;
		assert (m_zoids [index].num_children == num_children_best_case) ;
		assert (m_zoids [index].num_children <= total_num_subzoids) ;
	}
#ifndef NDEBUG
	if (sim_can_cut)
	{
		m_zoids [index].stime = space_cut_elapsed_time ;
	}
#endif
	time_type divide_and_conquer_time = 0 ;
	if (time_cut && sim_can_cut)
	{
		if (space_cut_elapsed_time < time_cut_elapsed_time)
		{
			//space cut is better
			divide_and_conquer_time = space_cut_elapsed_time ;
			//decision is already set for space cut.
		}
		else
		{
			//time cut is better
			divide_and_conquer_time = time_cut_elapsed_time ;
			decision = 1 ;
			m_zoids [index] = bak ; //restore the backup
			assert (m_zoids [index].num_children == 2) ;
		}
	}
	else if (time_cut)
	{
		//time cut is the only choice
		divide_and_conquer_time = time_cut_elapsed_time ;
		decision = 1 ;
		assert (m_zoids [index].num_children == 2) ;
	}
	else if (sim_can_cut)
	{
		//space cut is the only choice
		divide_and_conquer_time = space_cut_elapsed_time ;
		//decision is already set for space cut.
	}

	bool force_divide = false ;
	bool child_divides = false ;
#ifdef SUBSUMPTION_SPACE
	//int max_num_level_divide = -1 ;
	for (int i = 0 ; i < m_zoids [index].num_children ; i++)
	{
		unsigned long child_index = m_zoids [index].children [i] ;
		decision_type d = m_zoids [child_index].decision ;
		if (d & m_space_cut_mask || d & 1)
		{
			child_divides = true ;
			//test if child passed inversion test
			if ((int) m_zoids [child_index].num_level_divide >= 1) 
			{
				force_divide = true ;
				m_zoids [index].num_level_divide = 
							m_zoids [child_index].num_level_divide ;
			}
		}
	}
	/*{
		unsigned long child_index = m_zoids [index].children [i] ;
		decision_type d = m_zoids [child_index].decision ;
		if (d & m_space_cut_mask || d & 1) 
		{
			//find the # of consecutive descendants that divided
			max_num_level_divide = max(max_num_level_divide, 
							(int) m_zoids [child_index].num_level_divide) ;
		}
	}
	m_zoids [index].num_level_divide = min(max_num_level_divide + 1, 
										   DIVIDE_COUNTER) ;
	if (m_zoids [index].num_level_divide >= DIVIDE_COUNTER)
	{
		force_divide = true ;
	}*/
#endif
#ifndef SUBSUMPTION_TIME
	max_loop_time = 0 ;
#endif

    //base case
	//suppose loop_time(z) >= loop_time(z'), z' \in tree(z)
	//if divide_and_conquer_time(z) < max_{z' \in tree(z)} loop_time(z')
	//	then avoid computing loop_time(z).
	//else compute loop_time(z).
	time_type loop_time = LONG_MAX ;
	if ((divide_and_conquer && divide_and_conquer_time < max_loop_time)
		|| force_divide || lt > (m_initial_height + 1) / 2)
	{
		//do not compute loop_time.
		m_zoids [index].decision |= decision ;
		m_zoids [index].time = divide_and_conquer_time ;
	}
	else if (! divide_and_conquer)
	{
		//determine the looping time on the zoid
//		for (int i = 0 ; i < 2 ; i++)
		{
			//stopwatch_stop_and_start(ptr, t) ;
			stopwatch_start(ptr) ;
			f(t0, t1, grid);
			//stopwatch_stop_and_start(ptr, t) ;
			stopwatch_stop(ptr) ;
			stopwatch_get_elapsed_time(ptr, t) ;
			//stopwatch_get_elapsed_time_corrected(ptr, t) ;
			//time_type t2 = stopwatch_stop_and_start(&m_stopwatch) ;
			loop_time = min (t, loop_time) ;
		}
		assert (loop_time >= 0) ;
#ifndef NDEBUG
		m_zoids [index].ltime = loop_time ;
#endif
		//max_loop_time = max(loop_time, max_loop_time) ;
		//max_loop_time = loop_time ;
		//assert (loop_time >= max_loop_time) ;
		//set a flag to indicate that we looped on z.
		m_zoids [index].decision = (decision_type) 1 << 
					  (zoid_type::NUM_BITS_DECISION - 2) ;
		m_zoids [index].time = loop_time ;
	}
	else 
	{
		assert (divide_and_conquer && 
				divide_and_conquer_time >= max_loop_time) ;

		//time_type zoid_loop_time = 0 ; //the loop time of z
#if 0
//#ifdef MARCHUP
		//find the new max loop time in tree(z) instead of looping at z.
		//to do : can we use a max_loop_decision?
		//set the best decision found so far. 
		m_zoids [index].decision |= decision ;
		//loop_time will be the maximum loop time in tree(z). 
        trap_find_mlt_space_time_interior(t0, t1, grid, &(m_zoids [index]),
				necessary_time + projected_time1, f, loop_time, zoid_loop_time);
#else
		//determine the looping time on the zoid
#ifdef SUBSUMPTION_SPACE
		//for (int i = 0 ; i < 2 ; i++)
#endif
		{
			//stopwatch_stop_and_start(ptr, t) ;
			stopwatch_start(ptr) ;
			f(t0, t1, grid);
			//stopwatch_stop_and_start(ptr, t) ;
			stopwatch_stop(ptr) ;
			stopwatch_get_elapsed_time(ptr, t) ;
			//stopwatch_get_elapsed_time_corrected(ptr, t) ;
			//time_type t = stopwatch_stop_and_start(&m_stopwatch) ;
			loop_time = min(t, loop_time) ;
		}
		assert (loop_time >= 0) ;
		//set a flag to indicate that we looped on z.
		m_zoids [index].decision |= (decision_type) 1 << 
					  (zoid_type::NUM_BITS_DECISION - 2) ;
#endif
		
#ifndef NDEBUG
		//check if looping happened at z.
		if (m_zoids [index].decision & 
				1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2))
		{
			m_zoids [index].ltime = loop_time ;
		}
#endif
		//max_loop_time = max(loop_time, max_loop_time) ;
		max_loop_time = loop_time ;
		//if we looped at z, then compare divide and conquer time with
		//loop_time 
		if (m_zoids [index].decision & (decision_type) 1 <<
				(zoid_type::NUM_BITS_DECISION - 2))
		{ 
			if(divide_and_conquer_time < zoid_type::FUZZ * loop_time)
			{
				//choose divide and conquer
				m_zoids [index].decision |= decision ;
				m_zoids [index].time = divide_and_conquer_time ;
#ifdef SUBSUMPTION_SPACE
				if (child_divides)
				{
					m_zoids [index].num_level_divide = 
						(int) m_zoids [index].num_level_divide + 1 ;
				}
#endif
			}
			else
			{
				//choose loop.
				//set decision to loop.
				m_zoids [index].decision = (decision_type) 1 << 
						  (zoid_type::NUM_BITS_DECISION - 2) ;
				m_zoids [index].time = loop_time ;
				/*if (force_divide)
				{
					//an inversion occurs
					cout << "inversion" << endl ;
					for (int i = N_RANK-1; i >= 0; --i)
					{
						cout << " x0 [" << i << "] " << grid.x0 [i] 
						 << " x1 [" << i << "] " << grid.x1 [i] 
						<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
						<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
						<< " lt " << lt << endl ;
					}
					cout << "loop time " << zoid_loop_time <<
							"divide time " << necessary_time + 
							projected_time1 << endl ;
					//exit(0) ;
				}*/
			}
		}
		else
		{
			//we didn't loop at z and found a zoid z' in tree(z) such that
			//divide_and_conquer_time(z) < loop time(z')
			assert (divide_and_conquer_time < loop_time) ;
			m_zoids [index].decision |= decision ;
			m_zoids [index].time = divide_and_conquer_time ;
		}
	}
#ifdef SUBSUMPTION_TIME
	m_zoids [index].max_loop_time = max_loop_time ;
#endif
	
	child_time += m_zoids [index].time ;
	//start measuring linkage time
	//stopwatch_stop_and_start(ptr, t) ;
	stopwatch_start(ptr) ;
}

template <int N_RANK> template <typename F, typename BF>
inline void auto_tune<N_RANK>::symbolic_trap_space_time_cut_boundary(
	int t0, int t1,	grid_info<N_RANK> const & grid, 
	unsigned long parent_index, int child_index, time_type & linkage_time, 
	time_type & child_time,
	F const & f, BF const & bf, time_type & max_loop_time)
{
	time_type t = 0 ;
	stopwatch * ptr = &m_stopwatch ;
	//stopwatch_stop_and_start(ptr, t) ;
	stopwatch_stop(ptr) ;
	stopwatch_get_elapsed_time(ptr, t) ;
	linkage_time += t ;
	//linkage_time += stopwatch_stop_and_start(&m_stopwatch) ;
    const int lt = t1 - t0;

    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;

	int centroid = 0, width = 1 ; 
	int total_num_subzoids = 1 ;
	int num_subzoids [N_RANK] ;
	unsigned long key = 0 ;
	decision_type decision = 0 ;

	struct timespec start1, end1 ;
	int dim_touching_bdry = 0 ;
	int centroid_dim_touching_bdry = 0, width_dim_touching_bdry = 1 ;
    for (int i = N_RANK-1; i >= 0; --i) {
        unsigned long lb, tb;
        int thres ;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
        lb = (grid.x1[i] - grid.x0[i]);
		centroid = pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) * width + 
					centroid ;
		assert (centroid >= 0) ;
		width *= phys_length_ [i] ;
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
		/*cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
			<< " lt " << lt << endl ;*/
        thres = slope_[i] * lt ;
        /*if (lb == phys_length_[i] && grid.dx0[i] == 0 && grid.dx1[i] == 0) 
		{ 
			//set if initial cut on the dimension 
			decision |= 1 << i + 1 + 3 * N_RANK ;
		}*/
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
			dim_touching_bdry |= 1 << i ; 
			centroid_dim_touching_bdry = 
				pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) * 
				width_dim_touching_bdry + centroid_dim_touching_bdry ;
			assert (centroid_dim_touching_bdry >= 0) ;
			width_dim_touching_bdry *= phys_length_ [i] ;
		}
		else
		{
			limit = dx_recursive_[i] ;
		}
		num_subzoids [i] = 0 ;
		if (short_side >= 2 * thres) 
		{
			space_cut = true ;
			//set if a space cut can be done in dimension i
			decision |= 1 << (i + 1) ;
			num_subzoids [i] = 3 ;
			total_num_subzoids *= num_subzoids [i] ;
		}
        sim_can_cut |= space_cut ;

		key <<= num_bits_width ;
		key |= lb ;
		key <<= num_bits_width ;
		key |= tb ;

        call_boundary |= l_touch_boundary;
    }
    if (call_boundary)
	{
        l_dt_stop = dt_recursive_boundary_;
	}
    else
	{
        l_dt_stop = dt_recursive_;
	}
	
	unsigned long index ;
	bool projection_exists = false ;
	if (call_boundary)
	{
#ifdef TIME_INVARIANCE_BOUNDARY
		projection_exists = check_and_create_time_invariant_replica (key, 
			lt, centroid, centroid, index, l_father_grid,(1 << N_RANK) - 1);
#else
		//space-time invariance at boundary
		//you can use space invariance in some dimension
		projection_exists = check_and_create_time_invariant_replica (key, 
			lt, centroid_dim_touching_bdry, centroid, index, l_father_grid, 
			dim_touching_bdry) ;
#endif
	}
	else
	{
#ifdef TIME_INVARIANCE_INTERIOR
		projection_exists = check_and_create_time_invariant_replica (key, 
			lt, centroid, centroid, index, l_father_grid,(1 << N_RANK) - 1);
#else
		//space-time invariance at interior
		projection_exists = check_and_create_space_time_invariant_replica (
						key, lt, centroid, index, l_father_grid) ;
#endif
	}
	zoid_type & z = m_zoids [index] ;
	zoid_type & parent = m_zoids [parent_index] ;
	//add the zoid as a child of the parent
	parent.add_child(&z, child_index, index) ;
	assert (index < m_zoids.size()) ;
	if (projection_exists)
	{ 
		//Add the time of the zoid to the parent.
#ifdef SUBSUMPTION_TIME
		max_loop_time = z.max_loop_time ;
#endif
		//a zoid with the projection already exists. return
		child_time += z.time ;
		//start measuring linkage time
		//stopwatch_stop_and_start(ptr, t) ;
		stopwatch_start(ptr) ;
		return ;
	}
	z.decision |= call_boundary << 
				  (zoid<N_RANK>::NUM_BITS_DECISION - 1) ;

	bool divide_and_conquer = false ;
	bool time_cut = false ;
	time_type time_cut_elapsed_time = 0, space_cut_elapsed_time = LONG_MAX ;
#ifdef FIXED_TIME_CUT
	//cut in time only when space cut is not possible
	if (lt > l_dt_stop && ! sim_can_cut)
#else
	if (lt > l_dt_stop)  //time cut
#endif
	{
		divide_and_conquer = true ;
		time_cut = true ;
		m_zoids [index].set_capacity(max (2, total_num_subzoids)) ;
		m_zoids [index].resize_children(2) ;
	
		time_type time1 = 0, time2 = 0, ltime = 0, ctime = 0 ;
        // cut into time 
		//stopwatch_stop_and_start(ptr, t) ;
		stopwatch_start(ptr) ;
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
    	/*for (int i = N_RANK-1; i >= 0; --i) {
        	touch_boundary(i, lt, l_father_grid) ;
    	}*/
        if (call_boundary) {
            symbolic_trap_space_time_cut_boundary(t0, t0+halflt, l_son_grid,
					index, 0, ltime, ctime, f, bf, time1);
        } else {
            symbolic_trap_space_time_cut_interior(t0, t0+halflt, l_son_grid,
					index, 0, ltime, ctime, f, time1);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            symbolic_trap_space_time_cut_boundary(t0+halflt, t1, l_son_grid,
					index, 1, ltime, ctime, f, bf, time2);
        } else {
            symbolic_trap_space_time_cut_interior(t0+halflt, t1, l_son_grid,
					index, 1, ltime, ctime, f, time2);
        }
		//calculate the linkage time.
		//stopwatch_stop_and_start(ptr, t) ;
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t) ;
		time_cut_elapsed_time = t + ltime ;
		//time_cut_elapsed_time = stopwatch_stop_and_start(&m_stopwatch) 
		//						+ ltime ;
		time1 = max(time1, time2) ;
		max_loop_time = max(time1, max_loop_time) ;
		assert (ltime >= 0) ;
		assert (time_cut_elapsed_time >= 0) ;
		assert (m_zoids [index].num_children == 2) ;
		time_cut_elapsed_time += ctime ;
		/*for (int i = 0 ; i < m_zoids [index].num_children ; i++)
		{
			unsigned long child_index = m_zoids [index].children [i] ;
			//add the time of children.
			time_cut_elapsed_time += m_zoids [child_index].time ;
		}*/
		assert (time_cut_elapsed_time >= 0) ;
#ifndef NDEBUG
		m_zoids [index].ttime = time_cut_elapsed_time ;
#endif
    }
	zoid_type bak ;
	if (sim_can_cut)
	{
		assert (decision) ;
		divide_and_conquer = true ;
		if (time_cut)
		{
			//back up the time cut children data
			bak = m_zoids [index] ;
			assert (bak.num_children == 2) ;
		}
		m_zoids [index].set_capacity(total_num_subzoids) ;
	}
    if (sim_can_cut) 
	{
		zoid_type bak2 ;
		int num_cases = 1 << N_RANK ;
		int best_case = 0 ; 
#ifdef FIXED_SPACE_CUT
		//do a hyper space cut
		for (int i = (decision >> 1) ; i == (decision >> 1) ; i--)
#else
		for (int i = num_cases - 1 ; i > 0 ; i--)
#endif
		{
			int invalid_case = 0 ;
#ifdef FIXED_SPACE_CUT
			int num_children = total_num_subzoids ;
#else
			int num_children = 1 ;
			for (int j = 0 ; j < N_RANK && ! invalid_case ; j++)
			{
				int like_to_cut_dim_j = i & 1 << j ;
				int can_cut_dim_j = decision & 1 << (j + 1) ;
				//if we like to cut dim j but cannot cut dim j,
				//the case is invalid.
				invalid_case = like_to_cut_dim_j && ! can_cut_dim_j ;
				//find # of subzoids in dim j
				if (like_to_cut_dim_j && can_cut_dim_j)
				{
					num_children *= num_subzoids [j] ;
				}
			}
#endif
			if (invalid_case)
			{
				//invalid case.
				continue ;
			}
			assert (num_children <= total_num_subzoids) ;
			m_zoids [index].resize_children(num_children) ;
			//stopwatch_stop_and_start(ptr, t) ;
			stopwatch_start(ptr) ;
			/*for (int i = N_RANK-1; i >= 0; --i) {
				touch_boundary(i, lt, l_father_grid) ;
			}*/
			/* cut into space */
			time_type time = 0, ltime = 0, elapsed_time = 0, ctime = 0  ;
			if (call_boundary) 
			{
				symbolic_trap_space_cut_boundary(t0, t1, l_father_grid, 
					index, f, bf, num_subzoids, ltime, ctime, time, i) ;
			}
			else
			{
				symbolic_trap_space_cut_interior(t0, t1, l_father_grid, 
					index, f, num_subzoids, ltime, ctime, time, i) ;
			}
			//calculate the linkage time
			//stopwatch_stop_and_start(ptr, t) ;
			stopwatch_stop(ptr) ;
			stopwatch_get_elapsed_time(ptr, t) ;
			elapsed_time = t + ltime;
			//elapsed_time = stopwatch_stop_and_start(&m_stopwatch) + ltime;
			max_loop_time = max(time, max_loop_time) ;
			/*
			if (elapsed_time < 0)
			{
				cout << "elapsed_time " << elapsed_time << endl ;
				cout << "stopwatch_get_elapsed_time(&s2) " << stopwatch_get_elapsed_time(&s2) << " rtime " << rtime << endl ;
				for (int i = N_RANK - 1 ; i >= 0 ; i--)
				{	
					cout << " x0 [" << i << "] " << grid.x0 [i] 
					 << " x1 [" << i << "] " << grid.x1 [i] 
					<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
					<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
					<< " lt " << lt << endl ;
				}
			}*/
			assert (ltime >= 0) ;
			assert (elapsed_time >= 0) ;
			assert (m_zoids [index].num_children == num_children) ;
			assert (m_zoids [index].num_children <= total_num_subzoids) ;
			elapsed_time += ctime ;
			/*for (int i = 0 ; i < m_zoids [index].num_children ; i++)
			{
				unsigned long child_index = m_zoids [index].children [i] ;
				//add the time of children.
				elapsed_time += m_zoids [child_index].time ;
			}*/
			assert (elapsed_time >= 0) ;
			if (elapsed_time < space_cut_elapsed_time)
			{
				space_cut_elapsed_time = elapsed_time ;
				best_case = i ;
				//back up the zoid with its children.
				bak2 = m_zoids [index] ;
			}
		}
		assert (space_cut_elapsed_time >= 0) ;
		//set the decision with the best case found
		for (int j = 0 ; j < N_RANK ; j++)
		{
			int bit = best_case & 1 << j ;
			decision = (decision & ~(1 << (j + 1))) | ((bit != 0) << (j + 1)) ;
#ifndef NDEBUG
			int decision_bit = decision & 1 << (j + 1) ;
			assert (decision_bit == bit << 1) ;
#endif
		}
		//restore the back up.
		m_zoids [index] = bak2 ;
		assert (m_zoids [index].num_children <= total_num_subzoids) ;
	}

#ifndef NDEBUG
	if (sim_can_cut)
	{
		m_zoids [index].stime = space_cut_elapsed_time ;
	}
#endif
	time_type divide_and_conquer_time = 0 ;
	if (time_cut && sim_can_cut)
	{
		if (space_cut_elapsed_time < time_cut_elapsed_time)
		{
			//space cut is better
			divide_and_conquer_time = space_cut_elapsed_time ;
			//decision is already set for space cut.
		}
		else
		{
			//time cut is better
			divide_and_conquer_time = time_cut_elapsed_time ;
			decision = 1 ;
			m_zoids [index] = bak ; //restore the backup
			assert (m_zoids [index].num_children == 2) ;
		}
	}
	else if (time_cut)
	{
		//time cut is the only choice
		divide_and_conquer_time = time_cut_elapsed_time ;
		decision = 1 ;
		assert (m_zoids [index].num_children == 2) ;
	}
	else if (sim_can_cut)
	{
		//space cut is the only choice
		divide_and_conquer_time = space_cut_elapsed_time ;
		//decision is already set for space cut.
	}

	bool force_divide = false ;
	bool child_divides = false ;
#ifdef SUBSUMPTION_SPACE
	//int max_num_level_divide = -1 ;
	for (int i = 0 ; i < m_zoids [index].num_children ; i++)
	{
		unsigned long child_index = m_zoids [index].children [i] ;
		decision_type d = m_zoids [child_index].decision ;
		if (d & m_space_cut_mask || d & 1)
		{
			child_divides = true ;
			//test if child passed inversion test
			if ((int) m_zoids [child_index].num_level_divide >= 1) 
			{
				force_divide = true ;
				m_zoids [index].num_level_divide = 
							m_zoids [child_index].num_level_divide ;
			}
		}
	}
	/*{
		unsigned long child_index = m_zoids [index].children [i] ;
		decision_type d = m_zoids [child_index].decision ;
		if (d & m_space_cut_mask || d & 1) 
		{
			//find the # of consecutive descendants that divided
			max_num_level_divide = max(max_num_level_divide, 
							(int) m_zoids [child_index].num_level_divide) ;
		}
	}
	m_zoids [index].num_level_divide = min(max_num_level_divide + 1, 
										   DIVIDE_COUNTER) ;
	if (m_zoids [index].num_level_divide >= DIVIDE_COUNTER)
	{
		force_divide = true ;
	}*/
#endif
#ifndef SUBSUMPTION_TIME
	max_loop_time = 0 ;
#endif
	// base case
	//suppose loop_time(z) >= loop_time(z'), z' \in tree(z)
	//if divide_and_conquer_time(z) < max_{z' \in tree(z)} loop_time(z')
	//	then avoid computing loop_time(z).
	//else compute loop_time(z).
	time_type loop_time = LONG_MAX ;
	if ((divide_and_conquer && divide_and_conquer_time < max_loop_time)
		|| force_divide || lt > (m_initial_height + 1) / 2)
	{
		//do not compute loop_time.
		m_zoids [index].decision |= decision ;
		m_zoids [index].decision |= call_boundary <<
					(zoid_type::NUM_BITS_DECISION - 1) ;
		m_zoids [index].time = divide_and_conquer_time ;
	}
	else if (! divide_and_conquer)
	{
		//determine the looping time on the zoid
//		for (int i = 0 ; i < 2 ; i++)
		{
			//stopwatch_stop_and_start(ptr, t) ;
			stopwatch_start(ptr) ;
			if (call_boundary)
			{
				base_case_kernel_boundary(t0, t1, l_father_grid, bf) ;
			} 
			else 
			{ 
				f(t0, t1, l_father_grid) ;
			}
			//stopwatch_stop_and_start(ptr, t) ;
			stopwatch_stop(ptr) ;
			stopwatch_get_elapsed_time(ptr, t) ;
			//	stopwatch_get_elapsed_time_corrected(ptr, t) ;
			//time_type t = stopwatch_stop_and_start(&m_stopwatch) ;
			loop_time = min (loop_time, t) ;
		}
		//loop_time = tdiff2(&end1, &start1) ;
		assert (loop_time >= 0) ;
#ifndef NDEBUG
		m_zoids [index].ltime = loop_time ;
#endif
		//max_loop_time = max(loop_time, max_loop_time) ;
		//max_loop_time = loop_time ;
		//set a flag to indicate that we looped on z.
		m_zoids [index].decision = (decision_type) 1 << 
					  (zoid_type::NUM_BITS_DECISION - 2) ;
		m_zoids [index].decision |= call_boundary << 
					  (zoid_type::NUM_BITS_DECISION - 1) ;
		m_zoids [index].time = loop_time ;
	}
	else 
	{
		assert (divide_and_conquer && 
				divide_and_conquer_time >= max_loop_time);
		//time_type zoid_loop_time = 0 ; //the loop time of z
#if 0
//#ifdef MARCHUP
		//find the new max loop time in tree(z) instead of looping at z.
		//to do : can we use a max_loop_decision?
		//set the best decision found so far. 
		m_zoids [index].decision |= decision ;
		//loop_time will be the maximum loop time in tree(z). 
        trap_find_mlt_space_time_boundary(t0, t1, grid, &(m_zoids [index]),
			necessary_time + projected_time1, f, bf, loop_time, zoid_loop_time);
#else
		//determine the looping time on the zoid
#ifdef SUBSUMPTION_SPACE
		//for (int i = 0 ; i < 2 ; i++)
#endif
		{
			//stopwatch_stop_and_start(ptr, t) ;
			stopwatch_start(ptr) ;
			if (call_boundary)
			{
				base_case_kernel_boundary(t0, t1, l_father_grid, bf);
			} 
			else 
			{
				f(t0, t1, l_father_grid);
			}
			//stopwatch_stop_and_start(ptr, t) ;
			stopwatch_stop(ptr) ;
			stopwatch_get_elapsed_time(ptr, t) ;
			//stopwatch_get_elapsed_time_corrected(ptr, t) ;
			//time_type t = stopwatch_stop_and_start(&m_stopwatch) ;
			loop_time = min(t, loop_time) ;
		}
		assert (loop_time >= 0) ;
		//set a flag to indicate that we looped on z.
		m_zoids [index].decision |= (decision_type) 1 << 
					  (zoid_type::NUM_BITS_DECISION - 2) ;
#endif	
#ifndef NDEBUG
		//check if looping happened at z.
		if (m_zoids [index].decision & 
			1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2))
		{
			m_zoids [index].ltime = loop_time ;
		}
#endif
		//max_loop_time = max(loop_time, max_loop_time) ;
		max_loop_time = loop_time ;
		//if we looped at z, then compare divide and conquer time with
		//loop_time 
		if (m_zoids [index].decision & (decision_type) 1 <<
				(zoid_type::NUM_BITS_DECISION - 2))
		{ 
			if(divide_and_conquer_time < zoid_type::FUZZ * loop_time)
			{
				//choose divide and conquer
				m_zoids [index].decision |= decision ;
				m_zoids [index].time = divide_and_conquer_time ;
#ifdef SUBSUMPTION_SPACE
				if (child_divides)
				{
					m_zoids [index].num_level_divide = 
						(int) m_zoids [index].num_level_divide + 1 ;
				}
#endif
			}
			else
			{
				//choose loop.
				//set the decision to loop.
				m_zoids [index].decision = (decision_type) 1 << 
						  (zoid_type::NUM_BITS_DECISION - 2) ;
				m_zoids [index].time = loop_time ;
			}
		}
		else
		{
			//we didn't loop at z and found a zoid z' in tree(z) such that
			//divide_and_conquer_time(z) < loop time(z')
			assert(divide_and_conquer_time < loop_time) ;
			m_zoids [index].decision |= decision ;
			m_zoids [index].time = divide_and_conquer_time ;
		}
		m_zoids [index].decision |= call_boundary << 
						  (zoid_type::NUM_BITS_DECISION - 1) ;
	}
#ifdef SUBSUMPTION_TIME
	m_zoids [index].max_loop_time = max_loop_time ;
#endif

	child_time += m_zoids [index].time ;
	//start measuring linkage time
	//stopwatch_stop_and_start(ptr, t) ;
	stopwatch_start(ptr) ;
}

// trap space cuts. 
template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::trap_space_cut_interior_measure(
				int t0, int t1, grid_info<N_RANK> const & grid, 
				simple_zoid_type * projection_zoid, F const & f)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];
	const int lt = t1 - t0 ;

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }

	int child_index = 0 ;
	assert (projection_zoid) ;
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
					//assert(child_index < projection_zoid->num_children) ;
					unsigned long index = projection_zoid->children [child_index] ;
					simple_zoid_type * child = &(m_simple_zoids [index]) ;
					//zoid_type * child = &(m_zoids [index]) ;
					child_index++ ; //looks like a race
                    trap_space_time_cut_interior(l_son->t0, 
						l_son->t1, l_son->grid, child, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                pop_queue(curr_dep_pointer);
				//assert(child_index < projection_zoid->num_children) ;
				unsigned long index = projection_zoid->children [child_index] ;
				simple_zoid_type * child = &(m_simple_zoids [index]) ;
				//zoid_type * child = &(m_zoids [index]) ;
				child_index++ ;
                if (queue_len_[curr_dep_pointer] == 0)
				{
                    this->trap_space_time_cut_interior_measure(l_father->t0,
 							l_father->t1, l_father->grid, child, f);
				}
                else
				{
                    this->trap_space_time_cut_interior_measure(
						l_father->t0, l_father->t1, l_father->grid, child, f) ;
				}
#endif
            } else {
                // performing a space cut on dimension 'level' 
                pop_queue(curr_dep_pointer);
                const grid_info<N_RANK> l_father_grid = l_father->grid;
                const int t0 = l_father->t0, t1 = l_father->t1;
                const int level = l_father->level;
                const int thres = slope_[level] * lt;
                const int lt = (t1 - t0);
                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
                const bool cut_lb = (lb < tb);
				if ((projection_zoid->decision & 1 << (level + 1)) == 0) {
                    // if we can't cut into this dimension, just directly push 
                    // it into the circular queue 
                    push_queue(curr_dep_pointer, level-1, t0, t1, 
							l_father_grid) ;
                } else  {
                    // can_cut! 
					assert ((projection_zoid->decision & 1 << (level + 1)) != 0) ;
                    if (cut_lb) {
                        const int mid = (lb/2);
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        // push the middle triangular minizoid (gray) into 
                        // circular queue of (curr_dep) 
                        l_son_grid.x0[level] = l_start + mid - thres;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_start + mid + thres;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
                        // cilk_sync 
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        // push the left big trapezoid (black)
                        // into circular queue of (curr_dep + 1)
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start + mid - thres;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

                        // push the right big trapezoid (black)
                        // into circular queue of (curr_dep + 1)
                        l_son_grid.x0[level] = l_start + mid + thres;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                    } // end if (cut_lb) 
                    else {
                        // cut_tb 
                        const int mid = (tb/2);
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);
                        const int ul_start = (l_father_grid.x0[level] + l_father_grid.dx0[level] * lt);

                        // push left black sub-grid into circular queue of (curr_dep) 
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = ul_start + mid;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        // push right black sub-grid into circular queue of (curr_dep) 
                        l_son_grid.x0[level] = ul_start + mid;;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
                        // cilk_sync 
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        // push the middle gray triangular minizoid into 
                        // circular queue of (curr_dep + 1)
                        l_son_grid.x0[level] = ul_start + mid;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = ul_start + mid;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
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

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::trap_space_cut_interior(
				int t0, int t1, grid_info<N_RANK> const & grid, 
				simple_zoid_type * projection_zoid, F const & f)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];
	const int lt = t1 - t0 ;

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }

	int child_index = 0 ;
	assert (projection_zoid) ;
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
					//assert(child_index < projection_zoid->num_children) ;
					unsigned long index = projection_zoid->children [child_index] ;
					simple_zoid_type * child = &(m_simple_zoids [index]) ;
					//zoid_type * child = &(m_zoids [index]) ;
					child_index++ ; //looks like a race
                    trap_space_time_cut_interior(l_son->t0, 
						l_son->t1, l_son->grid, child, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                pop_queue(curr_dep_pointer);
				//assert(child_index < projection_zoid->num_children) ;
				unsigned long index = projection_zoid->children [child_index] ;
				simple_zoid_type * child = &(m_simple_zoids [index]) ;
				//zoid_type * child = &(m_zoids [index]) ;
				child_index++ ;
                if (queue_len_[curr_dep_pointer] == 0)
				{
                    this->trap_space_time_cut_interior(l_father->t0,
 							l_father->t1, l_father->grid, child, f);
				}
                else
				{
                    this->trap_space_time_cut_interior(
						l_father->t0, l_father->t1, l_father->grid, child, f) ;
				}
#endif
            } else {
                // performing a space cut on dimension 'level' 
                pop_queue(curr_dep_pointer);
                const grid_info<N_RANK> l_father_grid = l_father->grid;
                const int t0 = l_father->t0, t1 = l_father->t1;
                const int level = l_father->level;
                const int thres = slope_[level] * lt;
                const int lt = (t1 - t0);
                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
                const bool cut_lb = (lb < tb);
				if ((projection_zoid->decision & 1 << (level + 1)) == 0) {
                    // if we can't cut into this dimension, just directly push 
                    // it into the circular queue 
                    push_queue(curr_dep_pointer, level-1, t0, t1, 
							l_father_grid) ;
                } else  {
                    // can_cut! 
					assert ((projection_zoid->decision & 1 << (level + 1)) != 0) ;
                    if (cut_lb) {
                        const int mid = (lb/2);
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        // push the middle triangular minizoid (gray) into 
                        // circular queue of (curr_dep) 
                        l_son_grid.x0[level] = l_start + mid - thres;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_start + mid + thres;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
                        // cilk_sync 
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        // push the left big trapezoid (black)
                        // into circular queue of (curr_dep + 1)
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start + mid - thres;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

                        // push the right big trapezoid (black)
                        // into circular queue of (curr_dep + 1)
                        l_son_grid.x0[level] = l_start + mid + thres;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                    } // end if (cut_lb) 
                    else {
                        // cut_tb 
                        const int mid = (tb/2);
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);
                        const int ul_start = (l_father_grid.x0[level] + l_father_grid.dx0[level] * lt);

                        // push left black sub-grid into circular queue of (curr_dep) 
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = ul_start + mid;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        // push right black sub-grid into circular queue of (curr_dep) 
                        l_son_grid.x0[level] = ul_start + mid;;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
                        // cilk_sync 
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        // push the middle gray triangular minizoid into 
                        // circular queue of (curr_dep + 1)
                        l_son_grid.x0[level] = ul_start + mid;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = ul_start + mid;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
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

template <int N_RANK> template <typename F, typename BF>
inline void auto_tune<N_RANK>::trap_space_cut_boundary_measure(
	int t0, int t1,	grid_info<N_RANK> const & grid, 
	simple_zoid_type * projection_zoid, F const & f, BF const & bf)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];
    const int lt = t1 - t0 ;

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }
	
	int child_index = 0 ;
	assert (projection_zoid) ;
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
					//assert(child_index < projection_zoid->num_children) ;
					unsigned long index = projection_zoid->children [child_index] ;
					simple_zoid_type * child = &(m_simple_zoids [index]) ;
					//zoid_type * child = &(m_zoids [index]) ;
					child_index++ ;
                    trap_space_time_cut_boundary(l_son->t0, 
						l_son->t1, l_son->grid, child, f, bf);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                pop_queue(curr_dep_pointer);
				//assert(child_index < projection_zoid->num_children) ;
				unsigned long index = projection_zoid->children [child_index] ;
				simple_zoid_type * child = &(m_simple_zoids [index]) ;
				//zoid_type * child = &(m_zoids [index]) ;
				child_index++ ;
                if (queue_len_[curr_dep_pointer] == 0) {
                    this->trap_space_time_cut_boundary_measure(l_father->t0,
						l_father->t1, l_father->grid, child, f, bf);
                } else {
                    this->trap_space_time_cut_boundary_measure(
						l_father->t0, l_father->t1, l_father->grid, child, f, bf);
                }
#endif
            } else {
                // performing a space cut on dimension 'level' 
                pop_queue(curr_dep_pointer);
                grid_info<N_RANK> l_father_grid = l_father->grid;
                const int t0 = l_father->t0, t1 = l_father->t1;
                const int level = l_father->level;
                const int thres = slope_[level] * lt;
                const int lt = (t1 - t0);
                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
                const bool cut_lb = (lb < tb);
				if ((projection_zoid->decision & 1 << (level + 1)) == 0) {
                    // if we can't cut into this dimension, just directly push
                    // it into the circular queue
                    push_queue(curr_dep_pointer, level-1, t0, t1, 
							l_father_grid) ;
                } else {
					assert ((projection_zoid->decision & 1 << (level + 1))!= 0);
                    if (cut_lb) {
                        // if cutting lb, there's no initial cut! 
                        assert(lb != phys_length_[level] || l_father_grid.dx0[level] != 0 || l_father_grid.dx1[level] != 0);
                        const int mid = lb/2;
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        // push the middle gray minizoid
                        // into circular queue of (curr_dep) 
                        l_son_grid.x0[level] = l_start + mid - thres;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_start + mid + thres;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        // cilk_sync 
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        // push one sub-grid into circular queue of (curr_dep + 1)
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start + mid - thres;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

                        // push one sub-grid into circular queue of (curr_dep + 1)
                        l_son_grid.x0[level] = l_start + mid + thres;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                    } // end if (cut_lb) 
                    else { // cut_tb 
                        if (lb == phys_length_[level] && l_father_grid.dx0[level] == 0 && l_father_grid.dx1[level] == 0) { /* initial cut on the dimension */
						   /* initial cut on the dimension */
                            assert(l_father_grid.dx0[level] == 0);
                            assert(l_father_grid.dx1[level] == 0);
                            //const int mid = tb/2;
							//const int dx = slope_ [level] * lt ;
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);

                            l_son_grid.x0[level] = l_start ;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end ;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
							
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            // push triangle into circular queue of (curr_dep + 1)
                            l_son_grid.x0[level] = l_end ;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = l_end ;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
							/*
                            //draw a triangle with a vertex at midpoint of 
							//top base.
                            l_son_grid.x0[level] = mid - dx ;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = mid + dx ;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            // cilk_sync 
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            // push trapezoid into circular queue of (curr_dep + 1)
                            l_son_grid.x0[level] = mid + dx ;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = l_end + mid - dx ;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
							*/
                        } else { /* NOT the initial cut! */
                            const int mid = tb/2;
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);
                            const int ul_start = (l_father_grid.x0[level] + l_father_grid.dx0[level] * lt);
                            // push one sub-grid into circular queue of (curr_dep) 
                            l_son_grid.x0[level] = l_start;
                            l_son_grid.dx0[level] = l_father_grid.dx0[level];
                            l_son_grid.x1[level] = ul_start + mid;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            // push one sub-grid into circular queue of (curr_dep) 
                            l_son_grid.x0[level] = ul_start + mid;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end;
                            l_son_grid.dx1[level] = l_father_grid.dx1[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            // cilk_sync 
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            // push one sub-grid into circular queue of (curr_dep + 1)
                            l_son_grid.x0[level] = ul_start + mid;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = ul_start + mid;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                        }                    
                    } // end if (cut_tb) 
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
inline void auto_tune<N_RANK>::trap_space_cut_boundary(
	int t0, int t1,	grid_info<N_RANK> const & grid, 
	simple_zoid_type * projection_zoid, F const & f, BF const & bf)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];
    const int lt = t1 - t0 ;

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }
	
	int child_index = 0 ;
	assert (projection_zoid) ;
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
					//assert(child_index < projection_zoid->num_children) ;
					unsigned long index = projection_zoid->children [child_index] ;
					simple_zoid_type * child = &(m_simple_zoids [index]) ;
					//zoid_type * child = &(m_zoids [index]) ;
					child_index++ ;
                    trap_space_time_cut_boundary(l_son->t0, 
						l_son->t1, l_son->grid, child, f, bf);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                pop_queue(curr_dep_pointer);
				//assert(child_index < projection_zoid->num_children) ;
				unsigned long index = projection_zoid->children [child_index] ;
				simple_zoid_type * child = &(m_simple_zoids [index]) ;
				//zoid_type * child = &(m_zoids [index]) ;
				child_index++ ;
                if (queue_len_[curr_dep_pointer] == 0) {
                    this->trap_space_time_cut_boundary(l_father->t0,
						l_father->t1, l_father->grid, child, f, bf);
                } else {
                    this->trap_space_time_cut_boundary(
						l_father->t0, l_father->t1, l_father->grid, child, f, bf);
                }
#endif
            } else {
                // performing a space cut on dimension 'level' 
                pop_queue(curr_dep_pointer);
                grid_info<N_RANK> l_father_grid = l_father->grid;
                const int t0 = l_father->t0, t1 = l_father->t1;
                const int level = l_father->level;
                const int thres = slope_[level] * lt;
                const int lt = (t1 - t0);
                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
                const bool cut_lb = (lb < tb);
				if ((projection_zoid->decision & 1 << (level + 1)) == 0) {
                    // if we can't cut into this dimension, just directly push
                    // it into the circular queue
                    push_queue(curr_dep_pointer, level-1, t0, t1, 
							l_father_grid) ;
                } else {
					assert ((projection_zoid->decision & 1 << (level + 1))!= 0);
                    if (cut_lb) {
                        // if cutting lb, there's no initial cut! 
                        assert(lb != phys_length_[level] || l_father_grid.dx0[level] != 0 || l_father_grid.dx1[level] != 0);
                        const int mid = lb/2;
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        // push the middle gray minizoid
                        // into circular queue of (curr_dep) 
                        l_son_grid.x0[level] = l_start + mid - thres;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_start + mid + thres;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        // cilk_sync 
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        // push one sub-grid into circular queue of (curr_dep + 1)
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start + mid - thres;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

                        // push one sub-grid into circular queue of (curr_dep + 1)
                        l_son_grid.x0[level] = l_start + mid + thres;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                    } // end if (cut_lb) 
                    else { // cut_tb 
                        if (lb == phys_length_[level] && l_father_grid.dx0[level] == 0 && l_father_grid.dx1[level] == 0) { /* initial cut on the dimension */
						   /* initial cut on the dimension */
                            assert(l_father_grid.dx0[level] == 0);
                            assert(l_father_grid.dx1[level] == 0);
                            //const int mid = tb/2;
							//const int dx = slope_ [level] * lt ;
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);

                            l_son_grid.x0[level] = l_start ;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end ;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
							
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            // push triangle into circular queue of (curr_dep + 1)
                            l_son_grid.x0[level] = l_end ;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = l_end ;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
							/*
                            //draw a triangle with a vertex at midpoint of 
							//top base.
                            l_son_grid.x0[level] = mid - dx ;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = mid + dx ;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            // cilk_sync 
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            // push trapezoid into circular queue of (curr_dep + 1)
                            l_son_grid.x0[level] = mid + dx ;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = l_end + mid - dx ;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
							*/
                        } else { /* NOT the initial cut! */
                            const int mid = tb/2;
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);
                            const int ul_start = (l_father_grid.x0[level] + l_father_grid.dx0[level] * lt);
                            // push one sub-grid into circular queue of (curr_dep) 
                            l_son_grid.x0[level] = l_start;
                            l_son_grid.dx0[level] = l_father_grid.dx0[level];
                            l_son_grid.x1[level] = ul_start + mid;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            // push one sub-grid into circular queue of (curr_dep) 
                            l_son_grid.x0[level] = ul_start + mid;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end;
                            l_son_grid.dx1[level] = l_father_grid.dx1[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            // cilk_sync 
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            // push one sub-grid into circular queue of (curr_dep + 1)
                            l_son_grid.x0[level] = ul_start + mid;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = ul_start + mid;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                        }                    
                    } // end if (cut_tb) 
                } // end if (can_cut) 
			} // end if (performing a space cut) 
        } // end while (queue_len_[curr_dep] > 0) 
#if !USE_CILK_FOR
        cilk_sync;
#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}
template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::
trap_space_time_cut_interior_measure(int t0,
			int t1, grid_info<N_RANK> const & grid, 
			simple_zoid_type * projection_zoid, F const & f)
{
	//measure linkage time
	time_type t = 0 ;
	stopwatch * ptr = &m_stopwatch ;
	stopwatch_stop_and_start_sub(ptr, t) ;
	m_actual_time += t ;

    const int lt = t1 - t0;
    grid_info<N_RANK> l_son_grid;

	assert (projection_zoid) ;
#ifndef NDEBUG
    for (int i = N_RANK-1; i >= 0; --i) {
		grid_info <N_RANK> & grid2 = projection_zoid->info ;
		/*int x0 = pmod(grid.x0 [i], phys_length_ [i]) ;
		int x1 = pmod(grid.x1 [i], phys_length_ [i]) ;

		int x0_ = pmod(grid2.x0 [i], phys_length_ [i]) ;
		int x1_ = pmod(grid2.x1 [i], phys_length_ [i]) ;

		int x2 = pmod(grid.x0[i] + grid.dx0[i] * lt, phys_length_ [i]) ;
		int x3 = pmod(grid.x1[i] + grid.dx1[i] * lt, phys_length_ [i]) ;

		int x2_ = pmod(grid2.x0[i] + grid2.dx0[i] * lt, phys_length_ [i]) ;
		int x3_ = pmod(grid2.x1[i] + grid2.dx1[i] * lt, phys_length_ [i]) ;*/

		int x0 = grid.x0 [i] ;
		int x1 = grid.x1 [i] ;

		int x0_ = grid2.x0 [i] ;
		int x1_ = grid2.x1 [i] ;

		int x2 = grid.x0[i] + grid.dx0[i] * lt ;
		int x3 = grid.x1[i] + grid.dx1[i] * lt ;

		int x2_ = grid2.x0[i] + grid2.dx0[i] * lt ;
		int x3_ = grid2.x1[i] + grid2.dx1[i] * lt ;

		//if (x0 != x0_ || x1 != x1_ || x2 != x2_ || x3 != x3_)
		if (x1 - x0 != x1_ - x0_ || x3 - x2 != x3_ - x2_)
		{
			cout << "zoid and proj zoid differ " << endl ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
				cout << " x0 [" << i << "] " << grid2.x0 [i]
				 << " x1 [" << i << "] " << grid2.x1 [i]
				<< " x2 [" << i << "] " << grid2.x0[i] + grid2.dx0[i] * lt
				<< " x3 [" << i << "] " << grid2.x1[i] + grid2.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
			assert(0) ;
		}
	}
#endif

	//struct timespec start, end ;
	//clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start) ;
	//measure the time to lookup the zoid
	decision_type decision = projection_zoid->decision ;
	//clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end) ;
	//m_dag_lookup_time += tdiff2(&end, &start) ;

    if (decision & m_space_cut_mask) 
	{
        /* cut into space */
		//measure linkage time + time of children.
		stopwatch_stop_and_start_sub(ptr, t) ;
        //return trap_space_cut_interior(t0, t1, grid, projection_zoid, f) ;
        trap_space_cut_interior_measure(t0, t1, grid, projection_zoid, f) ;
		stopwatch_stop_and_start_sub(ptr, t) ;
		m_actual_time += t ;
		return ;
    } 
	else if (decision & 1)
	{
		assert (projection_zoid->num_children == 2) ;
        /* cut into time */
		assert (projection_zoid->children [0]) ;
		assert (projection_zoid->children [1]) ;
        //assert(lt > dt_recursive_);
		unsigned long index1 = projection_zoid->children [0] ;
		unsigned long index2 = projection_zoid->children [1] ;
		//measure linkage time + time of children.
		stopwatch_stop_and_start_sub(ptr, t) ;
        int halflt = lt / 2;
        l_son_grid = grid;
        trap_space_time_cut_interior_measure(t0, t0+halflt, 
									l_son_grid, &(m_simple_zoids [index1]), f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        //return trap_space_time_cut_interior(t0+halflt, t1, 
        trap_space_time_cut_interior_measure(t0+halflt, t1, 
							l_son_grid, &(m_simple_zoids [index2]), f);
		stopwatch_stop_and_start_sub(ptr, t) ;
		m_actual_time += t ;
		return ;
    }
	else
	{
#ifdef TIME_INVARIANCE_INTERIOR
		assert (decision == 
					1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ; // ||
				//projection_zoid->num_children == 0) ;
#else
		assert (decision == 
				3 << (zoid<N_RANK>::NUM_BITS_DECISION - 2) ||
				decision ==
				1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ; // ||
				//projection_zoid->num_children == 0) ;
#endif
		//loop
#ifndef NDEBUG
		struct timespec start, end;
		clock_gettime(CLOCK_MONOTONIC, &start) ;
#endif
		//measure loop time
		stopwatch_stop_and_start_sub(ptr, t) ;
		f(t0, t1, grid);
		stopwatch_stop_and_start_sub(ptr, t) ;
		m_actual_time += t ;
#ifndef NDEBUG
		clock_gettime(CLOCK_MONOTONIC, &end) ;
		double time = tdiff2(&end, &start) ;
		if (time > 100 * projection_zoid->time)
		{
			cout << "runtime " << time * 1e3 << " ms exceeds predicted time " 
				<< projection_zoid->time * 1e3 << " ms" << endl ;
			grid_info <N_RANK> & grid2 = projection_zoid->info ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
			cout << endl ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid2.x0 [i]
				 << " x1 [" << i << "] " << grid2.x1 [i]
				<< " x2 [" << i << "] " << grid2.x0[i] + grid2.dx0[i] * lt
				<< " x3 [" << i << "] " << grid2.x1[i] + grid2.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
		}
#endif
	}
	//measure linkage time
	stopwatch_stop_and_start_sub(ptr, t) ;
}

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::
trap_space_time_cut_interior(int t0,
			int t1, grid_info<N_RANK> const & grid, 
			simple_zoid_type * projection_zoid, F const & f)
{
    const int lt = t1 - t0;
    grid_info<N_RANK> l_son_grid;

	assert (projection_zoid) ;
#ifndef NDEBUG
    for (int i = N_RANK-1; i >= 0; --i) {
		grid_info <N_RANK> & grid2 = projection_zoid->info ;
		/*int x0 = pmod(grid.x0 [i], phys_length_ [i]) ;
		int x1 = pmod(grid.x1 [i], phys_length_ [i]) ;

		int x0_ = pmod(grid2.x0 [i], phys_length_ [i]) ;
		int x1_ = pmod(grid2.x1 [i], phys_length_ [i]) ;

		int x2 = pmod(grid.x0[i] + grid.dx0[i] * lt, phys_length_ [i]) ;
		int x3 = pmod(grid.x1[i] + grid.dx1[i] * lt, phys_length_ [i]) ;

		int x2_ = pmod(grid2.x0[i] + grid2.dx0[i] * lt, phys_length_ [i]) ;
		int x3_ = pmod(grid2.x1[i] + grid2.dx1[i] * lt, phys_length_ [i]) ;*/

		int x0 = grid.x0 [i] ;
		int x1 = grid.x1 [i] ;

		int x0_ = grid2.x0 [i] ;
		int x1_ = grid2.x1 [i] ;

		int x2 = grid.x0[i] + grid.dx0[i] * lt ;
		int x3 = grid.x1[i] + grid.dx1[i] * lt ;

		int x2_ = grid2.x0[i] + grid2.dx0[i] * lt ;
		int x3_ = grid2.x1[i] + grid2.dx1[i] * lt ;

		//if (x0 != x0_ || x1 != x1_ || x2 != x2_ || x3 != x3_)
		if (x1 - x0 != x1_ - x0_ || x3 - x2 != x3_ - x2_)
		{
			cout << "zoid and proj zoid differ " << endl ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
				cout << " x0 [" << i << "] " << grid2.x0 [i]
				 << " x1 [" << i << "] " << grid2.x1 [i]
				<< " x2 [" << i << "] " << grid2.x0[i] + grid2.dx0[i] * lt
				<< " x3 [" << i << "] " << grid2.x1[i] + grid2.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
			assert(0) ;
		}
	}
#endif
    if (projection_zoid->decision & m_space_cut_mask) 
	{
        /* cut into space */
        return trap_space_cut_interior(t0, t1, grid, projection_zoid, f) ;
    } 
	else if (projection_zoid->decision & 1)
	{
		assert (projection_zoid->num_children == 2) ;
        /* cut into time */
		assert (projection_zoid->children [0]) ;
		assert (projection_zoid->children [1]) ;
        //assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
		unsigned long index = projection_zoid->children [0] ;
        trap_space_time_cut_interior(t0, t0+halflt, 
									l_son_grid, &(m_simple_zoids [index]), f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
		index = projection_zoid->children [1] ;
        return trap_space_time_cut_interior(t0+halflt, t1, 
							l_son_grid, &(m_simple_zoids [index]), f);
    }
	else
	{
#ifdef WRITE_DAG
		//file_interior << lt ;
		for (int i = 0 ; i < N_RANK ; i++) 
		{
			unsigned long lb, tb;
			lb = (grid.x1[i] - grid.x0[i]);
			tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
			/*cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;*/
			//file_interior << "," << lb << "," << tb ;
			file_interior [i] << lt << " , " << max(lb, tb)  << endl ;
		}
		//file_interior << endl ;
#endif
#ifdef TIME_INVARIANCE_INTERIOR
		assert (projection_zoid->decision == 
					1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ; // ||
				//projection_zoid->num_children == 0) ;
#else
		assert (projection_zoid->decision == 
				3 << (zoid<N_RANK>::NUM_BITS_DECISION - 2) ||
				projection_zoid->decision ==
				1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ; // ||
				//projection_zoid->num_children == 0) ;
#endif
		//loop
#ifndef NDEBUG
		struct timespec start, end;
		clock_gettime(CLOCK_MONOTONIC, &start) ;
#endif
		f(t0, t1, grid);
#ifndef NDEBUG
		clock_gettime(CLOCK_MONOTONIC, &end) ;
		double time = tdiff2(&end, &start) ;
		if (time > 100 * projection_zoid->time)
		{
			cout << "runtime " << time * 1e3 << " ms exceeds predicted time " 
				<< projection_zoid->time * 1e3 << " ms" << endl ;
			grid_info <N_RANK> & grid2 = projection_zoid->info ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
			cout << endl ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid2.x0 [i]
				 << " x1 [" << i << "] " << grid2.x1 [i]
				<< " x2 [" << i << "] " << grid2.x0[i] + grid2.dx0[i] * lt
				<< " x3 [" << i << "] " << grid2.x1[i] + grid2.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
		}
#endif
		return ;
	}
}

template <int N_RANK> template <typename F, typename BF>
inline void auto_tune<N_RANK>::
trap_space_time_cut_boundary_measure(int t0,
				int t1,	grid_info<N_RANK> const & grid, 
				simple_zoid_type * projection_zoid, F const & f, BF const & bf)
{
	//measure linkage time
	time_type t = 0 ;
	stopwatch * ptr = &m_stopwatch ;
	stopwatch_stop_and_start_sub(ptr, t) ;
	m_actual_time += t ;

    const int lt = t1 - t0;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
	assert (projection_zoid) ;
#ifndef NDEBUG
	bool dim_touching_bdry [N_RANK] ;
#endif
	decision_type call_boundary = 0 ;
    for (int i = N_RANK-1; i >= 0; --i) {
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid) ;
#ifndef NDEBUG
		dim_touching_bdry [i] = l_touch_boundary ;
#endif
		call_boundary |= l_touch_boundary ;
    }
	
#ifndef NDEBUG
    for (int i = N_RANK-1; i >= 0; --i) {
		grid_info <N_RANK> & grid2 = projection_zoid->info ;
		int x0, x1, x2, x3, x0_, x1_, x2_, x3_ ;
		bool error = false ;
		if (dim_touching_bdry [i])
		{
			x0 = pmod(grid.x0 [i], phys_length_ [i]) ;
			x1 = pmod(grid.x1 [i], phys_length_ [i]) ;

			x0_ = pmod(grid2.x0 [i], phys_length_ [i]) ;
			x1_ = pmod(grid2.x1 [i], phys_length_ [i]) ;

			x2 = pmod(grid.x0[i] + grid.dx0[i] * lt, phys_length_ [i]) ;
			x3 = pmod(grid.x1[i] + grid.dx1[i] * lt, phys_length_ [i]) ;

			x2_ = pmod(grid2.x0[i] + grid2.dx0[i] * lt, phys_length_ [i]) ;
			x3_ = pmod(grid2.x1[i] + grid2.dx1[i] * lt, phys_length_ [i]) ;
			if (x0 != x0_ || x1 != x1_ || x2 != x2_ || x3 != x3_) 
			{
				error = true ;
			}
		}
		else
		{
			x0 = grid.x0 [i] ;
			x1 = grid.x1 [i] ;

			x0_ = grid2.x0 [i] ;
			x1_ = grid2.x1 [i] ;

			x2 = grid.x0[i] + grid.dx0[i] * lt ;
			x3 = grid.x1[i] + grid.dx1[i] * lt ;

			x2_ = grid2.x0[i] + grid2.dx0[i] * lt ;
			x3_ = grid2.x1[i] + grid2.dx1[i] * lt ;
			if (x1 - x0 != x1_ - x0_ || x3 - x2 != x3_ - x2_)
			{
				error = true ;
			}
		}
		//if (error && projection_zoid->num_children > 0)
		if (error)
		{
			cout << "zoid and proj zoid differ " << endl ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
			cout << endl ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid2.x0 [i]
				 << " x1 [" << i << "] " << grid2.x1 [i]
				<< " x2 [" << i << "] " << grid2.x0[i] + grid2.dx0[i] * lt
				<< " x3 [" << i << "] " << grid2.x1[i] + grid2.dx1[i] * lt
				<< " lt " << lt << endl ;
			}

			cout << "decision " << (int) projection_zoid->decision << endl ;
			assert(0) ;
		}
	}
#endif
	//struct timespec start, end ;
	//clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start) ;
	//measure the time to lookup the zoid
	decision_type decision = projection_zoid->decision ;
	//clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end) ;
	//m_dag_lookup_time += tdiff2(&end, &start) ;

	if (decision & m_space_cut_mask)
	{
		//cout << "space cut " << endl ;
		//cut into space 
		//measure linkage time + time of children.
		stopwatch_stop_and_start_sub(ptr, t) ;
		if (call_boundary) 
		{
			//return trap_space_cut_boundary(t0, t1, l_father_grid, 
			trap_space_cut_boundary_measure(t0, t1, l_father_grid, 
										projection_zoid, f, bf) ;
		}
		else
		{
			//return trap_space_cut_interior(t0, t1, l_father_grid, 
			trap_space_cut_interior_measure(t0, t1, l_father_grid, 
										projection_zoid, f) ;
		}
		stopwatch_stop_and_start_sub(ptr, t) ;
		m_actual_time += t ;
		return ;
	} 
	else if (decision & 1)
	{
		//cout << "time cut " << endl ;
		assert (projection_zoid->num_children == 2) ;
		// cut into time 
		assert (projection_zoid->children [0]) ;
		assert (projection_zoid->children [1]) ;
		
		unsigned long index1 = projection_zoid->children [0] ;
		unsigned long index2 = projection_zoid->children [1] ;
		//measure linkage time + time of children.
		stopwatch_stop_and_start_sub(ptr, t) ;
		int halflt = lt / 2;
		l_son_grid = l_father_grid;
		if (call_boundary) {
			trap_space_time_cut_boundary_measure(t0, t0+halflt, 
								l_son_grid, &(m_simple_zoids [index1]), f, bf);
		} else {
			trap_space_time_cut_interior_measure(t0, t0+halflt, 
								l_son_grid, &(m_simple_zoids [index1]), f);
		}

		for (int i = 0; i < N_RANK; ++i) {
			l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
			l_son_grid.dx0[i] = l_father_grid.dx0[i];
			l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
			l_son_grid.dx1[i] = l_father_grid.dx1[i];
		}
		if (call_boundary) {
			//return trap_space_time_cut_boundary(t0+halflt, t1, 
			trap_space_time_cut_boundary_measure(t0+halflt, t1, 
								l_son_grid, &(m_simple_zoids [index2]), f, bf);
		} else {
			//return trap_space_time_cut_interior(t0+halflt, t1, 
			trap_space_time_cut_interior_measure(t0+halflt, t1, 
								l_son_grid, &(m_simple_zoids [index2]), f);
		}
		stopwatch_stop_and_start_sub(ptr, t) ;
		m_actual_time += t ;
	}
	else
	{
		//loop
		assert (decision == 
				3 << (zoid<N_RANK>::NUM_BITS_DECISION - 2) ||
				decision ==
				1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ;
#ifndef NDEBUG
		struct timespec start, end;
		clock_gettime(CLOCK_MONOTONIC, &start) ;
#endif
		//measure loop time
		stopwatch_stop_and_start_sub(ptr, t) ;
		if (call_boundary) {
#ifdef TIME_INVARIANCE_BOUNDARY
			assert (decision == 
					3 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ; 
#endif
            base_case_kernel_boundary(t0, t1, l_father_grid, bf);
        } else { 
#ifdef TIME_INVARIANCE_INTERIOR
			assert (decision == 
					1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ;
#endif
            f(t0, t1, l_father_grid);
        }
		stopwatch_stop_and_start_sub(ptr, t) ;
		m_actual_time += t ;
#ifndef NDEBUG
		clock_gettime(CLOCK_MONOTONIC, &end) ;
		double time = tdiff2(&end, &start) ;
		if (time > 100 * projection_zoid->time &&
			projection_zoid->time > 0)
		{
			cout << "runtime " << time * 1e3 << " ms exceeds predicted time " 
				<< projection_zoid->time * 1e3 << " ms" << endl ;
			grid_info <N_RANK> & grid2 = projection_zoid->info ;
		
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
			cout << endl ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid2.x0 [i]
				 << " x1 [" << i << "] " << grid2.x1 [i]
				<< " x2 [" << i << "] " << grid2.x0[i] + grid2.dx0[i] * lt
				<< " x3 [" << i << "] " << grid2.x1[i] + grid2.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
		}
#endif
	}
	//measure linkage time
	stopwatch_stop_and_start_sub(ptr, t) ;
}

template <int N_RANK> template <typename F, typename BF>
inline void auto_tune<N_RANK>::
trap_space_time_cut_boundary(int t0,
				int t1,	grid_info<N_RANK> const & grid, 
				simple_zoid_type * projection_zoid, F const & f, BF const & bf)
{
    const int lt = t1 - t0;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;

	assert (projection_zoid) ;
#ifndef NDEBUG
	bool dim_touching_bdry [N_RANK] ;
#endif
	decision_type call_boundary = 0 ;
    for (int i = N_RANK-1; i >= 0; --i) {
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid) ;
#ifndef NDEBUG
		dim_touching_bdry [i] = l_touch_boundary ;
#endif
		call_boundary |= l_touch_boundary ;
    }
	//decision_type call_boundary = projection_zoid->decision >> 
	//				zoid<N_RANK>::NUM_BITS_DECISION - 1 ;
	//assert (cb == call_boundary) ;
	//decision_type call_boundary = cb ;
	
#ifndef NDEBUG
    for (int i = N_RANK-1; i >= 0; --i) {
		grid_info <N_RANK> & grid2 = projection_zoid->info ;
		int x0, x1, x2, x3, x0_, x1_, x2_, x3_ ;
		bool error = false ;
		if (dim_touching_bdry [i])
		{
			x0 = pmod(grid.x0 [i], phys_length_ [i]) ;
			x1 = pmod(grid.x1 [i], phys_length_ [i]) ;

			x0_ = pmod(grid2.x0 [i], phys_length_ [i]) ;
			x1_ = pmod(grid2.x1 [i], phys_length_ [i]) ;

			x2 = pmod(grid.x0[i] + grid.dx0[i] * lt, phys_length_ [i]) ;
			x3 = pmod(grid.x1[i] + grid.dx1[i] * lt, phys_length_ [i]) ;

			x2_ = pmod(grid2.x0[i] + grid2.dx0[i] * lt, phys_length_ [i]) ;
			x3_ = pmod(grid2.x1[i] + grid2.dx1[i] * lt, phys_length_ [i]) ;
			if (x0 != x0_ || x1 != x1_ || x2 != x2_ || x3 != x3_) 
			{
				error = true ;
			}
		}
		else
		{
			x0 = grid.x0 [i] ;
			x1 = grid.x1 [i] ;

			x0_ = grid2.x0 [i] ;
			x1_ = grid2.x1 [i] ;

			x2 = grid.x0[i] + grid.dx0[i] * lt ;
			x3 = grid.x1[i] + grid.dx1[i] * lt ;

			x2_ = grid2.x0[i] + grid2.dx0[i] * lt ;
			x3_ = grid2.x1[i] + grid2.dx1[i] * lt ;
			if (x1 - x0 != x1_ - x0_ || x3 - x2 != x3_ - x2_)
			{
				error = true ;
			}
		}
		//if (error && projection_zoid->num_children > 0)
		if (error)
		{
			cout << "zoid and proj zoid differ " << endl ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
			cout << endl ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid2.x0 [i]
				 << " x1 [" << i << "] " << grid2.x1 [i]
				<< " x2 [" << i << "] " << grid2.x0[i] + grid2.dx0[i] * lt
				<< " x3 [" << i << "] " << grid2.x1[i] + grid2.dx1[i] * lt
				<< " lt " << lt << endl ;
			}

			cout << "decision " << (int) projection_zoid->decision << endl ;
			assert(0) ;
		}
	}
#endif
	if (projection_zoid->decision & m_space_cut_mask)
	{
		//cout << "space cut " << endl ;
		//cut into space 
		if (call_boundary) 
		{
			return trap_space_cut_boundary(t0, t1, l_father_grid, 
										projection_zoid, f, bf) ;
		}
		else
		{
			return trap_space_cut_interior(t0, t1, l_father_grid, 
										projection_zoid, f) ;
		}
	} 
	else if (projection_zoid->decision & 1)
	{
		//cout << "time cut " << endl ;
		assert (projection_zoid->num_children == 2) ;
		// cut into time 
		assert (projection_zoid->children [0]) ;
		assert (projection_zoid->children [1]) ;
		
		int halflt = lt / 2;
		l_son_grid = l_father_grid;
		unsigned long index = projection_zoid->children [0] ;
		if (call_boundary) {
			trap_space_time_cut_boundary(t0, t0+halflt, 
								l_son_grid, &(m_simple_zoids [index]), f, bf);
		} else {
			trap_space_time_cut_interior(t0, t0+halflt, 
								l_son_grid, &(m_simple_zoids [index]), f);
		}

		for (int i = 0; i < N_RANK; ++i) {
			l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
			l_son_grid.dx0[i] = l_father_grid.dx0[i];
			l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
			l_son_grid.dx1[i] = l_father_grid.dx1[i];
		}
		index = projection_zoid->children [1] ;
		if (call_boundary) {
			return trap_space_time_cut_boundary(t0+halflt, t1, 
								l_son_grid, &(m_simple_zoids [index]), f, bf);
		} else {
			return trap_space_time_cut_interior(t0+halflt, t1, 
								l_son_grid, &(m_simple_zoids [index]), f);
		}
	}
	else
	{
#ifdef WRITE_DAG
		/*if (call_boundary)
		{
			file_boundary << lt ;
		}
		else
		{
			file_interior << lt ;
		}*/
		for (int i = 0 ; i < N_RANK ; i++) 
		{
			unsigned long lb, tb;
			lb = (grid.x1[i] - grid.x0[i]);
			tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
			/*cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;*/
			if (call_boundary)
			{
				//file_boundary << "," << lb << "," << tb ;
				file_boundary [i] << lt << " , " << max(lb, tb)  << endl ;
			}
			else
			{
				//file_interior << "," << lb << "," << tb ;
				file_interior [i] << lt << " , " << max(lb, tb)  << endl ;
			}
		}
		/*if (call_boundary)
		{
			file_boundary << endl ;
		}
		else
		{
			file_interior << endl ;
		}*/
#endif
		//cout << "decision " << (int) projection_zoid->decision << endl ;
		//loop
		assert (projection_zoid->decision == 
				3 << (zoid<N_RANK>::NUM_BITS_DECISION - 2) ||
				projection_zoid->decision ==
				1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ; // ||
				//projection_zoid->num_children == 0) ;
#ifndef NDEBUG
		struct timespec start, end;
		clock_gettime(CLOCK_MONOTONIC, &start) ;
#endif
		if (call_boundary) {
#ifdef TIME_INVARIANCE_BOUNDARY
			assert (projection_zoid->decision == 
					3 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ; // ||
					//projection_zoid->num_children == 0) ;
#endif
            base_case_kernel_boundary(t0, t1, l_father_grid, bf);
        } else { 
#ifdef TIME_INVARIANCE_INTERIOR
			assert (projection_zoid->decision == 
					1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ; // ||
					//projection_zoid->num_children == 0) ;
#endif
            f(t0, t1, l_father_grid);
        }
#ifndef NDEBUG
		clock_gettime(CLOCK_MONOTONIC, &end) ;
		double time = tdiff2(&end, &start) ;
		if (time > 100 * projection_zoid->time &&
			projection_zoid->time > 0)
		{
			cout << "runtime " << time * 1e3 << " ms exceeds predicted time " 
				<< projection_zoid->time * 1e3 << " ms" << endl ;
			grid_info <N_RANK> & grid2 = projection_zoid->info ;
		
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
			cout << endl ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid2.x0 [i]
				 << " x1 [" << i << "] " << grid2.x1 [i]
				<< " x2 [" << i << "] " << grid2.x0[i] + grid2.dx0[i] * lt
				<< " x3 [" << i << "] " << grid2.x1[i] + grid2.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
		}
#endif
	}
}

#undef dx_recursive_boundary_  
#undef dx_recursive_ 
#undef dt_recursive_boundary_ 
#undef dt_recursive_ 
#undef slope_ 
#undef touch_boundary 
#undef base_case_kernel_boundary
#undef base_case_kernel_boundary_rectangle

#endif
