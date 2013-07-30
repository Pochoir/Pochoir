/*
 **********************************************************************************
 *  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 *  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 *                           Charles E. Leiserson <cel@mit.edu>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   Suggestsions:                  yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 *********************************************************************************
 */

#ifndef POCHOIR_WALK_RECURSIVE_HPP
#define POCHOIR_WALK_RECURSIVE_HPP

#include "pochoir_common.hpp"
#include "pochoir_walk.hpp"

#define initial_cut(i) (lb[i] == phys_length_[i])
/* grid.x1[i] >= phys_grid_.x1[i] - stride_[i] - slope_[i] 
 * because we compute the kernel with range [a, b)
 */
template <int N_RANK>
inline bool Algorithm<N_RANK>::touch_boundary(int i, int lt, grid_info<N_RANK> & grid) 
{
    bool interior = false;
    if (grid.x0[i] >= uub_boundary[i] 
     && grid.x0[i] + grid.dx0[i] * lt >= uub_boundary[i]) {
#if (KLEIN == 0)
#if 1
        /* this is for NON klein bottle */
        interior = true;
        /* by this branch, we are assuming the shape is NOT a Klein bottle */
        grid.x0[i] -= phys_length_[i];
        grid.x1[i] -= phys_length_[i];
#else
        interior = false;
#endif
#else
        /* this is for klein bottle! */
#if 1
#if DEBUG
        fprintf(stderr, "Before klein_region: \n");
        print_grid(stderr, 0, lt, grid);
#endif
        interior = true;
//        fprintf(stderr, "call klein_region!\n");
        klein_region(grid, phys_grid_);
#if DEBUG
        fprintf(stderr, "After klein_region: \n");
        print_grid(stderr, 0, lt, grid);
#endif
#else
        interior = false;
#endif
#endif
    } else if (grid.x1[i] <= ulb_boundary[i] 
            && grid.x1[i] + grid.dx1[i] * lt <= ulb_boundary[i]
            && grid.x0[i] >= lub_boundary[i]
            && grid.x0[i] + grid.dx0[i] * lt >= lub_boundary[i]) {
        interior = true;
    } else {
        interior = false;
    }
    return !interior;
}

template <int N_RANK>
inline bool Algorithm<N_RANK>::within_boundary(int t0, int t1, grid_info<N_RANK> & grid)
{
    bool l_touch_boundary = false;
    int lt = t1 - t0;
    for (int i = 0; i < N_RANK; ++i) {
        l_touch_boundary = l_touch_boundary || touch_boundary(i, lt, grid);
    }
    return !l_touch_boundary;
}


/*
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
*/
/* ************************************************************************************** */
/* following are the procedures for obase with duality , always cutting based on shorter bar
 */

template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::shorter_duo_sim_obase_space_cut(int t0, int t1, grid_info<N_RANK> const grid, F const & f)
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
                    shorter_duo_sim_obase_bicut(l_son->t0, l_son->t1, l_son->grid, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0)
                    shorter_duo_sim_obase_bicut(l_father->t0, l_father->t1, l_father->grid, f);
                else
                    cilk_spawn shorter_duo_sim_obase_bicut(l_father->t0, l_father->t1, l_father->grid, f);
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
                if (!can_cut) {
                    // if we can't cut into this dimension, just directly push 
                    // it into the circular queue 
                    //
                    push_queue(curr_dep_pointer, level-1, t0, t1, l_father_grid);
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
                } // end if (can_cut) 
            } // end if (performing a space cut) 
        } // end while (queue_len_[curr_dep] > 0) 
#if !USE_CILK_FOR
        cilk_sync;
#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}

/* This is for boundary region space cut! , always cutting based on the shorter bar
 */

template <int N_RANK> template <typename F, typename BF>
inline void Algorithm<N_RANK>::shorter_duo_sim_obase_space_cut_p(int t0, int t1, grid_info<N_RANK> const grid, F const & f, BF const & bf)
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
                    shorter_duo_sim_obase_bicut_p(l_son->t0, l_son->t1, l_son->grid, f, bf);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0) {
                    shorter_duo_sim_obase_bicut_p(l_father->t0, l_father->t1, l_father->grid, f, bf);
                } else {
                    cilk_spawn shorter_duo_sim_obase_bicut_p(l_father->t0, l_father->t1, l_father->grid, f, bf);
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
                const bool can_cut = CAN_CUT_B ;
                if (!can_cut) {
                    // if we can't cut into this dimension, just directly push
                    // it into the circular queue
                    //
                    push_queue(curr_dep_pointer, level-1, t0, t1, l_father_grid);
                } else {
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
//                        if (lb == phys_length_[level] && l_father_grid.dx0[level] == 0 && l_father_grid.dx1[level] == 0) { // initial cut on the dimension 
//                            assert(l_father_grid.dx0[level] == 0);
//                            assert(l_father_grid.dx1[level] == 0);
//                            const int mid = tb/2;
//                            grid_info<N_RANK> l_son_grid = l_father_grid;
//                            const int l_start = (l_father_grid.x0[level]);
//                            const int l_end = (l_father_grid.x1[level]);
//                            const int ul_start = (l_father_grid.x0[level] + l_father_grid.dx0[level] * lt);
//                            // merge the big black trapezoids 
//                            l_son_grid.x0[level] = ul_start + mid;
//                            l_son_grid.dx0[level] = slope_[level];
//                            l_son_grid.x1[level] = l_end + (ul_start - l_start) + mid;
//                            l_son_grid.dx1[level] = -slope_[level];
//                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//
//                            // cilk_sync 
//                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
//                            // push middle minizoid into circular queue of (curr_dep + 1)
//                            l_son_grid.x0[level] = ul_start + mid;
//                            l_son_grid.dx0[level] = -slope_[level];
//                            l_son_grid.x1[level] = ul_start + mid;
//                            l_son_grid.dx1[level] = slope_[level];
//                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//                        } else { // NOT the initial cut! 
                        if (lb == phys_length_[level] && l_father_grid.dx0[level] == 0 && l_father_grid.dx1[level] == 0) { /* initial cut on the dimension */
						   /* initial cut on the dimension */
                            assert(l_father_grid.dx0[level] == 0);
                            assert(l_father_grid.dx1[level] == 0);
                            const int mid = tb/2;
							const int dx = slope_ [level] * lt ;
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);
                            //draw a triangle with a vertex at midpoint of 
							//top base.
                            l_son_grid.x0[level] = mid - dx ;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = mid + dx ;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            /* cilk_sync */
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            // push trapezoid into circular queue of (curr_dep + 1)
                            l_son_grid.x0[level] = mid + dx ;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = l_end + mid - dx ;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
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


/* This is the version for interior region cut! */
template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::shorter_duo_sim_obase_bicut(int t0, int t1, grid_info<N_RANK> const grid, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;
	//cout << " t1 " << t1 << " t0 " << t0 << endl ;
#if STAT
    int l_count_cut = 0;
    int l_bottom_total_area = 1;
    int l_top_total_area = 1;
    int l_total_points;
#endif

    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
		/*cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
			<< " lt " << lt << endl ;*/
        bool cut_lb = (lb < tb);
        thres = (slope_[i] * lt);
        sim_can_cut = SIM_CAN_CUT_I ;
        /* as long as there's one dimension can conduct a cut, we conduct a 
         * multi-dimensional cut!
         */
#if STAT
        l_count_cut = (l_can_cut ? l_count_cut+1 : l_count_cut);
#endif
    }

#if STAT
//    l_total_points = l_bottom_total_area * t1 / 3 - l_top_total_area * t0 / 3;
#endif
    if (sim_can_cut) {
        /* cut into space */
#if STAT
    // sim_count_cut[l_count_cut] = (l_count_cut > 0 ? sim_count_cut[l_count_cut] + 1 : sim_count_cut[l_count_cut]);
        ++sim_count_cut[l_count_cut];
#endif
        shorter_duo_sim_obase_space_cut(t0, t1, grid, f);
        return;
    // } else if (lt > dt_recursive_ && l_total_points > Z) {
    } else if (lt > dt_recursive_) {
        /* cut into time */
//        assert(dt_recursive_ >= r_t);
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        //int halflt = (lt + 1) / 2;
        l_son_grid = grid;
        shorter_duo_sim_obase_bicut(t0, t0+halflt, l_son_grid, f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        shorter_duo_sim_obase_bicut(t0+halflt, t1, l_son_grid, f);
        return;
    } else {
        // base case
#if DEBUG
        printf("call interior!\n");
        print_grid(stdout, t0, t1, grid);
        // fprintf(stderr, "l_total_points = %d\n", l_total_points);
#endif
#if STAT
        ++interior_region_count;
#endif
        f(t0, t1, grid);
//        base_case_kernel_interior(t0, t1, grid, f);
        return;
    }  
}

/* This is the version for boundary region cut! */
template <int N_RANK> template <typename F, typename BF>
inline void Algorithm<N_RANK>::shorter_duo_sim_obase_bicut_p(int t0, int t1, grid_info<N_RANK> const grid, F const & f, BF const & bf)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;

	//cout << " t1 " << t1 << " t0 " << t0 << endl ;
#if STAT
    int l_count_cut = 0;
    int l_bottom_total_area = 1;
    int l_top_total_area = 1;
    int l_total_points;
#endif

    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
		/*cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
			<< " lt " << lt << endl ; */
        thres = (slope_[i] * lt);
        /* l_father_grid may be mapped to a new region in touch_boundary() */
        /* for the initial cut, we exclude the begining and end point to minimize
         * the overhead on boundary
        */
        /* lb == phys_length_[i] indicates an initial cut! */
        bool cut_lb = (lb < tb);
        sim_can_cut = SIM_CAN_CUT_B ;
        call_boundary |= l_touch_boundary;
#if STAT
        l_count_cut = (l_can_cut ? l_count_cut + 1 : l_count_cut);
        l_bottom_total_area *= lb;
        l_top_total_area *= tb;
#endif
    }

#if STAT
//    l_total_points = l_bottom_total_area * t1 / 3 - l_top_total_area * t0 / 3;
#endif

    if (sim_can_cut) {
        /* cut into space */
        /* push the first l_father_grid that can be cut into the circular queue */
        /* boundary cuts! */
#if STAT
    // sim_count_cut[l_count_cut] = (l_count_cut > 0 ? sim_count_cut[l_count_cut] + 1 : sim_count_cut[l_count_cut]);
        ++sim_count_cut[l_count_cut];
#endif
        if (call_boundary) 
            shorter_duo_sim_obase_space_cut_p(t0, t1, l_father_grid, f, bf);
        else
            shorter_duo_sim_obase_space_cut(t0, t1, l_father_grid, f);
        return;
    } 

    if (call_boundary)
        l_dt_stop = dt_recursive_boundary_;
    else
        l_dt_stop = dt_recursive_;

    if (lt > l_dt_stop) {
        /* cut into time */
        int halflt = lt / 2;
        //int halflt = (lt + 1) / 2;
        l_son_grid = l_father_grid;
        if (call_boundary) {
            shorter_duo_sim_obase_bicut_p(t0, t0+halflt, l_son_grid, f, bf);
        } else {
            shorter_duo_sim_obase_bicut(t0, t0+halflt, l_son_grid, f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            shorter_duo_sim_obase_bicut_p(t0+halflt, t1, l_son_grid, f, bf);
        } else {
            shorter_duo_sim_obase_bicut(t0+halflt, t1, l_son_grid, f);
        }
        return;
    } 

    // if (l_total_area <= Z || base_cube_t) {
        /* for base_cube_t: -- prevent too small time cut! 
         *      (cut_lb && lb > dx_recursive_boundary_ && lb < 2 * thres)
         *  ||  (!cut_lb && tb > dx_recursive_boundary_ && lb < thres)
         */
        // base case
#if DEBUG
        printf("call boundary!\n");
        print_grid(stdout, t0, t1, l_father_grid);
#endif
#if STAT
        ++boundary_region_count;
        boundary_points_count += l_total_points;
#endif
        if (call_boundary) {
            base_case_kernel_boundary(t0, t1, l_father_grid, bf);
        } else {
            f(t0, t1, l_father_grid);
        }
        return;
}

#endif
