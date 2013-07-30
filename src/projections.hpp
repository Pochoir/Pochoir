/*
 * ============================================================================
 *       Filename:  projections.hpp
 *    Description:  Has routines 
 *					1. that count the # of projections in
 *					the default space/time cut scheme and the
 *					modified space/power of two time cut scheme.
 *					2. that implement the modified space/power of two
 *					time cuts.
 *					The code uses the time/space cut code framework in
 *					pochoir_walk_recursive.hpp
 *        Created:  10/02/2012
 *         Author:  Eka Palamadai, epn@mit.edu
 * ============================================================================
 */
#ifndef PROJECTIONS_HPP
#define PROJECTIONS_HPP

#include "pochoir_common.hpp"
#include "pochoir_walk.hpp"

#define initial_cut(i) (lb[i] == phys_length_[i])
/* grid.x1[i] >= phys_grid_.x1[i] - stride_[i] - slope_[i] 
 * because we compute the kernel with range [a, b)
 */
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
bool print_projections ;

#ifdef DEFAULT_SPACE_CUT 
template <int N_RANK>
inline void Algorithm<N_RANK>::space_cut_interior_projections(int t0, int t1, 
												grid_info<N_RANK> const grid,
												const int index)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }

    /* set up the initial grid */
    push_queue(0, N_RANK-1, t0, t1, grid);
    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
        const int curr_dep_pointer = (curr_dep & 0x1);
        while (queue_len_[curr_dep_pointer] > 0) {
            top_queue(curr_dep_pointer, l_father);
            if (l_father->level < 0) {
                /* spawn all the grids in circular_queue_[curr_dep][] */
#if USE_CILK_FOR 
                /* use cilk_for to spawn all the sub-grid */
// #pragma cilk_grainsize = 1
                for (int j = 0; j < queue_len_[curr_dep_pointer]; ++j) {
                    int i = pmod((queue_head_[curr_dep_pointer]+j), ALGOR_QUEUE_SIZE);
                    queue_info * l_son = &(circular_queue_[curr_dep_pointer][i]);
                    /* assert all the sub-grid has done N_RANK spatial cuts */
                    assert(l_son->level == -1);
                    space_time_cut_interior_projections(l_son->t0, l_son->t1, l_son->grid,
											index) ;
                } 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                /* use cilk_spawn to spawn all the sub-grid */
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0)
                    space_time_cut_interior_projections(l_father->t0, l_father->t1, 
											l_father->grid, index) ;
                else
                    space_time_cut_interior_projections(l_father->t0, l_father->t1, 
											l_father->grid, index) ;
#endif
            } else {
                /* performing a space cut on dimension 'level' */
                pop_queue(curr_dep_pointer);
                const grid_info<N_RANK> l_father_grid = l_father->grid;
                const int t0 = l_father->t0, t1 = l_father->t1;
                const int lt = (t1 - t0);
                const int level = l_father->level;
                const int thres = slope_[level] * lt;
                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
                const bool cut_lb = (lb < tb);
                bool can_cut = cut_lb ? (lb >= 2 * thres) : (tb >= 2 * thres);
				can_cut = can_cut && (lt > 1) ;
                if (!can_cut) {
                    /* if we can't cut into this dimension, just directly push 
                     * it into the circular queue 
                     */
                    push_queue(curr_dep_pointer, level-1, t0, t1, l_father_grid);
                } else {
                    /* can_cut! */
                    if (cut_lb) {
                        const int mid = (lb/2);
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        /* push the middle triangular minizoid (gray) into 
                         * circular queue of (curr_dep) 
                         */
                        l_son_grid.x0[level] = l_start + mid - thres;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_start + mid + thres;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        /* cilk_sync */
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        /* push the left big trapezoid (black)
                         * into circular queue of (curr_dep + 1)
                         */
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start + mid - thres;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

                        /* push the right big trapezoid (black)
                         * into circular queue of (curr_dep + 1)
                         */
                        l_son_grid.x0[level] = l_start + mid + thres;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

                    } /* end if (cut_lb) */
                    else {
                        /* cut_tb */
                        const int mid = (tb/2);
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);
                        const int ul_start = (l_father_grid.x0[level] + l_father_grid.dx0[level] * lt);

                        /* push left black sub-grid into circular queue of (curr_dep) */
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = ul_start + mid;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        /* push right black sub-grid into circular queue of (curr_dep) */
                        l_son_grid.x0[level] = ul_start + mid;;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        /* push the middle gray triangular minizoid into 
                         * circular queue of (curr_dep + 1)
                         */
                        l_son_grid.x0[level] = ul_start + mid;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = ul_start + mid;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                    } /* end else (cut_tb) */
                } /* end if (can_cut) */
            } /* end if (performing a space cut) */
        } /* end while (queue_len_[curr_dep] > 0) */
        assert(queue_len_[curr_dep_pointer] == 0);
    } /* end for (curr_dep < N_RANK+1) */
}

template <int N_RANK> 
inline void Algorithm<N_RANK>::space_cut_boundary_projections(int t0, int t1, 
												grid_info<N_RANK> const grid,
												const int index)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }

    /* set up the initial grid */
    push_queue(0, N_RANK-1, t0, t1, grid);
    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
        const int curr_dep_pointer = (curr_dep & 0x1);
        while (queue_len_[curr_dep_pointer] > 0) {
            top_queue(curr_dep_pointer, l_father);
            if (l_father->level < 0) {
                /* spawn all the grids in circular_queue_[curr_dep][] */
#if USE_CILK_FOR 
                /* use cilk_for to spawn all the sub-grid */
// #pragma cilk_grainsize = 1
                for (int j = 0; j < queue_len_[curr_dep_pointer]; ++j) {
                    int i = pmod((queue_head_[curr_dep_pointer]+j), ALGOR_QUEUE_SIZE);
                    queue_info * l_son = &(circular_queue_[curr_dep_pointer][i]);
                    /* assert all the sub-grid has done N_RANK spatial cuts */
                    assert(l_son->level == -1);
                    space_time_cut_boundary_projections(l_son->t0, l_son->t1, l_son->grid,
											index);
                } 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                /* use cilk_spawn to spawn all the sub-grid */
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0) {
                    space_time_cut_boundary_projections(l_father->t0, l_father->t1, 
											l_father->grid, index);
                } else {
                    space_time_cut_boundary_projections(l_father->t0, l_father->t1, 
											l_father->grid, index);
                }
#endif
            } else {
                /* performing a space cut on dimension 'level' */
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
        		bool can_cut = cut_lb ? (lb >= 2 * thres) : (tb >= 2 * thres) ;
				can_cut = can_cut && (lt > 1) ;
                if (!can_cut) {
                    /* if we can't cut into this dimension, just directly push
                     * it into the circular queue
                    */
                    push_queue(curr_dep_pointer, level-1, t0, t1, l_father_grid);
                } else {
                    /* can_cut */
                    if (cut_lb) {
                        /* if cutting lb, there's no initial cut! */
                        assert(lb != phys_length_[level] || l_father_grid.dx0[level] != 0 || l_father_grid.dx1[level] != 0);
                        const int mid = lb/2;
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        /* push the middle gray minizoid
                         * into circular queue of (curr_dep) 
                         */
                        l_son_grid.x0[level] = l_start + mid - thres;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_start + mid + thres;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        /* cilk_sync */
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        /* push one sub-grid into circular queue of (curr_dep + 1)*/
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start + mid - thres;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

                        /* push one sub-grid into circular queue of (curr_dep + 1)*/
                        l_son_grid.x0[level] = l_start + mid + thres;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                    } /* end if (cut_lb) */
                    else { /* cut_tb */
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
                            /* push one sub-grid into circular queue of (curr_dep) */
                            l_son_grid.x0[level] = l_start;
                            l_son_grid.dx0[level] = l_father_grid.dx0[level];
                            l_son_grid.x1[level] = ul_start + mid;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            /* push one sub-grid into circular queue of (curr_dep) */
                            l_son_grid.x0[level] = ul_start + mid;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end;
                            l_son_grid.dx1[level] = l_father_grid.dx1[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            /* cilk_sync */
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            /* push one sub-grid into circular queue of (curr_dep + 1)*/
                            l_son_grid.x0[level] = ul_start + mid;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = ul_start + mid;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                        }                    
                    } /* end if (cut_tb) */
                } /* end if (can_cut) */
            } /* end if (performing a space cut) */
        } /* end while (queue_len_[curr_dep] > 0) */
        assert(queue_len_[curr_dep_pointer] == 0);
    } /* end for (curr_dep < N_RANK+1) */
}

#else
//The new space cut algorithm for interior region
//A space cut in a trapezoid is done by connecting the ends of the shorter side 
//with the longer side
template <int N_RANK>
inline void Algorithm<N_RANK>::space_cut_interior_projections(int t0, int t1, 
												grid_info<N_RANK> const grid,
												const int index)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }

    /* set up the initial grid */
    push_queue(0, N_RANK-1, t0, t1, grid);
    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
        const int curr_dep_pointer = (curr_dep & 0x1);
        while (queue_len_[curr_dep_pointer] > 0) {
            top_queue(curr_dep_pointer, l_father);
            if (l_father->level < 0) {
                /* spawn all the grids in circular_queue_[curr_dep][] */
#if USE_CILK_FOR 
                /* use cilk_for to spawn all the sub-grid */
// #pragma cilk_grainsize = 1
                for (int j = 0; j < queue_len_[curr_dep_pointer]; ++j) {
                    int i = pmod((queue_head_[curr_dep_pointer]+j), ALGOR_QUEUE_SIZE);
                    queue_info * l_son = &(circular_queue_[curr_dep_pointer][i]);
                    /* assert all the sub-grid has done N_RANK spatial cuts */
                    assert(l_son->level == -1);
                    space_time_cut_interior_projections(l_son->t0, l_son->t1, l_son->grid,
											index) ;
                } 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                /* use cilk_spawn to spawn all the sub-grid */
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0)
                    space_time_cut_interior_projections(l_father->t0, l_father->t1, 
											l_father->grid, index) ;
                else
                    space_time_cut_interior_projections(l_father->t0, l_father->t1, 
											l_father->grid, index) ;
#endif
            } else {
                /* performing a space cut on dimension 'level' */
                pop_queue(curr_dep_pointer);
                const grid_info<N_RANK> l_father_grid = l_father->grid;
                const int t0 = l_father->t0, t1 = l_father->t1;
                const int lt = (t1 - t0);
                const int level = l_father->level;
                const int thres = slope_[level] * lt;
                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
                const bool cut_lb = (lb < tb);
                bool can_cut = cut_lb ? (lb >= 2 * thres) : (tb >= 2 * thres);
				can_cut = can_cut && (lt > 1) ;
                if (!can_cut) {
                    /* if we can't cut into this dimension, just directly push 
                     * it into the circular queue 
                     */
                    push_queue(curr_dep_pointer, level-1, t0, t1, l_father_grid);
                } else {
                    /* can_cut! */
                    if (cut_lb) {
                        const int mid = (lb/2);
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        l_son_grid.x0[level] = l_start ;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_end ;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start ;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                        l_son_grid.x0[level] = l_end ;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                    } /* end if (cut_lb) */
                    else {
                        /* cut_tb */
                        const int mid = (tb/2);
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

						const int dx0_h = l_father_grid.dx0[level] * lt;
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start + 2 * dx0_h ;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

						const int dx1_h = l_father_grid.dx1[level] * lt;
                        l_son_grid.x0[level] = l_end + 2 * dx1_h ;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        l_son_grid.x0[level] = l_start + 2 * dx0_h ;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end + 2 * dx1_h ;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                    } /* end else (cut_tb) */
                } /* end if (can_cut) */
            } /* end if (performing a space cut) */
        } /* end while (queue_len_[curr_dep] > 0) */
        assert(queue_len_[curr_dep_pointer] == 0);
    } /* end for (curr_dep < N_RANK+1) */
}

/* This is for boundary region space cut! , always cutting based on the shorter bar
 */
//new space cut algorithm
template <int N_RANK> 
inline void Algorithm<N_RANK>::space_cut_boundary_projections(int t0, int t1, 
												grid_info<N_RANK> const grid,
												const int index)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }

    /* set up the initial grid */
    push_queue(0, N_RANK-1, t0, t1, grid);
    for (int curr_dep = 0; curr_dep < N_RANK+1; ++curr_dep) {
        const int curr_dep_pointer = (curr_dep & 0x1);
        while (queue_len_[curr_dep_pointer] > 0) {
            top_queue(curr_dep_pointer, l_father);
            if (l_father->level < 0) {
                /* spawn all the grids in circular_queue_[curr_dep][] */
#if USE_CILK_FOR 
                /* use cilk_for to spawn all the sub-grid */
// #pragma cilk_grainsize = 1
                for (int j = 0; j < queue_len_[curr_dep_pointer]; ++j) {
                    int i = pmod((queue_head_[curr_dep_pointer]+j), ALGOR_QUEUE_SIZE);
                    queue_info * l_son = &(circular_queue_[curr_dep_pointer][i]);
                    /* assert all the sub-grid has done N_RANK spatial cuts */
                    assert(l_son->level == -1);
                    space_time_cut_boundary_projections(l_son->t0, l_son->t1, l_son->grid,
											index);
                } 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                /* use cilk_spawn to spawn all the sub-grid */
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0) {
                    space_time_cut_boundary_projections(l_father->t0, l_father->t1, 
											l_father->grid, index);
                } else {
                    space_time_cut_boundary_projections(l_father->t0, l_father->t1, 
											l_father->grid, index);
                }
#endif
            } else {
                /* performing a space cut on dimension 'level' */
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
        		bool can_cut = cut_lb ? (lb >= 2 * thres) : (tb >= 2 * thres) ;
				can_cut = can_cut && (lt > 1) ;
                if (!can_cut) {
                    /* if we can't cut into this dimension, just directly push
                     * it into the circular queue
                    */
                    push_queue(curr_dep_pointer, level-1, t0, t1, l_father_grid);
                } else {
                    /* can_cut */
                    if (cut_lb) {
                        /* if cutting lb, there's no initial cut! */
                        assert(lb != phys_length_[level] || l_father_grid.dx0[level] != 0 || l_father_grid.dx1[level] != 0);
                        const int mid = lb/2;
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        l_son_grid.x0[level] = l_start ;
                        l_son_grid.dx0[level] = slope_[level];
                        l_son_grid.x1[level] = l_end ;
                        l_son_grid.dx1[level] = -slope_[level];
                        push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
                        l_son_grid.x0[level] = l_start;
                        l_son_grid.dx0[level] = l_father_grid.dx0[level];
                        l_son_grid.x1[level] = l_start ;
                        l_son_grid.dx1[level] = slope_[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                        l_son_grid.x0[level] = l_end ;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level];
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                    } /* end if (cut_lb) */
                    else { /* cut_tb */
                        if (lb == phys_length_[level] && l_father_grid.dx0[level] == 0 && l_father_grid.dx1[level] == 0) { /* initial cut on the dimension */
                            assert(l_father_grid.dx0[level] == 0);
                            assert(l_father_grid.dx1[level] == 0);
                            const int mid = tb/2;
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);
                            const int ul_start = (l_father_grid.x0[level] + l_father_grid.dx0[level] * lt);
                            l_son_grid.x0[level] = l_start ;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end ;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            l_son_grid.x0[level] = l_end ;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = l_end ;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                        } else { /* NOT the initial cut! */
                            const int mid = tb/2;
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);
							const int dx0_h = l_father_grid.dx0[level] * lt;
                            
                            l_son_grid.x0[level] = l_start;
                            l_son_grid.dx0[level] = l_father_grid.dx0[level];
                            l_son_grid.x1[level] = l_start + 2 * dx0_h ;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

							const int dx1_h = l_father_grid.dx1[level] * lt;
                            l_son_grid.x0[level] = l_end + 2 * dx1_h ;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end;
                            l_son_grid.dx1[level] = l_father_grid.dx1[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            l_son_grid.x0[level] = l_start + 2 * dx0_h ;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = l_end + 2 * dx1_h ;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                        }                    
                    } /* end if (cut_tb) */
                } /* end if (can_cut) */
            } /* end if (performing a space cut) */
        } /* end while (queue_len_[curr_dep] > 0) */
        assert(queue_len_[curr_dep_pointer] == 0);
    } /* end for (curr_dep < N_RANK+1) */
}
#endif

/* This is the version for interior region cut! */
template <int N_RANK>
inline void Algorithm<N_RANK>::space_time_cut_interior_projections(int t0, 
												int t1, 
												grid_info<N_RANK> const grid,
												const int index)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;

	print_projection (t0, t1, grid, index) ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
        bool cut_lb = (lb < tb);
        thres = (slope_[i] * lt);
        sim_can_cut = sim_can_cut || (cut_lb ? (lb >= 2 * thres) : (tb >= 2 * thres)) ; 
		sim_can_cut = sim_can_cut && (lt > 1) ;
    }

    if (sim_can_cut) {
        /* cut into space */
        space_cut_interior_projections(t0, t1, grid, index);
        return;
    } else if (lt > 1) {
        /* cut into time */
        int halflt = lt / 2;
        l_son_grid = grid;
        space_time_cut_interior_projections(t0, t0+halflt, l_son_grid, index);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        space_time_cut_interior_projections(t0+halflt, t1, l_son_grid, index);
        return;
    } else {
        // base case
		//find homogeneity in the base case
        return;
    }  
}

/* This is the version for boundary region cut! */
template <int N_RANK> 
inline void Algorithm<N_RANK>::space_time_cut_boundary_projections(int t0, 
												int t1, 
												grid_info<N_RANK> const grid,
												const int index) 
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;

	print_projection (t0, t1, grid, index) ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
        thres = (slope_[i] * lt);
        bool cut_lb = (lb < tb);
        sim_can_cut = sim_can_cut || cut_lb ? (lb >= 2 * thres) : (tb >= 2 * thres) ;
		sim_can_cut = sim_can_cut && (lt > 1) ;
        call_boundary |= l_touch_boundary;
    }


    if (sim_can_cut) {
        /* cut into space */
        if (call_boundary) 
            space_cut_boundary_projections(t0, t1, l_father_grid, index);
        else
            space_cut_interior_projections(t0, t1, l_father_grid, index);
        return;
    } 

    if (lt > 1) {
        /* cut into time */
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
        if (call_boundary) {
            space_time_cut_boundary_projections(t0, t0+halflt, l_son_grid, index);
        } else {
            space_time_cut_interior_projections(t0, t0+halflt, l_son_grid, index);
        }
		//find projection of the bottom zoid and top zoid

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            space_time_cut_boundary_projections(t0+halflt, t1, l_son_grid, index);
        } else {
            space_time_cut_interior_projections(t0+halflt, t1, l_son_grid, index);
        }
        return;
    } 

        // base case
		//find homogeneity in the base case
		//call the function that traverses the projection in the base case
		//and finds the homogeneity of the base case zoid.
        return;
}

template <int N_RANK> 
inline void Algorithm<N_RANK>::print_projection(int t0, int t1, 
								grid_info<N_RANK> const & grid, 
								const int index)
{
}

template <> 
inline void Algorithm<1>::print_projection(int t0, int t1, 
								grid_info<1> const & grid,
								const int index)
{
	int dt = t1 - t0 ;
	assert (dt >= 1) ;
	unsigned long qx1 = grid.x0 [0] ;
	unsigned long qx2 = grid.x1 [0] ;
	unsigned long qx3 = grid.x0 [0] + grid.dx0 [0] * dt ; //(dt - 1) ;
	unsigned long qx4 = grid.x1 [0] + grid.dx1 [0] * dt ; //(dt - 1) ;

	assert (qx1 <= qx2) ;
	assert (qx3 <= qx4) ;
	
	if (qx1 == qx3 && qx2 == qx4)
	{
		//rectangle
		return ;
	}
	bool found = false ;
	int x, x2 ;
	//unsigned long begin_index ;
	unsigned long unwrapped_index ;
	unsigned long mid_point ;
	//find the bigger base
	if (qx2 - qx1 > qx4 - qx3)
	{
		x = qx2 - qx1 ;
		x2 = qx4 - qx3 ;
		unwrapped_index = qx1 ;
		//begin_index = pmod(qx1, phys_length_ [0]) ;

		mid_point = qx1 + ((qx2 - qx1) >> 1) ;
		mid_point = pmod(mid_point, phys_length_ [0]) ;
	}
	else
	{
		x = qx4 - qx3 ;
		x2 = qx2 - qx1 ;
		unwrapped_index = qx3 ;
		//begin_index = pmod(qx3, phys_length_ [0]) ;

		mid_point = qx3 + ((qx4 - qx3) >> 1) ;
		mid_point = pmod(mid_point, phys_length_ [0]) ;
	}
	if (x < 0)
	{
		cout << "error : projection length is less than zero " << endl ;
		return ;
	}
	if (x == 0)
	{
		return ;
	}
	//cout << begin_index << " " << x << endl ;
	//set<int> & s = m_1d [index] [begin_index] ;
	set<int> & s = m_1d [index] [mid_point] ;
	set<int>::iterator pos = s.lower_bound(x) ;
	if (*pos != x)
	{
		/*cout << " t0 " << t0  << " t1 " << t1  << " x0 " << qx1 
		<< " x1 " << qx2 << " x2 " << qx3 << " x3 " << qx4
		<< " dx0 " << grid.dx0 [0] << " dx1 " << grid.dx1 [0]
		<< endl  ; 
		if (mid_point >= phys_length_ [0])
		{
			cout << "mid_point " << mid_point << endl ;
		}*/
		s.insert(pos, x) ;
		num_projections_1d [index]++ ;
		m_1d_proj_length [x]++ ;
		m_1d_index_by_length [x].insert(unwrapped_index);
		//map <int, int> & m = m_1d_map [begin_index] ;
		map <int, int> & m = m_1d_map [mid_point] ;
		m.insert(std::make_pair(x, t1 - t0)) ; //insert(length, height) pair
		/*if (print_projections)
		{
			cout << "new proj [" << unwrapped_index << ", " << unwrapped_index + x << "]" << endl ; 
		}*/
	}
//	else
//	{
//		map <int, int> & m = m_1d_map [begin_index] ;
//		map <int, int>::iterator pos = m.find(x) ;
//		if (pos == m.end())
//		{
//			cout << "Error : Projection not found " << endl ; 
//		}
//		if (pos->second == t1 - t0)  //duplicate projection with same height
//		{
//			cout << "duplicate proj [" << unwrapped_index << ", " << unwrapped_index + x << "]" << endl ; 
//			cout << " t0 " << t0  << " t1 " << t1  << " x0 " << qx1 
//			<< " x1 " << qx2 << " x2 " << qx3 << " x3 " << qx4
//			<< endl  ; 
//		}
//	}
}

template <> 
inline void Algorithm<2>::print_projection(int t0, int t1, 
								grid_info<2> const & grid,
								const int index)
{
	int dt = t1 - t0 ;
		
	int x2 = grid.x0 [1] + grid.dx0 [1] * dt ; //(dt - 1) ;
	int y2 = grid.x0 [0] + grid.dx0 [0] * dt ; //(dt - 1) ;
	/*{
		cout << " (t0, t1) (" << t0  - 1 << "," << t1 - 1 << ") (x0, x1) " << 
		"(" << grid.x0 [1] << "," << grid.x1 [1] << ") " << " (x2, x3) " << 
		"(" << grid.x0 [1] + grid.dx0 [1] * (dt - 1) << "," << 
			grid.x1 [1] + grid.dx1 [1] * (dt - 1) << ") " << " (y0, y1) " <<
		"(" << grid.x0 [0] << "," << grid.x1 [0] << ")  " << " (y2, y3) " << 
		"(" << grid.x0 [0] + grid.dx0 [0] * (dt - 1) << "," << 
			grid.x1 [0] + grid.dx1 [0] * (dt - 1) << ") " << endl ; 
	}*/
	if (grid.x0 [1] <= x2 && grid.x0 [0] <= y2)
	{
		//bottom rectangle is the projection
		int x0 = pmod (grid.x0 [1], phys_length_ [1]) ;
		int y0 = pmod (grid.x0 [0], phys_length_ [0]) ;
		//choose (x0, y0) as the reference point
		int ref_point = y0 * phys_length_ [1] + x0 ;

		int x1 = pmod (grid.x1 [1], phys_length_ [1]) ;
		int y1 = pmod (grid.x1 [0], phys_length_ [0]) ;
		int corner = y1 * phys_length_ [1] + x1 ;
		if (ref_point == corner)
		{
			return ; //projection is of zero area
		}
		set <int> & s = m_2d_r [index] [ref_point] ;
		set <int>::iterator pos = s.lower_bound(corner) ;
		if (*pos != corner)
		{
			s.insert(pos, corner) ;
			num_projections_2d_r [index]++ ;
		}
	}
	else if (grid.x0 [1] > x2 && grid.x0 [0] > y2)
	{
		//top rectangle is the projection
		//choose (x2, y2) as the reference point
		int ref_point = pmod (y2, phys_length_ [0]) * phys_length_ [1] + 
						pmod (x2, phys_length_ [1]) ;
		int x3 = grid.x1 [1] + grid.dx1 [1] * dt ; // (dt - 1) ; 
		int y3 = grid.x1 [0] + grid.dx1 [0] * dt ; //(dt - 1) ;
		int corner = pmod (y3, phys_length_ [0]) * phys_length_ [1] + 
					pmod (x3, phys_length_ [1]) ;
		if (ref_point == corner)
		{
			return ; //projection is of zero area
		}
		set <int> & s = m_2d_r [index] [ref_point] ;
		set <int>::iterator pos = s.lower_bound(corner) ;
		if (*pos != corner)
		{
			s.insert(pos, corner) ;
			num_projections_2d_r [index]++ ;
		}
	}
	else
	{
		//projection is an octagon
		int ref_point ;
		proj_2d p ;
		p.x [0] = pmod (grid.x0 [1], phys_length_ [1]) ;
		p.y [0] = pmod (grid.x0 [0], phys_length_ [0]) ;
		p.x [1] = pmod (grid.x1 [1], phys_length_ [1]) ;
		p.y [1] = pmod (grid.x1 [0], phys_length_ [0]) ;
		p.x [3] = grid.x1 [1] + grid.dx1 [1] * dt ; // (dt - 1) ; 
		p.y [3] = grid.x1 [0] + grid.dx1 [0] * dt ; // (dt - 1) ;
		p.x [3] = pmod (p.x [3], phys_length_ [1]) ;
		p.y [3] = pmod (p.y [3], phys_length_ [0]) ;
		if (grid.x0 [1] <= x2)
		{
			//x dimension converges
			p.y [2] = pmod (y2, phys_length_ [0]) ;
			p.x [2] = pmod (x2, phys_length_ [1]) ;
			ref_point = p.y [0] * phys_length_ [1] + p.x [0] ;
			p.octagon_type = 0 ;
		}
		else
		{
			//x dimension diverges
			p.y [2] = pmod (y2, phys_length_ [0]) ;
			p.x [2] = pmod (x2, phys_length_ [1]) ;
			ref_point = p.y [2] * phys_length_ [1] + p.x [2] ;
			p.octagon_type = 1 ;
		}
		bool found = false ;
		vector<proj_2d > & vec = (m_2d_o [index]) [ref_point] ;
		for (int i = 0 ; i < vec.size() ; i++)
		{
			proj_2d & p2 = vec [i] ;
			if (p.octagon_type == p2.octagon_type)
			{
				if (p2.x [0] == p.x [0] && p2.x [1] == p.x [1] &&
					p2.x [2] == p.x [2] && p2.x [3] == p.x [3] &&
					p2.y [0] == p.y [0] && p2.y [1] == p.y [1] && 
					p2.y [2] == p.y [2] && p2.y [3] == p.y [3])
				{
					found = true ;
				}
			}
			else 
			{
				if (p2.x [2] == p.x [0] && p2.x [3] == p.x [1] &&
					p2.x [0] == p.x [2] && p2.x [1] == p.x [3] &&
					p2.y [2] == p.y [0] && p2.y [3] == p.y [1] && 
					p2.y [0] == p.y [2] && p2.y [1] == p.y [3])
				{
					found = true ;
				}
			}
			if (found)
			{
				break ;
			}
		}
		if (! found)
		{
			vec.push_back(p);
			num_projections_2d_o [index]++ ;
		}
	}
}

#ifdef DEFAULT_TIME_CUT
template <int N_RANK> 
inline void Algorithm<N_RANK>::compute_projections(int t0, int t1, 
				grid_info<N_RANK> const grid) 
{
	int T = t1 - t0 ;
	cout << "to " << t0 << " t1 " << t1 << endl ;
	int W = 0 ;  //max_width among all dimensions
	int slope ;
	for (int i = 0 ; i < N_RANK ; i++)
	{
		cout << "dim " << i << " length " << phys_length_ [i] << endl ;
		if (phys_length_ [i] > W)
		{
			W = phys_length_ [i] ;
			slope = slope_ [i] ;
		}		
	}
	if (N_RANK == 1)
	{
		m_1d [0].reserve (phys_length_ [0]) ;
		m_1d [0].resize (phys_length_ [0]) ;
		m_1d [1].reserve (phys_length_ [0]) ;
		m_1d [1].resize (phys_length_ [0]) ;

		m_1d_proj_length.resize (phys_length_ [0] + 1, 0) ;
		m_1d_index_by_length.resize (phys_length_ [0] + 1) ; 
		m_1d_map.resize (phys_length_ [0]) ;
		num_projections_1d [0] = 0 ;
		num_projections_1d [1] = 0 ;
	}
	if (N_RANK == 2)
	{
		m_2d_r [0].reserve (phys_length_ [0] * phys_length_ [1]) ;
		m_2d_r [0].resize (phys_length_ [0] * phys_length_ [1]) ;
		m_2d_o [0].reserve (phys_length_ [0] * phys_length_ [1]) ;
		m_2d_o [0].resize (phys_length_ [0] * phys_length_ [1]) ;

		m_2d_r [1].reserve (phys_length_ [0] * phys_length_ [1]) ;
		m_2d_r [1].resize (phys_length_ [0] * phys_length_ [1]) ;
		m_2d_o [1].reserve (phys_length_ [0] * phys_length_ [1]) ;
		m_2d_o [1].resize (phys_length_ [0] * phys_length_ [1]) ;

		num_projections_2d_r [0] = 0 ;
		num_projections_2d_r [1] = 0 ;
		num_projections_2d_o [0] = 0 ;
		num_projections_2d_o [1] = 0 ;
	}
	cout << "time shift " << time_shift << endl ;
	//Time cuts are needed if W < 2 * slope * T
	if (W < 2 * slope * T)
	{
		//compute 2^k where k = ceil(lg (2 * slope * T / W))
		//int k = 8 * sizeof(int) - __builtin_clz((T - 1) * 2 * slope / W) ; 
		int k = 8 * sizeof(int) - __builtin_clz((T * 2 * slope - 1) / W) ; 
		cout << "k " << k << endl ;
		int two_to_the_k = 1 << k ; 
		cout << "width " << W << " T " << T <<  " 2^k "  << two_to_the_k << endl ;
		cout << "slope " << slope << endl ;
		//h1 = floor (T/(2^k))
		int h1 = T / two_to_the_k ;
		if (W >= 2 * slope * T)
		{
			h1 = T ;	
		}
		assert (W >= 2 * slope * h1) ; 
		/*if (N_RANK == 1)
		{
			//simulate the upright triangle and inverted trapezoid separately
			int mid = W / 2 ;
			int dx = slope_ [0] * h1 ;
			grid_info<N_RANK> grid2 ;
			grid2.x0[0] = 0 ;
			grid2.dx0[0] = slope_[0];
			grid2.x1[0] = 2 * dx ;
			grid2.dx1[0] = -slope_[0];
			//space_time_cut_boundary_projections(t0, t0 + h1, grid2, 0) ;*/
			/*
			grid2.x0[0] = mid - dx ;
			grid2.dx0[0] = slope_[0];
			grid2.x1[0] = mid + dx ;
			grid2.dx1[0] = -slope_[0];
			space_time_cut_boundary_projections(t0, t0 + h1, grid2, 0) ;
			cout << " mid - dx " << mid - dx << " mid + dx " << mid + dx << endl;
			cout << "h1 " << h1 << endl ;	
			cout << "mid " << mid << " dx " << dx << " mid + dx " << mid+dx <<
					" W + mid - dx " << W + mid - dx << endl ;
			grid2.x0[0] = mid + dx ;
			grid2.dx0[0] = -slope_[0];
			grid2.x1[0] = W + mid - dx ;
			grid2.dx1[0] = slope_[0];*/
			/*grid2.x0[0] = dx ;
			grid2.dx0[0] = -slope_[0];
			grid2.x1[0] = W - dx ;
			grid2.dx1[0] = slope_[0];*/
			/*grid2.x0[0] = 0 ;
			grid2.dx0[0] = slope_[0];
			grid2.x1[0] = W ;
			grid2.dx1[0] = -slope_[0];
			//cout << "dx " << dx << " W - dx " << W - dx << endl ;
			space_time_cut_boundary_projections(t0, t0 + h1, grid2, 0) ;
		}*/
		space_time_cut_boundary_projections(t0, time_shift + h1, grid, 0) ;
		int h2 = (T + two_to_the_k - 1) / two_to_the_k ;
		if (h2 != h1)
		{
			//you can have a case where W >= 2 * slope * h1 but
			//W < 2 * slope * h2 
			//Example W = 33, T = 33, slope = 1 yields h1 = 16 and h2 = 17
			if (W < 2 * slope * h2)
			{
				cout << " h1 " << h1 << " h2 " << h2 / 2 
					 << " h3 " << (h2 + 1) / 2 << endl ;
				space_time_cut_boundary_projections(t0, time_shift + h2 / 2, grid, 0) ;
				space_time_cut_boundary_projections(t0, time_shift + (h2 + 1) / 2, grid, 0);
			}
			else
			{
				cout << " h1 " << h1 << " h2 " << h2 << endl ;
				space_time_cut_boundary_projections(t0, time_shift + h2, grid, 0) ;
			}
		}
		else
		{
			cout << " h1 " << h1 << " h2 " << h2 << endl ;
		}
	}
	else
	{
		space_time_cut_boundary_projections(t0, t1, grid, 0) ;
	}
}
#else
template <int N_RANK> 
inline void Algorithm<N_RANK>::compute_projections(int t0, int t1, 
				grid_info<N_RANK> const grid) 
{
	int T = t1 - t0 ;
	cout << "to " << t0 << " t1 " << t1 << endl ;
	int W = 0 ;  //max_width among all dimensions
	int slope ;
	for (int i = 0 ; i < N_RANK ; i++)
	{
		cout << "dim " << i << " length " << phys_length_ [i] << 
				" x0 [" << i << "] " << grid.x0[i] <<
				" x1 [" << i << "] " << grid.x1[i] <<
				" dx0 [" << i << "] " << grid.dx0[i] <<
				" dx1 [" << i << "] " << grid.dx1[i] <<
				endl ;
		if (phys_length_ [i] > W)
		{
			W = phys_length_ [i] ;
			slope = slope_ [i] ;
		}		
	}
	if (N_RANK == 1)
	{
		m_1d [0].reserve (phys_length_ [0]) ;
		m_1d [0].resize (phys_length_ [0]) ;
		m_1d [1].reserve (phys_length_ [0]) ;
		m_1d [1].resize (phys_length_ [0]) ;

		m_1d_index_by_length.resize (phys_length_ [0] + 1) ; 
		m_1d_proj_length.resize (phys_length_ [0] + 1, 0) ;
		m_1d_map.resize (phys_length_ [0]) ;
		num_projections_1d [0] = 0 ;
		num_projections_1d [1] = 0 ;
	}
	if (N_RANK == 2)
	{
		m_2d_r [0].reserve (phys_length_ [0] * phys_length_ [1]) ;
		m_2d_r [0].resize (phys_length_ [0] * phys_length_ [1]) ;
		m_2d_o [0].reserve (phys_length_ [0] * phys_length_ [1]) ;
		m_2d_o [0].resize (phys_length_ [0] * phys_length_ [1]) ;

		m_2d_r [1].reserve (phys_length_ [0] * phys_length_ [1]) ;
		m_2d_r [1].resize (phys_length_ [0] * phys_length_ [1]) ;
		m_2d_o [1].reserve (phys_length_ [0] * phys_length_ [1]) ;
		m_2d_o [1].resize (phys_length_ [0] * phys_length_ [1]) ;

		num_projections_2d_r [0] = 0 ;
		num_projections_2d_r [1] = 0 ;
		num_projections_2d_o [0] = 0 ;
		num_projections_2d_o [1] = 0 ;
	}
	cout << "time shift " << time_shift << endl ;
	//find index of most significant bit that is set
	int Wn = W / (2 * slope) ;
	int index_msb = 8 * sizeof(int) - __builtin_clz(Wn) - 1 ;
	int h1 = 1 << index_msb ;
	cout << "t0 " << t0 << " t1 " << t1 << " Wn " << Wn << " 2^floor(lg(Wn)) "
				 << h1 << endl ;
	if (W > 2 * slope * T)
	{
		h1 = T ;
	}
	cout << "t0 " << t0 << " t1 " << t1 << 
			" h1 " << h1 << " t0 + h1 " <<
			t0 + h1 << endl ;
	space_time_cut_boundary_projections(t0, t0 + h1, grid, 0) ;
	int r = T / h1 * h1 ;
	t0 += r ;
	int h2 = t1 - t0 ;
	cout << "t0 " << t0 << " t1 " << t1 << " h2 " << h2 << 
		" t0 + h2 " << t0 + h2 << endl << endl ;
	while (h2 > 0)
	{
		//find index of most significant bit that is set
		index_msb = 8 * sizeof(int) - __builtin_clz(h2) - 1 ;
		int h = 1 << index_msb ;
		cout << "t0 " << t0 << " t1 " << t1 << 
			" h " << h << " t0 + h " <<
			t0 + h << endl ;
		space_time_cut_boundary_projections(t0, t0 + h, grid, 0) ;
		t0 += h ;
		h2 = t1 - t0 ;
	}
}
#endif

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
		/*cout << "t0 " << t0 << " t1 " << t1 << 
			" h1 " << h1 << " t0 + h1 " <<
			t0 + h1 << endl ;*/
		space_time_cut_boundary(t0, t0 + h1, grid, f, bf) ;
		t0 += h1 ;
	}

	int h2 = t1 - t0 ;
	//time cuts happen only if height > dt_recursive_
	while (h2 > dt_recursive_)
	{
		//find index of most significant bit that is set
		index_msb = (sizeof(int) << 3) - __builtin_clz(h2) - 1 ;
		int h = 1 << index_msb ;
		/*cout << "t0 " << t0 << " t1 " << t1 << 
			" h " << h << " t0 + h " <<
			t0 + h << endl ;*/
		space_time_cut_boundary(t0, t0 + h, grid, f, bf) ;
		t0 += h ;
		h2 = t1 - t0 ;
	}
	while (h2 > 1)
	{
		//find index of most significant bit that is set
		index_msb = (sizeof(int) << 3) - __builtin_clz(h2) - 1 ;
		int h = 1 << index_msb ;
		/*cout << "t0 " << t0 << " t1 " << t1 << 
			" h " << h << " t0 + h " <<
			t0 + h << endl ;*/
		bool abnormal = false ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			num_triangles [i] = dx_recursive_ [i] / ((slope_ [i] * h) << 1) ;
		}
		abnormal_region_space_time_cut_boundary(t0, t0 + h, grid, f, bf) ;
		t0 += h ;
		h2 = t1 - t0 ;
	}
	if (h2 == 1)
	{
		cout << "h = 1 t0 " << t0 << " t1 " << t1 << 
			 " t0 + h " << t0 + h2 << endl ;
		//base_case_kernel_boundary(t0, t0 + h2, grid, bf);
		shorter_duo_sim_obase_bicut_p(t0, t0 + h2, grid, f, bf) ;
	}
}
#endif
