#ifndef PROJECTIONS_HPP
#define PROJECTIONS_HPP

#include "pochoir_common.hpp"
#include "pochoir_walk.hpp"

#define initial_cut(i) (lb[i] == phys_length_[i])
/* grid.x1[i] >= phys_grid_.x1[i] - stride_[i] - slope_[i] 
 * because we compute the kernel with range [a, b)
 */
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

/* ************************************************************************************** */
/* following are the procedures for obase with duality , always cutting based on shorter bar
 */
template <int N_RANK>
inline void Algorithm<N_RANK>::space_cut_interior(int t0, int t1, 
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
                    space_time_cut_interior(l_son->t0, l_son->t1, l_son->grid,
											index) ;
                } 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                /* use cilk_spawn to spawn all the sub-grid */
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0)
                    space_time_cut_interior(l_father->t0, l_father->t1, 
											l_father->grid, index) ;
                else
                    space_time_cut_interior(l_father->t0, l_father->t1, 
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
                //const bool can_cut = cut_lb ? (lb >= 2 * thres && lb > dx_recursive_[level]) : (tb >= 2 * thres && lb > dx_recursive_[level]);
                //const bool can_cut = cut_lb ? (lb >= 2 * thres && lb > 1) : (tb >= 2 * thres && tb > 1);
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

/* This is for boundary region space cut! , always cutting based on the shorter bar
 */
template <int N_RANK> 
inline void Algorithm<N_RANK>::space_cut_boundary(int t0, int t1, 
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
                    space_time_cut_boundary(l_son->t0, l_son->t1, l_son->grid,
											index);
                } 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                /* use cilk_spawn to spawn all the sub-grid */
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0) {
                    space_time_cut_boundary(l_father->t0, l_father->t1, 
											l_father->grid, index);
                } else {
                    space_time_cut_boundary(l_father->t0, l_father->t1, 
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
                //const bool can_cut = cut_lb ? (l_touch_boundary ? (lb >= 2 * thres && lb > dx_recursive_boundary_[level]) : (lb >= 2 * thres && lb > dx_recursive_[level])) : (l_touch_boundary ? (tb >= 2 * thres && lb > dx_recursive_boundary_[level]) : (tb >= 2 * thres && lb > dx_recursive_[level]));
        		//const bool can_cut = cut_lb ? (lb >= 2 * thres && lb > 1) : (tb >= 2 * thres && tb > 1) ;
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
                            assert(l_father_grid.dx0[level] == 0);
                            assert(l_father_grid.dx1[level] == 0);
                            const int mid = tb/2;
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);
                            const int ul_start = (l_father_grid.x0[level] + l_father_grid.dx0[level] * lt);
                            /* merge the big black trapezoids */
                            l_son_grid.x0[level] = ul_start + mid;
                            l_son_grid.dx0[level] = slope_[level];
                            l_son_grid.x1[level] = l_end + (ul_start - l_start) + mid;
                            l_son_grid.dx1[level] = -slope_[level];
                            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            /* cilk_sync */
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
                            /* push middle minizoid into circular queue of (curr_dep + 1)*/
                            l_son_grid.x0[level] = ul_start + mid;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = ul_start + mid;
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


/* This is the version for interior region cut! */
template <int N_RANK>
inline void Algorithm<N_RANK>::space_time_cut_interior(int t0, int t1, 
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
		//cout << "dx_recursive " << dx_recursive_[i] << "dx_recursive_boundary_ " << dx_recursive_boundary_[i] << endl ;
		//cout << "lb " << lb <<  " tb " << tb << " lt " << lt << endl ; 
        //sim_can_cut = sim_can_cut || (cut_lb ? (lb >= 2 * thres & lb > dx_recursive_[i]) : (tb >= 2 * thres & tb > dx_recursive_[i]));
        //sim_can_cut = sim_can_cut || (cut_lb ? (lb >= 2 * thres && lb > dx_recursive_[i]) : (tb >= 2 * thres && tb > dx_recursive_[i]));
        //sim_can_cut = sim_can_cut || (cut_lb ? (lb >= 2 * thres && lb > 1) : (tb >= 2 * thres && tb > 1)) ; 
        sim_can_cut = sim_can_cut || (cut_lb ? (lb >= 2 * thres) : (tb >= 2 * thres)) ; 
		sim_can_cut = sim_can_cut && (lt > 1) ;
        /* as long as there's one dimension can conduct a cut, we conduct a 
         * multi-dimensional cut!
         */
    }

    if (sim_can_cut) {
        /* cut into space */
        space_cut_interior(t0, t1, grid, index);
        return;
    // } else if (lt > dt_recursive_ && l_total_points > Z) {
    //} else if (lt > dt_recursive_) {
    } else if (lt > 1) {
        /* cut into time */
//        assert(dt_recursive_ >= r_t);
        //assert(lt > dt_recursive_);
        int halflt = lt / 2;
        //int halflt = (lt + 1) / 2 ;
        l_son_grid = grid;
        space_time_cut_interior(t0, t0+halflt, l_son_grid, index);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        space_time_cut_interior(t0+halflt, t1, l_son_grid, index);
        return;
    } else {
        // base case
		//find homogeneity in the base case
        //f(t0, t1, grid);
        return;
    }  
}

/* This is the version for boundary region cut! */
template <int N_RANK> 
inline void Algorithm<N_RANK>::space_time_cut_boundary(int t0, int t1, 
												grid_info<N_RANK> const grid,
												const int index) 
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;

	/*
		Assume lt = Omega(max(x_i)) where x_i is the width of the zoid in 
		dimension i. Then we keep cutting time in half until
		the condition max(x_i) >= 2 * slope [i] * lt' is satisfied.
		At this point, the n-dimensional cuboid satisfies the
		condition 2 * slope [i] * 2 * lt' > max(x_i) >= 2 * slope [i] * lt' 
		Besides, the following conditions hold from this point onwards.
	 */
	//If two zoids have the same height and same projection
	//at the bottom time step and top time step, then all
	//subzoids cut from the two zoids will be identical.
	//Above statement is a strong condition.

	/*If two zoids have the same projection
		if their heights are the same, all subzoids cut from the two zoids
		will be identical.
		if their heights are not the same, then the smaller zoid is a subzoid
		of the bigger zoid.
	 */

	//add code to count projections here
	//find if the projection of the zoid exists already.
	//If yes, we can skip traversing the zoid completely.

	print_projection (t0, t1, grid, index) ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
        thres = (slope_[i] * lt);
        /* l_father_grid may be mapped to a new region in touch_boundary() */
        /* for the initial cut, we exclude the begining and end point to minimize
         * the overhead on boundary
        */
        /* lb == phys_length_[i] indicates an initial cut! */
        bool cut_lb = (lb < tb);
		//cout << "dx_recursive " << dx_recursive_[i] << "dx_recursive_boundary_ " << dx_recursive_boundary_[i] << endl ;
		//cout << "boundary " << "lb " << lb <<  " tb " << tb << " lt " << lt << endl ; 
        //sim_can_cut = sim_can_cut || (cut_lb ? (l_touch_boundary ? (lb >= 2 * thres & lb > dx_recursive_boundary_[i]) : (lb >= 2 * thres & lb > dx_recursive_[i])) : (l_touch_boundary ? (tb >= 2 * thres & tb > dx_recursive_boundary_[i]) : (tb > 2 * thres & tb > dx_recursive_[i])));
        //sim_can_cut = sim_can_cut || (cut_lb ? (l_touch_boundary ? (lb >= 2 * thres && lb > dx_recursive_boundary_[i]) : (lb >= 2 * thres && lb > dx_recursive_[i])) : (l_touch_boundary ? (tb >= 2 * thres && tb > dx_recursive_boundary_[i]) : (tb > 2 * thres && tb > dx_recursive_[i])));
        //sim_can_cut = sim_can_cut || cut_lb ? (lb >= 2 * thres && lb > 1) : (tb >= 2 * thres && tb > 1) ;
        sim_can_cut = sim_can_cut || cut_lb ? (lb >= 2 * thres) : (tb >= 2 * thres) ;
		sim_can_cut = sim_can_cut && (lt > 1) ;
        call_boundary |= l_touch_boundary;
    }


    if (sim_can_cut) {
        /* cut into space */
        /* push the first l_father_grid that can be cut into the circular queue */
        /* boundary cuts! */
        if (call_boundary) 
            space_cut_boundary(t0, t1, l_father_grid, index);
        else
            space_cut_interior(t0, t1, l_father_grid, index);
        return;
    } 

    /*if (call_boundary)
        l_dt_stop = dt_recursive_boundary_;
    else
        l_dt_stop = dt_recursive_;
	*/
    if (lt > 1) {
        /* cut into time */
        int halflt = lt / 2;
        //int halflt = (lt + 1) / 2 ;
        l_son_grid = l_father_grid;
        if (call_boundary) {
            space_time_cut_boundary(t0, t0+halflt, l_son_grid, index);
        } else {
            space_time_cut_interior(t0, t0+halflt, l_son_grid, index);
        }
		//find projection of the bottom zoid and top zoid

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            space_time_cut_boundary(t0+halflt, t1, l_son_grid, index);
        } else {
            space_time_cut_interior(t0+halflt, t1, l_son_grid, index);
        }
        return;
    } 

        // base case
		//find homogeneity in the base case
		//call the function that traverses the projection in the base case
		//and finds the homogeneity of the base case zoid.
        /*if (call_boundary) {
            base_case_kernel_boundary(t0, t1, l_father_grid, bf);
        } else {
            f(t0, t1, l_father_grid);
        }*/
        return;
}

template <int N_RANK> 
inline void Algorithm<N_RANK>::print_projection(int t0, int t1, 
								grid_info<N_RANK> const & grid, 
								const int index)
{
}

/*
template <> 
inline void Algorithm<1>::print_projection(int t0, int t1, 
								grid_info<1> const & grid)
{
	int dt = t1 - t0 ;
	unsigned long qx1 = grid.x0 [0] ;
	unsigned long qx2 = grid.x1 [0] ;
	unsigned long qx3 = grid.x0 [0] + grid.dx0 [0] * dt ;
	unsigned long qx4 = grid.x1 [0] + grid.dx1 [0] * dt ;

	//cout << " t0 " << t0  - 1 << " t1 " << t1 - 1 << " x0 " << qx1 
	//	<< " x1 " << qx2 << " x2 " << qx3 << " x3 " << qx4
	//	<< endl  ; 

	assert (qx1 <= qx2) ;
	assert (qx3 <= qx4) ;

	multimap<unsigned long, unsigned long>::iterator pos ;
	bool found = false ;
	//find the bigger base
	if (qx2 - qx1 > qx4 - qx3)
	{
		int x = qx2 - qx1 ;
		if (x <= 0)
		{
			return ;
		}
		pos = m_1d.lower_bound(x) ;
		for ( ; pos != m_1d.upper_bound(x) ; pos++)
		{
			assert (pos->first == x) ;
			if (pos->second == (qx1 % phys_length_[0]))
			{
				found = true ;
				break ;
			}
		}
		if (! found)
		{
			//cout << "inserting " << qx2 - qx1 << "," << 
			//	qx1 % phys_length_[0] << endl ;
			m_1d.insert(make_pair(x, qx1 % phys_length_[0]));
		}
	}
	else
	{
		int x = qx4 - qx3 ;
		if (x <= 0)
		{
			return ;
		}
		pos = m_1d.lower_bound(x) ;
		for ( ; pos != m_1d.upper_bound(x) ; pos++)
		{
			assert (pos->first == x) ;
			if (pos->second == qx3 % (phys_length_[0]))
			{
				found = true ;
				break ;
			}
		}
		if (! found)
		{
			//cout << "inserting " << x << ", " << qx3 % phys_length_[0] 
			//	<< endl ;
			m_1d.insert(make_pair(x, qx3 % phys_length_[0])) ;
		}
	}
}
*/

template <> 
inline void Algorithm<1>::print_projection(int t0, int t1, 
								grid_info<1> const & grid,
								const int index)
{
	int dt = t1 - t0 ;
	assert (dt >= 1) ;
	unsigned long qx1 = grid.x0 [0] ;
	unsigned long qx2 = grid.x1 [0] ;
	unsigned long qx3 = grid.x0 [0] + grid.dx0 [0] * (dt - 1) ;
	unsigned long qx4 = grid.x1 [0] + grid.dx1 [0] * (dt - 1) ;

	/*cout << " t0 " << t0  << " t1 " << t1  << " x0 " << qx1 
		<< " x1 " << qx2 << " x2 " << qx3 << " x3 " << qx4
		<< " dx0 " << grid.dx0 [0] << " dx1 " << grid.dx1 [0]
		<< endl  ; */
	
	assert (qx1 <= qx2) ;
	assert (qx3 <= qx4) ;
	
	if (qx1 == qx3 && qx2 == qx4)
	{
		//rectangle
		return ;
	}
	bool found = false ;
	int x ;
	unsigned long begin_index ;
	//find the bigger base
	if (qx2 - qx1 > qx4 - qx3)
	{
		x = qx2 - qx1 ;
		begin_index = qx1 % phys_length_ [0] ;
	}
	else
	{
		x = qx4 - qx3 ;
		begin_index = qx3 % phys_length_ [0] ;
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
	set<unsigned long> & s = m_1d [index] [begin_index] ;
	set<unsigned long>::iterator pos = s.lower_bound(x) ;
	if (*pos != x)
	{
		s.insert(pos, x) ;
		num_projections++ ;
	}
}

template <> 
inline void Algorithm<2>::print_projection(int t0, int t1, 
								grid_info<2> const & grid,
								const int index)
{
	int dt = t1 - t0 ;
	{
	//x[0] and x[1] are the boundaries of the base of the zoid in x dimension
	//x[2] and x[3] are the boundaries of the top of the zoid in x dimension
	proj_2d p ;
	p.x [0] = grid.x0 [1] ;
	p.x [1] = grid.x1 [1] ;

	//y[0] and y[1] are the boundaries of the base of the zoid in y dimension
	//y[2] and y[3] are the boundaries of the top of the zoid in y dimension
	p.y [0] = grid.x0 [0] ;
	p.y [1] = grid.x1 [0] ;
	  
	//find the other end of the zoid
	p.x [2] = p.x [0] + grid.dx0 [1] * (dt - 1) ;
	p.x [3] = p.x [1] + grid.dx1 [1] * (dt - 1) ;
	p.y [2] = p.y [0] + grid.dx0 [0] * (dt - 1) ;
	p.y [3] = p.y [1] + grid.dx1 [0] * (dt - 1) ;

	assert (p.x [0] <= p.x [1] && p.y [0] <= p.y [1]) ;
	assert (p.x [2] <= p.x [3] && p.y [2] <= p.y [3]) ;
	
	if (p.x [0] == p. x[1] && p.x [2] == p.x [3] && p.x [1] == p.x [2] &&
		p.y [0] == p. y[1] && p.y [2] == p.y [3] && p.y [1] == p.y [2])
	{
		//rectangle with zero area
		/*cout << "Area 0 " << endl ;
		cout << " (t0, t1) (" << t0  - 1 << "," << t1 - 1 << ") (x0, x1) " << 
		"(" << p.x [0] << "," << p.x [1] << ") " << " (x2, x3) " << 
		"(" << p.x [2] << "," << p.x [3] << ") " << " (y0, y1) " <<
		"(" << p.y [0] << "," << p.y [1] << ")  " << " (y2, y3) " << 
		"(" << p.y [2] << "," << p.y [3] << ") " << endl ; */
		return ;
	}
	//update the co-ordinates with modulo width
	p.x[0] = p.x [0] % phys_length_ [1] ;
	p.x[1] = p.x [1] % phys_length_ [1] ;
	p.x[2] = p.x [2] % phys_length_ [1] ;
	p.x[3] = p.x [3] % phys_length_ [1] ;
	
	p.y[0] = p.y [0] % phys_length_ [0] ;
	p.y[1] = p.y [1] % phys_length_ [0] ;
	p.y[2] = p.y [2] % phys_length_ [0] ;
	p.y[3] = p.y [3] % phys_length_ [0] ;

	int ref_point ;
	//case1 : x dim converges
	if (p.x [0] <= p.x [2])
	{
		assert (p. x[3] <= p. x[1]) ;
		ref_point = p. y[0] * phys_length_ [1] + p. x [0] ;
		//y dim converges
		if (p.y [0] <= p.y [2]) 
		{
			assert (p.y [3] <= p.y [1]) ;
			p.type = 0 ;
			p.top_right = p. y[1] * phys_length_ [1] + p. x [1] ;
		}
		//y dim diverges
		else
		{
			assert (p.y [1] < p.y [3]) ;
			p.type = 1 ;
			p.octagon_type = 0 ;
		}
	}
	//case 2 : x dim diverges
	else 
	{
		assert (p.x [3] > p.x [1]) ;
		ref_point = p. y[2] * phys_length_ [1] + p. x [2] ;
		//y dim diverges
		if (p.y [0] >= p.y [2])
		{
			assert (p.y [3] >= p.y [1]) ;
			p.type = 0 ;
			p.top_right = p. y[3] * phys_length_ [1] + p. x [3] ;
		}
		//y dim converges
		else
		{
			assert (p.y [3] < p.y [1]) ;
			p.type = 1 ;
			p.octagon_type = 1 ;
		}
	}
	bool found = false ;

	std::vector<proj_2d > & vec = (m_2d [index]) [ref_point] ;
	for (int i = 0 ; i < vec.size() ; i++)
	{
		proj_2d & p2 = vec [i] ;
		if (p.type == p2.type) 
		{
			if (p.type == 0)
			{
				//projection is rectangle
				if (p.top_right == p2.top_right)
				{
					found = true ;
				}
			}
			else
			{
				//projection is octagon
				if (p.octagon_type == p2.octagon_type)
				{
					if (p2.x [0] == p.x [0] && p2.x [1] == p.x [1] &&
						p2.x [2] == p.x [2] && p2.x [3] == p.x[3] &&
						p2.y [0] == p.y [0] && p2.y [1] == p.y [1] && 
						p2.y [2] == p.y [2] && p2.y [3] == p.y [3])
					{
						found = true ;
					}
				}
				else 
				{
					if (p2.x [2] == p.x [0] && p2.x [3] == p.x [1] &&
						p2.x [0] == p.x [2] && p2.x [1] == p.x[3] &&
						p2.y [2] == p.y [0] && p2.y [3] == p.y [1] && 
						p2.y [0] == p.y [2] && p2.y [1] == p.y [3])
					{
						found = true ;
					}
				}
			}
		}
		if (found)
		{
			break ;
		}
	}
	
	if (! found)
	{
		//m_2d [index].insert(make_pair(ref_point, p));
		vec.push_back(p);
		num_projections++ ;
	}
	}
	
	/*
	{
	proj_2d p ;
	//choose (x0, y0) as the reference point
	int x0 = grid.x0 [1] % phys_length_ [1] ;
	int y0 = grid.x0 [0] % phys_length_ [0] ;
	int bottom_ref_point = y0 * phys_length_ [1] + x0 ;
	//choose (x2, y2) as the other reference point
	int x2 = x0 + grid.dx0 [1] * (dt - 1) ; 
	int y2 = y0 + grid.dx0 [0] * (dt - 1) ;
	int top_ref_point = y2 * phys_length_ [1] + x2 ;

	int x1 = grid.x1 [1] % phys_length_ [1] ;
	int y1 = grid.x1 [0] % phys_length_ [0] ;
	p.top_right[0] = y1 * phys_length_ [1] + x1 ;
	int x3 = x1 + grid.dx1 [1] * (dt - 1) ; 
	int y3 = y1 + grid.dx1 [0] * (dt - 1) ;
	p.top_right[1] = y3 * phys_length_ [1] + x3 ;

	set <unsigned long> & s = m_2d [index] [bottom_ref_point] ;
	set <unsigned long>::iterator pos = s.lower_bound(top_ref_point) ;
	proj_2d & p1 = *pos ;
	//if (*pos != top_ref_point)
	if (p1.top_right [0] != p.top_right [0] || 
		p1.top_right [1] != p.top_right [1])
	{
		set <unsigned long> & s1 = m_2d [index] [top_ref_point] ;
		set <unsigned long>::iterator pos1 = s.lower_bound(bottom_ref_point) ;
		proj_2d & p2 = *pos1 ;
		if (p2.top_right [0] != p.top_right[])
		//if (*pos1 != bottom_ref_point)
		{
			//projection not found
			//s.insert(pos, top_ref_point) ;
			s.insert(pos, make_pair(top_ref_point, p)) ;
			num_projections++ ;
		}
		else
		{
			//compare the top roght corners
			
		}
	}
	}*/
}

/*
template <int N_RANK> 
inline void Algorithm<N_RANK>::compute_projections(int t0, int t1, 
				grid_info<N_RANK> const grid) 
{
	int dt = t1 - t0 ;
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
	}
	if (N_RANK == 2)
	{
		m_2d [0].reserve (phys_length_ [0] * phys_length_ [1]) ;
		m_2d [0].resize (phys_length_ [0] * phys_length_ [1]) ;
	}
	num_projections = 0 ;
	cout << "time shift " << time_shift << endl ;
	//Time cuts are needed if W < 2 * slope * dt
	if (W < 2 * slope * dt)
	{
		//compue 2^k where k = ceil(lg (2 * slope * dt / W))
		int k = 8 * sizeof(int) - __builtin_clz((dt - 1) * 2 * slope / W) ; 
		cout << "k " << k << endl ;
		int two_to_the_k = 1 << k ; 
		cout << "width " << W << " dt " << dt <<  " 2^k "  << two_to_the_k << endl ;
		cout << "slope " << slope << endl ;
		//h1 = floor (dt/(2^k))
		int h1 = dt / two_to_the_k ;
		assert (W >= 2 * slope * h1) ; 
		space_time_cut_boundary(t0, time_shift + h1, grid, 0) ;
		//h2 = ceil (dt/(2^k))
		int h2 = (dt + two_to_the_k - 1) / two_to_the_k ;
		if (h2 != h1)
		{
			//you can have a case where W >= 2 * slope * h1 but
			//W < 2 * slope * h2 
			//Example W = 33, dt = 33, slope = 1 yields h1 = 16 and h2 = 17
			if (W < 2 * slope * h2)
			{
				cout << " h1 " << h1 << " h2 " << h2 / 2 
					 << " h3 " << (h2 + 1) / 2 << endl ;
				space_time_cut_boundary(t0, time_shift + h2 / 2, grid, 0) ;
				space_time_cut_boundary(t0, time_shift + (h2 + 1) / 2, grid, 0);
			}
			else
			{
				cout << " h1 " << h1 << " h2 " << h2 << endl ;
				space_time_cut_boundary(t0, time_shift + h2, grid, 0) ;
			}
		}
		else
		{
			cout << " h1 " << h1 << " h2 " << h2 << endl ;
		}
	}
	//space_time_cut_boundary(t0, t1, grid, 2) ;
}*/


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
	int Wn = W / (2 * slope) ;
	//Time cuts are needed if Wn < T
	if (Wn < T)
	{
		int h2 = T % Wn ;
		cout << "slope " << slope << endl ;
		space_time_cut_boundary(t0, time_shift + Wn, grid, 0) ;
		if (h2 > 0)
		{
			space_time_cut_boundary(t0, time_shift + h2, grid, 1) ;
		}
		cout << " Wn " << Wn << " T % Wn " << h2 << endl ;
	}
}
#endif
