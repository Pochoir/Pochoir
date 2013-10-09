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
#ifndef AUTO_TUNING_SAWZOID_HPP 
#define AUTO_TUNING_SAWZOID_HPP 

#include "auto_tuning_header.hpp"

#define dx_recursive_boundary_  (m_algo.dx_recursive_boundary_)
#define dx_recursive_ (m_algo.dx_recursive_)
#define dt_recursive_boundary_ (m_algo.dt_recursive_boundary_)
#define dt_recursive_ (m_algo.dt_recursive_)
#define slope_ m_algo.slope_
#define touch_boundary m_algo.touch_boundary
#define phys_length_ m_algo.phys_length_
#define base_case_kernel_boundary m_algo.base_case_kernel_boundary

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::symbolic_sawzoid_space_cut_interior(int t0,
			int t1, grid_info<N_RANK> const & grid, unsigned long parent_index, 
			F const & f, int * num_zoids, double & rcost, double & ncost)
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
                    symbolic_sawzoid_space_time_cut_interior(l_son->t0, 
						l_son->t1, l_son->grid, parent_index, child_index, 
						rcost, ncost, f);
					child_index++ ;
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				symbolic_sawzoid_space_time_cut_interior(l_father->t0, 
					l_father->t1, l_father->grid, parent_index, child_index, 
					rcost, ncost, f);
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
                if (num_zoids [level] == 0) {
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
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}

/* Boundary space cut. Uses sawzoid space cut.
 */
template <int N_RANK> template <typename F, typename BF> 
inline void auto_tune<N_RANK>::symbolic_sawzoid_space_cut_boundary(int t0,
		int t1, grid_info<N_RANK> const & grid, unsigned long parent_index, 
		F const & f, BF const & bf, int * num_zoids, double & rcost, 
		double & ncost)
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
                    symbolic_sawzoid_space_time_cut_boundary(l_son->t0, 
						l_son->t1, l_son->grid, parent_index, child_index, 
						rcost, ncost, f, bf);
					child_index++ ; //this can be a race.
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				symbolic_sawzoid_space_time_cut_boundary(l_father->t0, 
					l_father->t1, l_father->grid, parent_index, child_index,
					rcost, ncost, f, bf);
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
                if (num_zoids [level] == 0) {
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
							if (num_zoids [level] == 5)
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
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::symbolic_sawzoid_space_time_cut_interior(
		int t0, int t1, grid_info<N_RANK> const & grid, 
		unsigned long parent_index, int child_index, 
		double & rcost, double & ncost, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;

	int centroid = 0, width = 1 ; 
	int offset = 0, total_num_subzoids = 1 ;
	int num_subzoids [N_RANK] ;
	unsigned long key = 0 ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
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
		}
		else
		{
			short_side = tb ;
		}
		num_subzoids [i] = 0 ;
		if (short_side >= (thres << 1) && lb > dx_recursive_[i])
		{
			space_cut = true ;
			num_subzoids [i] = 3 ;
			if (short_side - (thres << 2) >= 0)
			{
				num_subzoids [i]= 5 ;
			}
			total_num_subzoids *= num_subzoids [i] ;
		}
        sim_can_cut |= space_cut ;
    }
	//cout << "centroid " << centroid << " key " << key << endl ;
	unsigned long index ;
	bool projection_exists = check_and_create_projection (key, lt, 
											centroid, index, grid) ;
	//cout << " index " << index << endl ;
	//cout << " child index " << child_index << endl ;
	zoid_type & z = m_zoids [index];
	zoid_type & parent = m_zoids [parent_index] ;
	//cout << "parent id " << parent.id << endl ;
	//cout << "child id " << z.id << endl ;
	//cout << "address of parent " << &parent << endl ;
	//add the zoid as a child of the parent
	parent.add_child(&z, child_index, index) ;
	assert (index < m_zoids.size()) ;
	
	if (projection_exists)
	{ 
		//Add the cost of the zoid to the parent.
		//rcost += 0 ; //set redundant cost to 0.
		ncost += z.cost ;  //set necessary cost 
		//a zoid with the projection already exists. return
		return ;
	}
	//print_dag() ;
	double divide_and_conquer_cost = 0 ;
	//struct timeval start, end ;
	struct timespec start, end;
	bool divide_and_conquer = false ;
    if (sim_can_cut) 
	{
		double rcost1 = 0, ncost1 = 0 ;
		divide_and_conquer = true ;
		//gettimeofday(&start, 0);
		clock_gettime(CLOCK_MONOTONIC, &start);
		z.resize_children(total_num_subzoids) ;
        /* cut into space */
		symbolic_sawzoid_space_cut_interior(t0, t1, grid, index, f, 
										num_subzoids, rcost1, ncost1) ;
		//gettimeofday(&end, 0);
		clock_gettime(CLOCK_MONOTONIC, &end);
		//double diff = tdiff(&end, &start) ;
		double diff = tdiff2(&end, &start) ;
		divide_and_conquer_cost = diff - rcost1 + ncost1 ;
		assert (divide_and_conquer_cost >= 0.) ;
		/*if (divide_and_conquer_cost < (double) 0.)
		{
			cout << "ht " << lt << "diff " << diff << "rcost " << rcost1 << " ncost " << ncost1 << endl ;
			
    		for (int i = N_RANK-1; i >= 0; --i) 
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
			assert (0);
		}*/
    } 
	else if (lt > dt_recursive_) 
	{
		double rcost1 = 0, ncost1 = 0 ;
		divide_and_conquer = true ;
		//gettimeofday(&start, 0);
		clock_gettime(CLOCK_MONOTONIC, &start);
		z.resize_children(2) ;
        /* cut into time */
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
        symbolic_sawzoid_space_time_cut_interior(t0, t0+halflt, l_son_grid, 
				index, 0, rcost1, ncost1, f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        symbolic_sawzoid_space_time_cut_interior(t0+halflt, t1, l_son_grid, 
				index, 1, rcost1, ncost1, f);

		//gettimeofday(&end, 0);
		clock_gettime(CLOCK_MONOTONIC, &end);
		//double diff = tdiff(&end, &start) ;
		double diff = tdiff2(&end, &start) ;
		divide_and_conquer_cost = diff - rcost1 + ncost1 ;
		
		assert (divide_and_conquer_cost >= 0.) ;
		/*if (divide_and_conquer_cost < (double) 0)
		{
			cout << "ht " << lt << " diff " << diff << " rcost " << rcost1 << " ncost " << ncost1 << endl ;
    		for (int i = N_RANK-1; i >= 0; --i) 
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
			assert (0);
		}*/
    } 
    // base case
	//determine the looping cost on the zoid
	double loop_cost = 0. ;
	//gettimeofday(&start, 0);
	clock_gettime(CLOCK_MONOTONIC, &start) ;
	f(t0, t1, grid);
	//gettimeofday(&end, 0);
	clock_gettime(CLOCK_MONOTONIC, &end) ;
	//loop_cost = tdiff(&end, &start) ;
	loop_cost = tdiff2(&end, &start) ;
	//store the decision for the zoid  and pass the redundant cost 
	//to the parent
	//cout << "ht " << t1 - t0 << "divide n conquer cost " << 
	//	divide_and_conquer_cost << " loop cost " << loop_cost << endl ;
	if (divide_and_conquer && divide_and_conquer_cost < loop_cost)
	{
		m_zoids [index].decision = 1 ;
		rcost += loop_cost ;
		m_zoids [index].cost = divide_and_conquer_cost ;
	}
	else
	{
		m_zoids [index].decision = 0 ;
		rcost += divide_and_conquer_cost ;
		m_zoids [index].cost = loop_cost ;
	}
	//set necessary cost to zero since parent's timer includes the same.
	//ncost += 0 ; 
}

template <int N_RANK> template <typename F, typename BF>
inline void auto_tune<N_RANK>::symbolic_sawzoid_space_time_cut_boundary(
		int t0, int t1,	grid_info<N_RANK> const & grid, 
		unsigned long parent_index, int child_index, 
		double & rcost, double & ncost, F const & f, BF const & bf)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;

	int centroid = 0, width = 1 ; 
	int offset = 0, total_num_subzoids = 1 ;
	int num_subzoids [N_RANK] ;
	unsigned long key = 0 ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
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
		int short_side ;
		bool space_cut = false ;
		if (lb < tb)
		{
			short_side = lb ;
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
			num_subzoids [i] = 3 ;
			if (short_side - (thres << 2) >= 0)
			{
				num_subzoids [i] = 5 ;
			}
			total_num_subzoids *= num_subzoids [i] ;
		}
        sim_can_cut |= space_cut ;
		unsigned long dim_key = (unsigned long) lb << num_bits_width | tb ;
		key = key << offset | dim_key ;
		offset += num_bits_dim ;
        call_boundary |= l_touch_boundary;
    }
	//cout << "centroid " << centroid << " key " << key << endl ;
	unsigned long index ;
	bool projection_exists = check_and_create_projection (key, lt, 
									centroid, index, l_father_grid) ;
	//cout << " index " << index << endl ;
	zoid_type & z = m_zoids [index] ;
	zoid_type & parent = m_zoids [parent_index] ;
	//cout << "parent id " << parent.id << endl ;
	//cout << "address of parent " << &parent << endl ;
	//add the zoid as a child of the parent
	parent.add_child(&z, child_index, index) ;
	assert (index < m_zoids.size()) ;
	if (projection_exists)
	{ 
		//Add the cost of the zoid to the parent.
		//rcost += 0 ; //set redundant cost to 0.
		ncost += z.cost ;  //set necessary cost 
		//a zoid with the projection already exists. return
		return ;
	}
	//print_dag() ;
    if (call_boundary)
	{
        l_dt_stop = dt_recursive_boundary_;
	}
    else
	{
        l_dt_stop = dt_recursive_;
	}
	double divide_and_conquer_cost = 0 ;
	//struct timeval start, end ;
	struct timespec start, end ;
	bool divide_and_conquer = false ;
    if (sim_can_cut) 
	{
		double rcost1 = 0, ncost1 = 0 ;
		divide_and_conquer = true ;
		//gettimeofday(&start, 0);
		clock_gettime(CLOCK_MONOTONIC, &start) ;
		//cout << "space cut " << endl ;
		z.resize_children(total_num_subzoids) ;
        //cut into space 
        if (call_boundary) 
		{
            symbolic_sawzoid_space_cut_boundary(t0, t1, l_father_grid, index, 
							 f, bf, num_subzoids, rcost1, ncost1);
        }
		else
		{
            symbolic_sawzoid_space_cut_interior(t0, t1, l_father_grid, index, 
							f, num_subzoids, rcost1, ncost1);
		}
		//gettimeofday(&end, 0);
		clock_gettime(CLOCK_MONOTONIC, &end) ;
		//double diff = tdiff(&end, &start) ;
		double diff = tdiff2(&end, &start) ;
		divide_and_conquer_cost = diff - rcost1 + ncost1 ;
		assert (divide_and_conquer_cost >= 0.) ;
		/*if (divide_and_conquer_cost < (double) 0)
		{
			cout << "ht " << lt << " diff " << diff << " rcost " << rcost1 << " ncost " << ncost1 << endl ;
    		for (int i = N_RANK-1; i >= 0; --i) 
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
			assert (0);
		}*/
    } 
	else if (lt > l_dt_stop)  //time cut
	{
		double rcost1 = 0, ncost1 = 0 ;
		divide_and_conquer = true ;
		//gettimeofday(&start, 0);
		clock_gettime(CLOCK_MONOTONIC, &start) ;
		//cout << "time cut " << endl ;
		z.resize_children(2) ;
        // cut into time 
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
        if (call_boundary) {
            symbolic_sawzoid_space_time_cut_boundary(t0, t0+halflt, l_son_grid, 										index, 0, rcost1, ncost1, f, bf);
        } else {
            symbolic_sawzoid_space_time_cut_interior(t0, t0+halflt, l_son_grid,
										index, 0, rcost1, ncost1, f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            symbolic_sawzoid_space_time_cut_boundary(t0+halflt, t1, l_son_grid,
										index, 1, rcost1, ncost1, f, bf);
        } else {
            symbolic_sawzoid_space_time_cut_interior(t0+halflt, t1, l_son_grid,
										index, 1, rcost1, ncost1, f);
        }
		//gettimeofday(&end, 0);
		clock_gettime(CLOCK_MONOTONIC, &end) ;
		//double diff = tdiff(&end, &start) ;
		double diff = tdiff2(&end, &start) ;
		divide_and_conquer_cost = diff - rcost1 + ncost1 ;
		assert (divide_and_conquer_cost >= 0.) ;
		/*if (divide_and_conquer_cost < (double) 0)
		{
			cout << "ht " << lt << " diff " << diff << " rcost " << rcost1 << " ncost " << ncost1 << endl ;
    		for (int i = N_RANK-1; i >= 0; --i) 
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
			}
			assert (0);
		}*/
    } 
	// base case
	//determine the looping cost on the zoid
	double loop_cost = 0. ;
	//gettimeofday(&start, 0);
	clock_gettime(CLOCK_MONOTONIC, &start) ;
	if (call_boundary) 
	{
		base_case_kernel_boundary(t0, t1, l_father_grid, bf);
	} 
	else 
	{ 
		f(t0, t1, l_father_grid);
	}
	//gettimeofday(&end, 0);
	clock_gettime(CLOCK_MONOTONIC, &end) ;
	//loop_cost = tdiff(&end, &start) ;
	loop_cost = tdiff2(&end, &start) ;
	//store the decision for the zoid  and pass the redundant cost 
	//to the parent
	//cout << "ht " << t1 - t0 << "divide n conquer cost " << 
	//	divide_and_conquer_cost << " loop cost " << loop_cost << endl ;
	if (divide_and_conquer && divide_and_conquer_cost < loop_cost)
	{
		m_zoids [index].decision = 1 ;
		rcost += loop_cost ;
		m_zoids [index].cost = divide_and_conquer_cost ;
	}
	else
	{
		m_zoids [index].decision = 0 ;
		rcost += divide_and_conquer_cost ;
		m_zoids [index].cost = loop_cost ;
	}
	//set necessary cost to zero since parent's timer includes the same.
	//ncost += 0 ; 
}

// sawzoid space cuts. 
template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::sawzoid_space_cut_interior(
				int t0, int t1, grid_info<N_RANK> const & grid, 
				zoid_type * projection_zoid, F const & f)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

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
					assert(child_index < projection_zoid->num_children) ;
					unsigned long index = projection_zoid->children [child_index] ;
					zoid_type * child = &(m_zoids [index]) ;
					child_index++ ; //looks like a race
                    sawzoid_space_time_cut_interior(l_son->t0, 
						l_son->t1, l_son->grid, child, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				assert(child_index < projection_zoid->num_children) ;
				unsigned long index = projection_zoid->children [child_index] ;
				zoid_type * child = &(m_zoids [index]) ;
				child_index++ ;
                if (queue_len_[curr_dep_pointer] == 0)
				{
                    sawzoid_space_time_cut_interior(l_father->t0,
 							l_father->t1, l_father->grid, child, f);
				}
                else
				{
                    cilk_spawn sawzoid_space_time_cut_interior(
						l_father->t0, l_father->t1, l_father->grid, child, f) ;
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
						const int cut_more = ((lb - (thres << 2)) >= 0) ;
						if (cut_more)
						//if (num_zoids [level] == 5)
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
						const int cut_more = ((tb - (thres << 2)) >= 0) ;
						if (cut_more)
						//if (num_zoids [level] == 5)
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
inline void auto_tune<N_RANK>::sawzoid_space_cut_boundary(
	int t0, int t1,	grid_info<N_RANK> const & grid, 
	zoid_type * projection_zoid, F const & f, BF const & bf)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

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
					assert(child_index < projection_zoid->num_children) ;
					unsigned long index = projection_zoid->children [child_index] ;
					zoid_type * child = &(m_zoids [index]) ;
					child_index++ ;
                    sawzoid_space_time_cut_boundary(l_son->t0, 
						l_son->t1, l_son->grid, child, f, bf);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				assert(child_index < projection_zoid->num_children) ;
				unsigned long index = projection_zoid->children [child_index] ;
				zoid_type * child = &(m_zoids [index]) ;
				child_index++ ;
                if (queue_len_[curr_dep_pointer] == 0) {
                    sawzoid_space_time_cut_boundary(l_father->t0,
						l_father->t1, l_father->grid, child, f, bf);
                } else {
                    cilk_spawn sawzoid_space_time_cut_boundary(
						l_father->t0, l_father->t1, l_father->grid, child, f, bf);
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

						const int cut_more = ((lb - (thres << 2)) >= 0) ;
						if (cut_more)
						//if (num_zoids [level] == 5)
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
							const int cut_more = ((lb - (thres << 2)) >= 0) ;
							if (cut_more)
							//if (num_zoids [level] == 5)
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
							const int cut_more = ((tb - (thres << 2)) >= 0) ;
							if (cut_more)
							//if (num_zoids [level] == 5)
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
inline void auto_tune<N_RANK>::
sawzoid_space_time_cut_interior(int t0,
				int t1, grid_info<N_RANK> const & grid, 
				zoid_type * projection_zoid, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;
	int num_subzoids [N_RANK] ;

	assert (projection_zoid) ;
	assert (projection_zoid->height == lt) ;
	if (projection_zoid->decision == 0)
	{
		//loop
		/*cout << "ht of base case " << t1 - t0 << endl ;
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
		}*/
		f(t0, t1, grid);
		return ;
	}
	
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);

#ifndef NDEBUG
		zoid_type * z = projection_zoid ;
		grid_info <N_RANK> & grid2 = z->info ;
		int x0 = pmod(grid.x0 [i], phys_length_ [i]) ;
		int x1 = pmod(grid.x1 [i], phys_length_ [i]) ;

		int x0_ = pmod(grid2.x0 [i], phys_length_ [i]) ;
		int x1_ = pmod(grid2.x1 [i], phys_length_ [i]) ;

		int x2 = pmod(grid.x0[i] + grid.dx0[i] * lt, phys_length_ [i]) ;
		int x3 = pmod(grid.x1[i] + grid.dx1[i] * lt, phys_length_ [i]) ;

		int x2_ = pmod(grid2.x0[i] + grid2.dx0[i] * lt, phys_length_ [i]) ;
		int x3_ = pmod(grid2.x1[i] + grid2.dx1[i] * lt, phys_length_ [i]) ;

		if (x0 != x0_ || x1 != x1_ || x2 != x2_ || x3 != x3_)
		{
			cout << "zoid and proj zoid differ " << endl ;
			cout << "\nzoid " << z->id << " height " << z->height <<
				//" num children " << z->children.size() <<
				" num children " << z->num_children <<
				" num_parents " << z->parents.size() << " geneity " ;
			//print_bits(&(z->geneity), sizeof(word_type) * 8);
			int h = z->height ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
				cout << " x0 [" << i << "] " << grid2.x0 [i]
				 << " x1 [" << i << "] " << grid2.x1 [i]
				<< " x2 [" << i << "] " << grid2.x0[i] + grid2.dx0[i] * h
				<< " x3 [" << i << "] " << grid2.x1[i] + grid2.dx1[i] * h
				<< " h " << h << endl ;
			}
			assert(0) ;
		}
#endif
        bool cut_lb = (lb < tb);
        thres = (slope_[i] * lt);
        sim_can_cut = SIM_CAN_CUT_I ;
    }
    if (sim_can_cut) 
	{
        /* cut into space */
        sawzoid_space_cut_interior(t0, t1, grid, projection_zoid,
												f) ;
    } 
	else if (lt > dt_recursive_) 
	{
        /* cut into time */
		assert (projection_zoid->children [0]) ;
		assert (projection_zoid->children [1]) ;
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
		unsigned long index = projection_zoid->children [0] ;
        sawzoid_space_time_cut_interior(t0, t0+halflt, 
									l_son_grid, &(m_zoids [index]), f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
		index = projection_zoid->children [1] ;
        sawzoid_space_time_cut_interior(t0+halflt, t1, 
								l_son_grid, &(m_zoids [index]), f);
    }
}

template <int N_RANK> template <typename F, typename BF>
inline void auto_tune<N_RANK>::
sawzoid_space_time_cut_boundary(int t0,
					int t1,	grid_info<N_RANK> const & grid, 
					zoid_type * projection_zoid, F const & f, BF const & bf)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;
	int num_subzoids [N_RANK] ;

	assert (projection_zoid) ;
	assert (projection_zoid->height == lt) ;

    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
#ifndef NDEBUG
		zoid_type * z = projection_zoid ;
		grid_info <N_RANK> & grid2 = z->info ;
		int x0 = pmod(grid.x0 [i], phys_length_ [i]) ;
		int x1 = pmod(grid.x1 [i], phys_length_ [i]) ;

		int x0_ = pmod(grid2.x0 [i], phys_length_ [i]) ;
		int x1_ = pmod(grid2.x1 [i], phys_length_ [i]) ;

		int x2 = pmod(grid.x0[i] + grid.dx0[i] * lt, phys_length_ [i]) ;
		int x3 = pmod(grid.x1[i] + grid.dx1[i] * lt, phys_length_ [i]) ;

		int x2_ = pmod(grid2.x0[i] + grid2.dx0[i] * lt, phys_length_ [i]) ;
		int x3_ = pmod(grid2.x1[i] + grid2.dx1[i] * lt, phys_length_ [i]) ;

		if (x0 != x0_ || x1 != x1_ || x2 != x2_ || x3 != x3_)
		{
			cout << "zoid and proj zoid differ " << endl ;
			cout << "\nzoid " << z->id << " height " << z->height <<
				//" num children " << z->children.size() <<
				" num children " << z->num_children <<
				" num_parents " << z->parents.size() << " geneity " ;
			//print_bits(&(z->geneity), sizeof(word_type) * 8);
			int h = z->height ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
				cout << " x0 [" << i << "] " << grid2.x0 [i]
				 << " x1 [" << i << "] " << grid2.x1 [i]
				<< " x2 [" << i << "] " << grid2.x0[i] + grid2.dx0[i] * h
				<< " x3 [" << i << "] " << grid2.x1[i] + grid2.dx1[i] * h
				<< " h " << h << endl ;
			}
			assert(0) ;
		}
#endif
        thres = slope_[i] * lt ;
		bool cut_lb = (lb < tb);
		sim_can_cut = SIM_CAN_CUT_B ;
        call_boundary |= l_touch_boundary;
    }
	
	if (projection_zoid->decision == 0)
	{
		/*cout << "ht of base case " << t1 - t0 << endl ;
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
				cout << " x0 [" << i << "] " << grid.x0 [i] 
				 << " x1 [" << i << "] " << grid.x1 [i] 
				<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
				<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
				<< " lt " << lt << endl ;
		}*/
		//loop
		if (call_boundary) {
            base_case_kernel_boundary(t0, t1, l_father_grid, bf);
        } else { 
            f(t0, t1, l_father_grid);
        }
		return ;
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
            sawzoid_space_cut_boundary(t0, t1, l_father_grid, 
										projection_zoid, f, bf) ;
        }
		else
		{
            sawzoid_space_cut_interior(t0, t1, l_father_grid, 
										projection_zoid, f) ;
		}
    } 
	else if (lt > l_dt_stop)  //time cut
	{
		//cout << "time cut " << endl ;
        // cut into time 
		assert (projection_zoid->children [0]) ;
		assert (projection_zoid->children [1]) ;
		
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
		unsigned long index = projection_zoid->children [0] ;
        if (call_boundary) {
            sawzoid_space_time_cut_boundary(t0, t0+halflt, 
								l_son_grid, &(m_zoids [index]), f, bf);
        } else {
            sawzoid_space_time_cut_interior(t0, t0+halflt, 
								l_son_grid, &(m_zoids [index]), f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
		index = projection_zoid->children [1] ;
        if (call_boundary) {
            sawzoid_space_time_cut_boundary(t0+halflt, t1, 
								l_son_grid, &(m_zoids [index]), f, bf);
        } else {
            sawzoid_space_time_cut_interior(t0+halflt, t1, 
								l_son_grid, &(m_zoids [index]), f);
        }
    } 
}
//The following code is not needed for now.

/* This is the version for boundary region cut! */
//template <int N_RANK> template <typename F, typename BF>
//inline void auto_tune<N_RANK>::abnormal_space_time_cut_boundary(int t0, int t1, grid_info<N_RANK> const grid, F const & f, BF const & bf)
//{
//    const int lt = t1 - t0;
//    bool sim_can_cut = false, call_boundary = false;
//    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
//    int l_dt_stop;
//
//#if STAT
//    int l_count_cut = 0;
//    int l_bottom_total_area = 1;
//    int l_top_total_area = 1;
//    int l_total_points;
//#endif
//
//    for (int i = N_RANK-1; i >= 0; --i) {
//        int lb, thres, tb;
//        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
//        lb = (grid.x1[i] - grid.x0[i]);
//        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
//        thres = (slope_[i] * lt);
//        
//        bool cut_lb = (lb < tb);
//		sim_can_cut = SIM_CAN_CUT_B ; 
//        call_boundary |= l_touch_boundary;
//#if STAT
//        l_count_cut = (l_can_cut ? l_count_cut + 1 : l_count_cut);
//        l_bottom_total_area *= lb;
//        l_top_total_area *= tb;
//#endif
//    }
//	if (N_RANK == 1)
//	{
//		sim_can_cut = sim_can_cut && lt > 1 ;
//	}
//
//    if (sim_can_cut) {
//        /* cut into space */
//#if STAT
//        ++sim_count_cut[l_count_cut];
//#endif
//        if (call_boundary) 
//		{
//            abnormal_space_cut_boundary(t0, t1, l_father_grid, f, bf) ;
//		}
//        else
//		{
//            abnormal_space_cut_interior(t0, t1, l_father_grid, f);
//		}
//        return;
//    } 
//
//    if (call_boundary)
//        l_dt_stop = dt_recursive_boundary_;
//    else
//        l_dt_stop = dt_recursive_;
//
//    if (lt > l_dt_stop) {
//        /* cut into time */
//        int halflt = lt / 2;
//        l_son_grid = l_father_grid;
//        if (call_boundary) {
//            abnormal_space_time_cut_boundary(t0, t0+halflt, l_son_grid, f, bf);
//        } else {
//            abnormal_space_time_cut_interior(t0, t0+halflt, l_son_grid, f);
//        }
//
//        for (int i = 0; i < N_RANK; ++i) {
//            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
//            l_son_grid.dx0[i] = l_father_grid.dx0[i];
//            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
//            l_son_grid.dx1[i] = l_father_grid.dx1[i];
//        }
//        if (call_boundary) {
//            abnormal_space_time_cut_boundary(t0+halflt, t1, l_son_grid, f, bf);
//        } else {
//            abnormal_space_time_cut_interior(t0+halflt, t1, l_son_grid, f);
//        }
//        return;
//    } 
//
//        // base case
//#if DEBUG
//        printf("call boundary!\n");
//        print_grid(stdout, t0, t1, l_father_grid);
//#endif
//#if STAT
//        ++boundary_region_count;
//        boundary_points_count += l_total_points;
//#endif
//        if (call_boundary) {
//            base_case_kernel_boundary(t0, t1, l_father_grid, bf);
//        } else {
//            f(t0, t1, l_father_grid);
//        }
//        return;
//}

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::symbolic_abnormal_space_time_cut_boundary(
		int t0, int t1,	grid_info<N_RANK> const & grid, 
		unsigned long parent_index, int child_index, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;

	int centroid = 0, width = 1 ; //, bottom_volume = 1, top_volume = 1 ;
	int offset = 0, num_subzoids = 1 ;
	unsigned long key = 0 ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
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
		int short_side ;
		bool space_cut = false ;
		if (lb < tb)
		{
			short_side = lb ;
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
		//if (short_side >= (thres << 1) && short_side > limit)
		if (short_side >= (thres << 1) && lb > limit)
		{
			space_cut = true ;
			int num_pieces = 3 ;
			if (short_side - (thres << 2) >= 0)
			{
				num_pieces = 5 ;
			}
			num_subzoids *= num_pieces ;
		}
        sim_can_cut |= space_cut ;

		unsigned long dim_key = (unsigned long) lb << num_bits_width | tb ;
		key = key << offset | dim_key ;
		offset += num_bits_dim ;
        call_boundary |= l_touch_boundary;
    }
	//cout << "centroid " << centroid << " key " << key << endl ;
	unsigned long index ;
	bool projection_exists = check_and_create_projection (key, lt, 
											centroid, index, l_father_grid) ;
											//centroid, &z, l_father_grid) ;
	//cout << " index " << index << endl ;
	zoid_type & z = m_zoids [index] ;
	zoid_type & parent = m_zoids [parent_index] ;
	//have some confusion about 
	//do we want to find the projection of the given grid? yes
	//does the grid correspond to parent?  seems No.
	//who are the children of the parent?  those from time cuts and space cuts
	//how to add the children to a parent? from here and from space cut routines.
	//cout << "parent id " << parent.id << endl ;
	//cout << "address of parent " << &parent << endl ;
	//add the zoid as a child of the parent
	//parent->add_child(z, child_index) ;
	parent.add_child(&z, child_index, index) ;
	//print_dag() ;
	
	if (projection_exists)
	{ 
		//update the geneity of the parent.
		parent.geneity |= z.geneity ;
		//a zoid with the projection already exists. return
		return ;
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
		z.resize_children(num_subzoids) ;
        //cut into space 
        if (call_boundary) 
		{
            symbolic_abnormal_space_cut_boundary(t0, t1, l_father_grid, index, f);
        }
		else
		{
            symbolic_abnormal_space_cut_interior(t0, t1, l_father_grid, index, f);
		}
    } 
	else if (lt > l_dt_stop)  //time cut
	{
		//cout << "time cut " << endl ;
		z.resize_children(2) ;
        // cut into time 
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
        if (call_boundary) {
            symbolic_abnormal_space_time_cut_boundary(t0, t0+halflt, l_son_grid, index, 0, f);
        } else {
            symbolic_abnormal_space_time_cut_interior(t0, t0+halflt, l_son_grid, index, 0, f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            symbolic_abnormal_space_time_cut_boundary(t0+halflt, t1, l_son_grid, index, 1, f);
        } else {
            symbolic_abnormal_space_time_cut_interior(t0+halflt, t1, l_son_grid, index, 1, f);
        }
    } 
	else
	{
		// base case
		//determine the geneity of the leaf
		compute_geneity(lt, l_father_grid, z.geneity, f) ;
		//print_bits(&(z->geneity), sizeof(word_type) * 8) ;
	}
	//update the geneity of the parent.
	m_zoids [parent_index].geneity |= m_zoids [index].geneity ;
	//parent->geneity |= z->geneity ;
}

/* This is the version for interior region cut! */
//template <int N_RANK> template <typename F>
//inline void auto_tune<N_RANK>::abnormal_space_time_cut_interior(int t0, int t1, grid_info<N_RANK> const grid, F const & f)
//{
//    const int lt = t1 - t0;
//    bool sim_can_cut = false;
//    grid_info<N_RANK> l_son_grid;
//#if STAT
//    int l_count_cut = 0;
//    int l_bottom_total_area = 1;
//    int l_top_total_area = 1;
//    int l_total_points;
//#endif
//
//    for (int i = N_RANK-1; i >= 0; --i) {
//        int lb, thres, tb;
//        lb = (grid.x1[i] - grid.x0[i]);
//        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
//        bool cut_lb = (lb < tb);
//        thres = (slope_[i] * lt);
//		sim_can_cut = SIM_CAN_CUT_I ;
//#if STAT
//        l_count_cut = (l_can_cut ? l_count_cut+1 : l_count_cut);
//#endif
//    }
//	if (N_RANK == 1)
//	{
//		sim_can_cut = sim_can_cut && lt > 1 ;
//	}
//
//    if (sim_can_cut) {
//        /* cut into space */
//#if STAT
//        ++sim_count_cut[l_count_cut];
//#endif
//        abnormal_space_cut_interior(t0, t1, grid, f);
//        return;
//    } else if (lt > dt_recursive_) {
//        /* cut into time */
//        assert(lt > dt_recursive_);
//        int halflt = lt / 2;
//        l_son_grid = grid;
//        abnormal_space_time_cut_interior(t0, t0+halflt, l_son_grid, f);
//
//        for (int i = 0; i < N_RANK; ++i) {
//            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
//            l_son_grid.dx0[i] = grid.dx0[i];
//            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
//            l_son_grid.dx1[i] = grid.dx1[i];
//        }
//        abnormal_space_time_cut_interior(t0+halflt, t1, l_son_grid, f);
//        return;
//    } else {
//        // base case
//#if DEBUG
//        printf("call interior!\n");
//        print_grid(stdout, t0, t1, grid);
//#endif
//#if STAT
//        ++interior_region_count;
//#endif
//        f(t0, t1, grid);
//        return;
//    }  
//}

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::symbolic_abnormal_space_time_cut_interior(
		int t0, int t1, grid_info<N_RANK> const & grid, 
		unsigned long parent_index, int child_index, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;

	int centroid = 0, width = 1 ; 
	int offset = 0, num_subzoids = 1 ;
	unsigned long key = 0 ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
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
		}
		else
		{
			short_side = tb ;
		}

		//if (short_side >= (thres << 1) && short_side > dx_recursive_[i])
		if (short_side >= (thres << 1) && lb > dx_recursive_[i])
		{
			space_cut = true ;
			int num_pieces = 3 ;
			if (short_side - (thres << 2) >= 0)
			{
				num_pieces = 5 ;
			}
			num_subzoids *= num_pieces ;
		}
        sim_can_cut |= space_cut ;
    }
	//cout << "centroid " << centroid << " key " << key << endl ;
	unsigned long index ;
	bool projection_exists = check_and_create_projection (key, lt, 
											centroid, index, grid) ;
	//cout << " index " << index << endl ;
	zoid_type & z = m_zoids [index];
	zoid_type & parent = m_zoids [parent_index] ;
	//cout << "parent id " << parent.id << endl ;
	//cout << "address of parent " << &parent << endl ;
	//add the zoid as a child of the parent
	parent.add_child(&z, child_index, index) ;
	//print_dag() ;
	
	if (projection_exists)
	{ 
		//update the geneity of the parent.
		parent.geneity |= z.geneity ;
		//a zoid with the projection already exists. return
		return ;
	}
    if (sim_can_cut) 
	{
		z.resize_children(num_subzoids) ;
        /* cut into space */
        symbolic_abnormal_space_cut_interior(t0, t1, grid, index, f);
    } 
	else if (lt > dt_recursive_) 
	{
		z.resize_children(2) ;
        /* cut into time */
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
        symbolic_abnormal_space_time_cut_interior(t0, t0+halflt, l_son_grid, index, 0, f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        symbolic_abnormal_space_time_cut_interior(t0+halflt, t1, l_son_grid, index, 1, f);
    } 
	else 
	{
        // base case
		//determine the geneity of the leaf
		compute_geneity(lt, grid, z.geneity, f) ;
		//print_bits(&(z->geneity), sizeof(word_type) * 8) ;
	}
	//update the geneity of the parent.
	m_zoids [parent_index].geneity |= m_zoids [index].geneity ;
	//parent->geneity |= z->geneity ;
}

//template <int N_RANK> template <typename F, typename BF>
//inline void auto_tune<N_RANK>::abnormal_space_cut_boundary(int t0, int t1, grid_info<N_RANK> const grid, F const & f, BF const & bf)
//{
//    queue_info *l_father;
//    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
//    int queue_head_[2], queue_tail_[2], queue_len_[2];
//
//    for (int i = 0; i < 2; ++i) {
//        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
//    }
//	
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
//                    abnormal_space_time_cut_boundary(l_son->t0, l_son->t1, l_son->grid, f, bf);
//                } // end cilk_for 
//                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
//                queue_len_[curr_dep_pointer] = 0;
//#else
//                // use cilk_spawn to spawn all the sub-grid 
//                pop_queue(curr_dep_pointer);
//                if (queue_len_[curr_dep_pointer] == 0) {
//                    abnormal_space_time_cut_boundary(l_father->t0, l_father->t1, l_father->grid, f, bf);
//                } else {
//                    cilk_spawn abnormal_space_time_cut_boundary(l_father->t0, l_father->t1, l_father->grid, f, bf);
//                }
//#endif
//            } else {
//                // performing a space cut on dimension 'level' 
//                pop_queue(curr_dep_pointer);
//                grid_info<N_RANK> l_father_grid = l_father->grid;
//                const int t0 = l_father->t0, t1 = l_father->t1;
//                const int lt = (t1 - t0);
//                const int level = l_father->level;
//                const int thres = slope_[level] * lt;
//                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
//                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
//                const bool cut_lb = (lb < tb);
//                const bool l_touch_boundary = touch_boundary(level, lt, l_father_grid);
//				const bool can_cut = CAN_CUT_B ;
//                if (!can_cut) {
//                    // if we can't cut into this dimension, just directly push
//                    // it into the circular queue
//                    //
//                    push_queue(curr_dep_pointer, level-1, t0, t1, l_father_grid);
//                } else  {
//                    /* can_cut */
//                    if (cut_lb) {
//                        grid_info<N_RANK> l_son_grid = l_father_grid;
//                        const int l_start = (l_father_grid.x0[level]);
//                        const int l_end = (l_father_grid.x1[level]);
//
//                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
//                        l_son_grid.x0[level] = l_start ;
//                        l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
//                        l_son_grid.x1[level] = l_start ;
//                        l_son_grid.dx1[level] = slope_[level] ;
//                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//                        l_son_grid.x0[level] = l_end ;
//                        l_son_grid.dx0[level] = -slope_[level];
//                        l_son_grid.x1[level] = l_end ;
//                        l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
//                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//						const int cut_more = ((lb - (thres << 2)) >= 0) ;
//						if (cut_more)
//						{
//							const int offset = (thres << 1) ;
//							l_son_grid.x0[level] = l_start ;
//	                        l_son_grid.dx0[level] = slope_[level];
//    	                    l_son_grid.x1[level] = l_start + offset;
//        	                l_son_grid.dx1[level] = -slope_[level] ;
//            	            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//							
//							l_son_grid.x0[level] = l_end - offset ;
//	                        l_son_grid.dx0[level] = slope_[level];
//    	                    l_son_grid.x1[level] = l_end ;
//        	                l_son_grid.dx1[level] = -slope_[level] ;
//            	            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//
//							l_son_grid.x0[level] = l_start + offset;
//	                        l_son_grid.dx0[level] = -slope_[level];
//    	                    l_son_grid.x1[level] = l_end - offset ;
//        	                l_son_grid.dx1[level] = slope_[level] ;
//            	            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//						}
//						else
//						{
//							l_son_grid.x0[level] = l_start ;
//	                        l_son_grid.dx0[level] = slope_[level];
//    	                    l_son_grid.x1[level] = l_end;
//        	                l_son_grid.dx1[level] = -slope_[level] ;
//            	            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//						}
//                    } /* end if (cut_lb) */
//                    else { /* cut_tb */
//                        if (lb == phys_length_[level] && l_father_grid.dx0[level] == 0 && l_father_grid.dx1[level] == 0) { /* initial cut on the dimension */
//                            grid_info<N_RANK> l_son_grid = l_father_grid;
//                            const int l_start = (l_father_grid.x0[level]);
//                            const int l_end = (l_father_grid.x1[level]);
//                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
//							const int cut_more = ((lb - (thres << 2)) >= 0) ;
//							if (cut_more)
//							{
//								const int offset = (thres << 1) ;
//                            	l_son_grid.x0[level] = l_start ;
//                            	l_son_grid.dx0[level] = slope_[level];
//                            	l_son_grid.x1[level] = l_start + offset ;
//                            	l_son_grid.dx1[level] = -slope_[level];
//                            	push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//
//                            	l_son_grid.x0[level] = l_end - offset ;
//                            	l_son_grid.dx0[level] = slope_[level];
//                            	l_son_grid.x1[level] = l_end ;
//                            	l_son_grid.dx1[level] = -slope_[level];
//                            	push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//
//                            	l_son_grid.x0[level] = l_start + offset ;
//                            	l_son_grid.dx0[level] = -slope_[level];
//                            	l_son_grid.x1[level] = l_end - offset ;
//                            	l_son_grid.dx1[level] = slope_[level];
//                            	push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//							}
//							else
//							{
//                            	l_son_grid.x0[level] = l_start ;
//                            	l_son_grid.dx0[level] = slope_[level];
//                            	l_son_grid.x1[level] = l_end ;
//                            	l_son_grid.dx1[level] = -slope_[level];
//                            	push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//							}
//                            l_son_grid.x0[level] = l_end ;
//                            l_son_grid.dx0[level] = -slope_[level];
//                            l_son_grid.x1[level] = l_end ;
//                            l_son_grid.dx1[level] = slope_[level];
//                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
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
//							push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//
//							l_son_grid.x0[level] = l_end - offset ;
//							l_son_grid.dx0[level] = slope_[level];
//							l_son_grid.x1[level] = l_end;
//							l_son_grid.dx1[level] = l_father_grid.dx1[level];
//							push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//
//							const int next_dep_pointer = (curr_dep + 1) & 0x1;
//							const int cut_more = ((tb - (thres << 2)) >= 0) ;
//							if (cut_more)
//							{
//								l_son_grid.x0[level] = l_start + offset ;
//								l_son_grid.dx0[level] = -slope_[level];
//								l_son_grid.x1[level] = l_start + offset;
//								l_son_grid.dx1[level] = slope_[level];
//								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//
//								l_son_grid.x0[level] = l_end - offset ;
//								l_son_grid.dx0[level] = -slope_[level];
//								l_son_grid.x1[level] = l_end - offset ;
//								l_son_grid.dx1[level] = slope_[level];
//								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//
//								l_son_grid.x0[level] = l_start + offset ;
//								l_son_grid.dx0[level] = slope_[level];
//								l_son_grid.x1[level] = l_end - offset ;
//								l_son_grid.dx1[level] = -slope_[level];
//								push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//							}
//							else
//							{
//								l_son_grid.x0[level] = l_start + offset ;
//								l_son_grid.dx0[level] = -slope_[level];
//								l_son_grid.x1[level] = l_end - offset ;
//								l_son_grid.dx1[level] = slope_[level];
//								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
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

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::symbolic_abnormal_space_cut_boundary(
		int t0, int t1, grid_info<N_RANK> const & grid, 
		unsigned long parent_index, F const & f)
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
                    symbolic_abnormal_space_time_cut_boundary(l_son->t0, 
						l_son->t1, l_son->grid, parent_index, child_index, f);
					child_index++ ; //this can be a race.
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				symbolic_abnormal_space_time_cut_boundary(l_father->t0, 
					l_father->t1, l_father->grid, parent_index, child_index, f);
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
                const bool l_touch_boundary = touch_boundary(level, lt, l_father_grid);
				const bool can_cut = CAN_CUT_B ;
                if (!can_cut) {
                    // if we can't cut into this dimension, just directly push
                    // it into the circular queue
                    //
                    push_queue(curr_dep_pointer, level-1, t0, t1, l_father_grid);
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
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                        l_son_grid.x0[level] = l_end ;
                        l_son_grid.dx0[level] = -slope_[level];
                        l_son_grid.x1[level] = l_end ;
                        l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
                        push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
						const int cut_more = ((lb - (thres << 2)) >= 0) ;
						if (cut_more)
						{
							const int offset = (thres << 1) ;
							l_son_grid.x0[level] = l_start ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_start + offset;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
							
							l_son_grid.x0[level] = l_end - offset ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_end ;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

							l_son_grid.x0[level] = l_start + offset;
	                        l_son_grid.dx0[level] = -slope_[level];
    	                    l_son_grid.x1[level] = l_end - offset ;
        	                l_son_grid.dx1[level] = slope_[level] ;
            	            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
						}
						else
						{
							l_son_grid.x0[level] = l_start ;
	                        l_son_grid.dx0[level] = slope_[level];
    	                    l_son_grid.x1[level] = l_end;
        	                l_son_grid.dx1[level] = -slope_[level] ;
            	            push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
						}
                    } /* end if (cut_lb) */
                    else { /* cut_tb */
                        if (lb == phys_length_[level] && l_father_grid.dx0[level] == 0 && l_father_grid.dx1[level] == 0) { /* initial cut on the dimension */
                            grid_info<N_RANK> l_son_grid = l_father_grid;
                            const int l_start = (l_father_grid.x0[level]);
                            const int l_end = (l_father_grid.x1[level]);
                            const int next_dep_pointer = (curr_dep + 1) & 0x1;
							const int cut_more = ((lb - (thres << 2)) >= 0) ;
							if (cut_more)
							{
								const int offset = (thres << 1) ;
                            	l_son_grid.x0[level] = l_start ;
                            	l_son_grid.dx0[level] = slope_[level];
                            	l_son_grid.x1[level] = l_start + offset ;
                            	l_son_grid.dx1[level] = -slope_[level];
                            	push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            	l_son_grid.x0[level] = l_end - offset ;
                            	l_son_grid.dx0[level] = slope_[level];
                            	l_son_grid.x1[level] = l_end ;
                            	l_son_grid.dx1[level] = -slope_[level];
                            	push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

                            	l_son_grid.x0[level] = l_start + offset ;
                            	l_son_grid.dx0[level] = -slope_[level];
                            	l_son_grid.x1[level] = l_end - offset ;
                            	l_son_grid.dx1[level] = slope_[level];
                            	push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
							}
							else
							{
                            	l_son_grid.x0[level] = l_start ;
                            	l_son_grid.dx0[level] = slope_[level];
                            	l_son_grid.x1[level] = l_end ;
                            	l_son_grid.dx1[level] = -slope_[level];
                            	push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
							}
                            l_son_grid.x0[level] = l_end ;
                            l_son_grid.dx0[level] = -slope_[level];
                            l_son_grid.x1[level] = l_end ;
                            l_son_grid.dx1[level] = slope_[level];
                            push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
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
							push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

							l_son_grid.x0[level] = l_end - offset ;
							l_son_grid.dx0[level] = slope_[level];
							l_son_grid.x1[level] = l_end;
							l_son_grid.dx1[level] = l_father_grid.dx1[level];
							push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

							const int next_dep_pointer = (curr_dep + 1) & 0x1;
							const int cut_more = ((tb - (thres << 2)) >= 0) ;
							if (cut_more)
							{
								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = -slope_[level];
								l_son_grid.x1[level] = l_start + offset;
								l_son_grid.dx1[level] = slope_[level];
								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

								l_son_grid.x0[level] = l_end - offset ;
								l_son_grid.dx0[level] = -slope_[level];
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = slope_[level];
								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = slope_[level];
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = -slope_[level];
								push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
							}
							else
							{
								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = -slope_[level];
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = slope_[level];
								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
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

//template <int N_RANK> template <typename F>
//inline void auto_tune<N_RANK>::abnormal_space_cut_interior(int t0, int t1, grid_info<N_RANK> const grid, F const & f)
//{
//    queue_info *l_father;
//    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
//    int queue_head_[2], queue_tail_[2], queue_len_[2];
//
//    for (int i = 0; i < 2; ++i) {
//        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
//    }
//
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
//                    abnormal_space_time_cut_interior(l_son->t0, l_son->t1, l_son->grid, f);
//                } // end cilk_for 
//                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
//                queue_len_[curr_dep_pointer] = 0;
//#else
//                // use cilk_spawn to spawn all the sub-grid 
//                pop_queue(curr_dep_pointer);
//                if (queue_len_[curr_dep_pointer] == 0)
//                    abnormal_space_time_cut_interior(l_father->t0, l_father->t1, l_father->grid, f);
//                else
//                    cilk_spawn abnormal_space_time_cut_interior(l_father->t0, l_father->t1, l_father->grid, f);
//#endif
//            } else {
//                // performing a space cut on dimension 'level' 
//                pop_queue(curr_dep_pointer);
//                const grid_info<N_RANK> l_father_grid = l_father->grid;
//                const int t0 = l_father->t0, t1 = l_father->t1;
//                const int lt = (t1 - t0);
//                const int level = l_father->level;
//                const int thres = slope_[level] * lt;
//                const int lb = (l_father_grid.x1[level] - l_father_grid.x0[level]);
//                const int tb = (l_father_grid.x1[level] + l_father_grid.dx1[level] * lt - l_father_grid.x0[level] - l_father_grid.dx0[level] * lt);
//                const bool cut_lb = (lb < tb);
//				const bool can_cut = CAN_CUT_I ;
//                if (!can_cut) {
//                    // if we can't cut into this dimension, just directly push 
//                    // it into the circular queue 
//                    //
//                    push_queue(curr_dep_pointer, level-1, t0, t1, l_father_grid);
//                } else  {
//                    /* can_cut! */
//                    if (cut_lb) {
//                        grid_info<N_RANK> l_son_grid = l_father_grid;
//                        const int l_start = (l_father_grid.x0[level]);
//                        const int l_end = (l_father_grid.x1[level]);
//
//                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
//						const int offset = (thres << 1) ;
//						//instead of creating triangle
//						//check if we can create a trapezoid which is a
//						//set of adjacent triangles. This is done to coarsen
//						//the base case in an abnormal region
//						if (tb - (offset * num_triangles [level] << 1) >= 0)
//						{
//							const int offset_2 = (num_triangles [level] - 1) * 
//													offset ;
//                        	l_son_grid.x0[level] = l_start ;
//                        	l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
//                        	l_son_grid.x1[level] = l_start + offset_2 ;
//                       		l_son_grid.dx1[level] = slope_[level] ;
//                        	push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//                        	l_son_grid.x0[level] = l_end - offset_2 ;
//                        	l_son_grid.dx0[level] = -slope_[level];
//                        	l_son_grid.x1[level] = l_end ;
//                        	l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
//                        	push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//							l_son_grid.x0[level] = l_start + offset_2 ;
//							l_son_grid.dx0[level] = slope_[level];
//							l_son_grid.x1[level] = l_end - offset_2 ;
//							l_son_grid.dx1[level] = -slope_[level] ;
//							push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//						}
//						else
//						{
//							l_son_grid.x0[level] = l_start ;
//							l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
//							l_son_grid.x1[level] = l_start ;
//							l_son_grid.dx1[level] = slope_[level] ;
//							push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//							l_son_grid.x0[level] = l_end ;
//							l_son_grid.dx0[level] = -slope_[level];
//							l_son_grid.x1[level] = l_end ;
//							l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
//							push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//							const int cut_more = ((lb - (thres << 2)) >= 0) ;
//							if (cut_more)
//							{
//								l_son_grid.x0[level] = l_start ;
//								l_son_grid.dx0[level] = slope_[level];
//								l_son_grid.x1[level] = l_start + offset ;
//								l_son_grid.dx1[level] = -slope_[level] ;
//								push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//								
//								l_son_grid.x0[level] = l_end - offset ;
//								l_son_grid.dx0[level] = slope_[level];
//								l_son_grid.x1[level] = l_end ;
//								l_son_grid.dx1[level] = -slope_[level] ;
//								push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//
//								l_son_grid.x0[level] = l_start + offset ;
//								l_son_grid.dx0[level] = -slope_[level] ;
//								l_son_grid.x1[level] = l_end - offset ;
//								l_son_grid.dx1[level] = slope_[level] ;
//								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//							}
//							else
//							{
//								l_son_grid.x0[level] = l_start ;
//								l_son_grid.dx0[level] = slope_[level];
//								l_son_grid.x1[level] = l_end ;
//								l_son_grid.dx1[level] = -slope_[level] ;
//								push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//							}
//						}
//                    } /* end if (cut_lb) */
//                    else {
//                        /* cut_tb */
//                        grid_info<N_RANK> l_son_grid = l_father_grid;
//                        const int l_start = (l_father_grid.x0[level]);
//                        const int l_end = (l_father_grid.x1[level]);
//						const int offset = (thres << 1) ;
//						const int next_dep_pointer = (curr_dep + 1) & 0x1;
//
//						if (lb - (offset * num_triangles [level] << 1)  >= 0)
//						{
//							const int offset_2 = num_triangles [level] * offset ;
//                        	l_son_grid.x0[level] = l_start ;
//                        	l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
//                        	l_son_grid.x1[level] = l_start + offset_2 ;
//                       		l_son_grid.dx1[level] = -slope_[level] ;
//                        	push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//                        	l_son_grid.x0[level] = l_end - offset_2 ;
//                        	l_son_grid.dx0[level] = slope_[level];
//                        	l_son_grid.x1[level] = l_end ;
//                        	l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
//                        	push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//							l_son_grid.x0[level] = l_start + offset_2 ;
//							l_son_grid.dx0[level] = -slope_[level];
//							l_son_grid.x1[level] = l_end - offset_2 ;
//							l_son_grid.dx1[level] = slope_[level] ;
//							push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//						}
//						else
//						{
//							l_son_grid.x0[level] = l_start;
//							l_son_grid.dx0[level] = l_father_grid.dx0[level];
//							l_son_grid.x1[level] = l_start + offset;
//							l_son_grid.dx1[level] = -slope_[level];
//							push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//
//							l_son_grid.x0[level] = l_end - offset ;
//							l_son_grid.dx0[level] = slope_[level];
//							l_son_grid.x1[level] = l_end;
//							l_son_grid.dx1[level] = l_father_grid.dx1[level];
//							push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//
//							const int cut_more = ((tb - (thres << 2)) >= 0) ;
//							if (cut_more)
//							{
//								l_son_grid.x0[level] = l_start + offset ;
//								l_son_grid.dx0[level] = -slope_[level];
//								l_son_grid.x1[level] = l_start + offset;
//								l_son_grid.dx1[level] = slope_[level];
//								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//
//								l_son_grid.x0[level] = l_end - offset ;
//								l_son_grid.dx0[level] = -slope_[level];
//								l_son_grid.x1[level] = l_end - offset ;
//								l_son_grid.dx1[level] = slope_[level];
//								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//
//								l_son_grid.x0[level] = l_start + offset ;
//								l_son_grid.dx0[level] = slope_[level];
//								l_son_grid.x1[level] = l_end - offset ;
//								l_son_grid.dx1[level] = -slope_[level];
//								push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
//							}
//							else
//							{
//								l_son_grid.x0[level] = l_start + offset ;
//								l_son_grid.dx0[level] = -slope_[level];
//								l_son_grid.x1[level] = l_end - offset ;
//								l_son_grid.dx1[level] = slope_[level];
//								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
//							}
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

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::symbolic_abnormal_space_cut_interior(
				int t0, int t1, grid_info<N_RANK> const & grid, 
				unsigned long parent_index, F const & f)
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
                    symbolic_abnormal_space_time_cut_interior(l_son->t0, 
						l_son->t1, l_son->grid, parent_index, child_index, f);
					child_index++ ; //this can be a race
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				symbolic_abnormal_space_time_cut_interior(l_father->t0, 
					l_father->t1, l_father->grid, parent_index, child_index, f);
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
				const bool can_cut = CAN_CUT_I ;
                if (!can_cut) {
                    // if we can't cut into this dimension, just directly push 
                    // it into the circular queue 
                    //
                    push_queue(curr_dep_pointer, level-1, t0, t1, l_father_grid);
                } else  {
                    /* can_cut! */
                    if (cut_lb) {
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);

                        const int next_dep_pointer = (curr_dep + 1) & 0x1;
						const int offset = (thres << 1) ;
						//instead of creating triangle
						//check if we can create a trapezoid which is a
						//set of adjacent triangles. This is done to coarsen
						//the base case in an abnormal region
						if (tb - (offset * num_triangles [level] << 1) >= 0)
						{
							const int offset_2 = (num_triangles [level] - 1) * 
													offset ;
                        	l_son_grid.x0[level] = l_start ;
                        	l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
                        	l_son_grid.x1[level] = l_start + offset_2 ;
                       		l_son_grid.dx1[level] = slope_[level] ;
                        	push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
                        	l_son_grid.x0[level] = l_end - offset_2 ;
                        	l_son_grid.dx0[level] = -slope_[level];
                        	l_son_grid.x1[level] = l_end ;
                        	l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
                        	push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
							l_son_grid.x0[level] = l_start + offset_2 ;
							l_son_grid.dx0[level] = slope_[level];
							l_son_grid.x1[level] = l_end - offset_2 ;
							l_son_grid.dx1[level] = -slope_[level] ;
							push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
						}
						else
						{
							l_son_grid.x0[level] = l_start ;
							l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
							l_son_grid.x1[level] = l_start ;
							l_son_grid.dx1[level] = slope_[level] ;
							push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
							l_son_grid.x0[level] = l_end ;
							l_son_grid.dx0[level] = -slope_[level];
							l_son_grid.x1[level] = l_end ;
							l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
							push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
							const int cut_more = ((lb - (thres << 2)) >= 0) ;
							if (cut_more)
							{
								l_son_grid.x0[level] = l_start ;
								l_son_grid.dx0[level] = slope_[level];
								l_son_grid.x1[level] = l_start + offset ;
								l_son_grid.dx1[level] = -slope_[level] ;
								push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
								
								l_son_grid.x0[level] = l_end - offset ;
								l_son_grid.dx0[level] = slope_[level];
								l_son_grid.x1[level] = l_end ;
								l_son_grid.dx1[level] = -slope_[level] ;
								push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = -slope_[level] ;
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = slope_[level] ;
								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
							}
							else
							{
								l_son_grid.x0[level] = l_start ;
								l_son_grid.dx0[level] = slope_[level];
								l_son_grid.x1[level] = l_end ;
								l_son_grid.dx1[level] = -slope_[level] ;
								push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
							}
						}
                    } /* end if (cut_lb) */
                    else {
                        /* cut_tb */
                        grid_info<N_RANK> l_son_grid = l_father_grid;
                        const int l_start = (l_father_grid.x0[level]);
                        const int l_end = (l_father_grid.x1[level]);
						const int offset = (thres << 1) ;
						const int next_dep_pointer = (curr_dep + 1) & 0x1;

						if (lb - (offset * num_triangles [level] << 1)  >= 0)
						{
							const int offset_2 = num_triangles [level] * offset ;
                        	l_son_grid.x0[level] = l_start ;
                        	l_son_grid.dx0[level] = l_father_grid.dx0[level] ;
                        	l_son_grid.x1[level] = l_start + offset_2 ;
                       		l_son_grid.dx1[level] = -slope_[level] ;
                        	push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
                        	l_son_grid.x0[level] = l_end - offset_2 ;
                        	l_son_grid.dx0[level] = slope_[level];
                        	l_son_grid.x1[level] = l_end ;
                        	l_son_grid.dx1[level] = l_father_grid.dx1[level] ;
                        	push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
							l_son_grid.x0[level] = l_start + offset_2 ;
							l_son_grid.dx0[level] = -slope_[level];
							l_son_grid.x1[level] = l_end - offset_2 ;
							l_son_grid.dx1[level] = slope_[level] ;
							push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
						}
						else
						{
							l_son_grid.x0[level] = l_start;
							l_son_grid.dx0[level] = l_father_grid.dx0[level];
							l_son_grid.x1[level] = l_start + offset;
							l_son_grid.dx1[level] = -slope_[level];
							push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

							l_son_grid.x0[level] = l_end - offset ;
							l_son_grid.dx0[level] = slope_[level];
							l_son_grid.x1[level] = l_end;
							l_son_grid.dx1[level] = l_father_grid.dx1[level];
							push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);

							const int cut_more = ((tb - (thres << 2)) >= 0) ;
							if (cut_more)
							{
								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = -slope_[level];
								l_son_grid.x1[level] = l_start + offset;
								l_son_grid.dx1[level] = slope_[level];
								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

								l_son_grid.x0[level] = l_end - offset ;
								l_son_grid.dx0[level] = -slope_[level];
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = slope_[level];
								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);

								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = slope_[level];
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = -slope_[level];
								push_queue(curr_dep_pointer, level-1, t0, t1, l_son_grid);
							}
							else
							{
								l_son_grid.x0[level] = l_start + offset ;
								l_son_grid.dx0[level] = -slope_[level];
								l_son_grid.x1[level] = l_end - offset ;
								l_son_grid.dx1[level] = slope_[level];
								push_queue(next_dep_pointer, level-1, t0, t1, l_son_grid);
							}
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
#undef dx_recursive_boundary_  
#undef dx_recursive_ 
#undef dt_recursive_boundary_ 
#undef dt_recursive_ 
#undef slope_ 
#undef touch_boundary 
#undef base_case_kernel_boundary

#endif
