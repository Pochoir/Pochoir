#ifndef GENEITY_PROBLEM_TRAP_HPP
#define GENEITY_PROBLEM_TRAP_HPP

#include "geneity_problem.hpp"

#define dx_recursive_boundary_  (m_algo.dx_recursive_boundary_)
#define dx_recursive_ (m_algo.dx_recursive_)
#define dt_recursive_boundary_ (m_algo.dt_recursive_boundary_)
#define dt_recursive_ (m_algo.dt_recursive_)
#define slope_ m_algo.slope_
#define touch_boundary m_algo.touch_boundary
#define phys_length_ m_algo.phys_length_
#define base_case_kernel_boundary m_algo.base_case_kernel_boundary

//grid - grid of the parent
//parent - zoid_type ptr of parent
//do space cuts on grid and call the symbolic_space_time_cut_interior routine
//The routine does a symbolic space cut for interior zoids. 
template <int N_RANK> template <typename F>
inline void geneity_problem<N_RANK>::symbolic_space_cut_interior(int t0, int t1,
				grid_info<N_RANK> const & grid, unsigned int parent_index, 
				F const & f)
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
                    symbolic_space_time_cut_interior(l_son->t0, l_son->t1, 
									l_son->grid, parent_index, child_index, f);
									//l_son->grid, parent, child_index, f);
					child_index++ ; //this can be a race
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                /*if (queue_len_[curr_dep_pointer] == 0)
				{
                    symbolic_space_time_cut_interior(l_father->t0, l_father->t1,
							 l_father->grid, parent_index, child_index, f);
									 //l_father->grid, parent, child_index, f);
				}
                else
				{
                    //cilk_spawn symbolic_space_time_cut_interior(l_father->t0, 
                    symbolic_space_time_cut_interior(l_father->t0, 
						l_father->t1, l_father->grid, parent_index, child_index, f);
				}*/
				symbolic_space_time_cut_interior(l_father->t0, 
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
#ifdef NDEBUG
        //cilk_sync;
#endif
#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}

//grid - grid of the parent
//parent - zoid_type ptr of parent
//do space cuts on grid and call the symbolic_space_time_cut_boundary routine
//The routine does a symbolic space cut for boundary zoids. 
template <int N_RANK> template <typename F>
inline void geneity_problem<N_RANK>::symbolic_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned int parent_index, F const & f)
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
                    symbolic_space_time_cut_boundary(l_son->t0, l_son->t1, 
									l_son->grid, parent_index, child_index, f);
										//l_son->grid, parent, child_index, f);
					child_index++ ; //this can be a race.
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                /*if (queue_len_[curr_dep_pointer] == 0) {
                    symbolic_space_time_cut_boundary(l_father->t0, l_father->t1,
								 l_father->grid, parent_index, child_index, f);
								 //l_father->grid, parent, child_index, f);
                } else {
                    //cilk_spawn symbolic_space_time_cut_boundary(l_father->t0, 
                    symbolic_space_time_cut_boundary(l_father->t0, 
						l_father->t1, l_father->grid, parent_index, child_index,f);
						//l_father->t1, l_father->grid, parent, child_index,f);
                }*/
				symbolic_space_time_cut_boundary(l_father->t0, 
					l_father->t1, l_father->grid, parent_index, child_index,f);
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
#ifdef NDEBUG
//        cilk_sync;
#endif
#endif
        assert(queue_len_[curr_dep_pointer] == 0);
    } // end for (curr_dep < N_RANK+1) 
}


//grid - grid of the child
//parent - zoid_type ptr of parent
//This routine does a symbolic space/time cut on interior zoids
template <int N_RANK> template <typename F>
inline void geneity_problem<N_RANK>::symbolic_space_time_cut_interior(int t0, 
		int t1, grid_info<N_RANK> const & grid, unsigned int parent_index, 
		int child_index, F const & f)
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
		//cout << "grid.x0[ " << i << "]" << grid.x0[i] << "phys_length_ [ " << i
		//	<< "]" << phys_length_ [i] << endl ;
		//cout << "pmod " << pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) << endl ;

		unsigned long dim_key = (unsigned long) lb << num_bits_width | tb ;
		key = key << offset | dim_key ;
		offset += num_bits_dim ;
		/*cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
			<< " lt " << lt << endl ;*/
        bool cut_lb = (lb < tb);
        thres = (slope_[i] * lt);
		bool space_cut = CAN_CUT_IN ;
        //sim_can_cut = SIM_CAN_CUT_I ;
        sim_can_cut |= space_cut ;
		if (space_cut)
		{
			num_subzoids *= 3 ;
		}
    }
	//cout << "centroid " << centroid << " key " << key << endl ;
	unsigned int index ;
	bool projection_exists = check_and_create_projection (key, lt, 
											centroid, index, grid) ;
	//cout << " index " << index << endl ;
	zoid_type & z = m_zoids [index] ;
	zoid_type & parent = m_zoids [parent_index] ;
	//cout << "parent id " << parent.id << endl ;
	//cout << "address of parent " << &parent << endl ;
	//add the zoid as a child of the parent
	m_children[parent.get_start() + child_index] = index ;
#ifndef NDEBUG
	parent.add_child(&z) ;
#endif
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
		//z.resize_children(num_subzoids) ;
        z.set_start(m_children.size()) ;
		z.set_num_children(num_subzoids) ;
		for (int i = 0 ; i < num_subzoids ; i++)
		{
			m_children.push_back(index) ;
		}
        /* cut into space */
        //symbolic_space_cut_interior(t0, t1, grid, z, f);
        symbolic_space_cut_interior(t0, t1, grid, index, f);
    } 
	else if (lt > dt_recursive_) 
	{
		//z.resize_children(2) ;
        z.set_start(m_children.size()) ;
		z.set_num_children(2) ;
		m_children.push_back(index) ;
		m_children.push_back(index) ;
        /* cut into time */
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
        //symbolic_space_time_cut_interior(t0, t0+halflt, l_son_grid, z, 0, f);
        symbolic_space_time_cut_interior(t0, t0+halflt, l_son_grid, index, 0, f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        //symbolic_space_time_cut_interior(t0+halflt, t1, l_son_grid, z, 1, f);
        symbolic_space_time_cut_interior(t0+halflt, t1, l_son_grid, index, 1, f);
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


//grid - grid of the child
//parent - zoid_type ptr of parent
//This routine does a symbolic space/time cut on boundary zoid
template <int N_RANK> template <typename F>
inline void geneity_problem<N_RANK>::symbolic_space_time_cut_boundary(int t0, 
		int t1,	grid_info<N_RANK> const & grid, unsigned int parent_index, 
		int child_index, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;

	int centroid = 0, width = 1 ;
	int offset = 0, num_subzoids = 1 ;
	unsigned long key = 0 ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
        lb = (grid.x1[i] - grid.x0[i]);
		centroid = pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) * width + 
					centroid ;
		//cout << "grid.x0[ " << i << "]" << grid.x0[i] << "phys_length_ [ " << i
		//	<< "]" << phys_length_ [i] << endl ;
		//cout << "pmod " << pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) << endl ;
		width = phys_length_ [i] ;
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
		/*cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
			<< " lt " << lt << endl ;*/
        thres = slope_[i] * lt ;
        bool cut_lb = (lb < tb);

		unsigned long dim_key = (unsigned long) lb << num_bits_width | tb ;
		key = key << offset | dim_key ;
		offset += num_bits_dim ;
		bool space_cut = CAN_CUT_BO ;
        //sim_can_cut = SIM_CAN_CUT_B ;
        sim_can_cut |= space_cut ;
		//if (sim_can_cut)
		if (space_cut)
		{
			//if initial space cut, there are 2 subzoids
			if (lb == phys_length_[i] && grid.dx0[i] == 0 && grid.dx1[i] == 0)
			{
				num_subzoids *= 2 ;
			}
			else
			{
				num_subzoids *= 3 ;
			}
		}
        call_boundary |= l_touch_boundary;
    }
	//cout << "centroid " << centroid << " key " << key << endl ;
	unsigned int index ;
	bool projection_exists = check_and_create_projection (key, lt, 
											centroid, index, l_father_grid) ;
	//cout << " index " << index << endl ;
	zoid_type & z = m_zoids [index] ;
	zoid_type & parent = m_zoids [parent_index] ;
	//cout << "parent id " << parent.id << endl ;
	//cout << "address of parent " << &parent << endl ;
	//add the zoid as a child of the parent
	m_children[parent.get_start() + child_index] = index ;
#ifndef NDEBUG
	parent.add_child(&z) ;
#endif
	//print_dag() ;
	
	if (projection_exists)
	{ 
		//update the geneity of the parent.
		parent.geneity |= z.geneity ;
		//a zoid with the projection already exists. return
		return ;
	}
	//projection doesn't exist.
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
		//z.resize_children(num_subzoids) ;
		z.set_start(m_children.size()) ;
		z.set_num_children(num_subzoids) ;
		for (int i = 0 ; i < num_subzoids ; i++)
		{
			m_children.push_back(index) ;
		}
        //cut into space 
        if (call_boundary) 
		{
            symbolic_space_cut_boundary(t0, t1, l_father_grid, index, f);
        }
		else
		{
            symbolic_space_cut_interior(t0, t1, l_father_grid, index, f);
		}
    } 
	else if (lt > l_dt_stop)  //time cut
	{
		//cout << "time cut " << endl ;
		//z.resize_children(2) ;
		z.set_start(m_children.size()) ;
		z.set_num_children(2) ;
		m_children.push_back(index) ;
		m_children.push_back(index) ;
        // cut into time 
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
        if (call_boundary) {
            symbolic_space_time_cut_boundary(t0, t0+halflt, l_son_grid, index, 
                                             0, f);
        } else {
            symbolic_space_time_cut_interior(t0, t0+halflt, l_son_grid, index, 
                                             0, f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            symbolic_space_time_cut_boundary(t0+halflt, t1, l_son_grid, index, 
                                             1, f);
        } else {
            symbolic_space_time_cut_interior(t0+halflt, t1, l_son_grid, index, 
                                             1, f);
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

//The routine does a homogeneous space cut for interior zoids. 
template <int N_RANK> template <typename F>
inline void geneity_problem<N_RANK>::homogeneous_space_cut_interior(int t0, 
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
					assert(child_index < projection_zoid->num_children) ;
#if 0
					unsigned int index = projection_zoid->children [child_index] ;
					zoid_type * child = &(m_zoids [index]) ;
					child_index++ ; //looks like a race
                    homogeneous_space_time_cut_interior(l_son->t0, l_son->t1, 
											l_son->grid, child, f);
#else
                    homogeneous_space_time_cut_interior(l_son->t0, l_son->t1, 
											l_son->grid, f);
#endif
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				//assert(child_index < projection_zoid->num_children) ;
//#ifndef NDEBUG
#if 0
				unsigned int index = projection_zoid->children [child_index] ;
				zoid_type * child = &(m_zoids [index]) ;
				child_index++ ; 
                if (queue_len_[curr_dep_pointer] == 0)
                    homogeneous_space_time_cut_interior(l_father->t0, 
							l_father->t1, l_father->grid, child, f);
                else
                    cilk_spawn homogeneous_space_time_cut_interior(l_father->t0, 							l_father->t1, l_father->grid, child, f);
#else
                if (queue_len_[curr_dep_pointer] == 0)
                    homogeneous_space_time_cut_interior(l_father->t0, 
							l_father->t1, l_father->grid, f);
                else
                    cilk_spawn homogeneous_space_time_cut_interior(l_father->t0, 	l_father->t1, l_father->grid, f);
#endif
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

//The routine does a homogeneous space cut for boundary zoids. 
template <int N_RANK> template <typename F>
inline void geneity_problem<N_RANK>::homogeneous_space_cut_boundary(int t0, 
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
					assert(child_index < projection_zoid->num_children) ;
#if 0
					unsigned int index = projection_zoid->children [child_index] ;
					zoid_type * child = &(m_zoids [index]) ;
					child_index++ ; 
                    homogeneous_space_time_cut_boundary(l_son->t0, l_son->t1, 
											l_son->grid, child, f);
#else
                    homogeneous_space_time_cut_boundary(l_son->t0, l_son->t1, 
											l_son->grid, f);
#endif
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				//assert(child_index < projection_zoid->num_children) ;
//#ifndef NDEBUG
#if 0
				unsigned int index = projection_zoid->children [child_index] ;
				zoid_type * child = &(m_zoids [index]) ;
				child_index++ ; 
                if (queue_len_[curr_dep_pointer] == 0) {
                    homogeneous_space_time_cut_boundary(l_father->t0, 
						l_father->t1, l_father->grid, child, f);
                } else {
                    cilk_spawn homogeneous_space_time_cut_boundary(l_father->t0, 			l_father->t1, l_father->grid, child, f) ;
                }
#else
                if (queue_len_[curr_dep_pointer] == 0) {
                    homogeneous_space_time_cut_boundary(l_father->t0, 
						l_father->t1, l_father->grid, f);
                } else {
                    cilk_spawn homogeneous_space_time_cut_boundary(l_father->t0, l_father->t1, l_father->grid, f) ;
                }
#endif
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


//This routine does a homogeneous space/time cut on interior zoids
template <int N_RANK> template <typename F>
inline void geneity_problem<N_RANK>::homogeneous_space_time_cut_interior(int t0, 
					int t1, grid_info<N_RANK> const & grid, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;

	//assert (projection_zoid) ;
	//assert (projection_zoid->height == lt) ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);

//#ifndef NDEBUG
#if 0
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
			cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
			<< " lt " << lt << endl ;
			cout << "\nzoid " << z->id << " height " << z->height <<
				" num children " << z->num_children <<
				//" num children " << z->children.size() <<
				" num_parents " << z->parents.size() << " geneity " ;
			print_bits(&(z->geneity), sizeof(word_type) * 8);
			int h = z->height ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
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
//#ifndef NDEBUG
#if 0
        homogeneous_space_cut_interior(t0, t1, grid, projection_zoid, f);
#else
        homogeneous_space_cut_interior(t0, t1, grid, f);
#endif
    } 
	else if (lt > dt_recursive_) 
	{
        /* cut into time */
		//assert (projection_zoid->children.size() == 2) ;
		//assert (projection_zoid->children [0]) ;
		//assert (projection_zoid->children [1]) ;
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
//#ifndef NDEBUG
#if 0
		unsigned int index = projection_zoid->children [0] ;
        homogeneous_space_time_cut_interior(t0, t0+halflt, l_son_grid, 
					&(m_zoids [index]), f);
					//projection_zoid->children [0], f);
#else
        homogeneous_space_time_cut_interior(t0, t0+halflt, l_son_grid, f);
#endif

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
//#ifndef NDEBUG
#if 0
		index = projection_zoid->children [1] ;
        homogeneous_space_time_cut_interior(t0+halflt, t1, l_son_grid, 
					&(m_zoids [index]), f);
					//projection_zoid->children [1], f);
#else
        homogeneous_space_time_cut_interior(t0+halflt, t1, l_son_grid, f);
#endif
    } 
	else 
	{
        // base case
		f(t0, t1, grid);
	}
}


//This routine does a homogeneous space/time cut on boundary zoid
template <int N_RANK> template <typename F> 
inline void geneity_problem<N_RANK>::homogeneous_space_time_cut_boundary(int t0, 
					int t1,	grid_info<N_RANK> const & grid, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;

	//assert (projection_zoid) ;
	//assert (projection_zoid->height == lt) ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int lb, thres, tb;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
//#ifndef NDEBUG
#if 0
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
			cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
			<< " lt " << lt << endl ;
			cout << "\nzoid " << z->id << " height " << z->height <<
				//" num children " << z->children.size() <<
				" num children " << z->num_children <<
				" num_parents " << z->parents.size() << " geneity " ;
			print_bits(&(z->geneity), sizeof(word_type) * 8);
			int h = z->height ;
			for (int i = N_RANK - 1 ; i >= 0 ; i--)
			{
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
//#ifndef NDEBUG
#if 0
        if (call_boundary) 
		{
            homogeneous_space_cut_boundary(t0, t1, l_father_grid, 
											projection_zoid, f);
        }
		else
		{
            homogeneous_space_cut_interior(t0, t1, l_father_grid, 
											projection_zoid, f);
		}
#else
        if (call_boundary) 
		{
            homogeneous_space_cut_boundary(t0, t1, l_father_grid, f);
        }
		else
		{
            homogeneous_space_cut_interior(t0, t1, l_father_grid, f);
		}
#endif
    } 
	else if (lt > l_dt_stop)  //time cut
	{
		//cout << "time cut " << endl ;
        // cut into time 
		//assert (projection_zoid->children.size() == 2) ;
		//assert (projection_zoid->children [0]) ;
		//assert (projection_zoid->children [1]) ;
		
        int halflt = lt / 2;
        l_son_grid = l_father_grid;

//#ifndef NDEBUG
#if 0
		unsigned int index = projection_zoid->children [0] ;
        if (call_boundary) {
            homogeneous_space_time_cut_boundary(t0, t0+halflt, l_son_grid, 
								&(m_zoids [index]), f);
								//projection_zoid->children [0], f);
        } else {
            homogeneous_space_time_cut_interior(t0, t0+halflt, l_son_grid, 
								&(m_zoids [index]), f);
								//projection_zoid->children [0], f);
        }
#else
        if (call_boundary) {
            homogeneous_space_time_cut_boundary(t0, t0+halflt, l_son_grid, f);
        } else {
            homogeneous_space_time_cut_interior(t0, t0+halflt, l_son_grid, f);
        }
#endif

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
//#ifndef NDEBUG
#if 0
		index = projection_zoid->children [1] ;
        if (call_boundary) {
            homogeneous_space_time_cut_boundary(t0+halflt, t1, l_son_grid, 
								&(m_zoids [index]), f);
								//projection_zoid->children [1], f);
        } else {
            homogeneous_space_time_cut_interior(t0+halflt, t1, l_son_grid, 
								&(m_zoids [index]), f);
								//projection_zoid->children [1], f);
        }
#else
        if (call_boundary) {
            homogeneous_space_time_cut_boundary(t0+halflt, t1, l_son_grid, f);
        } else {
            homogeneous_space_time_cut_interior(t0+halflt, t1, l_son_grid, f);
        }
#endif
    } 
	else
	{
		// base case
		if (call_boundary) {
            base_case_kernel_boundary(t0, t1, l_father_grid, f);
        } else { 
            f(t0, t1, l_father_grid);
        }
	}
}

//The routine does a heterogeneous space cut for interior zoids. 
template <int N_RANK> template <typename F>
inline void geneity_problem<N_RANK>::heterogeneous_space_cut_interior(int t0, 
				int t1, grid_info<N_RANK> const & grid, 
				zoid_type * projection_zoid, F const & f)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }

	assert (projection_zoid) ;
	int child_index = projection_zoid->start ;
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
					assert(child_index < projection_zoid->end) ;
					assert(child_index >= projection_zoid->start) ;
                    int pos = m_children [child_index] ;
					zoid_type * child = &(m_zoids [pos]) ;
					child_index++ ; //looks like a race
                    heterogeneous_space_time_cut_interior(l_son->t0, l_son->t1, 
											l_son->grid, child, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				assert(child_index < projection_zoid->end) ;
				assert(child_index >= projection_zoid->start) ;
                int pos = m_children [child_index] ;
				zoid_type * child = &(m_zoids [pos]) ;
				child_index++ ;
                if (queue_len_[curr_dep_pointer] == 0)
                    heterogeneous_space_time_cut_interior(l_father->t0, 
							l_father->t1, l_father->grid, child, f);
                else
                    cilk_spawn heterogeneous_space_time_cut_interior(
							l_father->t0, l_father->t1, l_father->grid, 
							child, f);
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

//The routine does a heterogeneous space cut for boundary zoids. 
template <int N_RANK> template <typename F, typename BF>
inline void geneity_problem<N_RANK>::heterogeneous_space_cut_boundary(int t0, 
			int t1,	grid_info<N_RANK> const & grid, 
			zoid_type * projection_zoid, F const & f, BF const & bf)
{
    queue_info *l_father;
    queue_info circular_queue_[2][ALGOR_QUEUE_SIZE];
    int queue_head_[2], queue_tail_[2], queue_len_[2];

    for (int i = 0; i < 2; ++i) {
        queue_head_[i] = queue_tail_[i] = queue_len_[i] = 0;
    }
	
	assert (projection_zoid) ;
	int child_index = projection_zoid->start ;
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
                    assert(child_index < projection_zoid->end) ;
					assert(child_index >= projection_zoid->start) ;
                    int pos = m_children [child_index] ;
					zoid_type * child = &(m_zoids [pos]) ;
					child_index++ ;
                    heterogeneous_space_time_cut_boundary(l_son->t0, l_son->t1, 
											l_son->grid, child, f, bf);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				assert(child_index < projection_zoid->end) ;
				assert(child_index >= projection_zoid->start) ;
				int pos = m_children [child_index] ;
                zoid_type * child = &(m_zoids [pos]) ;
				child_index++ ;
                if (queue_len_[curr_dep_pointer] == 0) {
                    heterogeneous_space_time_cut_boundary(l_father->t0, 
						l_father->t1, l_father->grid, child, f, bf);
                } else {
                    cilk_spawn heterogeneous_space_time_cut_boundary(
									l_father->t0, l_father->t1, l_father->grid,
									child, f, bf);
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


//This routine does a heterogeneous space/time cut on interior zoids
template <int N_RANK> template <typename F>
inline void geneity_problem<N_RANK>::heterogeneous_space_time_cut_interior(
				int t0, int t1, grid_info<N_RANK> const & grid, 
				zoid_type * projection_zoid, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;

	assert (projection_zoid) ;
	assert (projection_zoid->height == lt) ;
#ifdef GENEITY_TEST
	if (__builtin_popcount(projection_zoid->geneity) == 1)
	{
		//zoid is homogeneous
		int index = __builtin_ffs(projection_zoid->geneity) ;
		//cout << "zoid is homogeneous" << endl ;
		//print_bits(&(projection_zoid->geneity), sizeof(word_type) * 8);
//#ifndef NDEBUG
#if 0
		return homogeneous_space_time_cut_interior(t0, t1,	grid, 
								projection_zoid, (*m_clone_array) [index]) ; 
#else
		return homogeneous_space_time_cut_interior(t0, t1,	grid, 
								(*m_clone_array) [index]) ; 
#endif
	}
#endif
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
            int num_children = z->end - z->start ;
			cout << "zoid and proj zoid differ " << endl ;
			cout << "\nzoid " << z->id << " height " << z->height <<
				" num children " << num_children <<
				" num_parents " << z->parents.size() << " geneity " ;
			print_bits(&(z->geneity), sizeof(word_type) * 8);
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
			//assert(0) ;
		}
#endif
        bool cut_lb = (lb < tb);
        thres = (slope_[i] * lt);
        sim_can_cut = SIM_CAN_CUT_I ;
    }
    if (sim_can_cut) 
	{
        /* cut into space */
        heterogeneous_space_cut_interior(t0, t1, grid, projection_zoid, f);
    } 
	else if (lt > dt_recursive_) 
	{
        /* cut into time */
		assert (projection_zoid->end - projection_zoid->start == 2) ;
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
		int index = projection_zoid->start ;
		assert (index < projection_zoid->end) ;
		int pos = m_children [index] ;
		assert (pos < m_num_vertices) ;
        heterogeneous_space_time_cut_interior(t0, t0+halflt, l_son_grid, 
					&(m_zoids [pos]), f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        index++ ;
		assert (index < projection_zoid->end) ;
		pos = m_children [index] ;
		assert (pos < m_num_vertices) ;
        heterogeneous_space_time_cut_interior(t0+halflt, t1, l_son_grid, 
					&(m_zoids [pos]), f);
    }
	else 
	{
#ifdef GENEITY_TEST
        // base case
		f(t0, t1, grid);
#else
		if (__builtin_popcount(projection_zoid->geneity) == 1)
		{
			//zoid is homogeneous
			int index = __builtin_ffs(projection_zoid->geneity) ;
			//cout << "zoid is homogeneous" << endl ;
			(*m_clone_array) [index] (t0, t1, grid);
		}
		else
		{
			f(t0, t1, grid);
		}
#endif
	}
}


//This routine does a heterogeneous space/time cut on boundary zoid
template <int N_RANK> template <typename F, typename BF>
inline void geneity_problem<N_RANK>::heterogeneous_space_time_cut_boundary(
					int t0,
					int t1,	grid_info<N_RANK> const & grid, 
					zoid_type * projection_zoid, F const & f, BF const & bf)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;

	assert (projection_zoid) ;
#ifndef NDEBUG
	if (projection_zoid->height != lt)
	{
		cout << "heights differ : " << projection_zoid->height << "  " <<
				lt << endl ;
		assert (projection_zoid->height == lt) ;
	}
#endif
#ifdef GENEITY_TEST
	if (__builtin_popcount(projection_zoid->geneity) == 1)
	{
		//zoid is homogeneous
		int index = __builtin_ffs(projection_zoid->geneity) ;
		//cout << "zoid is homogeneous" << endl ;
		//print_bits(&(projection_zoid->geneity), sizeof(word_type) * 8);
//#ifndef NDEBUG
#if 0
		return homogeneous_space_time_cut_boundary(t0, t1,	grid, 
								projection_zoid, (*m_clone_array) [index]) ; 
#else
		return homogeneous_space_time_cut_boundary(t0, t1,	grid, 
								(*m_clone_array) [index]) ; 
#endif
	}
#endif
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
            int num_children = z->end - z->start ;
			cout << "zoid and proj zoid differ " << endl ;
			cout << "\nzoid " << z->id << " height " << z->height <<
				" num children " << num_children <<
				" num_parents " << z->parents.size() << " geneity " ;
			print_bits(&(z->geneity), sizeof(word_type) * 8);
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
            heterogeneous_space_cut_boundary(t0, t1, l_father_grid, 
											projection_zoid, f, bf);
        }
		else
		{
            heterogeneous_space_cut_interior(t0, t1, l_father_grid, 
											projection_zoid, f);
		}
    } 
	else if (lt > l_dt_stop)  //time cut
	{
		//cout << "time cut " << endl ;
        // cut into time 
        assert (projection_zoid->end - projection_zoid->start == 2) ;
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
		int index = projection_zoid->start ;
        assert (index < projection_zoid->end) ;
		int pos = m_children [index] ;
		assert (pos < m_num_vertices) ;
        if (call_boundary) {
            heterogeneous_space_time_cut_boundary(t0, t0+halflt, l_son_grid, 
								&(m_zoids [pos]), f, bf);
        } else {
            heterogeneous_space_time_cut_interior(t0, t0+halflt, l_son_grid, 
								&(m_zoids [pos]), f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
		index++ ;
        assert (index < projection_zoid->end) ;
		pos = m_children [index] ;
		assert (pos < m_num_vertices) ;
        if (call_boundary) {
            heterogeneous_space_time_cut_boundary(t0+halflt, t1, l_son_grid, 
								&(m_zoids [pos]), f, bf);
        } else {
            heterogeneous_space_time_cut_interior(t0+halflt, t1, l_son_grid, 
								&(m_zoids [pos]), f);
        }
    } 
	else
	{
#ifdef GENEITY_TEST
		// base case
		if (call_boundary) {
            base_case_kernel_boundary(t0, t1, l_father_grid, bf);
        } else { 
            f(t0, t1, l_father_grid);
        }
#else
		if (__builtin_popcount(projection_zoid->geneity) == 1)
		{
			//zoid is homogeneous
			int index = __builtin_ffs(projection_zoid->geneity) ;
			//cout << "zoid is homogeneous" << endl ;
			if (call_boundary) {
				base_case_kernel_boundary(t0, t1, l_father_grid, 
										(*m_clone_array) [index]);
			} else { 
				(*m_clone_array) [index] (t0, t1, l_father_grid);
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
#endif
