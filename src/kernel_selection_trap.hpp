#ifndef KERNEL_SELECTION_TRAP_HPP
#define KERNEL_SELECTION_TRAP_HPP

#include "kernel_selection.hpp"

#define dx_recursive_boundary_  (m_algo.dx_recursive_boundary_)
#define dx_recursive_ (m_algo.dx_recursive_)
#define dt_recursive_boundary_ (m_algo.dt_recursive_boundary_)
#define dt_recursive_ (m_algo.dt_recursive_)
#define slope_ m_algo.slope_
#define touch_boundary m_algo.touch_boundary
#define phys_length_ m_algo.phys_length_
#define base_case_kernel_boundary m_algo.base_case_kernel_boundary

//do space cuts on grid and call the symbolic_space_time_cut_interior routine
//The routine does a symbolic space cut for interior zoids. 
template <int N_RANK> template <typename F>
inline void kernel_selection<N_RANK>::symbolic_space_cut_interior(int t0, int t1, 
				grid_info<N_RANK> const & grid, F const & f)
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
                    symbolic_space_time_cut_interior(l_son->t0, l_son->t1, 
									l_son->grid, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				symbolic_space_time_cut_interior(l_father->t0, 
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

//do space cuts on grid and call the symbolic_space_time_cut_boundary routine
//The routine does a symbolic space cut for boundary zoids. 
template <int N_RANK> template <typename F>
inline void kernel_selection<N_RANK>::symbolic_space_cut_boundary(int t0, int t1, 
		grid_info<N_RANK> const & grid, F const & f)
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
                    symbolic_space_time_cut_boundary(l_son->t0, l_son->t1, 
									l_son->grid, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
				symbolic_space_time_cut_boundary(l_father->t0, 
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


//This routine does a symbolic space/time cut on interior zoids
template <int N_RANK> template <typename F>
inline void kernel_selection<N_RANK>::symbolic_space_time_cut_interior(int t0, 
		int t1, grid_info<N_RANK> const & grid, F const & f)
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
        symbolic_space_cut_interior(t0, t1, grid, f);
    } 
	else if (lt > dt_recursive_) 
	{
        /* cut into time */
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
        symbolic_space_time_cut_interior(t0, t0+halflt, l_son_grid, f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        symbolic_space_time_cut_interior(t0+halflt, t1, l_son_grid, f);
    } 
	else 
	{
		// base case
		//determine the geneity of the leaf
		word_type geneity = 0 ;
		compute_geneity(lt, grid, geneity, f) ;
		//print_bits(&(z->geneity), sizeof(word_type) * 8) ;
		int index = 0, width = 1 ; 
		int offset = 0 ;
		unsigned long key = 0 ;
	    for (int i = N_RANK-1; i >= 0; --i) 
		{
			int lb, tb;
			lb = (grid.x1[i] - grid.x0[i]);
			tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
			//index = pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) * width + 
			//			index ;
			index = pmod(grid.x0[i], phys_length_ [i]) * width + 
						index ;
			width = phys_length_ [i] ;
			unsigned long dim_key = (unsigned long) lb << num_bits_width | tb ;
			key = key << offset | dim_key ;
			//key = key << offset | lb ;
			offset += num_bits_dim ;
		}
		int kernel = -1 ;
		if (__builtin_popcount(geneity) == 1)
		{
			//zoid is homogeneous
			kernel = __builtin_ffs(geneity) ;
			//cout << "zoid is homogeneous" << endl ;
			//kernel_map [index][centroid] = kernel ;
		}
		else
		{
			//kernel_map [index][centroid] = m_clone_array->clones.size() - 1 ;
			kernel = m_clone_array->clones.size() - 1 ;
		}
		assert (kernel != -1) ;
		assert (kernel < m_clone_array->clones.size()) ;
		//int index = (lt + min_leaf_height - 1) / min_leaf_height ;
		//index = index - (lt == 1) ;
		height_leaf_kernel_table & table = 
							m_arr_height_leaf_kernel_table [index] ;
		//search the hash table with height as key.
		std::pair<hlk_table_iterator, hlk_table_iterator> p = 
											table.equal_range (lt) ;
		if (p.first != p.second)
		{
			//height lt exists
			//assert (p.second - p.first == 1) ; //only one height exists
			hlk_table_iterator start = p.first ;
			assert (start->first == lt) ;
			leaf_kernel_table & table2 = start->second ;
			std::pair<lk_table_iterator, lk_table_iterator> p2 = 
											table2.equal_range (key) ;
			if (p2.first != p2.second)
			{
				//(key,kernel) pair exists.
				//only one (key,kernel) pair may exist
				//assert (p2.second - p2.first == 1) ;
				assert (p2.first->first == key) ;
#ifndef NDEBUG
				zoid_type & z = p2.first->second ;
				if(kernel != (int) z.kernel)
				{
					cout << "Error. Two different leaves have same key" << 
							endl ;
					grid_info <N_RANK> g2 = z.info ;
					int h = z.height ;
					for (int i = N_RANK - 1 ; i >= 0 ; i--)
					{
						cout << " x0 [" << i << "] " << grid.x0 [i] 
						<< " x1 [" << i << "] " << grid.x1 [i] 
						<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
						<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
						<< " lt " << lt << endl ;
						cout << " x0 [" << i << "] " << g2.x0 [i]
						<< " x1 [" << i << "] " << g2.x1 [i]
						<< " x2 [" << i << "] " << g2.x0[i] + g2.dx0[i] * h
						<< " x3 [" << i << "] " << g2.x1[i] + g2.dx1[i] * h
						<< " h " << h << endl ;
					}
					cout << "kernel " << kernel << " z.kernel " <<
						(int) z.kernel << endl ;
					assert (kernel == (int) z.kernel) ;
				}
#else
				if(kernel != p2.first->second)
				{
					cout << "Error. Two different leaves have same key" << 
							endl ;
					cout << "kernel " << kernel << " p2.first->second " <<
						(int) p2.first->second << endl ;
				}
#endif
			}
			else
			{
#ifndef NDEBUG
				zoid_type z ;
				z.kernel = kernel ;
				z.info = grid ;
				z.height = lt ;
				//insert (key,zoid) pair
				table2.insert(std::pair<unsigned long, zoid_type>(key, z)) ;
#else
				//insert (key,kernel) pair
				table2.insert(std::pair<unsigned long, char>(key, (char) kernel)) ;
#endif
			}
		}
		else
		{
			//height lt doesn't exist.
			leaf_kernel_table table2 ;
#ifndef NDEBUG
			zoid_type z ;
			z.kernel = kernel ;
			z.info = grid ;
			z.height = lt ;
			//insert (key,zoid) pair
			table2.insert(std::pair<unsigned long, zoid_type>(key, z)) ;
#else
			//insert (key,kernel) pair
			table2.insert(std::pair<unsigned long, char>(key, (char) kernel)) ;
#endif
			//insert (height, leaf_kernel_table) pair
			table.insert(std::pair<int, leaf_kernel_table>(lt, table2)) ;
		}
	}
}


//This routine does a symbolic space/time cut on boundary zoid
template <int N_RANK> template <typename F>
inline void kernel_selection<N_RANK>::symbolic_space_time_cut_boundary(int t0, 
		int t1,	grid_info<N_RANK> const & grid, F const & f)
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
        //cut into space 
        if (call_boundary) 
		{
            symbolic_space_cut_boundary(t0, t1, l_father_grid, f);
        }
		else
		{
            symbolic_space_cut_interior(t0, t1, l_father_grid, f);
		}
    } 
	else if (lt > l_dt_stop)  //time cut
	{
        // cut into time 
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
        if (call_boundary) {
            symbolic_space_time_cut_boundary(t0, t0+halflt, l_son_grid, f);
        } else {
            symbolic_space_time_cut_interior(t0, t0+halflt, l_son_grid, f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            symbolic_space_time_cut_boundary(t0+halflt, t1, l_son_grid, f);
        } else {
            symbolic_space_time_cut_interior(t0+halflt, t1, l_son_grid, f);
        }
    } 
	else
	{
		// base case
		//determine the geneity of the leaf
		word_type geneity = 0 ;
		compute_geneity(lt, l_father_grid, geneity, f) ;
		//print_bits(&(z->geneity), sizeof(word_type) * 8) ;
		int index = 0, width = 1 ; 
		int offset = 0 ;
		unsigned long key = 0 ;
	    for (int i = N_RANK-1; i >= 0; --i) 
		{
			int lb, tb;
			lb = (grid.x1[i] - grid.x0[i]);
			tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
			//index = pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) * width + 
			//			index ;
			index = pmod(grid.x0[i], phys_length_ [i]) * width + 
						index ;
			width = phys_length_ [i] ;
			unsigned long dim_key = (unsigned long) lb << num_bits_width | tb ;
			key = key << offset | dim_key ;
			//key = key << offset | lb ;
			offset += num_bits_dim ;
		}
		int kernel = -1 ;
		if (__builtin_popcount(geneity) == 1)
		{
			//zoid is homogeneous
			kernel = __builtin_ffs(geneity) ;
			//cout << "zoid is homogeneous" << endl ;
			//kernel_map [index][centroid] = kernel ;
		}
		else
		{
			//kernel_map [index][centroid] = m_clone_array->clones.size() - 1 ;
			kernel = m_clone_array->clones.size() - 1 ;
		}
		assert (kernel != -1) ;
		assert (kernel < m_clone_array->clones.size()) ;
		//int index = (lt + min_leaf_height - 1) / min_leaf_height ;
		//index = index - (lt == 1) ;
		height_leaf_kernel_table & table = 
							m_arr_height_leaf_kernel_table [index] ;
		//search the hash table with height as key.
		std::pair<hlk_table_iterator, hlk_table_iterator> p = 
											table.equal_range (lt) ;
		if (p.first != p.second)
		{
			//height lt exists
			//assert (p.second - p.first == 1) ; //only one height exists
			hlk_table_iterator start = p.first ;
			assert (start->first == lt) ;
			leaf_kernel_table & table2 = start->second ;
			std::pair<lk_table_iterator, lk_table_iterator> p2 = 
											table2.equal_range (key) ;
			if (p2.first != p2.second)
			{
				//(key,kernel) pair exists.
				//only one (key,kernel) pair may exist
				//assert (p2.second - p2.first == 1) ;
				assert (p2.first->first == key) ;
#ifndef NDEBUG
				zoid_type & z = p2.first->second ;
				if(kernel != (int) z.kernel)
				{
					cout << "Error. Two different leaves have same key" << 
							endl ;
					cout << "index " << index << " key " << key << endl ;
					grid_info <N_RANK> g2 = z.info ;
					int h = z.height ;
					offset = 0 ;
					key = 0 ;
					width = 1 ;
					index = 0 ;
					for (int i = N_RANK - 1 ; i >= 0 ; i--)
					{
						cout << " x0 [" << i << "] " << grid.x0 [i] 
						<< " x1 [" << i << "] " << grid.x1 [i] 
						<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * lt
						<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * lt
						<< " lt " << lt << endl ;
						cout << " x0 [" << i << "] " << g2.x0 [i]
						<< " x1 [" << i << "] " << g2.x1 [i]
						<< " x2 [" << i << "] " << g2.x0[i] + g2.dx0[i] * h
						<< " x3 [" << i << "] " << g2.x1[i] + g2.dx1[i] * h
						<< " h " << h << endl ;

						int lb = (g2.x1[i] - g2.x0[i]);
						int tb = (g2.x1[i] + g2.dx1[i] * lt - g2.x0[i] - g2.dx0[i] * lt);
						cout << " pmod(g2.x0[i], phys_length_ [i]) " << 
								pmod(g2.x0[i], phys_length_ [i]) << endl ;
						//index = pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) * width + index ;
						index = pmod(g2.x0[i], phys_length_ [i]) * width + 
								index ;
						width = phys_length_ [i] ;
						cout << "index " << index << " width " << width << endl ;
						unsigned long dim_key = (unsigned long) lb << num_bits_width | tb ;
						key = key << offset | dim_key ;
						//key = key << offset | lb ;
						offset += num_bits_dim ;
					}
					cout << "index " << index << " key " << key << endl ;
					cout << "kernel " << kernel << " z.kernel " <<
						(int) z.kernel << endl ;
					assert (kernel == (int) z.kernel) ;
				}
#else
				if(kernel != p2.first->second)
				{
					cout << "Error. Two different leaves have same key" << 
							endl ;
					cout << "kernel " << kernel << " p2.first->second " <<
						(int) p2.first->second << endl ;
				}
#endif
			}
			else
			{
#ifndef NDEBUG
				zoid_type z ;
				z.kernel = kernel ;
				z.info = grid ;
				z.height = lt ;
				//insert (key,zoid) pair
				table2.insert(std::pair<unsigned long, zoid_type>(key, z)) ;
#else
				//insert (key,kernel) pair
				table2.insert(std::pair<unsigned long, char>(key, (char) kernel)) ;
#endif
			}
		}
		else
		{
			//height lt doesn't exist.
			leaf_kernel_table table2 ;
#ifndef NDEBUG
			zoid_type z ;
			z.kernel = kernel ;
			z.info = grid ;
			z.height = lt ;
			//insert (key,zoid) pair
			table2.insert(std::pair<unsigned long, zoid_type>(key, z)) ;
#else
			//insert (key,kernel) pair
			table2.insert(std::pair<unsigned long, char>(key, (char) kernel)) ;
#endif
			//insert (height, leaf_kernel_table) pair
			table.insert(std::pair<int, leaf_kernel_table>(lt, table2)) ;
		}
	}
}

//The routine does a heterogeneous space cut for interior zoids. 
template <int N_RANK> template <typename F>
inline void kernel_selection<N_RANK>::heterogeneous_space_cut_interior(int t0, 
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
                    heterogeneous_space_time_cut_interior(l_son->t0, l_son->t1, 
											l_son->grid, f);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0)
                    heterogeneous_space_time_cut_interior(l_father->t0, 
							l_father->t1, l_father->grid, f);
                else
                    cilk_spawn heterogeneous_space_time_cut_interior(
							l_father->t0, l_father->t1, l_father->grid, f);
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
inline void kernel_selection<N_RANK>::heterogeneous_space_cut_boundary(int t0, 
			int t1,	grid_info<N_RANK> const & grid, 
			F const & f, BF const & bf)
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
                    heterogeneous_space_time_cut_boundary(l_son->t0, l_son->t1, 
											l_son->grid, f, bf);
                } // end cilk_for 
                queue_head_[curr_dep_pointer] = queue_tail_[curr_dep_pointer] = 0;
                queue_len_[curr_dep_pointer] = 0;
#else
                // use cilk_spawn to spawn all the sub-grid 
                pop_queue(curr_dep_pointer);
                if (queue_len_[curr_dep_pointer] == 0) {
                    heterogeneous_space_time_cut_boundary(l_father->t0, 
						l_father->t1, l_father->grid, f, bf);
                } else {
                    cilk_spawn heterogeneous_space_time_cut_boundary(
									l_father->t0, l_father->t1, l_father->grid,
									f, bf);
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
inline void kernel_selection<N_RANK>::heterogeneous_space_time_cut_interior(int t0,
				int t1, grid_info<N_RANK> const & grid, F const & f)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;
	int p_lb [N_RANK], p_tb [N_RANK] ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int thres ;
        int lb = p_lb [i] = (grid.x1[i] - grid.x0[i]);
        int tb = p_tb [i] = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
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
        heterogeneous_space_cut_interior(t0, t1, grid, f);
    } 
	else if (lt > dt_recursive_) 
	{
        /* cut into time */
        assert(lt > dt_recursive_);
        int halflt = lt / 2;
        l_son_grid = grid;
        heterogeneous_space_time_cut_interior(t0, t0+halflt, l_son_grid, f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
        heterogeneous_space_time_cut_interior(t0+halflt, t1, l_son_grid, f);
    }
	else 
	{
		//continue from here
		// base case
		int index = 0, width = 1 ; 
		int offset = 0 ;
		unsigned long key = 0 ;
	    for (int i = N_RANK-1; i >= 0; --i) 
		{
			//index = pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) * width + 
			//			index ;
			index = pmod(grid.x0[i], phys_length_ [i]) * width + 
						index ;
			width = phys_length_ [i] ;
			unsigned long dim_key = (unsigned long) p_lb [i] << num_bits_width 
									| p_tb [i] ;
			key = key << offset | dim_key ;
			//key = key << offset | p_lb [i] ;
			offset += num_bits_dim ;
		}
		int kernel = -1 ;
		//int index = (lt + min_leaf_height - 1) / min_leaf_height ;
		//index = index - (lt == 1) ;
		height_leaf_kernel_table & table = 
							m_arr_height_leaf_kernel_table [index] ;
		//search the hash table with height as key.
		std::pair<hlk_table_iterator, hlk_table_iterator> p = 
											table.equal_range (lt) ;
		assert (p.first != p.second) ;
		hlk_table_iterator start = p.first ;
		assert (start->first == lt) ;
		leaf_kernel_table & table2 = start->second ;
		std::pair<lk_table_iterator, lk_table_iterator> p2 = 
										table2.equal_range (key) ;
		assert (p2.first != p2.second) ;
		assert (p2.first->first == key) ;
#ifndef NDEBUG
		kernel = p2.first->second.kernel ;
#else
		kernel = p2.first->second ;
#endif
		assert (kernel != -1) ;
		(*m_clone_array) [kernel] (t0, t1, grid);
	}
}


//This routine does a heterogeneous space/time cut on boundary zoid
template <int N_RANK> template <typename F, typename BF>
inline void kernel_selection<N_RANK>::heterogeneous_space_time_cut_boundary(int t0,
					int t1,	grid_info<N_RANK> const & grid, 
					F const & f, BF const & bf)
{
    const int lt = t1 - t0;
    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;
	
	int p_lb [N_RANK], p_tb [N_RANK] ;
    for (int i = N_RANK-1; i >= 0; --i) {
        int thres ;
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
        int lb = p_lb [i] = (grid.x1[i] - grid.x0[i]);
        int tb = p_tb [i] = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
        thres = slope_[i] * lt ;
        bool cut_lb = (lb < tb) ;

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
            heterogeneous_space_cut_boundary(t0, t1, l_father_grid, f, bf);
        }
		else
		{
            heterogeneous_space_cut_interior(t0, t1, l_father_grid, f);
		}
    } 
	else if (lt > l_dt_stop)  //time cut
	{
        // cut into time 
		
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
        if (call_boundary) {
            heterogeneous_space_time_cut_boundary(t0, t0+halflt, l_son_grid, 
								f, bf);
        } else {
            heterogeneous_space_time_cut_interior(t0, t0+halflt, l_son_grid, 
								f);
        }

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = l_father_grid.dx0[i];
            l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = l_father_grid.dx1[i];
        }
        if (call_boundary) {
            heterogeneous_space_time_cut_boundary(t0+halflt, t1, l_son_grid, 
								f, bf);
        } else {
            heterogeneous_space_time_cut_interior(t0+halflt, t1, l_son_grid, 
								f);
        }
    } 
	else
	{
		// base case
		int index = 0, width = 1 ; 
		int offset = 0 ;
		unsigned long key = 0 ;
	    for (int i = N_RANK-1; i >= 0; --i) 
		{
			//index = pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) * width + 
			//			index ;
			index = pmod(grid.x0[i], phys_length_ [i]) * width + 
						index ;
			width = phys_length_ [i] ;
			unsigned long dim_key = (unsigned long) p_lb [i] << num_bits_width 
									| p_tb [i] ;
			key = key << offset | dim_key ;
			//key = key << offset | p_lb [i] ;
			offset += num_bits_dim ;
		}
		int kernel = -1 ;
		//int index = (lt + min_leaf_height - 1) / min_leaf_height ;
		//index = index - (lt == 1) ;
		height_leaf_kernel_table & table = 
							m_arr_height_leaf_kernel_table [index] ;
		//search the hash table with height as key.
		std::pair<hlk_table_iterator, hlk_table_iterator> p = 
											table.equal_range (lt) ;
		assert (p.first != p.second) ;
		hlk_table_iterator start = p.first ;
		assert (start->first == lt) ;
		leaf_kernel_table & table2 = start->second ;
		std::pair<lk_table_iterator, lk_table_iterator> p2 = 
										table2.equal_range (key) ;
		assert (p2.first != p2.second) ;
		assert (p2.first->first == key) ;
#ifndef NDEBUG
		kernel = p2.first->second.kernel ;
#else
		kernel = p2.first->second ;
#endif
		assert (kernel != -1) ;
		
		//cout << "kernel index " << kernel << endl ;
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
