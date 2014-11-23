/*
 * ============================================================================
 *       Filename:  auto_tuning_trap_seq_space_cuts.hpp
 *    Description:  Has routines for autotuning trapezoidal divide-and-conquer
 *					algorithm, which performs non hyperspace cuts.
 *        Created:  10/02/2013
 *         Author:  Eka Palamadai, epn@mit.edu
 * ============================================================================
 */
#ifndef AUTO_TUNING_TRAP_SEQ_SPACE_CUTS_HPP
#define AUTO_TUNING_TRAP_SEQ_SPACE_CUTS_HPP

#include "auto_tuning_homogeneous_header.hpp"

#define dx_recursive_boundary_  (m_algo.dx_recursive_boundary_)
#define dx_recursive_ (m_algo.dx_recursive_)
#define dt_recursive_boundary_ (m_algo.dt_recursive_boundary_)
#define dt_recursive_ (m_algo.dt_recursive_)
#define slope_ m_algo.slope_
#define touch_boundary m_algo.touch_boundary
#define phys_length_ m_algo.phys_length_
#define base_case_kernel_boundary m_algo.base_case_kernel_boundary
#define base_case_kernel_boundary_tune m_algo.base_case_kernel_boundary_tune
#define uub_boundary m_algo.uub_boundary
#define ulb_boundary m_algo.ulb_boundary
#define lub_boundary m_algo.lub_boundary

#ifdef STOP_TUNING_EARLY
#define SPACE_CUT_INTERIOR(k_, x0_, x1_, dx0_, dx1_) \
if (elapsed_time + ltime + ctime < best_time) \
{ \
	stopwatch_stop(ptr) ; \
	stopwatch_get_elapsed_time(ptr, t) ; \
	elapsed_time += t ; \
	l_son_grid.x0[i] = x0_; \
	l_son_grid.dx0[i] = dx0_; \
	l_son_grid.x1[i] = x1_ ; \
	l_son_grid.dx1[i] = dx1_ ; \
	symbolic_trap_space_time_cut_interior(t0, t1, l_son_grid, \
			index, k_, ltime, ctime, f, best_time) ; \
	stopwatch_start(ptr) ; \
}

#define LOOP_INTERIOR \
{ \
	cilk_spawn loop_interior(t0, t1, grid, f, loop_time) ; \
	timespec ts ; \
	ts.tv_nsec = best_time ; \
	ts.tv_sec = 0 ; \
	if (nanosleep(&ts, 0) == 0) \
	{ \
		t1 = 0 ; \
	} \
	cilk_sync ; \
}

#else
#define SPACE_CUT_INTERIOR(k_, x0_, x1_, dx0_, dx1_) \
{ \
	l_son_grid.x0[i] = x0_; \
	l_son_grid.dx0[i] = dx0_; \
	l_son_grid.x1[i] = x1_ ; \
	l_son_grid.dx1[i] = dx1_ ; \
	symbolic_trap_space_time_cut_interior(t0, t1, l_son_grid, \
			index, k_, ltime, ctime, f, best_time) ; \
}

#define LOOP_INTERIOR loop_interior(t0, t1, grid, f, loop_time)
#endif

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::loop_interior(int t0, VOLATILE_INT t1, 
	grid_info<N_RANK> const & grid, F const & f, time_type & loop_time)
{
	stopwatch * ptr = &m_stopwatch ;
	time_type t = 0 ;
	stopwatch_start(ptr) ;
	f(t0, t1, grid) ;
	stopwatch_stop(ptr) ;
	stopwatch_get_elapsed_time(ptr, t) ;
	loop_time = t ;
	assert (loop_time >= 0) ;
}

template <int N_RANK> template <typename F, typename BF>
inline void auto_tune<N_RANK>::loop_boundary(int t0, VOLATILE_INT t1, 
	grid_info<N_RANK> const & grid, F const & f, BF const & bf, 
	time_type & loop_time, bool call_boundary)
{
	stopwatch * ptr = &m_stopwatch ;
	time_type t = 0 ;
	stopwatch_start(ptr) ;
	if (call_boundary)
	{
#ifdef STOP_TUNING_EARLY
		base_case_kernel_boundary_tune(t0, t1, grid, bf);
#else
		base_case_kernel_boundary(t0, t1, grid, bf);
#endif
	} 
	else 
	{
		f(t0, t1, grid);
	}
	stopwatch_stop(ptr) ;
	stopwatch_get_elapsed_time(ptr, t) ;
	loop_time = t ;
	assert (loop_time >= 0) ;
}

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::symbolic_trap_space_time_cut_interior(
		int t0, int t1, grid_info<N_RANK> const & grid, 
		unsigned long parent_index, int child_index, time_type & linkage_time,
		time_type & child_time, F const & f, time_type best_time)
{
	time_type t = 0 ;
	stopwatch * ptr = &m_stopwatch ;
    const int lt = t1 - t0 ;
	stopwatch_stop(ptr) ;
	stopwatch_get_elapsed_time(ptr, t) ;
	linkage_time += t ;

    bool sim_can_cut = false;
    grid_info<N_RANK> l_son_grid;

	int centroid = 0, width = 1 ; 
	int num_subzoids = 2 ;
	
	unsigned long key = 0 ;
	decision_type decision = 0 ;
	unsigned long num_grid_points = 1 ;
	//bool empty_zoid = false ;
	vector<int> space_cut_dims ;
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
			<< " lt " << lt << endl ; */
		
		assert (lb >= 0) ;
		assert (tb >= 0) ;
		assert (centroid >= 0) ;
		width *= phys_length_ [i] ;

		key <<= (num_bits_dim - 4) ;
		key |= lb ;
		key <<= 2 ;
		int dx0 = grid.dx0[i] ;
		if (dx0 > 0)
		{
			key |= 1 ;
		}
		else if (dx0 < 0)
		{
			key |= 2 ;
		}
		else
		{
			key |= 0 ; //this is a no-op
		}
		key <<= 2 ; //shift by 2 bits
		
		int dx1 = grid.dx1[i] ;
		if (dx1 > 0)
		{
			key |= 1 ;
		}
		else if (dx1 < 0)
		{
			key |= 2 ;
		}
		else
		{
			key |= 0 ; //this is a no-op
		}
		
        thres = slope_[i] * lt ;
		unsigned long short_side = min(lb, tb) ;
		if (short_side >= 2 * thres)
		{
			sim_can_cut = true ;
			num_subzoids = 3 ;
			space_cut_dims.push_back(i) ;
		}
		num_grid_points *= ((lb + tb) / 2) ;
		/*if (lt == 1 && lb == 0)
		{
			empty_zoid = true ;
		}*/
    }
#ifdef PERMUTE
	if (sim_can_cut)
	{
		int size = space_cut_dims.size() ;
		bool flag [size] ;
		int permuted_dims [size] ;
		for (int i = 0 ; i < size ; i++)
		{
			flag [i] = false ;
		}
		int end = size - 1 ;
		bool cut_contiguous_dim = false ;
		if (space_cut_dims [end] == 0)
		{
			flag [end] = true ; //mark index 'end' as taken
			permuted_dims [end] = 0 ;
			cut_contiguous_dim = true ;
			end = size - 2 ;
		}
		//permute the non-contiguous dimensions
		for (int i = end ; i >= 0 ; i--)
		{
			while (1)
			{
				int index = rand () % size ;
				if (flag [index] == false)
				{
					permuted_dims [i] = space_cut_dims [index] ;
					//mark index as taken
					flag [index] = true ;
					break ;
				}
			}
		}
		for (int i = size - 1 ; i >= 0 ; i--)
		{
			space_cut_dims [i] = permuted_dims [i] ;
		}
#ifndef NDEBUG
		for (int i = 0 ; i < size ; i++)
		{
			assert (flag [i]) ;
		}
		if (cut_contiguous_dim)
		{
			assert (space_cut_dims [size - 1] == 0) ; 
		}
#endif
	}
#endif
	num_grid_points *= lt ;
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
		//a zoid with the projection already exists. return
		//Add the projected time of the zoid
		child_time += z.time ;
		//start measuring linkage time
		stopwatch_start(ptr) ;
		return ;
	}
	/*if (empty_zoid)
	{
		sim_can_cut = false ;
	}*/
	bool time_cut = false ;
	bool divide_and_conquer = false ;
	time_type time_cut_elapsed_time = 0, space_cut_elapsed_time = LONG_MAX ; 
	time_type max_loop_time = 0 ;
#ifdef FIXED_TIME_CUT
	//cut in time only when space cut is not possible
	if (lt > dt_recursive_ && ! sim_can_cut) 
#else
	if (lt > dt_recursive_)
#endif
	{
		divide_and_conquer = true ;
		time_cut = true ;
		m_zoids [index].set_capacity(num_subzoids) ;
		m_zoids [index].resize_children(2) ;
	
		time_type ltime = 0, ctime = 0 ;
        /* cut into time */
		stopwatch_start(ptr) ;
        int halflt = lt / 2;
        l_son_grid = grid;
        symbolic_trap_space_time_cut_interior(t0, t0+halflt, l_son_grid, 
				index, 0, ltime, ctime, f, best_time) ;
		for (int i = 0; i < N_RANK; ++i) {
			l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
			l_son_grid.dx0[i] = grid.dx0[i];
			l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
			l_son_grid.dx1[i] = grid.dx1[i];
		}
#ifdef STOP_TUNING_EARLY
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t) ;
		time_cut_elapsed_time += t ;

		if (t + ltime + ctime < best_time)
		{
			stopwatch_start(ptr) ;
        	symbolic_trap_space_time_cut_interior(t0+halflt, t1, l_son_grid, 
				index, 1, ltime, ctime, f, best_time) ;
		}
#else
        symbolic_trap_space_time_cut_interior(t0+halflt, t1, l_son_grid, 
			index, 1, ltime, ctime, f, best_time) ;
#endif
		//measure the remaining function call overhead.
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t) ;
		//ltime is division cost + partial function call overhead.
		time_cut_elapsed_time += t + ltime + ctime ;
		assert (ltime >= 0) ;
		assert (ctime >= 0) ;
		assert (t >= 0) ;
		assert (m_zoids [index].num_children == 2) ;
		best_time = min (best_time, time_cut_elapsed_time) ;
#ifndef NDEBUG
		m_zoids [index].ttime = time_cut_elapsed_time ;
#endif
		
#if defined (MEASURE_COLD_MISS) || defined (SUBSUMPTION3) || \
	defined (SUBSUMPTION_TIME)
		time_type sum_loop_time = 0 ;
		time_type cache_penalty_children = 0 ;
		for (int i = 0 ; i < 2 ; i++)
		{
			unsigned long child_index = m_zoids [index].children [i] ;
#ifdef SUBSUMPTION_TIME
			max_loop_time = max(m_zoids [child_index].loop_time, max_loop_time);
#endif
#ifdef SUBSUMPTION3
			sum_loop_time += m_zoids [child_index].loop_time ;
#endif
#ifdef MEASURE_COLD_MISS
			cache_penalty_children += m_zoids [child_index].cache_penalty_time ;
#endif
		}
		assert (sum_loop_time >= 0) ;
		assert (cache_penalty_children >= 0) ;
#ifdef SUBSUMPTION3
		max_loop_time = max (max_loop_time, sum_loop_time) ;
#endif
#ifdef MEASURE_COLD_MISS
		m_zoids [index].cache_penalty_time = max(cache_penalty_children, 
										m_zoids [index].cache_penalty_time) ;
#endif
#endif
    }
	zoid_type bak ;
    if (sim_can_cut) 
	{
		divide_and_conquer = true ; 
		if (time_cut)
		{
			//back up the time cut children data
			bak = m_zoids [index] ;
			assert(bak.num_children == 2) ;
		}
		m_zoids [index].set_capacity(3) ;
		m_zoids [index].resize_children(3) ;
	}

	zoid_type bak2 ;
	int best_case = -1 ; 
    for (int j = 0 ; j < space_cut_dims.size() ; j++) {
		int i = space_cut_dims [j] ;
		time_type ltime = 0, ctime = 0 ;
		time_type elapsed_time = 0 ;
	
		//measure the divide time
		stopwatch_start(ptr) ;
        unsigned long lb, tb;
        int thres ;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
        thres = slope_[i] * lt ;
        bool cut_lb = (lb < tb);
		if (cut_lb) {
			assert(lb != phys_length_[i] || grid.dx0[i] != 0 || grid.dx1[i] != 0);
			const int mid = lb/2;
			l_son_grid = grid;
			const int l_start = grid.x0[i];
			const int l_end = grid.x1[i];

			SPACE_CUT_INTERIOR(0, l_start + mid - thres, l_start + mid + thres, slope_[i], -slope_[i]) ;
			SPACE_CUT_INTERIOR(1, l_start, l_start + mid - thres, grid.dx0[i], slope_[i]) ;
			SPACE_CUT_INTERIOR(2, l_start + mid + thres, l_end, -slope_[i], grid.dx1[i]) ;
		} // end if (cut_lb) 
		else { // cut_tb 
				const int mid = tb/2;
				l_son_grid = grid;
				const int l_start = grid.x0[i];
				const int l_end = grid.x1[i];
				const int ul_start = grid.x0[i] + grid.dx0[i] * lt;

				SPACE_CUT_INTERIOR(0, l_start, ul_start + mid, grid.dx0[i], -slope_[i]) ;
				SPACE_CUT_INTERIOR(1, ul_start + mid, l_end, slope_[i], grid.dx1[i]) ;
				SPACE_CUT_INTERIOR(2, ul_start + mid, ul_start + mid, -slope_[i], slope_[i]) ;
		} // end if (cut_tb) 
		
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t) ;
		assert (ltime >= 0) ;
		assert (t >= 0) ;
		assert (ctime >= 0) ;
		//ltime is division cost + partial function call overhead.
		elapsed_time += ltime + ctime + t ;
		if (elapsed_time < space_cut_elapsed_time)
		{
			space_cut_elapsed_time = elapsed_time ;
			best_case = i ;
			//back up the zoid with its children.
			bak2 = m_zoids [index] ;
			best_time = min(best_time, elapsed_time) ;
		}
		
		assert (m_zoids [index].num_children == 3) ;
#if defined (MEASURE_COLD_MISS) || defined (SUBSUMPTION3) || \
	defined (SUBSUMPTION_TIME)
		time_type sum_loop_time = 0 ;
		time_type cache_penalty_children = 0 ;
		for (int k = 0 ; k < 3 ; k++)
		{
			unsigned long child_index = m_zoids [index].children [k] ;
#ifdef SUBSUMPTION_TIME
			max_loop_time = max(m_zoids [child_index].loop_time, max_loop_time);
#endif
#ifdef SUBSUMPTION3
			sum_loop_time += m_zoids [child_index].loop_time ;
#endif
#ifdef MEASURE_COLD_MISS
			cache_penalty_children += m_zoids [child_index].cache_penalty_time ;
#endif
		}
		assert (sum_loop_time >= 0) ;
		assert (cache_penalty_children >= 0) ;
#ifdef SUBSUMPTION3
		max_loop_time = max (max_loop_time, sum_loop_time) ;
#endif
#ifdef MEASURE_COLD_MISS
		m_zoids [index].cache_penalty_time = max(cache_penalty_children, 
									m_zoids [index].cache_penalty_time) ;
#endif
#endif

#ifdef FIXED_SPACE_CUT
		//cut only in one dimension.
		break ;
#endif
    }
	if (sim_can_cut) 
	{
		assert (space_cut_elapsed_time >= 0) ;
		assert (best_case >= 0 && best_case < N_RANK) ;
		//set the decision with the best case found
		decision = (decision_type) 1 << (best_case + 1) ;
		//restore the back up.
		m_zoids [index] = bak2 ;
		assert (m_zoids [index].num_children == 3) ;
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
	for (int i = 0 ; i < m_zoids [index].num_children ; i++)
	{
		unsigned long child_index = m_zoids [index].children [i] ;
		decision_type d = m_zoids [child_index].decision ;
		if (d & m_space_cut_mask || d & 1)
		{
			child_divides = true ;
			//test if grand child divided
			if ((int) m_zoids [child_index].num_level_divide >= 1) 
			{
				force_divide = true ;
				m_zoids [index].num_level_divide = 
							m_zoids [child_index].num_level_divide ;
			}
		}
	}
#endif

#if defined (SUBSUMPTION_TIME) || defined (SUBSUMPTION3)
	if (divide_and_conquer && divide_and_conquer_time < max_loop_time)
	{
		force_divide = true ;
	}
#endif
	if (num_grid_points < 1000 || //empty_zoid || 
		divide_and_conquer_time < ptr->measure_time)
	{
		force_divide = false ;
	}

    //base case
	//The subsumption in time works as follows.
	//suppose loop_time(z) >= loop_time(z'), z' \in tree(z)
	//if divide_and_conquer_time(z) < max_{z' \in tree(z)} loop_time(z')
	//	then avoid computing loop_time(z).
	//else compute loop_time(z).
	//Other subsumption heuristics work similarly to avoid computing 
	//loop_time(z).
	time_type loop_time = LONG_MAX ;
	if (force_divide || lt > (m_initial_height + 1) / 2 )
	{
		//do not compute loop_time.
		m_zoids [index].decision = decision ;
		m_zoids [index].time = divide_and_conquer_time ;
	}
	else if (! divide_and_conquer)
	{
		//determine the looping time on the zoid
		time_type t1_, t2_ ;
		stopwatch_start(ptr) ;
		f(t0, t1, grid);
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t1_) ;

		stopwatch_start(ptr) ;
		f(t0, t1, grid);
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t2_) ;
		loop_time = min (t1_, t2_) ;
		assert (loop_time >= 0) ;
#ifdef MEASURE_COLD_MISS
		m_zoids [index].cache_penalty_time = max(t1_ - t2_, (time_type) 0) ;
#endif
#ifndef NDEBUG
		m_zoids [index].ltime = loop_time ;
#endif
		//set a flag to indicate that we looped on z.
		m_zoids [index].decision = (decision_type) 1 << 
					  (zoid_type::NUM_BITS_DECISION - 2) ;
		m_zoids [index].time = loop_time ;
	}
	else 
	{
		//determine the looping time on the zoid
		t = 0 ;
		//cout << "best time " << best_time << endl ;
		//stopwatch_start(ptr) ;
		LOOP_INTERIOR ;
		//stopwatch_stop(ptr) ;
		//stopwatch_get_elapsed_time(ptr, t) ;
		//loop_time = t ;
		//assert (loop_time >= 0) ;
#ifndef NDEBUG
		m_zoids [index].ltime = loop_time ;
#endif
		//max_loop_time of z is its loop time.
		max_loop_time = loop_time ;
		//compare divide and conquer time with loop_time 
		if(divide_and_conquer_time < zoid_type::FUZZ * loop_time)
		{
			//choose divide and conquer
			m_zoids [index].decision = decision ;
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
		}
	}
#if defined (SUBSUMPTION_TIME) || defined (SUBSUMPTION3)
	m_zoids [index].loop_time = max_loop_time ;
#endif
	
	child_time += m_zoids [index].time ;
	//start measuring linkage time
	stopwatch_start(ptr) ;
}

#define SPACE_CUT_BOUNDARY(k) \
if (call_boundary) \
{ \
	symbolic_trap_space_time_cut_boundary(t0, t1, l_son_grid, \
			index, k, ltime, ctime, f, bf, best_time) ; \
} \
else \
{ \
	symbolic_trap_space_time_cut_interior(t0, t1, l_son_grid, \
			index, k, ltime, ctime, f, best_time) ; \
}

#ifdef STOP_TUNING_EARLY
#define SPACE_CUT_TUNE_BOUNDARY(k_, x0_, x1_, dx0_, dx1_) \
if (elapsed_time + ltime + ctime < best_time) \
{ \
	stopwatch_stop(ptr) ; \
	stopwatch_get_elapsed_time(ptr, t) ; \
	elapsed_time += t ; \
	l_son_grid.x0[i] = x0_; \
	l_son_grid.dx0[i] = dx0_; \
	l_son_grid.x1[i] = x1_ ; \
	l_son_grid.dx1[i] = dx1_ ; \
	SPACE_CUT_BOUNDARY(k_) ; \
	stopwatch_start(ptr) ; \
}

#define LOOP_BOUNDARY \
{ \
	cilk_spawn loop_boundary(t0, t1, grid, f, bf, loop_time, call_boundary) ; \
	timespec ts ; \
	ts.tv_nsec = best_time ; \
	ts.tv_sec = 0 ; \
	if (nanosleep(&ts, 0) == 0) \
	{ \
		t1 = 0 ; \
	} \
	cilk_sync ; \
}
#else
#define SPACE_CUT_TUNE_BOUNDARY(k_, x0_, x1_, dx0_, dx1_) \
	l_son_grid.x0[i] = x0_; \
	l_son_grid.dx0[i] = dx0_; \
	l_son_grid.x1[i] = x1_ ; \
	l_son_grid.dx1[i] = dx1_ ; \
	SPACE_CUT_BOUNDARY(k_) ; 

#define LOOP_BOUNDARY loop_boundary(t0, t1, grid, f, bf, loop_time, call_boundary)
#endif

template <int N_RANK> template <typename F, typename BF>
inline void auto_tune<N_RANK>::symbolic_trap_space_time_cut_boundary(
	int t0, int t1,	grid_info<N_RANK> const & grid, 
	unsigned long parent_index, int child_index, time_type & linkage_time, 
	time_type & child_time, F const & f, BF const & bf, time_type best_time)
{
	time_type t = 0 ;
	stopwatch * ptr = &m_stopwatch ;
    const int lt = t1 - t0;
	stopwatch_stop(ptr) ;
	stopwatch_get_elapsed_time(ptr, t) ;
	linkage_time += t ;

    bool sim_can_cut = false, call_boundary = false;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
    int l_dt_stop;

	int centroid = 0, width = 1 ; 
	unsigned long key = 0 ;
	decision_type decision = 0 ;
	int num_subzoids = 2 ;

	int dim_touching_bdry = 0 ;
	int centroid_dim_touching_bdry = 0, width_dim_touching_bdry = 1 ;
	unsigned long num_grid_points = 1 ;
	//bool empty_zoid = false ;
	vector<int> space_cut_dims ;

	//measure time to adjust boundaries of zoid
	stopwatch_start(ptr) ;
    for (int i = N_RANK-1; i >= 0; --i) {
        bool l_touch_boundary = touch_boundary(i, lt, l_father_grid);
        call_boundary |= l_touch_boundary;
	}
	stopwatch_stop(ptr) ;
	stopwatch_get_elapsed_time(ptr, t) ;
	time_type bdry_time = t ;

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
			<< " lt " << lt << endl ; */
        thres = slope_[i] * lt ;
		if (l_touch_boundary)
		{
			dim_touching_bdry |= 1 << i ; 
			centroid_dim_touching_bdry = 
				pmod(grid.x0[i] + (lb >> 1), phys_length_ [i]) * 
				width_dim_touching_bdry + centroid_dim_touching_bdry ;
			assert (centroid_dim_touching_bdry >= 0) ;
			width_dim_touching_bdry *= phys_length_ [i] ;
		}
		unsigned long short_side = min (lb, tb) ;
		if (short_side >= 2 * thres) 
		{
        	sim_can_cut = true ;
			if (lb != phys_length_[i] || grid.dx0[i] != 0 || grid.dx1[i] != 0) 
			{
				num_subzoids = 3 ;
			}
			space_cut_dims.push_back(i) ;
		}

		key <<= (num_bits_dim - 4) ;
		key |= lb ;
		key <<= 2 ;
		int dx0 = grid.dx0[i] ;
		if (dx0 > 0)
		{
			key |= 1 ;
		}
		else if (dx0 < 0)
		{
			key |= 2 ;
		}
		else
		{
			key |= 0 ; //this is a no-op
		}
		key <<= 2 ; //shift by 2 bits
		
		int dx1 = grid.dx1[i] ;
		if (dx1 > 0)
		{
			key |= 1 ;
		}
		else if (dx1 < 0)
		{
			key |= 2 ;
		}
		else
		{
			key |= 0 ; //this is a no-op
		}

		num_grid_points *= ((lb + tb) / 2) ;
		/*if (lb == 0 && lt == 1)
		{
			empty_zoid = true ;
		}*/
    }
#ifdef PERMUTE
	if (sim_can_cut)
	{
		int size = space_cut_dims.size() ;
		bool flag [size] ;
		int permuted_dims [size] ;
		for (int i = 0 ; i < size ; i++)
		{
			flag [i] = false ;
		}
		int end = size - 1 ;
		bool cut_contiguous_dim = false ;
		if (space_cut_dims [end] == 0)
		{
			flag [end] = true ; //mark index 'end' as taken
			permuted_dims [end] = 0 ;
			cut_contiguous_dim = true ;
			end = size - 2 ;
		}
		//permute the non-contiguous dimensions
		for (int i = end ; i >= 0 ; i--)
		{
			while (1)
			{
				int index = rand () % size ;
				if (flag [index] == false)
				{
					permuted_dims [i] = space_cut_dims [index] ;
					//mark index as taken
					flag [index] = true ;
					break ;
				}
			}
		}
		for (int i = size - 1 ; i >= 0 ; i--)
		{
			space_cut_dims [i] = permuted_dims [i] ;
		}
#ifndef NDEBUG
		for (int i = 0 ; i < size ; i++)
		{
			assert (flag [i]) ;
		}
		if (cut_contiguous_dim)
		{
			assert (space_cut_dims [size - 1] == 0) ; 
		}
#endif
	}
#endif
	num_grid_points *= lt ;
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
			lt, centroid, centroid, index, l_father_grid, (1 << N_RANK) - 1);
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
			lt, centroid, centroid, index, l_father_grid, (1 << N_RANK) - 1);
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
		//a zoid with the projection already exists. return
		//Add the time of the zoid.
		child_time += z.time ;
		//start measuring linkage time
		stopwatch_start(ptr) ;
		return ;
	}
	//if empty internal zoid, do not divide any further
	/*if (empty_zoid && ! call_boundary)
	{
		sim_can_cut = false ;
	}*/
	bool divide_and_conquer = false ;
	bool time_cut = false ;
	time_type time_cut_elapsed_time = 0, space_cut_elapsed_time = LONG_MAX ;
	time_type max_loop_time = 0 ;
#ifdef FIXED_TIME_CUT
	//cut in time only when space cut is not possible
	if (lt > l_dt_stop && ! sim_can_cut)
#else
	if (lt > l_dt_stop)  //time cut
#endif
	{
		divide_and_conquer = true ;
		time_cut = true ;
		m_zoids [index].set_capacity(max (2, num_subzoids)) ;
		m_zoids [index].resize_children(2) ;
	
		time_type ltime = 0, ctime = 0 ;
        // cut into time 
		stopwatch_start(ptr) ;
        int halflt = lt / 2;
        l_son_grid = l_father_grid;
        if (call_boundary) {
            symbolic_trap_space_time_cut_boundary(t0, t0+halflt, l_son_grid,
					index, 0, ltime, ctime, f, bf, best_time);
        } else {
            symbolic_trap_space_time_cut_interior(t0, t0+halflt, l_son_grid,
					index, 0, ltime, ctime, f, best_time);
        }
		for (int i = 0; i < N_RANK; ++i) {
			l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
			l_son_grid.dx0[i] = l_father_grid.dx0[i];
			l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
			l_son_grid.dx1[i] = l_father_grid.dx1[i];
		}
#ifdef STOP_TUNING_EARLY
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t) ;
		time_cut_elapsed_time += t ;

		if (t + ltime + ctime < best_time)
		{
			stopwatch_start(ptr) ;
#endif
			if (call_boundary) {
				symbolic_trap_space_time_cut_boundary(t0+halflt, t1, l_son_grid,
						index, 1, ltime, ctime, f, bf, best_time);
			} else {
				symbolic_trap_space_time_cut_interior(t0+halflt, t1, l_son_grid,
						index, 1, ltime, ctime, f, best_time);
			}
#ifdef STOP_TUNING_EARLY
		}
#endif
		//measure the remaining function call overhead.
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t) ;
		//ltime is division cost + partial function call overhead.
		time_cut_elapsed_time += t + ltime + ctime ;
		assert (ltime >= 0) ;
		assert (ctime >= 0) ;
		assert (t >= 0) ;
		assert (m_zoids [index].num_children == 2) ;
		best_time = min (best_time, time_cut_elapsed_time) ;
#ifndef NDEBUG
		m_zoids [index].ttime = time_cut_elapsed_time ;
#endif
		
#if defined (MEASURE_COLD_MISS) || defined (SUBSUMPTION3) || \
	defined (SUBSUMPTION_TIME)
		time_type sum_loop_time = 0 ;
		time_type cache_penalty_children = 0 ;
		for (int i = 0 ; i < 2 ; i++)
		{
			unsigned long child_index = m_zoids [index].children [i] ;
#ifdef SUBSUMPTION_TIME
			max_loop_time = max(m_zoids [child_index].loop_time, max_loop_time);
#endif
#ifdef SUBSUMPTION3
			sum_loop_time += m_zoids [child_index].loop_time ;
#endif
#ifdef MEASURE_COLD_MISS
			cache_penalty_children += m_zoids [child_index].cache_penalty_time ;
#endif
		}
		assert (sum_loop_time >= 0) ;
		assert (cache_penalty_children >= 0) ;
#ifdef SUBSUMPTION3
		max_loop_time = max (max_loop_time, sum_loop_time) ;
#endif
#ifdef MEASURE_COLD_MISS
		m_zoids [index].cache_penalty_time = max(cache_penalty_children, 
										m_zoids [index].cache_penalty_time) ;
#endif
#endif
    }
	zoid_type bak ;
	if (sim_can_cut)
	{
		divide_and_conquer = true ;
		if (time_cut)
		{
			//back up the time cut children data
			bak = m_zoids [index] ;
			assert (bak.num_children == 2) ;
		}
		m_zoids [index].set_capacity(num_subzoids) ;
	}
	
	zoid_type bak2 ;
	int best_case = -1 ; 
    for (int j = 0 ; j < space_cut_dims.size() ; j++) {
		int i = space_cut_dims [j] ;
		time_type ltime = 0, ctime = 0  ;
		time_type elapsed_time = 0 ;
		if (grid.x1[i] - grid.x0[i] == phys_length_[i] && 
			l_father_grid.dx0[i] == 0 && l_father_grid.dx1[i] == 0)
		{
			m_zoids [index].resize_children(2) ;
		}
		else
		{
			m_zoids [index].resize_children(3) ;
		}
	
		//measure the divide time
		stopwatch_start(ptr) ;
        unsigned long lb, tb;
        int thres ;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
        thres = slope_[i] * lt ;
        bool cut_lb = (lb < tb);
		if (cut_lb) {
			assert(lb != phys_length_[i] || l_father_grid.dx0[i] != 0 || l_father_grid.dx1[i] != 0);
			const int mid = lb/2;
			l_son_grid = l_father_grid;
			const int l_start = l_father_grid.x0[i];
			const int l_end = l_father_grid.x1[i];

			//process the middle gray zoid
			SPACE_CUT_TUNE_BOUNDARY(0, l_start + mid - thres, l_start + mid + thres, slope_[i], -slope_[i]) ;
			SPACE_CUT_TUNE_BOUNDARY(1, l_start, l_start + mid - thres, l_father_grid.dx0[i], slope_[i]) ;
			SPACE_CUT_TUNE_BOUNDARY(2, l_start + mid + thres, l_end, -slope_[i], l_father_grid.dx1[i]) ;
		} // end if (cut_lb) 
		else { // cut_tb 
			if (lb == phys_length_[i] && l_father_grid.dx0[i] == 0 && l_father_grid.dx1[i] == 0) {
			   /* initial cut on the dimension */
				const int mid = tb/2;
				l_son_grid = l_father_grid;
				const int l_start = l_father_grid.x0[i];
				const int l_end = l_father_grid.x1[i];

				SPACE_CUT_TUNE_BOUNDARY(0, l_start, l_end, slope_[i], -slope_[i]);
				SPACE_CUT_TUNE_BOUNDARY(1, l_end, l_end, -slope_[i], slope_[i]);
			} else { /* NOT the initial cut! */
				const int mid = tb/2;
				l_son_grid = l_father_grid;
				const int l_start = l_father_grid.x0[i];
				const int l_end = l_father_grid.x1[i];
				const int ul_start = l_father_grid.x0[i] + l_father_grid.dx0[i] * lt;
				SPACE_CUT_TUNE_BOUNDARY(0, l_start, ul_start + mid, l_father_grid.dx0[i], -slope_[i]) ;
				SPACE_CUT_TUNE_BOUNDARY(1, ul_start + mid, l_end, slope_[i], l_father_grid.dx1[i]) ;
				SPACE_CUT_TUNE_BOUNDARY(2, ul_start + mid, ul_start + mid, -slope_[i], slope_[i]) ;
			}
		} // end if (cut_tb) 
	
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t) ;
		assert (ltime >= 0) ;
		assert (t >= 0) ;
		assert (ctime >= 0) ;
		//ltime is division cost + partial function call overhead.
		elapsed_time += t + ltime + ctime ;
		if (elapsed_time < space_cut_elapsed_time)
		{
			space_cut_elapsed_time = elapsed_time ;
			best_case = i ;
			//back up the zoid with its children.
			bak2 = m_zoids [index] ;
			best_time = min(best_time, elapsed_time) ;
		}
		
		assert (m_zoids [index].num_children <= num_subzoids) ;
#if defined (MEASURE_COLD_MISS) || defined (SUBSUMPTION3) || \
	defined (SUBSUMPTION_TIME)
		time_type sum_loop_time = 0 ;
		time_type cache_penalty_children = 0 ;
		for (int k = 0 ; k < m_zoids [index].num_children ; k++)
		{
			unsigned long child_index = m_zoids [index].children [k] ;
#ifdef SUBSUMPTION_TIME
			max_loop_time = max(m_zoids [child_index].loop_time, max_loop_time);
#endif
#ifdef SUBSUMPTION3
			sum_loop_time += m_zoids [child_index].loop_time ;
#endif
#ifdef MEASURE_COLD_MISS
			cache_penalty_children += m_zoids [child_index].cache_penalty_time ;
#endif
		}
		assert (sum_loop_time >= 0) ;
		assert (cache_penalty_children >= 0) ;
#ifdef SUBSUMPTION3
		max_loop_time = max (max_loop_time, sum_loop_time) ;
#endif
#ifdef MEASURE_COLD_MISS
		m_zoids [index].cache_penalty_time = max(cache_penalty_children, 
									m_zoids [index].cache_penalty_time) ;
#endif
#endif

#ifdef FIXED_SPACE_CUT
		//cut only in one dimension.
		break ;
#endif
    }
	if (sim_can_cut)
	{
		assert (space_cut_elapsed_time >= 0) ;
		assert (best_case >= 0 && best_case < N_RANK) ;
		//set the decision with the best case found
		decision = (decision_type) 1 << (best_case + 1) ;
		//restore the back up.
		m_zoids [index] = bak2 ;
		assert (m_zoids [index].num_children <= num_subzoids) ;
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
	divide_and_conquer_time += bdry_time ;

	bool force_divide = false ;
	bool child_divides = false ;
#ifdef SUBSUMPTION_SPACE
	for (int i = 0 ; i < m_zoids [index].num_children ; i++)
	{
		unsigned long child_index = m_zoids [index].children [i] ;
		decision_type d = m_zoids [child_index].decision ;
		if (d & m_space_cut_mask || d & 1)
		{
			child_divides = true ;
			//test if grand child divided
			if ((int) m_zoids [child_index].num_level_divide >= 1) 
			{
				force_divide = true ;
				m_zoids [index].num_level_divide = 
							m_zoids [child_index].num_level_divide ;
			}
		}
	}
#endif

#if defined (SUBSUMPTION_TIME) || defined (SUBSUMPTION3)
	if (divide_and_conquer && divide_and_conquer_time < max_loop_time)
	{
		force_divide = true ;
	}
#endif
	if (num_grid_points < 1000 || //empty_zoid || 
		divide_and_conquer_time < ptr->measure_time)
	{
		force_divide = false ;
	}

	// base case
	//The subsumption in time works as follows.
	//suppose loop_time(z) >= loop_time(z'), z' \in tree(z)
	//if divide_and_conquer_time(z) < max_{z' \in tree(z)} loop_time(z')
	//	then avoid computing loop_time(z).
	//else compute loop_time(z).
	//Other subsumption heuristics work similarly to avoid computing 
	//loop_time(z).
	time_type loop_time = LONG_MAX ;
	if (force_divide || lt > (m_initial_height + 1) / 2 )
	{
		//do not compute loop_time.
		m_zoids [index].decision = decision ;
		m_zoids [index].decision |= call_boundary <<
					(zoid_type::NUM_BITS_DECISION - 1) ;
		m_zoids [index].time = divide_and_conquer_time ;
	}
	else if (! divide_and_conquer)
	{
		//determine the looping time on the zoid
		time_type t1_, t2_ ;
		stopwatch_start(ptr) ;
		if (call_boundary)
		{
			base_case_kernel_boundary(t0, t1, l_father_grid, bf) ;
		} 
		else 
		{ 
			f(t0, t1, l_father_grid) ;
		}
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t1_) ;
		
		stopwatch_start(ptr) ;
		if (call_boundary)
		{
			base_case_kernel_boundary(t0, t1, l_father_grid, bf) ;
		} 
		else 
		{ 
			f(t0, t1, l_father_grid) ;
		}
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t2_) ;
		loop_time = min (t1_, t2_) + bdry_time ;
#ifdef MEASURE_COLD_MISS
		m_zoids [index].cache_penalty_time = max(t1_ - t2_, (time_type) 0) ;
#endif
	
		assert (loop_time >= 0) ;
#ifndef NDEBUG
		m_zoids [index].ltime = loop_time ;
#endif
		//set a flag to indicate that we looped on z.
		m_zoids [index].decision = (decision_type) 1 << 
					  (zoid_type::NUM_BITS_DECISION - 2) ;
		m_zoids [index].decision |= call_boundary << 
					  (zoid_type::NUM_BITS_DECISION - 1) ;
		m_zoids [index].time = loop_time ;
	}
	else 
	{
		//determine the looping time on the zoid
		/*
		stopwatch_start(ptr) ;
		if (call_boundary)
		{
			base_case_kernel_boundary(t0, t1, l_father_grid, bf);
		} 
		else 
		{
			f(t0, t1, l_father_grid);
		}
		stopwatch_stop(ptr) ;
		stopwatch_get_elapsed_time(ptr, t) ;
		loop_time = t + bdry_time ;
		assert (loop_time >= 0) ;
		*/
		LOOP_BOUNDARY ;
		loop_time += bdry_time ;
#ifndef NDEBUG
		m_zoids [index].ltime = loop_time ;
#endif
		//max_loop_time of z is its loop time.
		max_loop_time = loop_time ;
		//compare divide and conquer time with loop_time 
		if(divide_and_conquer_time < zoid_type::FUZZ * loop_time)
		{
			//choose divide and conquer
			m_zoids [index].decision = decision ;
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
		m_zoids [index].decision |= (decision_type) call_boundary << 
						  (zoid_type::NUM_BITS_DECISION - 1) ;
	}
#if defined (SUBSUMPTION_TIME) || defined (SUBSUMPTION3)
	m_zoids [index].max_loop_time = max_loop_time ;
#endif

	child_time += m_zoids [index].time ;
	//start measuring linkage time
	stopwatch_start(ptr) ;
}

template <int N_RANK> template <typename F>
inline void auto_tune<N_RANK>::
trap_space_time_cut_interior(int t0,
			int t1, grid_info<N_RANK> const & grid, 
			unsigned long zoid_index, F const & f)
{
    const int lt = t1 - t0;
    grid_info<N_RANK> l_son_grid;
	simple_zoid_type * projection_zoid = &(m_simple_zoids [zoid_index]) ;
	assert (projection_zoid) ;
#ifndef NDEBUG
    for (int i = N_RANK-1; i >= 0; --i) {
		grid_info <N_RANK> & grid2 = projection_zoid->info ;

		int x0 = grid.x0 [i] ;
		int x1 = grid.x1 [i] ;

		int x0_ = grid2.x0 [i] ;
		int x1_ = grid2.x1 [i] ;

		int x2 = grid.x0[i] + grid.dx0[i] * lt ;
		int x3 = grid.x1[i] + grid.dx1[i] * lt ;

		int x2_ = grid2.x0[i] + grid2.dx0[i] * lt ;
		int x3_ = grid2.x1[i] + grid2.dx1[i] * lt ;

		if (x1 - x0 != x1_ - x0_ || x3 - x2 != x3_ - x2_ ||
			grid.dx0 [i] != grid2.dx0 [i] || grid.dx1[i] != grid2.dx1[i])
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
		decision_type dim_to_cut = projection_zoid->decision & m_space_cut_mask;
		assert (__builtin_popcount((unsigned int) dim_to_cut) == 1) ;
		int i = __builtin_ctz((int) dim_to_cut) - 1 ;
		//cout << i << " " ; 
#ifndef NDEBUG
		decision_type d = (decision_type) 1 << (i + 1) ;
		assert (d & projection_zoid->decision) ;
#endif
		assert (i >= 0 && i < N_RANK) ;
        int lb, thres, tb;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
        bool cut_lb = (lb < tb);
        thres = (slope_[i] * lt);
	
		if (cut_lb) {
			// if cutting lb, there's no initial cut! 
			assert(lb != phys_length_[i] || grid.dx0[i] != 0 || grid.dx1[i] != 0);
			const int mid = lb/2;
			l_son_grid = grid;
			const int l_start = grid.x0[i];
			const int l_end = grid.x1[i];

			//process the middle gray zoid
			l_son_grid.x0[i] = l_start + mid - thres;
			l_son_grid.dx0[i] = slope_[i];
			l_son_grid.x1[i] = l_start + mid + thres;
			l_son_grid.dx1[i] = -slope_[i];
			unsigned long index = projection_zoid->children [0] ;
        	trap_space_time_cut_interior(t0, t1, l_son_grid, index, f);

			l_son_grid.x0[i] = l_start;
			l_son_grid.dx0[i] = grid.dx0[i];
			l_son_grid.x1[i] = l_start + mid - thres;
			l_son_grid.dx1[i] = slope_[i];
			index = projection_zoid->children [1] ;
        	trap_space_time_cut_interior(t0, t1, l_son_grid, index, f);

			l_son_grid.x0[i] = l_start + mid + thres;
			l_son_grid.dx0[i] = -slope_[i];
			l_son_grid.x1[i] = l_end;
			l_son_grid.dx1[i] = grid.dx1[i];
			index = projection_zoid->children [2] ;
        	trap_space_time_cut_interior(t0, t1, l_son_grid, index, f);
		} // end if (cut_lb) 
		else { // cut_tb 
			const int mid = tb/2;
			l_son_grid = grid;

			const int l_start = (grid.x0[i]);
			const int l_end = (grid.x1[i]);
			const int ul_start = (grid.x0[i] + grid.dx0[i] * lt);
			l_son_grid.x0[i] = l_start;
			l_son_grid.dx0[i] = grid.dx0[i];
			l_son_grid.x1[i] = ul_start + mid;
			l_son_grid.dx1[i] = -slope_[i];
			unsigned long index = projection_zoid->children [0] ;
        	trap_space_time_cut_interior(t0, t1, l_son_grid, index, f);

			l_son_grid.x0[i] = ul_start + mid;
			l_son_grid.dx0[i] = slope_[i];
			l_son_grid.x1[i] = l_end;
			l_son_grid.dx1[i] = grid.dx1[i];
			index = projection_zoid->children [1] ;
        	trap_space_time_cut_interior(t0, t1, l_son_grid, index, f);

			l_son_grid.x0[i] = ul_start + mid;
			l_son_grid.dx0[i] = -slope_[i];
			l_son_grid.x1[i] = ul_start + mid;
			l_son_grid.dx1[i] = slope_[i];
			index = projection_zoid->children [2] ;
        	trap_space_time_cut_interior(t0, t1, l_son_grid, index, f);
		} // end if (cut_tb) 
    }
	else if (projection_zoid->decision & 1)
	{
		assert (projection_zoid->num_children == 2) ;
        /* cut into time */
		assert (projection_zoid->children [0]) ;
		assert (projection_zoid->children [1]) ;
        int halflt = lt / 2;
        l_son_grid = grid;
		unsigned long index = projection_zoid->children [0] ;
        trap_space_time_cut_interior(t0, t0+halflt, l_son_grid, index, f);

        for (int i = 0; i < N_RANK; ++i) {
            l_son_grid.x0[i] = grid.x0[i] + grid.dx0[i] * halflt;
            l_son_grid.dx0[i] = grid.dx0[i];
            l_son_grid.x1[i] = grid.x1[i] + grid.dx1[i] * halflt;
            l_son_grid.dx1[i] = grid.dx1[i];
        }
		index = projection_zoid->children [1] ;
        trap_space_time_cut_interior(t0+halflt, t1, l_son_grid, index, f);
    }
	else
	{
#ifdef WRITE_DAG
		bool empty_zoid = false ;
		for (int i = 0 ; lt == 1 && i < N_RANK ; i++) 
		{
			unsigned long lb = (grid.x1[i] - grid.x0[i]);
			if (lb == 0)
			{
				empty_zoid = true ;
			}
		}
		for (int i = 0 ; !empty_zoid && i < N_RANK ; i++) 
		{
			unsigned long lb, tb;
			lb = (grid.x1[i] - grid.x0[i]);
			tb = (grid.x1[i]+ grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
			file_interior [i] << lt << " , " << max(lb, tb)  << endl ;
		}
#endif
#ifdef TIME_INVARIANCE_INTERIOR
		assert (projection_zoid->decision == 
				1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ;
#else
		assert (projection_zoid->decision == 
				3 << (zoid<N_RANK>::NUM_BITS_DECISION - 2) ||
				projection_zoid->decision ==
				1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ;
#endif
		//loop
#ifdef MEASURE_STATISTICS
		stopwatch * ptr = &m_stopwatch ;
		stopwatch_start(ptr) ;
#endif
		f(t0, t1, grid);
#ifdef MEASURE_STATISTICS
		stopwatch_stop(ptr) ;
		time_type t = 0 ;
		stopwatch_get_elapsed_time(ptr, t) ;
		assert (t >= 0) ;
		double time = stopwatch_time_to_double(t) ;
		zoid_statistics <N_RANK> & z = m_statistics [zoid_index] ;
		z.total += time ;
		z.count++ ;
		z.variance += (time * time) ;
		z.min = min(t, z.min) ;
		z.max = max(t, z.max) ;
		z.height = lt ;
#ifndef NDEBUG
		grid_info <N_RANK> & grid1 = z.info ;
		grid_info <N_RANK> & grid2 = projection_zoid->info ;
		
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
			assert (grid1.x0 [i] == grid2.x0 [i]) ;
			assert (grid1.x1 [i] == grid2.x1 [i]) ;
			assert (grid1.dx0 [i] == grid2.dx0 [i]) ;
			assert (grid1.dx1 [i] == grid2.dx1 [i]) ;
		}
#endif
#endif
	}
}

#define SPACE_CUT_BOUNDARY2(k) \
if (call_boundary) \
{ \
	trap_space_time_cut_boundary(t0, t1, l_son_grid, projection_zoid->children [k], f, bf); \
} \
else \
{ \
	trap_space_time_cut_interior(t0, t1, l_son_grid, projection_zoid->children [k], f); \
}

template <int N_RANK> template <typename F, typename BF>
inline void auto_tune<N_RANK>::trap_space_time_cut_boundary(int t0, int t1,	
				grid_info<N_RANK> const & grid, unsigned long zoid_index, 
				F const & f, BF const & bf)
{
    const int lt = t1 - t0;
    grid_info<N_RANK> l_father_grid = grid, l_son_grid;
	simple_zoid_type * projection_zoid = &(m_simple_zoids [zoid_index]) ;
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
			x2 = pmod(grid.x0[i] + grid.dx0[i] * lt, phys_length_ [i]) ;
			x3 = pmod(grid.x1[i] + grid.dx1[i] * lt, phys_length_ [i]) ;

			x0_ = pmod(grid2.x0 [i], phys_length_ [i]) ;
			x1_ = pmod(grid2.x1 [i], phys_length_ [i]) ;
			x2_ = pmod(grid2.x0[i] + grid2.dx0[i] * lt, phys_length_ [i]) ;
			x3_ = pmod(grid2.x1[i] + grid2.dx1[i] * lt, phys_length_ [i]) ;
			if (x0 != x0_ || x1 != x1_ || x2 != x2_ || x3 != x3_ ||
				grid.dx0 [i] != grid2.dx0 [i] || grid.dx1[i] != grid2.dx1[i])
			{
				error = true ;
			}
		}
		else
		{
			x0 = grid.x0 [i] ;
			x1 = grid.x1 [i] ;
			x2 = grid.x0[i] + grid.dx0[i] * lt ;
			x3 = grid.x1[i] + grid.dx1[i] * lt ;

			x0_ = grid2.x0 [i] ;
			x1_ = grid2.x1 [i] ;
			x2_ = grid2.x0[i] + grid2.dx0[i] * lt ;
			x3_ = grid2.x1[i] + grid2.dx1[i] * lt ;
			if (x1 - x0 != x1_ - x0_ || x3 - x2 != x3_ - x2_ ||
				grid.dx0 [i] != grid2.dx0 [i] || grid.dx1[i] != grid2.dx1[i])
			{
				error = true ;
			}
		}
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
		decision_type dim_to_cut = projection_zoid->decision & m_space_cut_mask;
		assert (__builtin_popcount((unsigned int) dim_to_cut) == 1) ;
		int i = __builtin_ctz((int) dim_to_cut) - 1 ;
#ifndef NDEBUG
		decision_type d = (decision_type) 1 << (i + 1) ;
		assert (d & projection_zoid->decision) ;
#endif
		assert (i >= 0 && i < N_RANK) ;
        int lb, thres, tb;
        lb = (grid.x1[i] - grid.x0[i]);
        tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
        thres = (slope_[i] * lt);
        bool cut_lb = (lb < tb);

		// can_cut 
		if (cut_lb) {
			// if cutting lb, there's no initial cut! 
			assert(lb != phys_length_[i] || l_father_grid.dx0[i] != 0 || l_father_grid.dx1[i] != 0);
			const int mid = lb/2;
			l_son_grid = l_father_grid;
			const int l_start = l_father_grid.x0[i];
			const int l_end = l_father_grid.x1[i];

			//process the middle gray zoid
			l_son_grid.x0[i] = l_start + mid - thres;
			l_son_grid.dx0[i] = slope_[i];
			l_son_grid.x1[i] = l_start + mid + thres;
			l_son_grid.dx1[i] = -slope_[i];
			SPACE_CUT_BOUNDARY2(0) ;

			l_son_grid.x0[i] = l_start;
			l_son_grid.dx0[i] = l_father_grid.dx0[i];
			l_son_grid.x1[i] = l_start + mid - thres;
			l_son_grid.dx1[i] = slope_[i];
			SPACE_CUT_BOUNDARY2(1) ;

			l_son_grid.x0[i] = l_start + mid + thres;
			l_son_grid.dx0[i] = -slope_[i];
			l_son_grid.x1[i] = l_end;
			l_son_grid.dx1[i] = l_father_grid.dx1[i];
			SPACE_CUT_BOUNDARY2(2) ;
		} // end if (cut_lb) 
		else { // cut_tb 
			if (lb == phys_length_[i] && l_father_grid.dx0[i] == 0 && l_father_grid.dx1[i] == 0) {
			   /* initial cut on the dimension */
				const int mid = tb/2;
				l_son_grid = l_father_grid;
				const int l_start = l_father_grid.x0[i];
				const int l_end = l_father_grid.x1[i];
				
				l_son_grid.x0[i] = l_start ;
				l_son_grid.dx0[i] = slope_[i];
				l_son_grid.x1[i] = l_end ;
				l_son_grid.dx1[i] = -slope_[i];
				SPACE_CUT_BOUNDARY2(0) ;

				l_son_grid.x0[i] = l_end ;
				l_son_grid.dx0[i] = -slope_[i];
				l_son_grid.x1[i] = l_end ;
				l_son_grid.dx1[i] = slope_[i];
				SPACE_CUT_BOUNDARY2(1) ;
			} else { /* NOT the initial cut! */
				const int mid = tb/2;
				l_son_grid = l_father_grid;
				const int l_start = l_father_grid.x0[i];
				const int l_end = l_father_grid.x1[i];
				const int ul_start = l_father_grid.x0[i] + l_father_grid.dx0[i] * lt;
				l_son_grid.x0[i] = l_start;
				l_son_grid.dx0[i] = l_father_grid.dx0[i];
				l_son_grid.x1[i] = ul_start + mid;
				l_son_grid.dx1[i] = -slope_[i];
				SPACE_CUT_BOUNDARY2(0) ;

				l_son_grid.x0[i] = ul_start + mid;
				l_son_grid.dx0[i] = slope_[i];
				l_son_grid.x1[i] = l_end;
				l_son_grid.dx1[i] = l_father_grid.dx1[i];
				SPACE_CUT_BOUNDARY2(1) ;

				l_son_grid.x0[i] = ul_start + mid;
				l_son_grid.dx0[i] = -slope_[i];
				l_son_grid.x1[i] = ul_start + mid;
				l_son_grid.dx1[i] = slope_[i];
				SPACE_CUT_BOUNDARY2(2) ;
			}
		} // end if (cut_tb) 
	}
	else if (projection_zoid->decision & 1)
	{
		assert (projection_zoid->num_children == 2) ;
		// cut into time 
		assert (projection_zoid->children [0]) ;
		assert (projection_zoid->children [1]) ;
		
		int halflt = lt / 2;
		l_son_grid = l_father_grid;
		unsigned long index = projection_zoid->children [0] ;
		if (call_boundary) {
			trap_space_time_cut_boundary(t0, t0+halflt, l_son_grid, index,f,bf);
		} else {
			trap_space_time_cut_interior(t0, t0+halflt, l_son_grid, index, f);
		}

		for (int i = 0; i < N_RANK; ++i) {
			l_son_grid.x0[i] = l_father_grid.x0[i] + l_father_grid.dx0[i] * halflt;
			l_son_grid.dx0[i] = l_father_grid.dx0[i];
			l_son_grid.x1[i] = l_father_grid.x1[i] + l_father_grid.dx1[i] * halflt;
			l_son_grid.dx1[i] = l_father_grid.dx1[i];
		}
		index = projection_zoid->children [1] ;
		if (call_boundary) {
			trap_space_time_cut_boundary(t0+halflt, t1, l_son_grid, index,f,bf);
		} else {
			trap_space_time_cut_interior(t0+halflt, t1, l_son_grid, index, f);
		}
	}
	else
	{
#ifdef WRITE_DAG
		bool empty_zoid = false ;
		for (int i = 0 ; lt == 1 && i < N_RANK ; i++) 
		{
			unsigned long lb = (grid.x1[i] - grid.x0[i]);
			if (lb == 0)
			{
				empty_zoid = true ;
			}
		}
		for (int i = 0 ; !empty_zoid && i < N_RANK ; i++) 
		{
			unsigned long lb, tb;
			lb = (grid.x1[i] - grid.x0[i]);
			tb = (grid.x1[i] + grid.dx1[i] * lt - grid.x0[i] - grid.dx0[i] * lt);
			if (call_boundary)
			{
				file_boundary [i] << lt << " , " << max(lb, tb)  << endl ;
			}
			else
			{
				file_interior [i] << lt << " , " << max(lb, tb)  << endl ;
			}
		}
#endif
		//loop
		assert (projection_zoid->decision == 
				3 << (zoid<N_RANK>::NUM_BITS_DECISION - 2) ||
				projection_zoid->decision ==
				1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ;
#ifdef MEASURE_STATISTICS
		stopwatch * ptr = &m_stopwatch ;
		stopwatch_start(ptr) ;
#endif
		if (call_boundary) {
#ifdef TIME_INVARIANCE_BOUNDARY
			assert (projection_zoid->decision == 
					3 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ;
#endif
            base_case_kernel_boundary(t0, t1, l_father_grid, bf);
        } else { 
#ifdef TIME_INVARIANCE_INTERIOR
			assert (projection_zoid->decision == 
					1 << (zoid<N_RANK>::NUM_BITS_DECISION - 2)) ;
#endif
            f(t0, t1, l_father_grid);
        }
#ifdef MEASURE_STATISTICS
		stopwatch_stop(ptr) ;
		time_type t = 0 ;
		stopwatch_get_elapsed_time(ptr, t) ;
		assert (t >= 0) ;
		double time = stopwatch_time_to_double(t) ;
		zoid_statistics <N_RANK> & z = m_statistics [zoid_index] ;
		z.total += time ;
		z.count++ ;
		z.variance += (time * time) ;
		z.min = min(t, z.min) ;
		z.max = max(t, z.max) ;
		z.height = lt ;
		z.boundary = (call_boundary ? 1 : 0) ;
#ifndef NDEBUG
		grid_info <N_RANK> & grid1 = z.info ;
		grid_info <N_RANK> & grid2 = projection_zoid->info ;
		
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
			assert (grid1.x0 [i] == grid2.x0 [i]) ;
			assert (grid1.x1 [i] == grid2.x1 [i]) ;
			assert (grid1.dx0 [i] == grid2.dx0 [i]) ;
			assert (grid1.dx1 [i] == grid2.dx1 [i]) ;
		}
#endif
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
#undef base_case_kernel_boundary_tune

#endif
