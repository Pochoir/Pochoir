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
 *        Created:  03/25/2013
 *         Author:  Eka Palamadai, epn@mit.edu
 * ============================================================================
 */
#ifndef KERNEL_SELECTION_SAWZOID_HEADER_HPP
#define KERNEL_SELECTION_SAWZOID_HEADER_HPP

#include "rbq_common.h"
#include <deque>
#include <unordered_map>
#include <climits> 
using namespace std ;

typedef unsigned int word_type ;


template <int N_RANK>
class kernel_selection_sawzoid
{
private:
#ifndef NDEBUG
template <int N_RANK>
class zoid
{
	grid_info <N_RANK> info ;
	int height ;
	char kernel ;
	int t0, t1 ;
};
	typedef zoid <N_RANK> zoid_type ;
	//hash table of (key, zoid) pair
	typedef unordered_map<unsigned long, zoid_type> leaf_kernel_table ;
    typedef typename unordered_map<unsigned long, zoid_type>::iterator
                  lk_table_iterator ;
	//vector of zoids 
	vector <zoid_type> m_kernel_map ;
#else
	//hash table of (key, kernel) pair
	typedef unordered_map<unsigned long, char> leaf_kernel_table ;
    typedef typename unordered_map<unsigned long, char>::iterator
                    lk_table_iterator ;
	//vector of kernels
	vector <char> m_kernel_map ;
#endif
	//typedef unordered_map<int, leaf_kernel_table> height_leaf_kernel_table;
    //typedef typename unordered_map<int, leaf_kernel_table>::iterator
    //                hlk_table_iterator ;
	//vector <height_leaf_kernel_table> m_arr_height_leaf_kernel_table ;
	vector <leaf_kernel_table> m_arr_leaf_kernel_table ;
	int num_bits_dim ; //# of bits for bottom and top widths in a dimension
	int num_bits_width ; //# of bits to store width

	Algorithm <N_RANK> & m_algo ; // a reference to Algorithm
	pochoir_clone_array <N_RANK> * m_clone_array ; 
	typedef typename Algorithm<N_RANK>::queue_info queue_info ;
	//int m_extended_length [N_RANK] ;
	//int m_num_leaves [N_RANK] ;

	void initialize(grid_info<N_RANK> const & grid, int t0, int t1)
	{
		cout << "t0 " << t0 << " t1 " << t1 << endl ;
		int T = t1 - t0 ;
		unsigned long volume = 1 ; // extended_volume = 1 ;
		int total_num_leaves = 1 ;
		/*for (int i = 0 ; i < N_RANK ; i++)
		{
			int W = m_algo.phys_length_ [i] ;
			int slope = m_algo.slope_ [i] ;
		
			int Wn = W / (slope << 1) ;
			//find index of most significant bit that is set
			int index_msb = (sizeof(int) << 3) - __builtin_clz(Wn) - 1 ;
			//h = 2^floor(lg(Wn)). The zoid with height h undergoes a space cut.
			int h = 1 << index_msb ;
			if (T < h)
			{
				index_msb = (sizeof(int) << 3) - __builtin_clz(T) - 1 ;
				h = 1 << index_msb ;
			}
			//cout << " length [ " << i << " ] " << W + h << endl ; 
			//m_extended_length [i] = W + h ;
			//extended_volume *= (W + h) ; 
			volume *= W ;
			//m_num_leaves [N_RANK] = W / (2 * slope * dt_recursive_) * 2 ;
			//total_num_leaves *= m_num_leaves [N_RANK] ;
		}*/
		for (int i = 0 ; i < N_RANK ; i++)
		{
			volume *= (grid.x1[i] - grid.x0[i]) ;		
		}

		//m_kernel_map.reserve(2 * N_RANK * extended_volume) ; 
		//m_kernel_map.resize(2 * N_RANK * extended_volume) ; 
		//m_kernel_map.reserve(2 * extended_volume) ; 
		//m_kernel_map.resize(2 * extended_volume) ; 
		m_kernel_map.reserve(2 * volume) ; 
		m_kernel_map.resize(2 * volume) ; 
		//m_kernel_map.reserve(2 * total_num_leaves) ; 
		//m_kernel_map.resize(2 * total_num_leaves) ; 

		//m_arr_leaf_kernel_table.reserve(extended_volume) ;	
		//m_arr_leaf_kernel_table.resize(extended_volume) ;
		m_arr_leaf_kernel_table.reserve(volume) ;	
		m_arr_leaf_kernel_table.resize(volume) ;

		//for (int i = 0 ; i < 2 * N_RANK * extended_volume ; i++)
		//for (int i = 0 ; i < 2 * extended_volume ; i++)
		for (int i = 0 ; i < 2 * volume ; i++)
		{
#ifndef NDEBUG
			m_kernel_map [i].kernel = (char) 0 ;
#else
			m_kernel_map [i] = (char) 0 ;
#endif
		}

		num_bits_dim = sizeof (unsigned long) * 8 / N_RANK ;
		num_bits_width = sizeof (unsigned long) * 8 / (2 * N_RANK) ;
	}

	template <typename F>
	void compute_geneity(int h, grid_info<N_RANK> const & grid, 
						word_type & geneity, F const & f) ;


	void set_clone_array(pochoir_clone_array <N_RANK> * clone_array)
	{
		m_clone_array = clone_array ;
	}


	template <typename F>
	inline void symbolic_modified_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F>
	inline void symbolic_modified_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F>
	inline void symbolic_modified_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F>
	inline void symbolic_modified_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F, typename BF>
	inline void heterogeneous_modified_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, 
		F const & f, BF const & bf) ;

	template <typename F>
	inline void heterogeneous_modified_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F, typename BF>
	inline void heterogeneous_modified_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, 
		F const & f, BF const & bf) ;

	template <typename F>
	inline void heterogeneous_modified_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, 
		F const & f) ;

	public :

	kernel_selection_sawzoid(Algorithm<N_RANK> & alg, 
				grid_info<N_RANK> const & grid, bool power_of_two):m_algo(alg)
	{
		//initialize(grid, power_of_two) ;
	}


	~kernel_selection_sawzoid()
	{
		for (int i = 0 ; i < m_arr_leaf_kernel_table.size() ; i++)
		{
			m_arr_leaf_kernel_table [i].clear() ;
		}
		m_arr_leaf_kernel_table.clear() ;
		m_kernel_map.clear() ;
	}

	template <typename F, typename BF, typename P>
    inline void do_power_of_two_time_cut(int t0, int t1,
        grid_info<N_RANK> const & grid, F const & f, BF const & bf, P const & p)
	{
		initialize(grid, t0, t1) ;
		int T = t1 - t0 ;
		int W = 0 ;  //max_width among all dimensions
		int slope ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			if (m_algo.phys_length_ [i] > W)
			{
				W = m_algo.phys_length_ [i] ;
				slope = m_algo.slope_ [i] ;
			}
		}
		//find index of most significant bit that is set
		int Wn = W / (slope << 1) ;
		int index_msb = (sizeof(int) << 3) - __builtin_clz(Wn) - 1 ;
		//h1 = 2^floor(lg(Wn)). The zoid with height h1 undergoes a space cut.
		int h1 = 1 << index_msb ;
		if (T < h1)
		{
			index_msb = (sizeof(int) << 3) - __builtin_clz(T) - 1 ;
			h1 = 1 << index_msb ;
		}
		cout << "h1 " << h1 << endl ;
		struct timespec start, end;
		clock_gettime(CLOCK_MONOTONIC, &start) ;
		symbolic_modified_space_time_cut_boundary(t0, t0 + h1, grid, p) ;
		int offset = t0 + T / h1 * h1 ;
		int h2 = t1 - offset ;
		while (h2 >= m_algo.dt_recursive_)
		{
			//find index of most significant bit that is set
			index_msb = (sizeof(int) << 3) - __builtin_clz(h2) - 1 ;
			//h = 2^floor(lg(h2)).
			int h = 1 << index_msb ;
			cout <<	" offset " << offset << " offset + h " <<
				offset + h << " h " << h << endl ;
			symbolic_modified_space_time_cut_boundary(offset, offset + h, grid,
														 p) ;
			offset += h ;
			h2 = t1 - offset ;
		}
		clock_gettime(CLOCK_MONOTONIC, &end) ;
		double ks_time = tdiff2(&end, &start) ;
		std::cout << "KS : consumed time :" << 1.0e3 * ks_time
				<< "ms" << std::endl;
		int m = T / h1 ;
		cout << " m " << m << endl ;
		clock_gettime(CLOCK_MONOTONIC, &start) ;
		for (int i = 0 ; i < m ; i++)
		{
			/*cout << "t0 " << t0 << " t1 " << t1 << 
				" h1 " << h1 << " t0 + h1 " <<
				t0 + h1 << endl ;*/
			heterogeneous_modified_space_time_cut_boundary(t0, t0 + h1, grid, 
				f, bf) ;
			t0 += h1 ;
		}

		h2 = t1 - t0 ;
		//time cuts happen only if height > dt_recursive_
		while (h2 >= m_algo.dt_recursive_)
		{
			//find index of most significant bit that is set
			index_msb = (sizeof(int) << 3) - __builtin_clz(h2) - 1 ;
			int h = 1 << index_msb ;
			cout << "t0 " << t0 << " t1 " << t1 << 
				" h " << h << " t0 + h " <<
				t0 + h << endl ;
			heterogeneous_modified_space_time_cut_boundary(t0, t0 + h, grid, 
				f, bf) ;
			t0 += h ;
			h2 = t1 - t0 ;
		}
		while (h2 > 1)
		{
			//find index of most significant bit that is set
			index_msb = (sizeof(int) << 3) - __builtin_clz(h2) - 1 ;
			int h = 1 << index_msb ;
			cout << "t0 " << t0 << " t1 " << t1 << 
				" h " << h << " t0 + h " <<
				t0 + h << endl ;
			/*for (int i = 0 ; i < N_RANK ; i++)
			{
				int num_triangles = m_algo.dx_recursive_ [i] / ((m_algo.slope_ [i] * h) << 1) ;
				m_algo.num_triangles [i] = max(1, num_triangles) ;
				cout << "num_triangles [ " << i << " ] " << m_algo.num_triangles [i] << endl ;
			}*/
			//if (h > 2)
			//{
				m_algo.space_time_cut_boundary(t0, t0 + h, grid, f, bf) ;
			//}
			/*else
			{
				//seems to be a race bug for h <= 2
				m_algo.abnormal_region_space_time_cut_boundary(t0, t0 + h, grid, f, bf) ;
			}*/
			t0 += h ;
			h2 = t1 - t0 ;
		}
		if (h2 == 1)
		{
			cout << "h = 1 t0 " << t0 << " t1 " << t1 << 
				 " t0 + h " << t0 + h2 << endl ;
			//base_case_kernel_boundary(t0, t0 + h2, grid, bf);
			m_algo.shorter_duo_sim_obase_bicut_p(t0, t0 + h2, grid, f, bf) ;
		}
		clock_gettime(CLOCK_MONOTONIC, &end) ;
		std::cout << "Compute time :" << 1.0e3 * tdiff2(&end, &start)
				<< "ms" << std::endl;
	}
} ;

template<>
template <typename F>
void kernel_selection_sawzoid<1>::compute_geneity(int h, 
					grid_info<1> const & grid, 
					word_type & geneity, F const & f)
{
	int lb = grid.x1[0] - grid.x0[0] ;
	int x0_top = grid.x0[0] + grid.dx0[0] * (h - 1) ;
	int x1_top = grid.x1[0] + grid.dx1[0] * (h - 1) ;
	int tb = x1_top - x0_top ;
	//cout << "comp geneity" << endl ;
	if (lb > tb)
	{
		for (int x0 = grid.x0 [0] ; x0 < grid.x1 [0] ; x0++)
		{
			int index = f(pmod(x0, m_algo.phys_length_ [0])) ;
			set_bit(geneity, index, 1) ;
		}
	}
	else
	{
		for (int x0 = x0_top ; x0 < x1_top ; x0++)
		{
			int index = f(pmod(x0, m_algo.phys_length_ [0])) ;
			set_bit(geneity, index, 1) ;
		}
		
	}
	//m_heterogeneity.insert(geneity) ;
}

template<>
template <typename F>
void kernel_selection_sawzoid<2>::compute_geneity(int h, 
					grid_info<2> const & grid, 
					word_type & geneity, F const & f)
{
	//cout << "comp geneity " << endl ;
	int x0 = grid.x0 [1] ;
	int x1 = grid.x1 [1] ;
	int x2 = grid.x0 [1] + grid.dx0 [1] * (h - 1) ;
	int x3 = grid.x1 [1] + grid.dx1 [1] * (h - 1) ;

	int y0 = grid.x0 [0] ;
	int y1 = grid.x1 [0] ;
	int y2 = grid.x0 [0] + grid.dx0 [0] * (h - 1) ; 
	int y3 = grid.x1 [0] + grid.dx1 [0] * (h - 1) ; 
	//cout << "x0 " << x0 << " x1 " << x1 << " x2 " << x2 << " x3 " << x3 << endl ;
	//cout << "y0 " << y0 << " y1 " << y1 << " y2 " << y2 << " y3 " << y3 << endl ;

	if (x0 <= x2 && y0 <= y2)
	{
		//cout << "bottom rectangle " << endl ;
		//bottom rectangle is the projection
		for (int x = x0 ; x < x1 ; x++)
		{
			for (int y = y0 ; y < y1 ; y++)
			{
				int index = f(pmod(x, m_algo.phys_length_ [1]), 
								pmod(y, m_algo.phys_length_ [0])) ;
				set_bit(geneity, index, 1) ;
			}
		}
	}
	else if (x0 > x2 && y0 > y2)
	{
		//cout << "top rectangle " << endl ;
		//top rectangle is the projection
		for (int x = x2 ; x < x3 ; x++)
		{
			for (int y = y2 ; y < y3 ; y++)
			{
				int index = f(pmod(x, m_algo.phys_length_ [1]), 
								pmod(y, m_algo.phys_length_ [0])) ;
				set_bit(geneity, index, 1) ;
			}
		}
	}
	else
	{
		//cout << "octagon " << endl ;
		//projection is an octagon
		if (x0 <= x2)
		{
			//cout << "x dimension converges and y dimension diverges" << endl ;
			//x dimension converges and y dimension diverges
			for (int x = x0, start = y0, end = y1 ; x < x2 ; x++)
			{
				for (int y = start ; y < end ; y++)
				{
					int index = f(pmod(x, m_algo.phys_length_ [1]),
                                pmod(y, m_algo.phys_length_ [0])) ;
					set_bit(geneity, index, 1) ;
				}
				start += grid.dx0 [0];
				end += grid.dx1 [0] ;
			}
			for (int x = x2 ; x < x3 ; x++)
			{
				for (int y = y2 ; y < y3 ; y++)
				{
					int index = f(pmod(x, m_algo.phys_length_ [1]),
                                pmod(y, m_algo.phys_length_ [0])) ;
					set_bit(geneity, index, 1) ;
				}
			}
			for (int x = x3, start = y2, end = y3 ; x < x1 ; x++)
			{
				for (int y = start ; y < end ; y++)
				{
					int index = f(pmod(x, m_algo.phys_length_ [1]),
                                pmod(y, m_algo.phys_length_ [0])) ;
					set_bit(geneity, index, 1) ;
				}
				start -= grid.dx0 [0];
				end -= grid.dx1 [0] ;
			}
		}
		else
		{
			//cout << "x dimension diverges and y dimension converges " << endl ;
			//cout << "x2 " << x2 << " x0 " << x0 << endl ;
			//x dimension diverges and y dimension converges
			for (int x = x2, start = y2, end = y3 ; x < x0 ; x++)
			{
				//cout << "start " << start << " end " << end << endl ;
				for (int y = start ; y < end ; y++)
				{
					int index = f(pmod(x, m_algo.phys_length_ [1]),
                                pmod(y, m_algo.phys_length_ [0])) ;
					//cout << pmod(x, m_algo.phys_length_ [1]) << ", " <<
					//	pmod(y, m_algo.phys_length_ [0]) << endl ;
					set_bit(geneity, index, 1) ;
				}
				start -= grid.dx0 [0];
				end -= grid.dx1 [0] ;
			}
			for (int x = x0 ; x < x1 ; x++)
			{
				for (int y = y0 ; y < y1 ; y++)
				{
					int index = f(pmod(x, m_algo.phys_length_ [1]),
                                pmod(y, m_algo.phys_length_ [0])) ;
					set_bit(geneity, index, 1) ;
				}
			}
			for (int x = x1, start = y0, end = y1 ; x < x3 ; x++)
			{
				//cout << " loop 2 start " << start << " end " << end << endl ;
				for (int y = start ; y < end ; y++)
				{
					int index = f(pmod(x, m_algo.phys_length_ [1]),
                                pmod(y, m_algo.phys_length_ [0])) ;
					//cout << pmod(x, m_algo.phys_length_ [1]) << ", " <<
					//	pmod(y, m_algo.phys_length_ [0]) << endl ;
					set_bit(geneity, index, 1) ;
				}
				start += grid.dx0 [0];
				end += grid.dx1 [0] ;
			}
		}
	}
	//m_heterogeneity.insert(geneity) ;
}

#endif
