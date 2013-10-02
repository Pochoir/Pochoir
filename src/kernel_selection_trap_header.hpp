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
#ifndef KERNEL_SELECTION_TRAP_HEADER_HPP
#define KERNEL_SELECTION_TRAP_HEADER_HPP

#include "rbq_common.h"
#include <deque>
#include <unordered_map>
#include <climits> 
using namespace std ;

typedef unsigned int word_type ;


template <int N_RANK>
class kernel_selection_trap
{
private:
#ifndef NDEBUG
template <int N_RANK>
class zoid
{
	grid_info <N_RANK> info ;
	int height ;
	char kernel ;
};
	typedef zoid <N_RANK> zoid_type ;
	//hash table of (key, zoid) pair
	typedef unordered_map<unsigned long, zoid_type> leaf_kernel_table ;
    typedef typename unordered_map<unsigned long, zoid_type>::iterator
                    lk_table_iterator ;
#else
	//hash table of (key, kernel) pair
	typedef unordered_map<unsigned long, char> leaf_kernel_table ;
    typedef typename unordered_map<unsigned long, char>::iterator
                    lk_table_iterator ;
#endif
	//typedef unordered_map<int, leaf_kernel_table> height_leaf_kernel_table;
    //typedef typename unordered_map<int, leaf_kernel_table>::iterator
              //      hlk_table_iterator ;
	//vector of geneity of leaves.
	//vector <vector <int> > kernel_map ;
	//vector <height_leaf_kernel_table> m_arr_height_leaf_kernel_table ;
	vector <leaf_kernel_table> m_arr_leaf_kernel_table ;
	int num_bits_dim ; //# of bits for bottom and top widths in a dimension
	int num_bits_width ; //# of bits to store width

	//int min_leaf_height ;
	Algorithm <N_RANK> & m_algo ; // a reference to Algorithm
	pochoir_clone_array <N_RANK> * m_clone_array ; 
	typedef typename Algorithm<N_RANK>::queue_info queue_info ;
	//int m_extended_length [N_RANK] ;

	void initialize(grid_info<N_RANK> const & grid) 
	{
		unsigned long volume = 1 ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			int W = (grid.x1[i] - grid.x0[i]) ;
			volume *= W ;		
		}

		m_arr_leaf_kernel_table.reserve(volume) ; 
		m_arr_leaf_kernel_table.resize(volume) ; 
		//m_arr_height_leaf_kernel_table.reserve(volume) ; 
		//m_arr_height_leaf_kernel_table.resize(volume) ; 

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
	inline void symbolic_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F>
	inline void symbolic_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F>
	inline void symbolic_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F>
	inline void symbolic_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, F const & f) ;


	template <typename F, typename BF>
	inline void heterogeneous_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, F const & f, BF const & bf) ;

	template <typename F>
	inline void heterogeneous_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F, typename BF>
	inline void heterogeneous_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, F const & f, BF const & bf) ;

	template <typename F>
	inline void heterogeneous_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, F const & f) ;


	public :

	kernel_selection_trap(Algorithm<N_RANK> & alg, grid_info<N_RANK> const & grid,
				  bool power_of_two):m_algo(alg)
	{
		initialize(grid) ;
	}


	~kernel_selection_trap()
	{
		/*for (int i = 0 ; i < m_arr_height_leaf_kernel_table.size() ; i++)
		{
			height_leaf_kernel_table & table = 
				m_arr_height_leaf_kernel_table [i] ;
			for (int j = 0 ; j < table.size() ; j++)
			{
				table [j].clear() ;
			}
			table.clear() ;
		}
		m_arr_height_leaf_kernel_table.clear() ;*/
		for (int i = 0 ; i < m_arr_leaf_kernel_table.size() ; i++)
		{
			m_arr_leaf_kernel_table [i].clear() ;
		}
		m_arr_leaf_kernel_table.clear() ;
	}
	

	template <typename F, typename BF, typename P>
    inline void do_default_space_time_cuts(int t0, int t1,
        grid_info<N_RANK> const & grid, F const & f, BF const & bf, P const & p)
	{
		int T = t1 - t0 ;
		cout << "t0 " << t0 << " t1 " << t1 << endl ;
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
		if (W >= 2 * slope * T)
		{
			symbolic_space_time_cut_boundary(t0, t1, grid, p) ;
			heterogeneous_space_time_cut_boundary(t0, t1, grid, f, bf) ;
		}
		else
		{
			//choose h1 to be the normalized width
			int h1 = W / (2 * slope) ;
			int h2 = T - T / h1 * h1 ;
			cout << "h1 " << h1 << " h2 " << h2 << endl ;
			symbolic_space_time_cut_boundary(t0, t0 + h1, grid, p) ;
			if (h2 > 0)
			{
				symbolic_space_time_cut_boundary(t0, t0 + h2, grid, p) ;
			}

#ifndef NDEBUG
			//print contents
			/*for (int i = 0 ; i < m_arr_height_leaf_kernel_table.size() ; i++)
			{
				height_leaf_kernel_table & table = 
					m_arr_height_leaf_kernel_table [i] ;
				for (hlk_table_iterator it = table.begin() ; it != table.end() ;
																it++)
				{
					int h = it->first ;
					leaf_kernel_table & table2 = it->second ;
					for (lk_table_iterator it2 = table2.begin() ; 
												it2 != table2.end(); it2++)
					{
						int kernel = it2->second ;
						cout << " index " << i << " height " << h << " key " <<
								it2->first << " kernel " << kernel << endl ;
					}
				}
			}*/
#endif
			int m = T / h1 ;
			cout << " m " << m << endl ;
			for (int i = 0 ; i < m ; i++)
			{
				heterogeneous_space_time_cut_boundary(t0, t0 + h1, grid, f, bf);
				t0 += h1 ;
			}
			if (h2 > 0)
			{
				//cout << "t0 " << t0 << endl ;
				heterogeneous_space_time_cut_boundary(t0, t0 + h2, grid, f, bf);
			}
		}
	}
} ;

template<>
template <typename F>
void kernel_selection_trap<1>::compute_geneity(int h, grid_info<1> const & grid, 
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
void kernel_selection_trap<2>::compute_geneity(int h, grid_info<2> const & grid, 
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
