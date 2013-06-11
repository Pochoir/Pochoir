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
#ifndef HETEROGENEITY_HPP
#define HETEROGENEITY_HPP

#include "rbq_common.h"
#include <deque>
#include <unordered_map>
using namespace std ;

typedef unsigned int word_type ;
/*template <int N_RANK>
class zoid_structure
{
	int height ;
	zoid_strcuture * children ;
	int width [N_RANK] ;
	int converging [N_RANK] ;
} ;*/

template <int N_RANK>
class zoid
{
	public :

	void resize_children(int size)
	{
		assert (size >= 0) ;
		children.resize(size, 0) ;	
	}

	void add_child(zoid * child, int child_index)
	{
		assert (child) ;
		//don't add the zoid as its own child
		if (this != child)
		{
			child->add_parent() ;
			//children.push_back(child) ;
			children [child_index] = child ;
#ifndef NDEBUG
			if (N_RANK == 1)
			{
				assert (child_index < 3) ;
			}
			child->add_parent(this->id) ;
#endif
		}
	}
	
	zoid() ;

	void add_parent()
	{
		num_parents++ ;
	}

	void remove_parent()
	{
		num_parents-- ;
	}

	int get_parent_count()
	{
		return num_parents ;
	}

	int add_parent(unsigned long parent_id)
	{
#ifndef NDEBUG
		parents.push_back(parent_id) ;
#endif
	}
	//destructor for zoid
	//need to delete children as well.
	//but same child may be referred by more than one parent.
	//can't delete a child pointer more than once
	~zoid()
	{
		assert (num_parents == 0) ;
		//cout << "zoid : destructor " << endl ;
#ifndef NDEBUG
		if (num_parents > 0)
		{
			cout << "zoid : error -- num_parents > 0 " << endl ;
		}
#endif
		for (int i = 0 ; i < children.size() ; i++)
		{
			if (children [i])
			{
				children [i]->remove_parent() ;
				//if child has no more parents, delete it.
				if (children [i]->get_parent_count() == 0)
				{
					delete children [i] ;
				}
			}
		}	
	}
	
	private :
	word_type geneity ;
	int height ;
	int num_parents ; //# of parents for the zoid
	//zoid_structure * structure ;
	vector<zoid <N_RANK> *> children ;
#ifndef NDEBUG
	grid_info <N_RANK> info ;
	unsigned long id ; //id of the zoid.
	vector<unsigned long> parents ;
#endif
} ;

template<>
zoid<1>::zoid()
{
	num_parents = 0 ;
	geneity = 0 ;
	children.reserve(2) ;
	//children.reserve(3) ;
	//children.resize(3,0) ;
}

template<>
zoid<2>::zoid()
{
	num_parents = 0 ;
	geneity = 0 ;
	children.reserve(2) ;
	//children.reserve(9) ;
	//children.resize(9,0) ;
}
/*
template <>
class zoid_sort_criterion
{
	public :
	bool operator() (const zoid<1> & z1, const zoid<1> & z2) const
	{
		grid_info <N_RANK> const & g1 = z1.info ;
		grid_info <N_RANK> const & g2 = z2.info ;
		int x0 = g1.x0 [0] ;
		int x1 = g1.x1 [0] ; 
		int x2 = g1.x0 [0] + g1.dx0 [0] * (dt - 1) ;
	}
} ; */

template <int N_RANK>
class heterogeneity
{
private:
	typedef zoid <N_RANK> zoid_type ;
	//typedef multimap<unsigned long, zoid_type *> mmap ;
	//typedef typename multimap<unsigned long, zoid_type *>::iterator mmap_iterator ;
	typedef unordered_multimap<unsigned long, zoid_type *> hash_table ;
	typedef typename unordered_multimap<unsigned long, zoid_type *>::iterator 
					hash_table_iterator ;

	void initialize(grid_info<N_RANK> const & grid)
	{
		unsigned long volume = 1 ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			volume *= (grid.x1[i] - grid.x0[i]) ;		
		}
		m_projections.reserve(volume) ;
		m_projections.resize(volume) ;
	}

	template <typename F>
	inline void build_heterogeneity_dag(int t0, int t1, 
						grid_info<N_RANK> const & grid, F const & f, int index)
	{
		assert (m_head [index] == (zoid_type *) 0) ;
		assert (m_projections.size()) ;
		//create a dummy head
		zoid_type dummy_head ;
		cout << "t0 " << t0 << " t1 " << t1 << endl ; 
		symbolic_space_time_cut_boundary(t0, t1, grid, &dummy_head, 0, f) ;
		//the first child of dummy_head is the actual head
		m_head [index] = dummy_head.children [0] ; 
		m_head [index]->remove_parent() ; //remove the dummy parent.
#ifndef NDEBUG
		//m_head->id = m_num_vertices++ ;
#endif
		//remove m_head [index] from children of dummy
		dummy_head.children [0] = 0 ; 
		vector<zoid_type *> ().swap(dummy_head.children) ;	
		cout << "# of parents of head " << m_head [index]->get_parent_count() 
			<< endl ;
	} 

	template <typename F>
	void compute_geneity(int h, grid_info<N_RANK> const & grid, 
						word_type & geneity, F const & f) ;
	/*{
		//for (int i = 0 ; i < N_RANK ; i++)
		int lb = grid.x1[0] - grid.x0[0] ;
		int x0_top = grid.x0[0] + grid.dx0[0] * (h - 1) ;
		int x1_top = grid.x1[0] + grid.dx1[0] * (h - 1) ;
		int tb = x1_top - x0_top ;
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
		m_heterogeneity.insert(geneity) ;
	}*/

	inline void destroy_heterogeneity_dag()
	{
		delete m_head [0] ;
		delete m_head [1] ;
		m_head [0] = 0 ;
		m_head [1] = 0 ;
		m_num_vertices = 0 ;
		m_num_projections = 0 ;
		for (int i = 0 ; i < m_projections.size() ; i++)
		{
			m_projections [i].clear() ;		//clear the projections
		}
		//vector<mmap>().swap(m_projections) ; //empty the projections vector.
		vector<hash_table>().swap(m_projections) ; //empty the projections vector.
		m_heterogeneity.clear() ;
	}
	
	//key is the bottom volume + top volume.
	inline bool check_and_create_projection (unsigned long key, 
					int height, int centroid, zoid_type ** zoid, 
					grid_info <N_RANK> const & grid)
	{
		assert (m_projections.size()) ;
		hash_table & h = m_projections [centroid] ;
		std::pair<hash_table_iterator, hash_table_iterator> p = 
													h.equal_range (key) ;
		
		//hash_table iterator has two elements, first and second.
		//atmost one zoid can exist at a centroid with a given 
		//top volume + bottom volume
		for (hash_table_iterator start = p.first ; start != p.second ; start++)
		//if (p.first != p.second)
		{
			//hash_table_iterator start = p.first ;
			//hash_table_iterator end = p.second ;
			assert (start->first == key) ;
			zoid_type * z = start->second ;
			//assert (z->height == height) ;
			if (z->height == height) 
			{
				*zoid = z ;
				return true ;
			}
		}
		//else
		//{
		zoid_type * z = new zoid_type() ;
		z->height = height ;
#ifndef NDEBUG
		z->info = grid ;
		z->id = m_num_vertices++ ;
		m_num_projections++ ;
		assert (m_num_vertices == m_num_projections) ;
		/*cout << "inserting zoid " << z->id << " key " << key << endl ;
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
			cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * height
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * height
			<< " h " << height << endl ; 
		}*/
#endif
		*zoid = z ;
		h.insert(std::pair<unsigned long, zoid_type *>(key, z)) ;
		//cout << "created zoid" << endl ;
		return false ;
		//}
	}

	void set_clone_array(pochoir_clone_array <N_RANK> * clone_array)
	{
		m_clone_array = clone_array ;
	}

	/*inline bool check_and_create_projection (int volume, int height, 
					int centroid, zoid_type ** zoid, 
					grid_info <N_RANK> const & grid)
	{
		assert (m_projections.size()) ;
		mmap & m = m_projections [centroid] ;
		std::pair<mmap_iterator, mmap_iterator> p = m.equal_range (volume) ;
		
		mmap_iterator start = p.first ;
		//multimap iterator has two elements, first and second.
		for (mmap_iterator end = p.second ; start != end ; start++)
		{
			assert (start->first == volume) ;
			zoid_type * z = start->second ;
			if (z->height == height)
			{
				*zoid = z ;
				return true ;
			}
		}*/
		/*multimap<int, zoid<N_RANK>* >::iterator pos = m.lower_bound (volume) ;
		
		//verify later if a given volume has atmost two zoids
		if (pos->first == volume && pos->second->height == height)
		{
			*zoid = pos->second ;
			return true ;
		}
		pos++ ;
		//atmost two zoids exist for a given projection
		if (pos->first == volume && pos->second->height == height)
		{
			*zoid = pos->second ;
			return true ;
		}*/
		/*zoid_type * z = new zoid_type() ;
		z->height = height ;
		z->info = grid ;
#ifndef NDEBUG
		z->id = m_num_vertices++ ;
		m_num_projections++ ;
		assert (m_num_vertices == m_num_projections) ;
#endif
		*zoid = z ;
		cout << "inserting zoid " << z->id << " volume " << volume << endl ;
		m.insert(start, std::pair<unsigned long, zoid_type *>(volume, z)) ;
		cout << "created zoid" << endl ;
		return false ;
	}*/

#ifndef NDEBUG
	void print_dag()
	{
		cout << "# vertices " << m_num_vertices << " # projections " <<
				m_num_projections << endl ;
		//do a BFS of the dag and print each zoid's info	
		vector <unsigned long> color ;
		color.reserve(m_num_vertices) ;
		color.resize(m_num_vertices) ;
		for (int j = 0 ; j < 2 ; j++)
		{
			if (m_head [j] == 0)
			{
				continue ;
			}
			cout << "head " << j << endl ;
			for (int i = 0 ; i < color.size() ; i++)
			{
				color [i] = 0 ;
			}
			
			color [0] = 1 ;
			deque<zoid_type *> que ;
			que.push_front(m_head [j]) ;

			while (!que.empty())
			{
				zoid_type * z = que.front() ;
				cout << "\nzoid " << z->id << " height " << z->height <<
					" num children " << z->children.size() << 
					" num_parents " << z->num_parents << " geneity " ;
				print_bits(&(z->geneity), sizeof(word_type) * 8);
				grid_info <N_RANK> & grid = z->info ;
				int h = z->height ;
				for (int i = N_RANK - 1 ; i >= 0 ; i--)
				{
					cout << " x0 [" << i << "] " << grid.x0 [i] 
					 << " x1 [" << i << "] " << grid.x1 [i] 
					<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * h
					<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * h
					<< " h " << h << endl ; 
				}
				if (z->geneity == 0)
				{
					int bottom_volume = 1, top_volume = 1 ;
					for (int i = 0 ; i < N_RANK ; i++)
					{
						int lb = grid.x1[i] - grid.x0[i] ;
						int x0_top = grid.x0[i] + grid.dx0[i] * (h - 1) ;
						int x1_top = grid.x1[i] + grid.dx1[i] * (h - 1) ;
						int tb = x1_top - x0_top ;
						bottom_volume *= lb ;
						top_volume *= tb ;
					}
					if (bottom_volume != 0 || top_volume != 0)
					{
						//zoid is not empty. It must have a geneity
						cout << "Error : geneity is 0 " << endl ;
						assert (z->geneity) ;
					}
				}
				if (z->num_parents > 1)
				{
					vector<unsigned long> & v = z->parents ;
					cout << "parents " << endl ;
					for (int i = 0 ; i < v.size() ; i++)
					{
						cout << v [i] << " " ;
					}
					cout << endl ;
				}
				que.pop_front() ;
				cout << "# of children " << z->children.size() << endl ;
				for (int i = 0 ; i < z->children.size() ; i++)
				{
					zoid_type * child = z->children[i] ;
					if (child && color [child->id] == 0)
					{
						color [child->id] = 1 ;
						que.push_back(child) ;
					}
				}
			}
		}
	}

	void print_heterogeneity()
	{
		cout << "# heterogeneity " << m_heterogeneity.size() << endl ;
		set<word_type>::iterator begin = m_heterogeneity.begin() ;
		for (set<word_type>::iterator end = m_heterogeneity.end() ; 
			 begin != end ; begin++)
		{
			print_bits(&(*begin), sizeof(word_type) * 8);
		}
	}
#endif

	template <typename F>
	inline void symbolic_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, zoid_type * parent, 
		int child_index, F const & f) ;

	template <typename F>
	inline void symbolic_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, zoid_type * parent, 
		int child_index, F const & f) ;

	template <typename F>
	inline void symbolic_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, zoid_type * parent, F const & f) ;

	template <typename F>
	inline void symbolic_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, zoid_type * parent, F const & f) ;

	template <typename F, typename BF>
	inline void heterogeneous_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, 
		F const & f, BF const & bf) ;

	template <typename F>
	inline void heterogeneous_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, F const & f) ;

	template <typename F, typename BF>
	inline void heterogeneous_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, 
		F const & f, BF const & bf) ;

	template <typename F>
	inline void heterogeneous_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, F const & f) ;

	template <typename F>
	inline void homogeneous_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, 
		F const & f) ;

	template <typename F>
	inline void homogeneous_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, F const & f) ;

	template <typename F>
	inline void homogeneous_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, 
		F const & f) ;

	template <typename F>
	inline void homogeneous_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, F const & f) ;

	//vector<mmap> m_projections ;
	vector<hash_table> m_projections ;
	zoid_type * m_head [2] ; // the start nodes of the dag
	Algorithm <N_RANK> & m_algo ; // a reference to Algorithm
	unsigned long m_num_vertices ; //# of zoids in the dag
	unsigned long m_num_projections ; //# of projections created
	const int NUM_BITS_IN_INT = 8 * sizeof(int) ;
	typedef typename Algorithm<N_RANK>::queue_info queue_info ;
	set<word_type> m_heterogeneity ;
	pochoir_clone_array <N_RANK> * m_clone_array ; 
	int num_bits_width ; //# of bits to store width
	int num_bits_dim ; //# of bits for bottom and top widths in a dimension

	public :

	heterogeneity(Algorithm<N_RANK> & alg, grid_info<N_RANK> const & grid):
										m_algo(alg)
	{
		m_head [0] = (zoid_type *) 0 ;
		m_head [1] = (zoid_type *) 0 ;
		m_num_vertices = 0 ;
		m_num_projections = 0 ;
		initialize(grid) ;
		if (N_RANK == 1)
		{
			num_bits_width = sizeof (unsigned long) * 8 / 2 ;
			num_bits_dim = sizeof (unsigned long) * 8 ; 
		}
		else if (N_RANK == 2)
		{
			num_bits_width = sizeof (unsigned long) * 8 / 4 ;
			num_bits_dim = sizeof (unsigned long) * 8 / 2 ; 
		}
		else if (N_RANK == 3)
		{
			num_bits_width = sizeof (unsigned long) * 8 / 6 ;
			num_bits_dim = sizeof (unsigned long) * 8 / 3 ; 
		}	
	}


	~heterogeneity()
	{
		//delete all zoids and clear the projections
		destroy_heterogeneity_dag() ;
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
		struct timeval start, end;
		double dag_time = 0. ;
		if (W >= 2 * slope * T)
		{
			//max width is >= 2 * sigma * h implies the zoid is ready for space cuts
			gettimeofday(&start, 0);
			build_heterogeneity_dag(t0, t1, grid, p, 0) ;
			gettimeofday(&end, 0);
			dag_time = tdiff(&end, &start) ;
			//print_dag() ;
			//print_heterogeneity() ;
			cout << "done building dag"  << endl ;
			heterogeneous_space_time_cut_boundary(t0, t1, grid, m_head [0], f, 
													bf) ;
		}
		else
		{
			//the zoid will be cut in time until it becomes normal
			//compute 2^k where k = ceil(lg (2 * slope * T / W))
			//int k = 8 * sizeof(int) - __builtin_clz((T - 1) * 2 * slope / W) ;
			//int k = 8 * sizeof(int) - __builtin_clz((T * 2 * slope - 1) / W) ;
			//int k = 8 * sizeof(int) - __builtin_clz(T * 2 * slope / W - 1) ; 
			//cout << "k " << k << endl ;
			//int two_to_the_k = 1 << k ; 
			//cout << "width " << W << " T " << T <<  " 2^k "  << two_to_the_k << endl ;
			//h1 = floor (T / two_to_the_k)
			//int h1 = T / two_to_the_k ;
			cout << "slope " << slope << endl ;
			//choose h1 to be the normalized width
			int h1 = W / (2 * slope) ;
			int h2 = T - T / h1 * h1 ;
			cout << "h1 " << h1 << " h2 " << h2 << endl ;
			struct timeval start, end ;
			gettimeofday(&start, 0);
			build_heterogeneity_dag(t0, t0 + h1, grid, p, 0) ;
			gettimeofday(&end, 0);
			dag_time = tdiff(&end, &start) ;
			cout << " t0 + T / h1 * h1  " << t0 + T / h1 * h1 << endl ;
			if (h2 > 0)
			{
				gettimeofday(&start, 0);
				build_heterogeneity_dag(t0 + T / h1 * h1, t1, grid, p, 1) ;
				gettimeofday(&end, 0);
				dag_time += tdiff(&end, &start) ;
			}
			std::cout << "DAG : consumed time :" << 1.0e3 * dag_time
					<< "ms" << std::endl;
			int m = T / h1 ;
			for (int i = 0 ; i < m ; i++)
			{
				//cout << "t0 " << t0 << endl ;
				heterogeneous_space_time_cut_boundary(t0, t0 + h1,  
							grid, m_head [0], f, bf) ;
				t0 += h1 ;
			}
			if (h2 > 0)
			{
				cout << "t0 " << t0 << endl ;
				/*gettimeofday(&start, 0);
				build_heterogeneity_dag(t0, t0 + h2, grid, p, 1) ;
				gettimeofday(&end, 0);
				dag_time += tdiff(&end, &start) ;*/
				//print_dag() ;
				//print_heterogeneity() ;
				heterogeneous_space_time_cut_boundary(t0, t0 + h2,  
							grid, m_head [1], f, bf) ;
			}
		}
	}
} ;

template<>
template <typename F>
void heterogeneity<1>::compute_geneity(int h, grid_info<1> const & grid, 
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
	m_heterogeneity.insert(geneity) ;
}

template<>
template <typename F>
void heterogeneity<2>::compute_geneity(int h, grid_info<2> const & grid, 
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
	m_heterogeneity.insert(geneity) ;
}

#endif
