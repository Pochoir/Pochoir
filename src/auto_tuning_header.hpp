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
#ifndef AUTO_TUNING_HEADER_HPP
#define AUTO_TUNING_HEADER_HPP

#include "rbq_common.h"
#include <deque>
#include <unordered_map>
#include <climits> 
using namespace std ;

typedef unsigned int word_type ;

template <int N_RANK>
class zoid
{
	public :

	void resize_children(int size)
	{
		//cout << "resize child for zoid " << id << " # children " << size << endl ;
		assert (size) ;
		children = new unsigned long [size];
		num_children = size ;
		for (int i = 0 ; i < size ; i++)
		{
			children [i] = 0 ;
		}
		//cout << "resize done " << endl ;
	}

	void add_child(zoid * child, int pos, unsigned long index)
	{
		//cout << "adding child for zoid " << id << endl ; 
		assert (num_children) ;
		assert (pos < num_children) ;
		assert (child) ;
		
		//don't add the zoid as its own child
		if (this != child)
		{
			children [pos] = index ;
#ifndef NDEBUG
			child->add_parent(this->id) ;
#endif
		}
	}
	
	zoid() 
	{
		//cout << "zoid : constructor " << endl ;
		decision = 0 ; //0 for loop
		children = 0 ;
		num_children = 0 ;
		cost = 0 ;
#ifndef NDEBUG
		id = ULONG_MAX ;
#endif
	};
	
	zoid & operator = (const zoid & z)
	{
		//cout << "zoid : assignment op for zoid " << z.id << endl ;
		if (this != &z)
		{
			delete [] children ;
			decision = z.decision ;
			cost = z.cost ;
			num_children = z.num_children ;
			children = 0 ;
			height = z.height ;
			if (num_children > 0)
			{
				children = new unsigned long [num_children];
				for (int i = 0 ; i < num_children ; i++)
				{
					children [i] = z.children [i] ;
				}
			}
#ifndef NDEBUG
			id = z.id ;
			info = z.info ;
			for (int i = 0 ; i < z.parents.size() ; i++)
			{
				parents.push_back(z.parents [i]) ;
			}
#endif
		}
		return *this ;
	}
	
	zoid(const zoid & z)
	{
		decision = z.decision ;
		cost = z.cost ;
		num_children = z.num_children ;
		//cout << "zoid : copy const for zoid " << z.id << " # children" << 
		//		num_children << endl ;
		children = 0 ;
		height = z.height ;
		if (num_children > 0)
		{
			children = new unsigned long [num_children];
			for (int i = 0 ; i < num_children ; i++)
			{
				children [i] = z.children [i] ;
			}
		}
#ifndef NDEBUG
		id = z.id ;
		info = z.info ;
		for (int i = 0 ; i < z.parents.size() ; i++)
		{
			parents.push_back(z.parents [i]) ;
		}
#endif
	}


	int add_parent(unsigned long parent_id)
	{
#ifndef NDEBUG
		parents.push_back(parent_id) ;
#endif
	}
	//destructor for zoid
	~zoid()
	{
		//cout << "zoid : destructor for zoid " << id << endl ;
		num_children = 0 ;
		decision = 0 ; // 0 for looping
		cost = 0 ;
		delete [] children ;
		children = 0 ;
		//cout << "zoid : end destructor for zoid " << id << endl ;
	}
	
	private :
	char decision ;
	int height ;
	unsigned long * children ;  
	int num_children ;
	double cost ;
#ifndef NDEBUG
	grid_info <N_RANK> info ;
	unsigned long id ; //id of the zoid.
	vector<unsigned long> parents ;
#endif
} ;


template <int N_RANK>
class auto_tune
{
private:
	typedef zoid <N_RANK> zoid_type ;
	typedef unordered_multimap<unsigned long, unsigned long> hash_table ;
	typedef typename unordered_multimap<unsigned long, unsigned long>::iterator 
					hash_table_iterator ;

	void initialize(grid_info<N_RANK> const & grid, bool power_of_two)
	{
		unsigned long volume = 1 ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			volume *= (grid.x1[i] - grid.x0[i]) ;		
		}
		m_projections.reserve(volume) ;
		m_projections.resize(volume) ;
		cout << "volume " << volume << endl ;
		//unsigned long P = volume / (1 << N_RANK) ; 
		//int lgP = 8 * sizeof(unsigned long) - __builtin_clzl(P - 1) ;
		//cout << "Expected # of projections P " << P  << " lg P " << lgP << endl ;
		if (power_of_two)
		{
			m_zoids.reserve(volume / 16) ;
			cout << "Expected # of projections P " << volume / 16 << endl ;
		}
		else
		{
			m_zoids.reserve(volume / 8) ;
			cout << "Expected # of projections P " << volume / 8 << endl ;
		}
	}

	template <typename F>
	inline void build_auto_tune_dag_trap(int t0, int t1, 
						grid_info<N_RANK> const & grid, F const & f, int index)
	{
		assert (m_head [index] == ULONG_MAX) ;
		assert (m_projections.size()) ;
		//create a dummy head
		m_zoids.push_back(zoid_type ()) ;
		zoid_type & dummy_head = m_zoids [m_num_vertices] ;
#ifndef NDEBUG
		dummy_head.id = m_num_vertices ;
#endif
		dummy_head.resize_children(1) ;
		m_num_vertices++ ;
		assert (m_num_vertices == m_zoids.size()) ;
		m_head [index] = m_num_vertices ;
		cout << "t0 " << t0 << " t1 " << t1 << endl ;
		symbolic_space_time_cut_boundary(t0, t1, grid, m_num_vertices - 1, 0, f) ;
	} 

	template <typename F, typename BF>
	inline void build_auto_tune_dag_sawzoid(int t0, int t1, 
						grid_info<N_RANK> const & grid, F const & f, 
						BF const & bf, int index)
	{
		assert (m_head [index] == ULONG_MAX) ;
		assert (m_projections.size()) ;
		//create a dummy head
		m_zoids.push_back(zoid_type ()) ;
		zoid_type & dummy_head = m_zoids [m_num_vertices] ;
#ifndef NDEBUG
		dummy_head.id = m_num_vertices ;
#endif
		dummy_head.resize_children(1) ;
		m_num_vertices++ ;
		assert (m_num_vertices == m_zoids.size()) ;
		m_head [index] = m_num_vertices ;
		cout << "t0 " << t0 << " t1 " << t1 << endl ;
		double rcost = 0, ncost = 0 ;
		symbolic_sawzoid_space_time_cut_boundary(t0, t1, grid, 
						m_num_vertices - 1, 0, rcost, ncost, f, bf) ;
	}

	template <typename F>
	void compute_geneity(int h, grid_info<N_RANK> const & grid, 
						word_type & geneity, F const & f) ;

	inline void clear_projections()
	{
		//cout << "begin clear proj" << endl ;
		//cout << " size of array of hash table " << m_projections.size() << endl ;
		for (int i = 0 ; i < m_projections.size() ; i++)
		{
			//cout << "size of hash table " << i << " : " 
			//	<< m_projections [i].size() << endl ;
			m_projections [i].clear() ;		//clear the projections
		}
		//cout << "done clearing contents " << endl ;
		m_projections.clear() ;
		//empty the projections vector.
		vector<hash_table>().swap(m_projections) ; 
		//cout << "end clear proj" << endl ;
	}

	inline void print_statistics(int T)
	{
		cout << " Triangles " << endl ;
		cout << "length , # of diff projs of length " << endl ;
		int num_triangle_lengths = 0 ;
		int num_trap_lengths = 0 ;
		for (int i = 0 ; i < m_1d_count_proj_length_triangle.size() ; i++)
		{
			if (m_1d_count_proj_length_triangle [i] >= T)
			{
				num_triangle_lengths++ ; 
				cout << i << " : " << m_1d_count_proj_length_triangle [i] << endl ;
			}
		}
		cout << " # of triangle lengths > " << T << " : " << num_triangle_lengths << endl ;
		cout << endl << endl << " Trapezoids " << endl ;
		cout << "length , # of diff projs of length " << endl ;
		for (int i = 0 ; i < m_1d_count_proj_length_trapezoid.size() ; i++)
		{
			if (m_1d_count_proj_length_trapezoid [i] >= T)
			{
				num_trap_lengths++ ;
				cout << i << " : " << m_1d_count_proj_length_trapezoid [i] << endl ;
			}
		}
		cout << " # of trapezoid lengths > " << T << " : " << num_trap_lengths << endl ;
	}
	/*{
		cout << " Triangles " << endl ;
		cout << "length , # of diff projs of length " << endl ;
		for (int i = 0 ; i < m_1d_count_proj_length_triangle.size() ; i++)
		{
			if (m_1d_count_proj_length_triangle [i] > 0)
			{
				cout << i << " : " << m_1d_count_proj_length_triangle [i] << endl ;
			}
		}
		for (int i = 0 ;  i < m_1d_index_by_length_triangle.size() ; i++)
		{
			set<unsigned long> & s = m_1d_index_by_length_triangle [i] ;
			if (s.size() > 0)
			{
				cout << endl << " length " << i
				<< " # of cells " << s.size() << endl ;
			}
			for (set<unsigned long> ::iterator it = s.begin() ; it != s.end() ;
													++it)
			{
				cout << *it << "," ;
			}
		}

		cout << endl << endl << " Trapezoids " << endl ;
		cout << "length , # of diff projs of length " << endl ;
		for (int i = 0 ; i < m_1d_count_proj_length_trapezoid.size() ; i++)
		{
			if (m_1d_count_proj_length_trapezoid [i] > 0)
			{
				cout << i << " : " << m_1d_count_proj_length_trapezoid [i] << endl ;
			}
		}
		for (int i = 0 ;  i < m_1d_index_by_length_trapezoid.size() ; i++)
		{
			set<unsigned long> & s = m_1d_index_by_length_trapezoid [i] ;
			if (s.size() > 0)
			{
				cout << endl << " length " << i
				<< " # of cells " << s.size() << endl ;
			}
			for (set<unsigned long> ::iterator it = s.begin() ; it != s.end() ;
													++it)
			{
				cout << *it << "," ;
			}
		}
	
	}*/

	inline void destroy_auto_tune_dag()
	{
		//delete m_head [0] ;
		//delete m_head [1] ;
		//m_head [0] = 0 ;
		//m_head [1] = 0 ;
		m_head.clear() ;
		m_num_vertices = 0 ;
		//m_num_projections = 0 ;
		clear_projections() ;
		m_heterogeneity.clear() ;
		m_zoids.clear() ;
		vector<zoid_type>().swap(m_zoids) ; //empty the zoids vector
	}
	
	//key is the bottom volume + top volume.
	inline bool check_and_create_projection (unsigned long key, 
					int height, int centroid, unsigned long & index, 
					//int height, int centroid, zoid_type ** zoid, 
					grid_info <N_RANK> const & grid)
	{
		assert (m_projections.size()) ;
		hash_table & h = m_projections [centroid] ;
		//cout << "centroid : "  << centroid << endl ;
		assert (centroid < m_projections.size()) ;
		//cout << "searching hashtable" << endl ;
		//cout << "size hashtable " << h.size() << endl ;
		std::pair<hash_table_iterator, hash_table_iterator> p = 
													h.equal_range (key) ;
		
		//hash_table iterator has two elements, first and second.
		for (hash_table_iterator start = p.first ; start != p.second ; start++)
		{
			assert (start->first == key) ;
			assert (start->second < m_num_vertices) ;
			zoid_type * z = &(m_zoids [start->second]) ;
			if (z->height == height) 
			{
				index = start->second ;
				//cout << "found entry" << endl ;
				return true ;
			}
		}
		//cout << "not found entry" << endl ;
		//cout << "pushing zoid" << endl ;
		if (m_num_vertices > m_zoids.capacity())
		{
			cout << "# nodes of DAG " << m_num_vertices << " exceeds capacity " 				<< m_zoids.capacity() << endl ;
		}
		m_zoids.push_back(zoid_type ()) ;
		//cout << "pushed zoid" << endl ;
		zoid_type & z = m_zoids [m_num_vertices] ;
		z.height = height ;
		//assert (m_num_vertices == m_num_projections) ;
#ifndef NDEBUG
		z.info = grid ;
		z.id = m_num_vertices ;
		//m_num_projections ;
		//assert (m_num_vertices == m_num_projections) ;
		cout << "inserting zoid " << z.id << " key " << key << endl ;
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
			cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * height
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * height
			<< " h " << height << endl ; 
		}
#endif
		//*zoid = z ;
		//h.insert(std::pair<unsigned long, zoid_type *>(key, z)) ;
		h.insert(std::pair<unsigned long, unsigned long>(key, m_num_vertices)) ;
		//cout << "inserted key" << endl ;
		index = m_num_vertices ;
		//cout << "created zoid " << m_zoids [index].id << endl ;
		m_num_vertices++ ;
		assert (m_num_vertices == m_zoids.size()) ;
		
		return false ;
	}

	/*void set_clone_array(pochoir_clone_array <N_RANK> * clone_array)
	{
		m_clone_array = clone_array ;
	}*/

	void dfs(unsigned long node, vector <zoid_type> & temp_zoids,
			 vector<unsigned long> & color, unsigned long & num_vertices)
	{
		color [node] = num_vertices ; //color node gray
		zoid_type & z = m_zoids [node] ;
		if (__builtin_popcount(z.geneity) == 1)
		{
			//do not store node's children
			//cout << "push back zoid " << z.id << endl ;
			temp_zoids.push_back(zoid_type()) ;
			zoid_type & z1 = temp_zoids [num_vertices] ;
			z1.geneity = z.geneity ;
			//z1.num_children = 0 ; //set in constructor
			z1.height = z.height ;
			num_vertices++ ;
			assert (num_vertices == temp_zoids.size()) ;
#ifndef NDEBUG
			z1.info = z.info ;
			z1.id = z.id ;
#endif
		}
		else
		{
			temp_zoids.push_back(z) ; //copy the zoid z
			unsigned long index = num_vertices ; //index into the vector
			num_vertices++ ;
			assert (num_vertices == temp_zoids.size()) ;
			for (int i = 0 ; i < z.num_children ; i++)
			{
				zoid_type & z1 = temp_zoids [index] ;	
				if (color [z.children [i]] == ULONG_MAX) //node is white
				{
					z1.children [i] = num_vertices ;
					dfs(z.children [i], temp_zoids, color, num_vertices) ;
				}
				else
				{
					//node is already visited.
					//assign the child's index 
					z1.children [i] = color [z.children [i]] ;
				}
			}
		}
	}

	//remove children of all homogeneous nodes.
	void compress_dag()
	{
		vector <unsigned long> color ;
		color.reserve(m_num_vertices) ;
		color.resize(m_num_vertices) ;

		vector<zoid_type> temp_zoids ;
		vector<unsigned long> head ;
		unsigned long num_vertices = 0 ;
		for (unsigned long j = 0 ; j < m_num_vertices ; j++)
		{
			color [j] = ULONG_MAX ; //color node white
		}

		for (int j = 0 ; j < m_head.size() ; j++)
		{
			if (color [m_head [j]] == ULONG_MAX)
			{
				head.push_back(num_vertices) ;
				dfs(m_head [j], temp_zoids, color, num_vertices) ;
			}
		}
		//swap the DAG and compressed DAG
		m_zoids.swap(temp_zoids) ;
		m_head.swap(head) ;
		m_num_vertices = num_vertices ;
	}

#ifndef NDEBUG
	void print_dag()
	{
		cout << "# vertices " << m_num_vertices << endl ;
		//cout << "# vertices " << m_num_vertices << " # projections " <<
		//		m_num_projections << endl ;
		//do a BFS of the dag and print each zoid's info	
		vector <unsigned long> color ;
		color.reserve(m_num_vertices) ;
		color.resize(m_num_vertices) ;
		for (int j = 0 ; j < m_head.size() ; j++)
		{
			//if (m_head [j] == 0)
			if (m_head [j] == ULONG_MAX)
			{
				continue ;
			}
			cout << "head " << j << endl ;
			for (int i = 0 ; i < color.size() ; i++)
			{
				color [i] = 0 ;
			}
			
			color [0] = 1 ;
			//deque<zoid_type *> que ;
			deque<unsigned long> que ;
			que.push_front(m_head [j]) ;

			while (!que.empty())
			{
				//zoid_type * z = que.front() ;
				unsigned long index = que.front() ;
				assert (index < m_num_vertices) ;
				zoid_type * z = &(m_zoids[index]) ;
				cout << "\nzoid " << z->id << " height " << z->height <<
					//" num children " << z->children.size() << 
					" num children " << z->num_children << 
					" num_parents " << z->parents.size() << 
					" cost " << z->cost << endl ;
					//" num_parents " << z->parents.size() << " geneity " ;
				//print_bits(&(z->geneity), sizeof(word_type) * 8);
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
				/*if (z->geneity == 0)
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
				}*/
				//if (z->num_parents > 1)
				vector<unsigned long> & v = z->parents ;
				cout << "parents " << endl ;
				for (int i = 0 ; i < v.size() ; i++)
				{
					cout << v [i] << " " ;
				}
				cout << endl ;

				que.pop_front() ;
				//cout << "# of children " << z->children.size() << endl ;
				//for (int i = 0 ; i < z->children.size() ; i++)
				for (int i = 0 ; i < z->num_children ; i++)
				{
					//zoid_type * child = z->children[i] ;
					unsigned long index = z->children[i] ;
					cout << "child index " << index << endl ;
					assert (index < m_num_vertices) ;
					zoid_type * child = &(m_zoids [index]);
					if (index && child && color [child->id] == 0)
					{
						assert (child->id < m_num_vertices) ;
						color [child->id] = 1 ;
						que.push_back(index) ;
					}
				}
			}
		}
	}

	void print_heterogeneity()
	{
		set<word_type>::iterator begin = m_heterogeneity.begin() ;
		for (set<word_type>::iterator end = m_heterogeneity.end() ; begin != end ;
																	begin++)
		{
			print_bits(&(*begin), sizeof(word_type) * 8);
		}
	}
#endif

	template <typename F>
	inline void symbolic_abnormal_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, unsigned long,
		int child_index, F const & f) ;

	template <typename F>
	inline void symbolic_abnormal_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, unsigned long,
		int child_index, F const & f) ;

	template <typename F>
	inline void symbolic_abnormal_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned long, F const & f) ;

	template <typename F>
	inline void symbolic_abnormal_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned long, F const & f) ;

	template <typename F, typename BF>
	inline void symbolic_sawzoid_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, unsigned long,
		int child_index, double &, double &, F const & f, BF const & bf) ;

	template <typename F>
	inline void symbolic_sawzoid_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, unsigned long,
		int child_index, double &, double &, F const & f) ;

	template <typename F, typename BF>
	inline void symbolic_sawzoid_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned long, F const & f, 
		BF const & bf, int *, double &, double &) ;

	template <typename F>
	inline void symbolic_sawzoid_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned long, F const & f, int *,
		double &, double &) ;

#if 0
	template <typename F>
	inline void symbolic_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, unsigned long,
		int child_index, F const & f) ;

	template <typename F>
	inline void symbolic_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, unsigned long,
		int child_index, F const & f) ;

	template <typename F>
	inline void symbolic_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned long, F const & f) ;

	template <typename F>
	inline void symbolic_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned long, F const & f) ;
#endif

	template <typename F, typename BF>
	inline void sawzoid_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, 
		F const & f, BF const & bf) ;

	template <typename F>
	inline void sawzoid_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, F const & f) ;

	template <typename F, typename BF>
	inline void sawzoid_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, 
		F const & f, BF const & bf) ;

	template <typename F>
	inline void sawzoid_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, 
		F const & f) ;

#if 0
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
#endif

	vector<zoid_type> m_zoids ; //the array of all nodes in the DAG
	vector<hash_table> m_projections ; //the array of hashtable of <key, zoid index>
	//zoid_type * m_head [2] ; // the start nodes of the dag
	//unsigned long m_head [2] ; // the indices of start nodes in the dag
	vector<unsigned long> m_head ; // the indices of start nodes in the dag
	Algorithm <N_RANK> & m_algo ; // a reference to Algorithm
	unsigned long m_num_vertices ; //# of zoids in the dag
	//unsigned long m_num_projections ; //# of projections created
	const int NUM_BITS_IN_INT = 8 * sizeof(int) ;
	typedef typename Algorithm<N_RANK>::queue_info queue_info ;
	set<word_type> m_heterogeneity ;
	//pochoir_clone_array <N_RANK> * m_clone_array ; 
	int num_bits_width ; //# of bits to store width
	int num_bits_dim ; //# of bits for bottom and top widths in a dimension
#ifdef COUNT_PROJECTIONS
	// keeps a count of each projection length in 1D
	vector<unsigned long> m_1d_count_proj_length_triangle ; 
	vector<unsigned long> m_1d_count_proj_length_trapezoid ; 
	//set of starting points for each projection length.
    vector<set<unsigned long> > m_1d_index_by_length_triangle ;
    vector<set<unsigned long> > m_1d_index_by_length_trapezoid ;
	unsigned long m_num_triangles ;
	unsigned long m_num_trapezoids ;
#endif

	public :

	auto_tune(Algorithm<N_RANK> & alg, grid_info<N_RANK> const & grid,
				  bool power_of_two):m_algo(alg)
	{
		m_head.reserve(2) ;
		//m_head.push_back (ULONG_MAX) ;
		//m_head.push_back (ULONG_MAX) ;
		//m_head [0] = ULONG_MAX ;
		//m_head [1] = ULONG_MAX ;
		m_num_vertices = 0 ;
		//m_num_projections = 0 ;
		initialize(grid, power_of_two) ;
		num_bits_dim = sizeof (unsigned long) * 8 / N_RANK ;
		num_bits_width = sizeof (unsigned long) * 8 / (2 * N_RANK) ;
		/*if (N_RANK == 1)
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
		}*/	
	}


	~auto_tune()
	{
		//delete all zoids and clear the projections
		destroy_auto_tune_dag() ;
	}

	template <typename F, typename BF>
    inline void do_power_of_two_time_cut(int t0, int t1,
        grid_info<N_RANK> const & grid, F const & f, BF const & bf)
	{
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
		struct timeval start, end;
		double dag_time = 0. ;

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
		m_head.push_back (ULONG_MAX) ;
		gettimeofday(&start, 0);
		build_auto_tune_dag_sawzoid(t0, t0 + h1, grid, f, bf, 0) ;
		if (m_zoids [m_head [0]].decision != 0)
		{
			cout << " decided to divide and conquer " << endl ;
		}

		int offset = t0 + T / h1 * h1 ;
		int h2 = t1 - offset ;
		int index = 1 ;
		while (h2 >= 1)
		//while (h2 > m_algo.dt_recursive_)
		{
			//find index of most significant bit that is set
			index_msb = (sizeof(int) << 3) - __builtin_clz(h2) - 1 ;
			int h = 1 << index_msb ;
			//cout << "t0 " << t0 << " t1 " << t1 << 
			cout <<	" offset " << offset << " offset + h " <<
				offset + h << " h " << h << endl ;
			m_head.push_back (ULONG_MAX) ;
			build_auto_tune_dag_sawzoid(offset, offset + h, grid,
											f, bf, index) ;
			offset += h ;
			h2 = t1 - offset ;
			index++ ;
		}
#ifndef NDEBUG
		print_dag() ;
#endif
		//do not build dag if h2 = 1. Use the white clone.
		gettimeofday(&end, 0);
		dag_time = tdiff(&end, &start) ;
		cout << "# vertices " << m_num_vertices << endl ;
		cout << "DAG capacity " << m_zoids.capacity() << endl ;
		std::cout << "DAG : consumed time :" << 1.0e3 * dag_time
				<< "ms" << std::endl;
		clear_projections() ;
#ifdef GENEITY_TEST
		//compress_dag () ;
		cout << "# vertices after compression" << m_num_vertices << endl ;
		cout << "DAG capacity after compression" << m_zoids.capacity() << endl ;
#endif
		int m = T / h1 ;
		for (int i = 0 ; i < m ; i++)
		{
			cout << "t0 " << t0 << " t1 " << t1 << 
				" h1 " << h1 << " t0 + h1 " <<
				t0 + h1 << endl ;
			sawzoid_space_time_cut_boundary(t0, t0 + h1, grid, 
				&(m_zoids [m_head [0]]), f, bf) ;
			t0 += h1 ;
		}

		h2 = t1 - t0 ;
		index = 1 ;
		//time cuts happen only if height > dt_recursive_
		//while (h2 > m_algo.dt_recursive_)
		while (h2 >= 1)
		{
			//find index of most significant bit that is set
			index_msb = (sizeof(int) << 3) - __builtin_clz(h2) - 1 ;
			int h = 1 << index_msb ;
			cout << "t0 " << t0 << " t1 " << t1 << 
				" h " << h << " t0 + h " <<
				t0 + h << endl ;
			sawzoid_space_time_cut_boundary(t0, t0 + h, grid, 
				&(m_zoids [m_head [index]]), f, bf) ;
			t0 += h ;
			h2 = t1 - t0 ;
			index++ ;
		}
		/*while (h2 > 1)
		{
			//find index of most significant bit that is set
			index_msb = (sizeof(int) << 3) - __builtin_clz(h2) - 1 ;
			int h = 1 << index_msb ;
			cout << "t0 " << t0 << " t1 " << t1 << 
				" h " << h << " t0 + h " <<
				t0 + h << endl ;
			bool abnormal = false ;
			for (int i = 0 ; i < N_RANK ; i++)
			{
				int num_triangles = m_algo.dx_recursive_ [i] / ((m_algo.slope_ [i] * h) << 1) ;
				m_algo.num_triangles [i] = max(1, num_triangles) ;
				cout << "num_triangles [ " << i << " ] " << m_algo.num_triangles [i] << endl ;
			}
			m_algo.abnormal_region_space_time_cut_boundary(t0, t0 + h, grid, f, bf) ;
			t0 += h ;
			h2 = t1 - t0 ;
		}
		if (h2 == 1)
		{
			cout << "h = 1 t0 " << t0 << " t1 " << t1 << 
				 " t0 + h " << t0 + h2 << endl ;
			//base_case_kernel_boundary(t0, t0 + h2, grid, bf);
			m_algo.shorter_duo_sim_obase_bicut_p(t0, t0 + h2, grid, f, bf) ;
		}*/
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
			m_head.push_back (ULONG_MAX) ;
#ifdef COUNT_PROJECTIONS 
			if (N_RANK == 1)
			{
				grid_info<N_RANK> grid2 ;
				int dx = slope * T ;
				grid2.x0[0] = dx ;
            	grid2.dx0[0] = -slope ;
           		grid2.x1[0] = W - dx ;
            	grid2.dx1[0] = slope ;
				cout << "Trap  x0 " << dx << " x1 " << W - dx << " dx0 " <<
					-slope << " dx1 " << slope << endl ;
				build_auto_tune_dag_trap(t0, t1, grid2, p, 0) ;
				cout << "# triangles " << m_num_triangles <<
					" # of trapezoids " << m_num_trapezoids <<
					" total " << m_num_triangles + m_num_trapezoids <<
					endl ;
			}
#else
			build_auto_tune_dag_trap(t0, t1, grid, p, 0) ;
#endif
			gettimeofday(&end, 0);
			dag_time = tdiff(&end, &start) ;
			//print_dag() ;
			//print_heterogeneity() ;
			cout << "done building dag"  << endl ;
			cout << "# vertices " << m_num_vertices << endl ;
			cout << "DAG capacity " << m_zoids.capacity() << endl ;
			std::cout << "DAG : consumed time :" << 1.0e3 * dag_time
					<< "ms" << std::endl;
			clear_projections() ;
#ifdef GENEITY_TEST
			compress_dag () ;
			cout << "# vertices after compression" << m_num_vertices << endl ;
			cout << "DAG capacity after compression" << m_zoids.capacity() << endl ;
#endif
#ifdef COUNT_PROJECTIONS 
			print_statistics(T) ;
#else
			heterogeneous_space_time_cut_boundary(t0, t1, grid, 
											&(m_zoids [m_head [0]]), f, bf) ;
#endif
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
			struct timeval start, end;
			gettimeofday(&start, 0);
			m_head.push_back (ULONG_MAX) ;
			build_auto_tune_dag_trap(t0, t0 + h1, grid, p, 0) ;
			gettimeofday(&end, 0);
			dag_time = tdiff(&end, &start) ;
			cout << " t0 + T / h1 * h1  " << t0 + T / h1 * h1 << endl ;
			if (h2 > 0)
			{
				gettimeofday(&start, 0);
				m_head.push_back (ULONG_MAX) ;
				build_auto_tune_dag_trap(t0 + T / h1 * h1, t1, grid, p, 1) ;
				gettimeofday(&end, 0);
				dag_time += tdiff(&end, &start) ;
			}
			cout << "# vertices " << m_num_vertices << endl ;
			cout << "DAG capacity " << m_zoids.capacity() << endl ;
			std::cout << "DAG : consumed time :" << 1.0e3 * dag_time
					<< "ms" << std::endl;
			clear_projections() ;
			//print_dag() ;
#ifdef GENEITY_TEST
			compress_dag () ;
			cout << "# vertices after compression" << m_num_vertices << endl ;
			cout << "DAG capacity after compression" << m_zoids.capacity() << endl ;
#endif
			//print_dag() ;
			int m = T / h1 ;
			for (int i = 0 ; i < m ; i++)
			{
				//cout << "t0 " << t0 << endl ;
				heterogeneous_space_time_cut_boundary(t0, t0 + h1,  
							grid, &(m_zoids [m_head [0]]), f, bf) ;
				t0 += h1 ;
			}
			if (h2 > 0)
			{
				cout << "t0 " << t0 << endl ;
				/*gettimeofday(&start, 0);
				build_auto_tune_dag_trap(t0, t0 + h2, grid, p, 1) ;
				gettimeofday(&end, 0);
				dag_time += tdiff(&end, &start) ;*/
				//print_dag() ;
				//print_heterogeneity() ;
				heterogeneous_space_time_cut_boundary(t0, t0 + h2,  
							grid, &(m_zoids [m_head [1]]), f, bf) ;
			}
		}
	}
} ;

template<>
template <typename F>
void auto_tune<1>::compute_geneity(int h, grid_info<1> const & grid, 
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
void auto_tune<2>::compute_geneity(int h, grid_info<2> const & grid, 
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
