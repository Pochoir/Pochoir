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
#ifndef GENEITY_PROBLEM_HPP
#define GENEITY_PROBLEM_HPP

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

	/*void resize_children(int size)
	{
		//cout << "resize child for zoid " << id << " # children " << size << endl ;
		assert (size) ;
		//children.resize(size, 0) ;	
		children = new unsigned long [size];
		num_children = size ;
		for (int i = 0 ; i < size ; i++)
		{
			children [i] = 0 ;
		}
		//cout << "resize done " << endl ;
	}*/

	void set_start(int s)
	{
		//set start and end indices of children
		assert (s >= 0) ;
		start = s ;
		end = s ;
	}

	int get_start()
	{
		return start ;
	}

	void set_num_children(int n)
	{
		assert (n >= 0) ;
		end = start + n ;
		assert (start <= end) ;
	}

#ifndef NDEBUG
	void add_child(zoid * child)
	{
		//cout << "adding child for zoid " << id << endl ; 
		assert (child) ;
		
		//don't add the zoid as its own child
		if (this != child)
		{
			child->add_parent(this->id) ;
		}
	}
#endif
	
	zoid() 
	{
		//cout << "zoid : constructor " << endl ;
		geneity = 0 ;
        height = 0 ;
		start = end = 0 ;
		decision = 0 ;
#ifndef NDEBUG
		id = ULONG_MAX ;
#endif
	};
	
	zoid & operator = (const zoid & z)
	{
		//cout << "zoid : assignment op for zoid " << z.id << endl ;
		if (this != &z)
		{
			geneity = z.geneity ;
			height = z.height ;
			start = z.start ;
			end = z.end ;
			decision = z.decision ;
			assert (start <= end) ;
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
		geneity = z.geneity ;
		//cout << "zoid : copy const for zoid " << z.id << " # children" << 
		//		num_children << endl ;
		height = z.height ;
		start = z.start ;
		end = z.end ;
		decision = z.decision ;
		assert (start <= end) ;
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
		start = end = 0 ;
		geneity = 0 ;
		decision = 0 ;
	}
	
	const int NUM_BITS_DECISION = sizeof(unsigned short) * 8 ;
	private :
	word_type geneity ;
	int height ;
	unsigned int start ;
	unsigned int end ;
	unsigned short decision ;
#ifndef NDEBUG
	grid_info <N_RANK> info ;
	unsigned long id ; //id of the zoid.
	vector<unsigned long> parents ;
#endif
} ;

// a compact representation of zoid
class simple_zoid
{
	word_type geneity ;
	int start; //start index of child
	int end ; //end index of child
} ;

template <int N_RANK>
class geneity_problem
{
private:
	typedef zoid <N_RANK> zoid_type ;
	//typedef multimap<unsigned long, zoid_type *> mmap ;
	//typedef typename multimap<unsigned long, zoid_type *>::iterator mmap_iterator ;
	//typedef unordered_multimap<unsigned long, zoid_type *> hash_table ;
	typedef unordered_multimap<unsigned long, unsigned int> hash_table ;
	//typedef typename unordered_multimap<unsigned long, zoid_type *>::iterator 
	typedef typename unordered_multimap<unsigned long, unsigned int>::iterator 
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
#ifdef COUNT_PROJECTIONS
		if (N_RANK == 1)
		{
			m_1d_count_proj_length_triangle.reserve(volume + 1) ;
			m_1d_count_proj_length_triangle.resize(volume + 1) ;
			m_1d_count_proj_length_trapezoid.reserve(volume + 1) ;
			m_1d_count_proj_length_trapezoid.resize(volume + 1) ;

			m_1d_index_by_length_triangle.reserve(volume + 1) ;
			m_1d_index_by_length_triangle.resize(volume + 1) ;
			m_1d_index_by_length_trapezoid.reserve(volume + 1) ;
			m_1d_index_by_length_trapezoid.resize(volume + 1) ;

			m_num_triangles = 0 ;
			m_num_trapezoids = 0 ;
		}
#endif
		//unsigned long P = volume / (1 << N_RANK) ; 
		//int lgP = 8 * sizeof(unsigned long) - __builtin_clzl(P - 1) ;
		//cout << "Expected # of projections P " << P  << " lg P " << lgP << endl ;
		if (power_of_two)
		{
			m_zoids.reserve(volume / 16) ;
            m_children.reserve(volume / 16) ;
			cout << "Expected # of projections P " << volume / 16 << endl ;
		}
		else
		{
			m_zoids.reserve(volume / 8) ;
			m_children.reserve(volume / 8) ;
            cout << "Expected # of projections P " << volume / 8 << endl ;
		}
	}

	template <typename F>
	inline void build_heterogeneity_dag(int t0, int t1, 
						grid_info<N_RANK> const & grid, F const & f, int index)
	{
		//assert (m_head [index] == (zoid_type *) 0) ;
		assert (m_head [index] == UINT_MAX) ;
		assert (m_projections.size()) ;
		//create a dummy head
		m_zoids.push_back(zoid_type ()) ;
		zoid_type & dummy_head = m_zoids [m_num_vertices] ;
#ifndef NDEBUG
		dummy_head.id = m_num_vertices ;
#endif
		//dummy_head.resize_children(1) ;
        dummy_head.set_start(0) ;
        dummy_head.set_num_children(1) ;
        m_children.push_back(0) ;
        
		m_num_vertices++ ;
		assert (m_num_vertices == m_zoids.size()) ;
		m_head [index] = m_num_vertices ;
		//assert (m_num_vertices == m_num_projections) ;
		cout << "t0 " << t0 << " t1 " << t1 << endl ;
		symbolic_space_time_cut_boundary(t0, t1, grid, m_num_vertices-1, 0, f) ;
	} 

	template <typename F>
	inline void build_heterogeneity_dag_modified(int t0, int t1, 
						grid_info<N_RANK> const & grid, F const & f, int index)
	{
		//assert (m_head [index] == (zoid_type *) 0) ;
		assert (m_head [index] == UINT_MAX) ;
		assert (m_projections.size()) ;
		//create a dummy head
		m_zoids.push_back(zoid_type ()) ;
		zoid_type & dummy_head = m_zoids [m_num_vertices] ;
#ifndef NDEBUG
		dummy_head.id = m_num_vertices ;
#endif
		//dummy_head.resize_children(1) ;
		dummy_head.set_start(0) ;
        dummy_head.set_num_children(1) ;
        m_children.push_back(0) ;
        
		m_num_vertices++ ;
		assert (m_num_vertices == m_zoids.size()) ;
		m_head [index] = m_num_vertices ;
		//assert (m_num_vertices == m_num_projections) ;
		cout << "t0 " << t0 << " t1 " << t1 << endl ;
		symbolic_modified_space_time_cut_boundary(t0, t1, grid, 
										m_num_vertices - 1, 0, f) ;
	}

	template <typename F>
	void compute_geneity(int h, grid_info<N_RANK> const & grid, 
						word_type & geneity, F const & f) ;

	inline void clear_projections()
	{
		for (int i = 0 ; i < m_projections.size() ; i++)
		{
			m_projections [i].clear() ;		//clear the projections
		}
		m_projections.clear() ;
		//empty the projections vector.
		vector<hash_table>().swap(m_projections) ; 
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

	inline void destroy_heterogeneity_dag()
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
					int height, int centroid, unsigned int & index, 
					grid_info <N_RANK> const & grid)
	{
		assert (m_projections.size()) ;
		hash_table & h = m_projections [centroid] ;
		//cout << "searching hashtable" << endl ;
		//cout << "size hashtable " << h.size() << endl ;
		std::pair<hash_table_iterator, hash_table_iterator> p = 
													h.equal_range (key) ;
		
		//hash_table iterator has two elements, first and second.
		//atmost one zoid can exist at a centroid with a given 
		//top volume + bottom volume
		for (hash_table_iterator start = p.first ; start != p.second ; start++)
		//if (p.first != p.second)
		{
			assert (start->first == key) ;
			assert (start->second < m_num_vertices) ;
			zoid_type * z = &(m_zoids [start->second]) ;
			//zoid_type * z = start->second ;
			//assert (z->height == height) ;
			if (z->height == height) 
			{
				//*zoid = z ;
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
#ifdef COUNT_PROJECTIONS
		if (N_RANK == 1)
		{
			unsigned long lb = grid.x1 [0] - grid.x0 [0] ;
			unsigned long tb = grid.x1[0] + grid.dx1[0] * height - 
								(grid.x0[0] + grid.dx0[0] * height) ; 	
			unsigned long length ;
			unsigned long index ;
			if (lb > tb)
			{
				length = lb ;
				index = grid.x0 [0] ;
			}
			else
			{
				length = tb ;
				index = grid.x0[0] + grid.dx0[0] * height ;
			}
			if (length > 0)
			{
				if (lb == 0 || tb == 0)
				{
					set <unsigned long> & s = 
						m_1d_index_by_length_triangle [length] ;
					if (s.find(index) == s.end())
					{
						//insert if index is not already in the set
						m_num_triangles++ ;
						m_1d_count_proj_length_triangle [length]++ ;
						m_1d_index_by_length_triangle [length].insert(index) ;
					}
				}			
				else
				{
					set <unsigned long> & s = 
						m_1d_index_by_length_trapezoid [length] ;
					if (s.find(index) == s.end())
					{
						//insert if index is not already in the set
						m_num_trapezoids++ ;
						m_1d_count_proj_length_trapezoid [length]++ ;
						m_1d_index_by_length_trapezoid [length].insert(index) ;
					}
				}
			}
		}
#endif
		//*zoid = z ;
		//h.insert(std::pair<unsigned long, zoid_type *>(key, z)) ;
		h.insert(std::pair<unsigned long, unsigned int>(key, m_num_vertices)) ;
		//cout << "inserted key" << endl ;
		index = m_num_vertices ;
		//cout << "created zoid " << m_zoids [index].id << endl ;
		m_num_vertices++ ;
		assert (m_num_vertices == m_zoids.size()) ;
		
		/*for (int i = 0 ; i < m_num_vertices ; i++)
		{
			cout << " zoid " << m_zoids [i].id << " address " << & (m_zoids [i]) 
				<< " # children " << m_zoids [i].num_children << endl ;
		}*/
		return false ;
		//}
	}

	void set_clone_array(pochoir_clone_array <N_RANK> * clone_array)
	{
		m_clone_array = clone_array ;
	}

	void dfs(unsigned int node, vector <zoid_type> & temp_zoids,
			 vector<unsigned int> & color, unsigned int & num_vertices,
             vector<unsigned int> & children)
	{
        //node is the index into m_zoids array
		color [node] = num_vertices ; //color node gray
		zoid_type & z = m_zoids [node] ;
		temp_zoids.push_back(z) ; //copy the zoid z
		zoid_type & z1 = temp_zoids [num_vertices] ;
		z1.start = children.size() ;
#ifndef NDEBUG
            z1.id = num_vertices ;
#endif
		num_vertices++ ;
		assert (num_vertices == temp_zoids.size()) ;
		if (__builtin_popcount(z.geneity) == 1)
		{
			//do not store node's children
            z1.set_num_children(0) ;
		}
		else
		{
			unsigned int index = num_vertices - 1 ;//store index of z1
            int num_children = z.end - z.start ;
			z1.set_num_children(num_children) ;

            for (int i = 0 ; i < num_children ; i++)
            {
                children.push_back(index);
            }
			             		
            for (unsigned int i = z.start ; i < z.end ; i++)
			{
            	zoid_type & z2 = temp_zoids [index] ;	
				unsigned int child_index = i - z.start ;
				unsigned int pos = m_children [i] ;
                assert (pos < m_num_vertices) ;
				if (color [pos] == UINT_MAX) //node is white
				{
					assert (z2.start + child_index < children.size()) ;
                    children [z2.start + child_index] = num_vertices ;
                    dfs(pos, temp_zoids, color, num_vertices, children) ;
				}
				else
				{
					//node is already visited.
					//assign the child's index 
					assert (z2.start + child_index < children.size()) ;
                    children [z2.start + child_index] = color [pos] ;
				}
			}
		}
	}

	//remove children of all homogeneous nodes.
	void compress_dag()
	{
		vector <unsigned int> color ;
		color.reserve(m_num_vertices) ;
		color.resize(m_num_vertices) ;
        
        vector <unsigned int> children ;
		vector<zoid_type> temp_zoids ;
		vector<unsigned int> head ;
		unsigned int num_vertices = 0 ;
		for (unsigned int j = 0 ; j < m_num_vertices ; j++)
		{
			color [j] = UINT_MAX ; //color node white
		}

		for (int j = 0 ; j < m_head.size() ; j++)
		{
			if (color [m_head [j]] == UINT_MAX)
			{
				head.push_back(num_vertices) ;
				dfs(m_head [j], temp_zoids, color, num_vertices, children) ;
			}
		}
		//swap the DAG and compressed DAG
		m_zoids.swap(temp_zoids) ;
		m_head.swap(head) ;
        m_children.swap(children) ;
		m_num_vertices = num_vertices ;
	}

	/*void create_simple_zoids()
	{
		assert (m_num_vertices == m_zoids.size()) ;
		m_simple_zoids.reserve (m_num_vertices) ;
		m_simple_zoids.resize (m_num_vertices) ;
		for (int i = 0 ; i < m_zoids.size() ; i++)
		{
			simple_zoid & dest = m_simple_zoids [i] ;
			zoid_type & source = m_zoids [i] ;
			dest.geneity = source.geneity ;
			if (source.num_children == 0)
			{
				dest.start = 0 ;
				dest.end = 0 ;
			}
			else
			{
				dest.start = source.children [0] ;
				dest.end = source.children[source.num_children - 1];
			}
		}
		//clear the contents of m_zoids.
		m_zoids.clear() ;
		vector<zoid_type> ().swap(m_zoids) ;
	}*/

#ifndef NDEBUG
	void print_dag()
	{
		cout << "# vertices " << m_num_vertices << endl ;
		//cout << "# vertices " << m_num_vertices << " # projections " <<
		//		m_num_projections << endl ;
		//do a BFS of the dag and print each zoid's info	
		vector <unsigned int> color ;
		color.reserve(m_num_vertices) ;
		color.resize(m_num_vertices) ;
		for (int j = 0 ; j < m_head.size() ; j++)
		{
			if (m_head [j] == UINT_MAX)
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
			deque<unsigned int> que ;
			que.push_front(m_head [j]) ;

			while (!que.empty())
			{
				//zoid_type * z = que.front() ;
				unsigned int index = que.front() ;
				assert (index < m_num_vertices) ;
				zoid_type * z = &(m_zoids[index]) ;
                int num_children = z->end - z->start ;
				cout << "\nzoid " << z->id << " height " << z->height <<
					//" num children " << z->children.size() << 
					" num children " << num_children << 
					" num_parents " << z->parents.size() << " geneity " ;
					//" num_parents " << z->num_parents << " geneity " ;
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
                for (int i = z->start ; i < z->end ; i++)
				{
                    unsigned int index = m_children [i] ;
					cout << "child index " << index << endl ;
					assert (index < m_num_vertices) ;
					zoid_type * child = &(m_zoids [index]);
					assert (child->id < m_num_vertices) ;
                    if (color [child->id] == 0)    
					{
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
		grid_info<N_RANK> const & grid, unsigned int,
		int child_index, F const & f) ;

	template <typename F>
	inline void symbolic_abnormal_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, unsigned int,
		int child_index, F const & f) ;

	template <typename F>
	inline void symbolic_abnormal_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned int, F const & f) ;

	template <typename F>
	inline void symbolic_abnormal_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned int, F const & f) ;

	template <typename F>
	inline void symbolic_modified_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, unsigned int,
		int child_index, F const & f) ;

	template <typename F>
	inline void symbolic_modified_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, unsigned int,
		int child_index, F const & f) ;

	template <typename F>
	inline void symbolic_modified_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned int, F const & f, int *) ;

	template <typename F>
	inline void symbolic_modified_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned int, F const & f, int *) ;

	template <typename F>
	inline void symbolic_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, unsigned int,
		int child_index, F const & f) ;

	template <typename F>
	inline void symbolic_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, unsigned int,
		int child_index, F const & f) ;

	template <typename F>
	inline void symbolic_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned int, F const & f) ;

	template <typename F>
	inline void symbolic_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned int, F const & f) ;

	template <typename F, typename BF>
	inline void heterogeneous_modified_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, 
		F const & f, BF const & bf) ;

	template <typename F>
	inline void heterogeneous_modified_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, F const & f) ;

	template <typename F, typename BF>
	inline void heterogeneous_modified_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, 
		F const & f, BF const & bf) ;
		//F const & f, BF const & bf, int *) ;

	template <typename F>
	inline void heterogeneous_modified_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, zoid_type * projection_zoid, 
		F const & f) ;
		//F const & f, int *) ;

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
	inline void homogeneous_modified_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F>
	inline void homogeneous_modified_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F>
	inline void homogeneous_modified_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, F const & f) ;
		//grid_info<N_RANK> const & grid, F const & f, int *) ;

	template <typename F>
	inline void homogeneous_modified_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, F const & f) ;
		//grid_info<N_RANK> const & grid, F const & f, int *) ;

	template <typename F>
	inline void homogeneous_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F>
	inline void homogeneous_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F>
	inline void homogeneous_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, F const & f) ;

	template <typename F>
	inline void homogeneous_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, F const & f) ;

	vector<zoid_type> m_zoids ; //the array of all nodes in the DAG
    // the array of indices of children of nodes in the DAG
    vector<unsigned int> m_children ; 
	//vector<simple_zoid> m_simple_zoids ; //the array of all nodes in the DAG
	vector<hash_table> m_projections ; //the array of hashtable of <key, zoid index>
	vector<unsigned int> m_head ; // the indices of start nodes in the dag
	Algorithm <N_RANK> & m_algo ; // a reference to Algorithm
	unsigned int m_num_vertices ; //# of zoids in the dag
	const int NUM_BITS_IN_INT = 8 * sizeof(int) ;
	typedef typename Algorithm<N_RANK>::queue_info queue_info ;
	set<word_type> m_heterogeneity ;
	pochoir_clone_array <N_RANK> * m_clone_array ; 
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

	geneity_problem(Algorithm<N_RANK> & alg, grid_info<N_RANK> const & grid,
				  bool power_of_two):m_algo(alg)
	{
		m_head.reserve(2) ;
		m_num_vertices = 0 ;
		initialize(grid, power_of_two) ;
		num_bits_dim = sizeof (unsigned long) * 8 / N_RANK ;
		num_bits_width = sizeof (unsigned long) * 8 / (2 * N_RANK) ;
	}


	~geneity_problem()
	{
		//delete all zoids and clear the projections
		destroy_heterogeneity_dag() ;
	}

	template <typename F, typename BF, typename P>
    inline void do_power_of_two_time_cut(int t0, int t1,
        grid_info<N_RANK> const & grid, F const & f, BF const & bf, P const & p)
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
		m_head.push_back (UINT_MAX) ;
		gettimeofday(&start, 0);
		build_heterogeneity_dag_modified(t0, t0 + h1, grid, p, 0) ;
		int offset = t0 + T / h1 * h1 ;
		int h2 = t1 - offset ;
		int index = 1 ;
		//while (h2 > 1)
		while (h2 > m_algo.dt_recursive_)
		{
			//find index of most significant bit that is set
			index_msb = (sizeof(int) << 3) - __builtin_clz(h2) - 1 ;
			int h = 1 << index_msb ;
			//cout << "t0 " << t0 << " t1 " << t1 << 
			cout <<	" offset " << offset << " offset + h " <<
				offset + h << " h " << h << endl ;
			m_head.push_back (UINT_MAX) ;
			build_heterogeneity_dag_modified(offset, offset + h, grid,
											p, index) ;
			offset += h ;
			h2 = t1 - offset ;
			index++ ;
		}
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
		//create_simple_zoids() ;
#endif
		int m = T / h1 ;
		for (int i = 0 ; i < m ; i++)
		{
			/*cout << "t0 " << t0 << " t1 " << t1 << 
				" h1 " << h1 << " t0 + h1 " <<
				t0 + h1 << endl ;*/
			heterogeneous_modified_space_time_cut_boundary(t0, t0 + h1, grid, 
				//&(m_simple_zoids [m_head [0]]), f, bf) ;
				&(m_zoids [m_head [0]]), f, bf) ;
			t0 += h1 ;
		}

		h2 = t1 - t0 ;
		index = 1 ;
		//time cuts happen only if height > dt_recursive_
		while (h2 > m_algo.dt_recursive_)
		{
			//find index of most significant bit that is set
			index_msb = (sizeof(int) << 3) - __builtin_clz(h2) - 1 ;
			int h = 1 << index_msb ;
			cout << "t0 " << t0 << " t1 " << t1 << 
				" h " << h << " t0 + h " <<
				t0 + h << endl ;
			heterogeneous_modified_space_time_cut_boundary(t0, t0 + h, grid, 
				//&(m_simple_zoids [m_head [index]]), f, bf) ;
				&(m_zoids [m_head [index]]), f, bf) ;
			t0 += h ;
			h2 = t1 - t0 ;
			index++ ;
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
			}
			m_algo.abnormal_region_space_time_cut_boundary(t0, t0 + h, grid, f, bf) ;		*/
			if (h > 2)
			{
				//continue from here
				m_algo.space_time_cut_boundary(t0, t0 + h, grid, f, bf) ;
			}
			else
			{
				m_algo.abnormal_region_space_time_cut_boundary(t0, t0 + h, grid, f, bf) ;	
			}
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
			m_head.push_back (UINT_MAX) ;
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
				build_heterogeneity_dag(t0, t1, grid2, p, 0) ;
				cout << "# triangles " << m_num_triangles <<
					" # of trapezoids " << m_num_trapezoids <<
					" total " << m_num_triangles + m_num_trapezoids <<
					endl ;
			}
#else
			build_heterogeneity_dag(t0, t1, grid, p, 0) ;
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
			//compress_dag () ;
			cout << "# vertices after compression" << m_num_vertices << endl ;
			cout << "DAG capacity after compression" << m_zoids.capacity() << endl ;
			//create_simple_zoids() ;
#endif
			//print_dag() ;
#ifdef COUNT_PROJECTIONS 
			print_statistics(T) ;
#else
			heterogeneous_space_time_cut_boundary(t0, t1, grid, 
                                    //&(m_simple_zoids [m_head [0]]), f, bf) ;
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
			m_head.push_back (UINT_MAX) ;
			build_heterogeneity_dag(t0, t0 + h1, grid, p, 0) ;
			gettimeofday(&end, 0);
			dag_time = tdiff(&end, &start) ;
			cout << " t0 + T / h1 * h1  " << t0 + T / h1 * h1 << endl ;
			if (h2 > 0)
			{
				gettimeofday(&start, 0);
				m_head.push_back (UINT_MAX) ;
				build_heterogeneity_dag(t0 + T / h1 * h1, t1, grid, p, 1) ;
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
			//compress_dag () ;
			cout << "# vertices after compression" << m_num_vertices << endl ;
			cout << "DAG capacity after compression" << m_zoids.capacity() << endl ;
			//create_simple_zoids() ;
#endif
			//print_dag() ;
			int m = T / h1 ;
			for (int i = 0 ; i < m ; i++)
			{
				//cout << "t0 " << t0 << endl ;
				heterogeneous_space_time_cut_boundary(t0, t0 + h1,  
							//grid, &(m_simple_zoids [m_head [0]]), f, bf) ;
							grid, &(m_zoids [m_head [0]]), f, bf) ;
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
							//grid, &(m_simple_zoids [m_head [1]]), f, bf) ;
							grid, &(m_zoids [m_head [1]]), f, bf) ;
			}
		}
	}
} ;

template<>
template <typename F>
void geneity_problem<1>::compute_geneity(int h, grid_info<1> const & grid, 
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
void geneity_problem<2>::compute_geneity(int h, grid_info<2> const & grid, 
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
