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
#ifndef AUTO_TUNING_HOMOGENEOUS_HEADER_HPP
#define AUTO_TUNING_HOMOGENEOUS_HEADER_HPP

#include "rbq_common.h"
#include <deque>
#include <unordered_map>
#include <climits> 
#include <time.h>
#include <cstring>
#include <fstream>
using namespace std ;

typedef unsigned int word_type ;
template <int N_RANK> class auto_tune ;

template <int N_RANK>
class zoid
{
	friend class auto_tune<N_RANK> ;
	public :
	typedef unsigned char decision_type ;
	inline void set_capacity(int size)
	{
		assert (size >= 0) ;
		if (capacity < size)
        {
            capacity = size ;
            delete [] children ;
            //if # of children increases, create and initialize a new children 
            //array
            children = new unsigned long [size];
            for (int i = 0 ; i < size ; i++)
            {
                children [i] = 0 ;
            }
        }
		assert (num_children <= capacity) ;
	}
	
	inline void resize_children(int size)
	{
		//cout << "resize child for zoid " << id << " # children " << size << endl ;
		assert (size >= 0) ;
		set_capacity(size) ;
		/*if (capacity < size)
		{
			capacity = size ;
			delete [] children ;
			//if # of children increases, create and initialize a new children 
			//array
			children = new unsigned long [size];
			for (int i = 0 ; i < size ; i++)
			{
				children [i] = 0 ;
			}
		}*/
		num_children = size ;
		assert (num_children <= capacity) ;
		//cout << "resize done " << endl ;
	}

	void add_child(zoid * child, int pos, unsigned long index)
	{
		//cout << "adding child for zoid " << id << " index " << index << 
		//		" pos " << pos << endl ; 
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
		capacity = 0 ;
		time = 0 ;
#ifdef SUBSUMPTION_SPACE
		num_level_divide = 0 ;
#endif
#ifdef SUBSUMPTION_TIME
		max_loop_time = 0 ;
#endif
		//height = 0 ;
#ifndef NDEBUG
		cache_penalty_time = 0 ;
		stime = 0 ;
		ttime = 0 ;
		ltime = 0 ;
		id = ULONG_MAX ;
#endif
	};
	
	//does a shallow copy of the contents of z.
	//does not copy children. 
	void shallow_copy(const zoid & z)
	{
		decision = z.decision ;
		time = z.time ;
#ifdef SUBSUMPTION_SPACE
		num_level_divide = z.num_level_divide ;
#endif
#ifdef SUBSUMPTION_TIME
		max_loop_time = z.max_loop_time ;
#endif
#ifndef NDEBUG
		cache_penalty_time = z.cache_penalty_time ;
		stime = z.stime ;
		ttime = z.ttime ;
		ltime = z.ltime ;
#endif
		//num_children = z.num_children ;
		num_children = 0 ;
		capacity = 0 ;
		//height = z.height ;
		//assert (num_children <= capacity) ;
#ifndef NDEBUG
		id = z.id ;
		info = z.info ;
		for (int i = 0 ; i < z.parents.size() ; i++)
		{
			parents.push_back(z.parents [i]) ;
		}
#endif
	}

	zoid & operator = (const zoid & z)
	{
		//cout << "zoid : assignment op for zoid " << z.id << endl ;
		if (this != &z)
		{
			decision = z.decision ;
			time = z.time ;
#ifdef SUBSUMPTION_SPACE
			num_level_divide = z.num_level_divide ;
#endif
#ifdef SUBSUMPTION_TIME
			max_loop_time = z.max_loop_time ;
#endif
#ifndef NDEBUG
			cache_penalty_time = z.cache_penalty_time ;
			stime = z.stime ;
			ttime = z.ttime ;
			ltime = z.ltime ;
#endif
			assert (z.num_children <= z.capacity) ;
			num_children = z.num_children ;
			//height = z.height ;
			/*if (num_children > 0)
			{
				if (capacity < z.capacity)
				{
					capacity = z.capacity ;
					delete [] children ;
					children = new unsigned long [capacity] ;
				}
				assert (children) ;
			}*/
			//resize the children array if necessary
			if (capacity < num_children)
			{
				capacity = num_children ;
				delete [] children ;
				children = new unsigned long [capacity] ;
			}
			assert (num_children <= capacity) ;
			for (int i = 0 ; i < num_children ; i++)
			{
				children [i] = z.children [i] ;
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
		time = z.time ;
#ifdef SUBSUMPTION_SPACE
		num_level_divide = z.num_level_divide ;
#endif
#ifdef SUBSUMPTION_TIME
		max_loop_time = z.max_loop_time ;
#endif
#ifndef NDEBUG
		cache_penalty_time = z.cache_penalty_time ;
		stime = z.stime ;
		ttime = z.ttime ;
		ltime = z.ltime ;
#endif
		num_children = z.num_children ;
		capacity = num_children ;
		//capacity = z.capacity ;
		assert (num_children <= capacity) ;
		//cout << "zoid : copy const for zoid " << z.id << " # children" << 
		//		num_children << endl ;
		children = 0 ;
		//height = z.height ;
		if (capacity > 0)
		{
			children = new unsigned long [capacity] ;
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


	void add_parent(unsigned long parent_id)
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
		capacity = 0 ;
		decision = 0 ; // 0 for looping
		time = 0 ;
#ifdef SUBSUMPTION_SPACE
		num_level_divide = 0 ;
#endif
#ifdef SUBSUMPTION_TIME
		max_loop_time = 0 ;
#endif
		//height = 0 ;
#ifndef NDEBUG
		cache_penalty_time = 0 ;
		stime = 0 ;
		ttime = 0 ;
		ltime = 0 ;
#endif
		delete [] children ;
		children = 0 ;
		//cout << "zoid : end destructor for zoid " << id << endl ;
	}
	
	static const int NUM_BITS_DECISION ;
	static const double FUZZ ;
	private :
	decision_type decision ;
	unsigned long * children ;  
	unsigned short capacity ;
	unsigned short num_children ;
	double time ;
#ifdef SUBSUMPTION_SPACE
	unsigned char num_level_divide ; //# of levels of consecutive divisions
#endif
#ifdef SUBSUMPTION_TIME
	double max_loop_time ;
#endif
#ifndef NDEBUG
	int height ;
	double cache_penalty_time ;
	grid_info <N_RANK> info ;
	unsigned long id ; //id of the zoid.
	vector<unsigned long> parents ;
	double stime ;
	double ttime ;
	double ltime ;
#endif
} ;

template <int N_RANK>
double const zoid<N_RANK>::FUZZ = 1 ;

template <int N_RANK>
int const zoid<N_RANK>::NUM_BITS_DECISION = sizeof(decision_type) * 8 ;
// a compact representation of zoid
template <int N_RANK>
class simple_zoid
{
	friend class auto_tune<N_RANK> ;
public :
	simple_zoid()
	{
		decision = 0 ;
		children = 0 ;
	}

	~simple_zoid()
	{
		decision = 0 ;
		delete [] children ;
		children = 0 ;
	}

	void resize_children(int size)
	{
		assert (size >= 0) ;
		assert (children == 0) ;
		if (size > 0)
		{
			children = new unsigned long [size];
		}
	}

	void resize_and_copy_children(int size, unsigned long * src)
	{
		//assert (size >= 0) ;
		//assert (children == 0) ;
		//children = new unsigned long [size];
		resize_children(size) ;
		for (int i = 0 ; i < size ; i++)
		{
			children [i] = src [i] ;
		}
	}
	typedef typename zoid<N_RANK>::decision_type decision_type ;
	//const int NUM_BITS_DECISION = zoid<N_RANK>::NUM_BITS_DECISION ;
private :
    decision_type decision ;
	unsigned long * children ;
#ifndef NDEBUG
	grid_info <N_RANK> info ;
#endif
} ;


template <int N_RANK>
class auto_tune
{
private:
	typedef zoid <N_RANK> zoid_type ;
	typedef simple_zoid <N_RANK> simple_zoid_type ;
	typedef unordered_multimap<unsigned long, unsigned long> hash_table ;
	typedef unordered_multimap<unsigned long, hash_table> two_level_hash_table ;
	typedef typename unordered_multimap<unsigned long, unsigned long>::iterator 
					hash_table_iterator ;
	typedef typename unordered_multimap<unsigned long, hash_table>::iterator 
					two_level_hash_table_iterator ;
	typedef typename zoid<N_RANK>::decision_type decision_type ;

	/*void flush_cache()
	{
		const int size = 20*1024*1024; // Allocate 20M. 
		//char *c = (char *)malloc(size);
		//for (int i = 0; i < 0xffff; i++)
		//cilk_for (int j = 0; j < size; j++)
		//	c[j] = (rand() % 1024)*j;
				//c[j] = i*j;
		//free (c) ;
		int r = rand() % 1024 ;
		memset(m_cache, r, size) ;
	}*/

	void create_simple_zoids()
	{
		cout << "num bits decision " << sizeof(decision_type) * 8 << endl ;
		assert (m_num_vertices == m_zoids.size()) ;
		m_simple_zoids.reserve (m_num_vertices) ;
		m_simple_zoids.resize (m_num_vertices) ;
		for (int i = 0 ; i < m_num_vertices ; i++)
		{
			simple_zoid_type & dest = m_simple_zoids [i] ;
			zoid_type & src = m_zoids [i] ;
			dest.decision = src.decision ;
			if (src.num_children > 0)
			{
				dest.resize_and_copy_children(src.num_children, src.children) ;
			}
#ifndef NDEBUG
			dest.info = src.info ;
#endif
		}
		//clear the contents of m_zoids.
		m_zoids.clear() ;
		vector<zoid_type> ().swap(m_zoids) ;
	}

	void fill_height_bucket(int h, int index, int level, int max_level)
	{
		if (h > 1)
		{
			//cout << "level " << level << " h " << h << endl ;
			assert (level < max_level) ;
			int floor_h = h / 2 ;
			int ceil_h = (h + 1) / 2 ;
			m_height_bucket [index][2 * level] = floor_h ;
			m_height_bucket [index][2 * level + 1] = ceil_h ;
			if (floor_h & 1)
			{
				fill_height_bucket(floor_h, index, level + 1, max_level) ;
			}
			else
			{
				fill_height_bucket(ceil_h, index, level + 1, max_level) ;
			}
		}
	}

	void initialize(grid_info<N_RANK> const & grid, int h1, int h2, 
					bool power_of_two)
	{
		assert (h1 > 0) ;
		cout << "FUZZ " << zoid_type::FUZZ << endl ;
	
		int max_level = log2(h1) + 1 ;
		cout << "max_level " << max_level << endl ;
		m_height_bucket [0].resize(2 * max_level) ;
		m_height_bucket [0][0] = h1 ;
		m_height_bucket [0][1] = h1 ;
		fill_height_bucket(h1, 0, 1, max_level) ;
	
		for (int i = 0 ; i < max_level ; i++)
		{
			cout << m_height_bucket [0] [2 * i] << " " 
				<< m_height_bucket [0] [2 * i + 1] << endl ;
		}
	
		if (h2 > 0 && h2 != h1)
		{
			max_level = log2(h2) + 1 ;
			cout << "max_level " << max_level << endl ;
			m_height_bucket [1].resize(2 * max_level) ;
			m_height_bucket [1][0] = h2 ;
			m_height_bucket [1][1] = h2 ;
			fill_height_bucket(h2, 1, 1, max_level) ;
		
			for (int i = 0 ; i < max_level ; i++)
			{
				cout << m_height_bucket [1] [2 * i] << " " 
					<< m_height_bucket [1] [2 * i + 1] << endl ;
			}
		}
		unsigned long volume = 1 ;
		//m_space_cut_mask = 0 ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			volume *= (grid.x1[i] - grid.x0[i]) ;		
			//m_space_cut_mask |= 1 << i + 1 ;
#ifdef WRITE_DAG
			file_interior << " N" << i+1 << " " << grid.x1[i] - grid.x0[i] ;	
			file_boundary << " N" << i+1 << " " << grid.x1[i] - grid.x0[i] ;	
#endif
		}
#ifdef WRITE_DAG
#ifdef TIME_INVARIANCE_INTERIOR
		file_interior << "\nTime invariant interior " << endl ;
#else
		file_interior << "\nSpace-time invariant interior " << endl ;
#endif
#ifdef TIME_INVARIANCE_BOUNDARY
		file_boundary << "\nTime invariant boundary " << endl ;
#else
		file_boundary << "\nSpace-time invariant boundary " << endl ;
#endif
#endif
		//cout << "space cut mask " << m_space_cut_mask << endl ;
		max_level = log2(h1) + 1 ;
		if (h2 > 0 && h1 != h2)
		{
			max_level += log2(h2) + 1 ;
		}
		m_projections_interior.reserve(2 * max_level) ;
		m_projections_interior.resize(2 * max_level) ; 
		int two_to_the_d = 1 << N_RANK ;
		/*for (int i = 0 ; i < two_to_the_d ; i++)
		{
			m_projections_interior [i].reserve(2 * max_level) ;
			m_projections_interior [i].resize(2 * max_level) ; 
		}*/
		for (int i = 0 ; i < two_to_the_d ; i++)
		{
			m_projections_boundary [i].reserve(2 * max_level) ;
			m_projections_boundary [i].resize(2 * max_level) ; 
		}
		cout << "volume " << volume << endl ;

		/*m_array = malloc (volume * m_type_size) ;
		if (! m_array)
		{
			cout << "auto tune :Malloc Failed " << endl ;
		}*/
		
		m_zoids.reserve(volume) ;
		/*if (power_of_two)
		{
			m_zoids.reserve(volume / 16) ;
			cout << "Expected # of projections P " << volume / 16 << endl ;
		}
		else
		{
			m_zoids.reserve(volume / 8) ;
			cout << "Expected # of projections P " << volume / 8 << endl ;
		}*/
	}

	inline void copy_data(void * dest, void * src, unsigned long length)
	{
		memcpy(dest, src, length * m_type_size) ;
	}

	template <typename F, typename BF>
	inline void build_auto_tune_dag_trap(int t0, int t1, 
						grid_info<N_RANK> const & grid, F const & f, 
						BF const & bf, int index)
	{
		assert (m_head [index] == ULONG_MAX) ;
		assert (m_projections_interior.size()) ;
		//assert (m_projections_boundary.size()) ;
		//create a dummy head
		m_zoids.push_back(zoid_type ()) ;
		zoid_type & dummy_head = m_zoids [m_num_vertices] ;
#ifndef NDEBUG
		dummy_head.id = m_num_vertices ;
#endif
		dummy_head.resize_children(1) ;
		unsigned long index_head = m_num_vertices ;
		m_num_vertices++ ;
		assert (m_num_vertices == m_zoids.size()) ;
		//m_head [index] = m_num_vertices ;
		cout << "t0 " << t0 << " t1 " << t1 << endl ;
		double rtime = 0, ntime = 0, max_loop_time = 0 ;
		symbolic_trap_space_time_cut_boundary(t0, t1, grid, 
					m_num_vertices - 1, 0, rtime, ntime, f, bf, max_loop_time) ;
		m_head [index] = m_zoids [index_head].children[0] ;
		m_zoids [index_head].resize_children (0) ;
#ifndef NDEBUG
		//remove the dummy parent of the zoid at m_head [index]
		m_zoids [m_head [index]].parents.pop_back() ;
		cout << "index " << index << " m_head [index] " << m_head [index] <<
		endl ;
		cout << " decision of head [" << index << " ] " << 
			(int) m_zoids [m_head [index]].decision 
			<< " time for space cut " << m_zoids [m_head [index]].stime << 
			" time for time cut " << m_zoids [m_head [index]].ttime << 
			" time to loop " << m_zoids [m_head [index]].ltime << endl ;
#endif
		cout << " decision of head [" << index << "] : " << 
			(int) m_zoids [m_head [index]].decision 
			<< " time " << m_zoids [m_head [index]].time * 1.0e3 << "ms" <<
			endl ;
	}

	template <typename F, typename BF>
	inline void build_auto_tune_dag_sawzoid(int t0, int t1, 
						grid_info<N_RANK> const & grid, F const & f, 
						BF const & bf, int index)
	{
		assert (m_head [index] == ULONG_MAX) ;
		assert (m_projections_interior.size()) ;
		//assert (m_projections_boundary.size()) ;
		//create a dummy head
		m_zoids.push_back(zoid_type ()) ;
		zoid_type & dummy_head = m_zoids [m_num_vertices] ;
#ifndef NDEBUG
		dummy_head.id = m_num_vertices ;
#endif
		dummy_head.resize_children(1) ;
		unsigned long index_head = m_num_vertices ;
		m_num_vertices++ ;
		assert (m_num_vertices == m_zoids.size()) ;
		//m_head [index] = m_num_vertices ;
		cout << "t0 " << t0 << " t1 " << t1 << endl ;
		double rtime = 0, ntime = 0, max_loop_time = 0 ;
		symbolic_sawzoid_space_time_cut_boundary(t0, t1, grid, 
					m_num_vertices - 1, 0, rtime, ntime, f, bf, max_loop_time);
		m_head [index] = m_zoids [index_head].children[0] ;
		m_zoids [index_head].resize_children (0) ;
#ifndef NDEBUG
		//remove the dummy parent of the zoid at m_head [index]
		m_zoids [m_head [index]].parents.pop_back() ;
		cout << "index " << index << " m_head [index] " << m_head [index] <<
		endl ;
		cout << " decision of head [" << index << " ] " << 
			(int) m_zoids [m_head [index]].decision 
			<< " time for space cut " << m_zoids [m_head [index]].stime << 
			" time for time cut " << m_zoids [m_head [index]].ttime << 
			" time to loop " << m_zoids [m_head [index]].ltime << endl ;
#endif
		cout << " decision of head [" << index << "] : " << 
			(int) m_zoids [m_head [index]].decision 
			<< " time " << m_zoids [m_head [index]].time * 1.0e3 << "ms" <<
			endl ;
	}

	inline void clear_projections()
	{
		int two_to_the_d = 1 << N_RANK ;
		/*
		for (int i = 0 ; i < two_to_the_d ; i++)
		{
			vector<hash_table> & p = m_projections_interior [i] ; 
			for (int i = 0 ; i < p.size() ; i++)
			{
				p [i].clear() ;	//clear the contents of hash table
			}
			p.clear() ; //clear the contents of vector
			//empty the vector.
			vector<hash_table>().swap(p) ; 
		}*/
		for (int i = 0 ; i < m_projections_interior.size() ; i++)
		{
			m_projections_interior [i].clear() ;//clear the contents of hash table
		}
		//clear the contents of vector
		m_projections_interior.clear() ;
		//empty the vector.
		vector<hash_table>().swap(m_projections_interior) ;

		/*for (int i = 0 ; i < m_projections_boundary.size() ; i++)
		{
			m_projections_boundary [i].clear() ;//clear the projections
		}
		m_projections_boundary.clear() ;
		//empty the projections vector.
		vector<hash_table>().swap(m_projections_boundary) ; */
		for (int k = 0 ; k < two_to_the_d ; k++)
		{
			vector <two_level_hash_table> & v = m_projections_boundary [k] ;
			for (int i = 0 ; i < v.size() ; i++)
			{
				two_level_hash_table & th = v [i] ;
				for (two_level_hash_table_iterator start = th.begin() ; 
							start != th.end() ; start++)
				{
					hash_table & h = start->second ;
					h.clear() ; //clear the contents of hash table
				}
				th.clear() ; //clear the contents of 2 level hash table.
			}
			v.clear() ;
			//empty the vector.
			vector<two_level_hash_table>().swap(v) ;
		}
	}


	inline void destroy_auto_tune_dag()
	{
		m_head.clear() ;
		m_num_vertices = 0 ;
		clear_projections() ;
		m_heterogeneity.clear() ;
		m_zoids.clear() ;
		vector<zoid_type>().swap(m_zoids) ; //empty the zoids vector
		//free (m_array) ;
	}
	
	inline bool check_and_create_time_invariant_replica(unsigned long const key,
					int const height, int const centroid, unsigned long & index,
					grid_info <N_RANK> const & grid, int dim_key)
	{

		//assert (m_projections_boundary.size()) ;
		assert (dim_key < (1 << N_RANK)) ;
		vector <two_level_hash_table> & projections_boundary = 
							m_projections_boundary [dim_key] ;
							//m_projections_boundary ;
		assert (projections_boundary.size()) ;
		bool found = false ;
		//try checking the 1st bucket
		int h1 = m_height_bucket [0][0] ;
		int k = log2((double) h1 / height) ;
		int two_to_the_k = 1 << k ; //k = 2^floor(log2 (h / h_k))
		if (height == h1 / two_to_the_k)
		{
			//height = floor (h1/2^k)
			//k is the level of the tree
			assert (height == m_height_bucket [0] [2 * k]) ;
			found = true ;
			k = 2 * k ; //index
		}
		else
		{
			//height may be ceil (h1/2^k)
			k = ceil(log2((double) h1 / height)) ;
			two_to_the_k = 1 << k ;
			if (height == (h1 + two_to_the_k - 1) / two_to_the_k)
			{
				assert (height == m_height_bucket [0] [2 * k + 1]) ;
				found = true ;
				k = 2 * k + 1 ;
			}
		}
		if (! found && ! m_height_bucket [1].empty())
		{
			//try checking the 2nd bucket
			int h2 = m_height_bucket [1][0] ;
			k = log2((double) h2 / height) ;
			two_to_the_k = 1 << k ; //k = 2^floor(log2 (h / h_k))
			int offset = 2 * (log2(h1) + 1) ;
			if (height == h2 / two_to_the_k)
			{
				//height = floor (h2/2^k)
				//k is the level of the tree
				assert (height == m_height_bucket [1] [2 * k]) ;
				found = true ;
				k = 2 * k + offset ;
			}
			else
			{
				//height may be ceil (h2/2^k)
				k = ceil(log2((double) h2 / height)) ;
				two_to_the_k = 1 << k ;
				if (height == (h2 + two_to_the_k - 1) / two_to_the_k)
				{
					assert (height == m_height_bucket [1] [2 * k + 1]) ;
					found = true ;
					k = 2 * k + 1 + offset ;
				}
			}
		}
		if (! found)
		{
			cout << "error in hash function. Height " << height << " not found "
				<< endl ;
			assert (found) ;
		}
		
		two_level_hash_table & th = projections_boundary [k] ;
//#ifndef NDEBUG
#if 0
		cout << "centroid : "  << centroid << " key " << key << endl ;
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
			cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * height
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * height
			<< " h " << height << endl ; 
		}
#endif
		std::pair<two_level_hash_table_iterator, 
				two_level_hash_table_iterator> p = th.equal_range (centroid) ;
		
		//hash_table iterator has two elements, first and second.
		for (two_level_hash_table_iterator start = p.first ; 
				start != p.second ; start++)
		{
			assert (start->first == centroid) ;
			hash_table & h = start->second ;
			std::pair<hash_table_iterator, hash_table_iterator> p1 = 
							h.equal_range (key) ;
			for (hash_table_iterator start1 = p1.first ; 
							start1 != p1.second ; start1++)
			{
				assert (start1->first == key) ;
				assert (start1->second < m_num_vertices) ;
				zoid_type * z = &(m_zoids [start1->second]) ;
				//assert (z->height == height) ;
				index = start1->second ;
#ifndef NDEBUG
				grid_info <N_RANK> grid2 = z->info ;
				int h = height ;
				for (int i = N_RANK - 1 ; i >= 0 ; i--)
				{
					bool error = false ;
					int x0 = grid.x0 [i], x1 = grid.x1 [i] ;
					int x2 = grid.x0 [i] + grid.dx0 [i] * h ;
					int x3 = grid.x1 [i] + grid.dx1 [i] * h ;

					int x0_ = grid2.x0 [i], x1_ = grid2.x1 [i] ;
					int x2_ = grid2.x0 [i] + grid2.dx0 [i] * h ; 
					int x3_ = grid2.x1 [i] + grid2.dx1 [i] * h ;
					if (dim_key == ((1 << N_RANK) - 1))
					{
						//use time invariance
						if (pmod(x0, m_algo.phys_length_ [i]) != 
							pmod(x0_, m_algo.phys_length_ [i]) ||
							pmod(x1, m_algo.phys_length_ [i]) != 
							pmod(x1_, m_algo.phys_length_ [i]) ||
							pmod(x2, m_algo.phys_length_ [i]) != 
							pmod(x2_, m_algo.phys_length_ [i]) ||
							pmod(x3, m_algo.phys_length_ [i]) != 
							pmod(x3_, m_algo.phys_length_ [i]))
						{
							error = true ;
						}
					}
					else
					{
						//use space-time invariance
						if (x1 - x0 != x1_ - x0_ || x3 - x2 != x3_ - x2_)
						{
							error = true ;
						}
					}
					if (error)
					{
						cout << "2 diff zoids hash to same key " << endl ;
						cout << "diff dim " << i << endl ;
						cout << "dim key " << dim_key << endl ;
						cout << "centroid " << centroid << endl ;
						cout << "key " << key << endl ;
						cout << " grid " << endl ;
						for (int j = N_RANK - 1 ; j >= 0 ; j--)
						{
							cout << " x0 [" << j << "] " << grid.x0 [j] 
							<< " x1 [" << j << "] " << grid.x1 [j] 
							<< " x2 [" << j << "] " << grid.x0[j] + grid.dx0[j] * h
							<< " x3 [" << j << "] " << grid.x1[j] + grid.dx1[j] * h
							<< " h " << h << endl ; 
						}
						cout << " grid 2 at index " << index << endl ;
						for (int j = N_RANK - 1 ; j >= 0 ; j--)
						{
							cout << " x0 [" << j << "] " << grid2.x0 [j] 
							<< " x1 [" << j << "] " << grid2.x1 [j] 
							<< " x2 [" << j << "] " << grid2.x0[j] + grid2.dx0[j] * h
							<< " x3 [" << j << "] " << grid2.x1[j] + grid2.dx1[j] * h
							<< " h " << h << endl ; 
						}
						assert (0) ;
					}
				}	
#endif
				return true ;
			}
			//key doesn't exist in hashtable
			m_zoids.push_back(zoid_type ()) ;
			zoid_type & z = m_zoids [m_num_vertices] ;
			//z.height = height ;
#ifndef NDEBUG
			z.info = grid ;
			z.id = m_num_vertices ;
#endif
			h.insert(std::pair<unsigned long, unsigned long>(key, 
												m_num_vertices)) ;
			index = m_num_vertices ;
			m_num_vertices++ ;
			assert (m_num_vertices == m_zoids.size()) ;
			return false ;
		}
		//centroid doesn't exist in the 2 level hash table
		m_zoids.push_back(zoid_type ()) ;
		zoid_type & z = m_zoids [m_num_vertices] ;
		//z.height = height ;
#ifndef NDEBUG
		z.info = grid ;
		z.id = m_num_vertices ;
#endif
		hash_table h ;
		h.insert(std::pair<unsigned long, unsigned long>(key, m_num_vertices)) ;
		index = m_num_vertices ;
		m_num_vertices++ ;
		assert (m_num_vertices == m_zoids.size()) ;
		projections_boundary [k].insert(
				std::pair<unsigned long, hash_table>(centroid, h)) ;
		
		return false ;
	}

	/*
	inline bool check_and_create_time_invariant_replica(unsigned long const key,
					int const height, int const centroid, unsigned long & index,
					grid_info <N_RANK> const & grid)
	{
		assert (m_projections_boundary.size()) ;
		assert (centroid < m_projections_boundary.size()) ;
		hash_table & h = m_projections_boundary [centroid] ;
//#ifndef NDEBUG
#if 0
		cout << "centroid : "  << centroid << " key " << key << endl ;
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
			cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * height
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * height
			<< " h " << height << endl ; 
		}
#endif
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
#ifndef NDEBUG
				grid_info <N_RANK> grid2 = z->info ;
				int h = height ;
				for (int i = N_RANK - 1 ; i >= 0 ; i--)
				{
					int x2 = grid.x0 [i] + grid.dx0 [i] * h ;
					int x3 = grid.x1 [i] + grid.dx1 [i] * h ;
					int x2_ = grid2.x0 [i] + grid2.dx0 [i] * h ; 
					int x3_ = grid2.x1 [i] + grid2.dx1 [i] * h ;
					if (pmod(grid.x0 [i], m_algo.phys_length_ [i]) != 
						pmod(grid2.x0 [i], m_algo.phys_length_ [i]) ||
						pmod(grid.x1 [i], m_algo.phys_length_ [i]) != 
						pmod(grid2.x1 [i], m_algo.phys_length_ [i]) ||
						pmod(x2, m_algo.phys_length_ [i]) != 
						pmod(x2_, m_algo.phys_length_ [i]) ||
						pmod(x3, m_algo.phys_length_ [i]) != 
						pmod(x3_, m_algo.phys_length_ [i]))
					{
						cout << "2 diff zoids hash to same key " << endl ;
						cout << "diff dim " << i << endl ;
						cout << "centroid " << centroid << endl ;
						cout << "key " << key << endl ;
						cout << " grid " << endl ;
						for (int j = N_RANK - 1 ; j >= 0 ; j--)
						{
							cout << " x0 [" << j << "] " << grid.x0 [j] 
							<< " x1 [" << j << "] " << grid.x1 [j] 
							<< " x2 [" << j << "] " << grid.x0[j] + grid.dx0[j] * h
							<< " x3 [" << j << "] " << grid.x1[j] + grid.dx1[j] * h
							<< " h " << h << endl ; 
						}
						cout << " grid 2 at index " << index << endl ;
						for (int j = N_RANK - 1 ; j >= 0 ; j--)
						{
							cout << " x0 [" << j << "] " << grid2.x0 [j] 
							<< " x1 [" << j << "] " << grid2.x1 [j] 
							<< " x2 [" << j << "] " << grid2.x0[j] + grid2.dx0[j] * h
							<< " x3 [" << j << "] " << grid2.x1[j] + grid2.dx1[j] * h
							<< " h " << h << endl ; 
						}
						assert (0) ;
					}
				}	
#endif
				return true ;
			}
		}
		if (m_num_vertices > m_zoids.capacity())
		{
			cout << "# nodes of DAG " << m_num_vertices << " exceeds capacity " 				<< m_zoids.capacity() << endl ;
		}
		m_zoids.push_back(zoid_type ()) ;
		zoid_type & z = m_zoids [m_num_vertices] ;
		z.height = height ;
#ifndef NDEBUG
		z.info = grid ;
		z.id = m_num_vertices ;
#endif
		m_projections_boundary [centroid].insert(std::pair<unsigned long, unsigned long>(key, m_num_vertices)) ;
		index = m_num_vertices ;
		m_num_vertices++ ;
		assert (m_num_vertices == m_zoids.size()) ;
		
		return false ;
	}*/

	inline bool check_and_create_space_time_invariant_replica(
					unsigned long const key, 
					int const height, unsigned long & index,
					grid_info <N_RANK> const & grid)
	{
		//assert (dim_key < (1 << N_RANK) - 1) ;
		vector <hash_table> & projections_interior = 
							m_projections_interior ;
							//m_projections_interior [dim_key] ;
		assert (projections_interior.size()) ;
		//int k = height - 1 ;  //index with the height for now.
		bool found = false ;
		//try checking the 1st bucket
		int h1 = m_height_bucket [0][0] ;
		int k = log2((double) h1 / height) ;
		int two_to_the_k = 1 << k ; //k = 2^floor(log2 (h / h_k))
		if (height == h1 / two_to_the_k)
		{
			//height = floor (h1/2^k)
			//k is the level of the tree
			assert (height == m_height_bucket [0] [2 * k]) ;
			found = true ;
			k = 2 * k ; //index
		}
		else
		{
			//height may be ceil (h1/2^k)
			k = ceil(log2((double) h1 / height)) ;
			two_to_the_k = 1 << k ;
			//cout << "k " << k << " 2^k " << two_to_the_k << endl ;
			//cout << "h1 " << h1 << endl ;
			//cout << "(h1 + 2^k - 1) / 2^k " << (h1 + two_to_the_k - 1) / two_to_the_k << endl ;
			if (height == (h1 + two_to_the_k - 1) / two_to_the_k)
			{
				assert (height == m_height_bucket [0] [2 * k + 1]) ;
				found = true ;
				k = 2 * k + 1 ;
			}
		}
		if (! found && ! m_height_bucket [1].empty())
		{
			//try checking the 2nd bucket
			int h2 = m_height_bucket [1][0] ;
			k = log2((double) h2 / height) ;
			two_to_the_k = 1 << k ; //k = 2^floor(log2 (h / h_k))
			int offset = 2 * (log2(h1) + 1) ;
			if (height == h2 / two_to_the_k)
			{
				//height = floor (h2/2^k)
				//k is the level of the tree
				assert (height == m_height_bucket [1] [2 * k]) ;
				found = true ;
				k = 2 * k + offset ;
			}
			else
			{
				//height may be ceil (h2/2^k)
				k = ceil(log2((double) h2 / height)) ;
				two_to_the_k = 1 << k ;
				if (height == (h2 + two_to_the_k - 1) / two_to_the_k)
				{
					assert (height == m_height_bucket [1] [2 * k + 1]) ;
					found = true ;
					k = 2 * k + 1 + offset ;
				}
			}
		}
		if (! found)
		{
			cout << "error in hash function. Height " << height << " not found "
				<< endl ;
			assert (found) ;
		}
		
		hash_table & h = projections_interior [k] ;
#if 0
		cout << "centroid : "  << centroid << " key " << key << endl ;
		cout << " height " << height << endl ;
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
			cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * height
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * height
			<< " h " << height << endl ; 
		}
#endif
		std::pair<hash_table_iterator, hash_table_iterator> p = 
													h.equal_range (key) ;
		
		//hash_table iterator has two elements, first and second.
		for (hash_table_iterator start = p.first ; start != p.second ; start++)
		{
			assert (start->first == key) ;
			assert (start->second < m_num_vertices) ;
#ifndef NDEBUG
			zoid_type * z = &(m_zoids [start->second]) ;
			//assert (z->height == height) ;
#endif
			index = start->second ;
			return true ;
		}
		if (m_num_vertices > m_zoids.capacity())
		{
			cout << "# nodes of DAG " << m_num_vertices << " exceeds capacity " 				<< m_zoids.capacity() << endl ;
		}
		m_zoids.push_back(zoid_type ()) ;
		zoid_type & z = m_zoids [m_num_vertices] ;
		//z.height = height ;
#ifndef NDEBUG
		z.info = grid ;
		z.id = m_num_vertices ;
		/*cout << "inserting zoid " << z.id << " key " << key << endl ;
		for (int i = N_RANK - 1 ; i >= 0 ; i--)
		{
			cout << " x0 [" << i << "] " << grid.x0 [i] 
			 << " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * height
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * height
			<< " h " << height << endl ; 
		}*/
#endif
		projections_interior [k].insert(std::pair<unsigned long, unsigned long>(key, m_num_vertices)) ;
		index = m_num_vertices ;
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
		//if leaf do not recurse further
		if (z.decision == 1 << (zoid_type::NUM_BITS_DECISION - 2) || 
			z.decision == 3 << (zoid_type::NUM_BITS_DECISION - 2))
		{
			//do not store node's children
			//cout << "push back zoid " << z.id << endl ;
			temp_zoids.push_back(zoid_type()) ;
			temp_zoids [num_vertices].shallow_copy(z) ;
			temp_zoids [num_vertices].resize_children(0) ;
			num_vertices++ ;
			assert (num_vertices == temp_zoids.size()) ;
		}
		else
		{
			temp_zoids.push_back(z) ; //copy the zoid z
			unsigned long index = num_vertices ; //index into the vector
			num_vertices++ ;
			assert (num_vertices == temp_zoids.size()) ;
			assert (z.num_children <= z.capacity) ;
			for (int i = 0 ; i < z.num_children ; i++)
			{
				zoid_type & z1 = temp_zoids [index] ;	
				//assert (z.children [i] > 0) ;
				assert (z.children [i] < m_num_vertices) ;
				if (color [z.children [i]] == ULONG_MAX) //node is white
				{
					z1.children [i] = num_vertices ;
					dfs(z.children [i], temp_zoids, color, num_vertices) ;
				}
				else
				{
					//node is already visited.
					assert (color [z.children [i]] < num_vertices) ;
					//assign the child's index 
					z1.children [i] = color [z.children [i]] ;
				}
			}
		}
	}

	//compress array by removing nodes that are not part of the final DAG
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
		//set color [0] = 0
		//m_zoids [0] is a dummy node.
		//a child of a zoid may also have index 0 if the
		//number of children was over allocated.
		//To avoid dfs into the dummy node, we set color [0] = 0.
		color [0] = 0 ;
		for (int j = 0 ; j < m_head.size() ; j++)
		{
			assert (m_head [j] < m_num_vertices) ;
			if (color [m_head [j]] == ULONG_MAX)
			{
				head.push_back(num_vertices) ;
				dfs(m_head [j], temp_zoids, color, num_vertices) ;
			}
			else
			{
				//node m_head [j] was already visited.
				assert (color [m_head [j]] < num_vertices) ;
				head.push_back(color [m_head [j]]) ;
			}
		}
#ifndef NDEBUG
		//set the id/parents.
		for (unsigned long j = 0 ; j < num_vertices ; j++)
		{
			zoid_type & z = temp_zoids [j] ;
			assert (z.id < m_num_vertices) ;
			assert (color [z.id] < num_vertices) ;
			z.id = color [z.id] ;
			for (int i = 0 ; i < z.parents.size() ; i++)
			{
				assert (z.parents [i] < m_num_vertices) ;
				/*if (color [z.parents [i]] >= num_vertices)
				{
					cout << "zoid " << z.id << " at index " << j <<
						": parent " << z.parents [i] <<
						" was not colored !" << endl ;
					cout << " color " << color [z.parents [i]] << endl ;
					assert (color [z.parents [i]] < num_vertices) ;
				}*/
				z.parents [i] = color [z.parents [i]] ;
			}
		}
#endif
		//swap the DAG and compressed DAG
		m_zoids.swap(temp_zoids) ;
		m_head.swap(head) ;
		m_num_vertices = num_vertices ;
	}

	bool read_dag_from_file(grid_info<N_RANK> const & grid, int T, int h1,
							double & time)
	{
		string name = m_problem_name ;
		if (name.size() == 0)
		{
			name = "auto_tune_dag" ;
		}
		char tmp [100] ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			sprintf(tmp, "_%d", grid.x1[i] - grid.x0[i]) ;
			name += tmp ;
		}
		sprintf(tmp, "_%d", T) ;
		name += tmp ;
#ifndef NDEBUG
		name += "_debug" ;
#endif
#ifdef TRAP
		name += "_trap" ;
#else
		name += "_sawzoid" ;
#endif
		ifstream dag(name.c_str()) ;
		if (! dag)
		{
			//file doesn't exist.
			cout << "dag file " << name << " doesn't exist" << endl ; 
			dag.close() ;
			return false ;
		}
		cout << " using dag file " << name << endl ;
		int h = 0 ;
		//read height
		dag >> h ;
		//cout << h << " " ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			int N = 0 ;
			//read grid dimension
			dag >> N ;
			//cout << N << " " ;
			if (N != grid.x1[i] - grid.x0[i])
			{
				return false ;
			}
		}
		//cout << endl ;
		string str ;
		//parse '# of head nodes "
		dag >> str ; dag >> str ; dag >> str ; dag >> str  ;
		int head_size ;
		//read size of head array
		dag >> head_size ;
		//cout << "# of head nodes " << head_size << endl ;
		//parse "head nodes"
		dag >> str ; dag >> str ; 
		//cout << "head nodes" << endl ;
		//vector <unsigned long> head ;
		m_head.resize(head_size) ;
		double t ;
		//read head [0]
		dag >> m_head [0] ;
		//read predicted time
		dag >> t ;
		//cout << m_head [0] <<  " " << t << endl ;
		time = t * (int) (T / h1) ;
		for (int i = 1 ; i < head_size ; i++)
		{
			dag >> m_head [i] ;
			dag >> t ;
			//cout << m_head [i] <<  " " << t << endl ;
			time += t ;
		}
		//parse # nodes
		dag >> str ; dag >> str ;
		//int num_nodes;
		//read num_nodes
		dag >> m_num_vertices ;
		//cout << "# nodes " << m_num_vertices << endl ;
		//parse nodes
		dag >> str ;
		//cout << "nodes" << endl ;
		//vector<simple_zoid_type> simple_zoids ; 
		m_simple_zoids.resize(m_num_vertices) ;
		for (int i = 0 ; i < m_num_vertices ; i++)
		{
			int decision ;
			int height, num_children ;
			dag >> decision ;
			//dag >> height ;
			dag >> num_children ;
			simple_zoid_type & z = m_simple_zoids [i] ;
			//z.set_decision(decision) ;
			z.decision = decision ;
			z.resize_children(num_children) ;
			//z.num_children = num_children ;
			//cout << decision << " " << height << " " << num_children << " " ;
			for (int i = 0 ; i < num_children ; i++)
			{
				int child ;
				dag >> child ;
				//z.set_child(i, child) ;
				z.children [i] = child ;
				//cout << child << " " ;
			}
			//cout << endl ;
#ifndef NDEBUG
			//read grid info
			grid_info<N_RANK> grid ;
			for (int i = 0 ; i < N_RANK ; i++)
			{
				dag >> grid.x0 [i] ; 
				dag >> grid.x1 [i] ;
				dag >> grid.dx0 [i] ; 
				dag >> grid.dx1 [i] ;
				//cout << grid.x0 [i] << " " << grid.x1 [i] << " " << 
				//	" " << grid.dx0 [i] << " " << grid.dx1 [i] << endl ;
			}
			//z.set_grid_info(grid) ;
			z.info = grid ;
#endif
		}
		dag.close() ;
		return true ;
	}

	void write_dag_to_file(grid_info<N_RANK> const & grid, int T)
	{
		if (m_problem_name.size() == 0)
		{
			cout << "auto tune : writing dag to file. problem name unspecified "
				<< endl ;
			m_problem_name = "auto_tune_dag" ;
		}
		char tmp [100] ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			sprintf(tmp, "_%d", grid.x1[i] - grid.x0[i]) ;
			m_problem_name += tmp ;
		}
		sprintf(tmp, "_%d", T) ;
		m_problem_name += tmp ;
#ifndef NDEBUG
		m_problem_name += "_debug" ;
#endif
#ifdef TRAP
		m_problem_name += "_trap" ;
#else
		m_problem_name += "_sawzoid" ;
#endif
		ofstream dag ;
		dag.open(m_problem_name.c_str()) ;
		//write height
		dag << T << " " ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			//write grid size
			dag << grid.x1[i] - grid.x0[i] << " " ;
		}
		dag << endl ;
		dag << "# of head nodes " << m_head.size() << endl ;
		dag << "head nodes" << endl ;
		for (int i = 0 ; i < m_head.size() ; i++)
		{
			dag << m_head [i] << " " << m_zoids [m_head [i]].time << endl ;
		}
		dag << "# nodes " << m_zoids.size() << endl ;
		dag << "nodes" << endl ;
		for (int i = 0 ; i < m_zoids.size() ; i++)
		{
			zoid_type & z = m_zoids [i] ;
			//dag << (int) z.decision << " " << z.height << " " << 
			dag << (int) z.decision << " " << //z.height << " " << 
				z.num_children << " " ;
			//cout << "z.num_children " << z.num_children << endl ; 
			for (int j = 0 ; j < z.num_children ; j++)
			{
				dag << z.children [j] << " " ;
			}
			dag << endl ;
#ifndef NDEBUG
			//write grid info
			grid_info<N_RANK> & grid = z.info ;
			for (int i = 0 ; i < N_RANK ; i++)
			{
				dag << grid.x0 [i] << " " << grid.x1 [i] << " " << 
					" " << grid.dx0 [i] << " " << grid.dx1 [i] << endl ;
			}
#endif
		}
		dag.close() ;
	}
#ifndef NDEBUG
	void print_dag()
	{
		cout << "# vertices " << m_num_vertices << endl ;
		//cout << "# vertices " << m_num_vertices << " # projections " <<
		//		m_num_projections_interior << endl ;
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
			cout << "head index " << m_head [j] << endl ;
			for (int i = 0 ; i < color.size() ; i++)
			{
				color [i] = 0 ;
			}
			assert (m_head [j] < m_num_vertices) ;
			assert (m_zoids [m_head [j]].id < m_num_vertices) ;
			//color [0] = 1 ; //looks like a bug. check it.
			//deque<zoid_type *> que ;
			deque<unsigned long> que ;
			color [m_zoids [m_head [j]].id] = 1 ;
			que.push_front(m_head [j]) ;

			while (!que.empty())
			{
				//zoid_type * z = que.front() ;
				unsigned long index = que.front() ;
				assert (index < m_num_vertices) ;
				zoid_type * z = &(m_zoids[index]) ;
				//cout << "\nid " << z->id << " h " << z->height <<
				cout << "\nid " << z->id << // " h " << z->height <<
					//" num children " << z->children.size() << 
					" # children " << z->num_children << endl ;
					//" num_parents " << z->parents.size() << 
				cout << " decision " << z->decision << 
					" time " << z->time << 
					" ptime " << z->cache_penalty_time << endl ;
				cout << " stime " << z->stime  <<
					" ttime " << z->ttime  <<
					" ltime " << z->ltime << endl ;
					//" num_parents " << z->parents.size() << " geneity " ;
				grid_info <N_RANK> & grid = z->info ;
				int h = z->height ;
				for (int i = N_RANK - 1 ; i >= 0 ; i--)
				{
					cout << " x0 [" << i << "] " << grid.x0 [i] 
					<< " x1 [" << i << "] " << grid.x1 [i] 
					<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * h
					<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * h
					<< endl ; 
				}
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
					//if (index && child && color [child->id] == 0)
					if (color [child->id] == 0)
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
	template <typename F, typename BF>
	inline void symbolic_trap_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, unsigned long,
		int child_index, double &, double &, F const & f, BF const & bf,
		double &) ;

	template <typename F>
	inline void symbolic_trap_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, unsigned long,
		int child_index, double &, double &, F const & f, double &) ;

	template <typename F, typename BF>
	inline void symbolic_trap_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned long, F const & f, 
		BF const & bf, int *, double &, double &, double &, int) ;

	template <typename F>
	inline void symbolic_trap_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned long, F const & f, int *,
		double &, double &, double &, int) ;

	template <typename F, typename BF>
	inline void trap_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, simple_zoid_type * projection_zoid, 
		F const & f, BF const & bf) ;

	template <typename F>
	inline void trap_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, simple_zoid_type * projection_zoid, 
		F const & f) ;

	template <typename F, typename BF>
	inline void trap_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, simple_zoid_type * projection_zoid, 
		F const & f, BF const & bf) ;

	template <typename F>
	inline void trap_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, simple_zoid_type * projection_zoid, 
		F const & f) ;

	template <typename F, typename BF>
	inline void symbolic_sawzoid_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, unsigned long,
		int child_index, double &, double &, F const & f, BF const & bf,
		double &) ;

	template <typename F, typename BF>
	inline void sawzoid_find_mlt_space_time_boundary(int t0, int t1, 
		grid_info<N_RANK> const & grid, zoid_type * zoid,
		double const root_dnc_time, F const & f, BF const & bf,
		double & max_loop_time, double & zoid_loop_time) ;

	template <typename F>
	inline void symbolic_sawzoid_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, unsigned long,
		int child_index, double &, double &, F const & f, double &) ;

	template <typename F>
	inline void sawzoid_find_mlt_space_time_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, zoid_type * zoid,
		double const root_dnc_time, F const & f,
		double & max_loop_time, double & zoid_loop_time) ;

	template <typename F, typename BF>
	inline void symbolic_sawzoid_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned long, F const & f, 
		BF const & bf, int *, double &, double &, double &, int) ;

	template <typename F, typename BF>
	inline void symbolic_sawzoid_space_cut_boundary_span(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned long, F const & f, 
		BF const & bf, int *, double &, double &, double &, double &, int) ;

	template <typename F, typename BF>
	inline void sawzoid_find_mlt_space_boundary(
		int t0, int t1, grid_info<N_RANK> const & grid, 
		zoid_type * projection_zoid, double const root_dnc_time, 
		F const & f, BF const & bf, double & max_loop_time, double &) ;

	template <typename F>
	inline void symbolic_sawzoid_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned long, F const & f, int *,
		double &, double &, double &, int) ;

	template <typename F>
	inline void symbolic_sawzoid_space_cut_interior_span(int t0, int t1,
		grid_info<N_RANK> const & grid, unsigned long, F const & f, int *,
		double &, double &, double &, double &, int) ;

	template <typename F>
	inline void sawzoid_find_mlt_space_interior(
		int t0, int t1, grid_info<N_RANK> const & grid, 
		zoid_type * projection_zoid, double const root_dnc_time, 
		F const & f, double & max_loop_time, double &) ;

	template <typename F, typename BF>
	inline void sawzoid_space_time_cut_boundary(int t0, int t1,  
		grid_info<N_RANK> const & grid, simple_zoid_type * projection_zoid, 
		F const & f, BF const & bf) ;

	template <typename F>
	inline void sawzoid_space_time_cut_interior(int t0, int t1, 
		grid_info<N_RANK> const & grid, simple_zoid_type * projection_zoid, 
		F const & f) ;

	template <typename F, typename BF>
	inline void sawzoid_space_cut_boundary(int t0, int t1,
		grid_info<N_RANK> const & grid, simple_zoid_type * projection_zoid, 
		F const & f, BF const & bf) ;

	template <typename F>
	inline void sawzoid_space_cut_interior(int t0, int t1,
		grid_info<N_RANK> const & grid, simple_zoid_type * projection_zoid, 
		F const & f) ;

	//vector<double> m_array ;
	//void * m_array ;
	decision_type m_space_cut_mask ;
	vector<zoid_type> m_zoids ; //the array of all nodes in the DAG
	vector<simple_zoid_type> m_simple_zoids ; //a compact array of nodes in the DAG
	//the array of hashtable of <key, zoid index> for interior region
	vector<hash_table> m_projections_interior ; 
	//the array of hashtable of <key, zoid index> for boundary region
	//vector<hash_table> m_projections_boundary ; 
	vector<two_level_hash_table> m_projections_boundary [1 << N_RANK] ; 
	vector<unsigned long> m_head ; // the indices of start nodes in the dag
	Algorithm <N_RANK> & m_algo ; // a reference to Algorithm
	unsigned long m_num_vertices ; //# of zoids in the dag
	const int NUM_BITS_IN_INT = 8 * sizeof(int) ;
	typedef typename Algorithm<N_RANK>::queue_info queue_info ;
	set<word_type> m_heterogeneity ;
	int num_bits_width ; //# of bits to store width
	int num_bits_dim ; //# of bits for bottom and top widths in a dimension
	int m_initial_height ; //initial height of the zoid
	string m_problem_name ; //name of the problem we solve.
	ofstream file_interior ; //write height, widths of a zoid in the interior
	ofstream file_boundary ; //write height, widths of a zoid at the boundary
	int m_type_size ; //size of type of date that is backed up in tuning

	vector<int> m_height_bucket [2] ; //array of vectors of height buckets.
	const int DIVIDE_COUNTER = 2 ;

	inline void sawzoid_space_cut_interior_core
		(int const, int const, int const, int const, 
		grid_info<N_RANK> const &, int, int, queue_info (*)[ALGOR_QUEUE_SIZE], 
		int *, int *, int *, int const, int const) ;

	inline void sawzoid_space_cut_boundary_core
		(int const, int const, int const, int const, 
		grid_info<N_RANK> const &, int, int, queue_info (*)[ALGOR_QUEUE_SIZE], 
		int *, int *, int *, int const, int const) ;
	public :

	auto_tune(Algorithm<N_RANK> & alg, grid_info<N_RANK> const & grid,
				  bool power_of_two, char * name, int T,
				  int type_size):m_algo(alg)
	{
		if (name != 0)
		{
			m_problem_name = name ;
		}
		else
		{
			m_problem_name = "" ;
		}
		//m_array = 0 ;
		m_type_size = type_size ;
		m_head.reserve(2) ;
		m_num_vertices = 0 ;
		//initialize(grid, T, power_of_two) ;
		num_bits_dim = sizeof (unsigned long) * 8 / N_RANK ;
		num_bits_width = sizeof (unsigned long) * 8 / (2 * N_RANK) ;

		m_space_cut_mask = 0 ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			m_space_cut_mask |= 1 << (i + 1) ;
		}
#ifdef WRITE_DAG
		time_t now = time(0);
		tm* localtm = localtime(&now);
		char time [100] ;
		strftime(time, 100, "_interior_%m%d%Y_%H%M%S", localtm) ;
		char interior_file [500], boundary_file [500] ;
		strcpy(interior_file, name) ;
		strcat(interior_file, time) ;
		cout << "interior_file " << interior_file << endl ;

		strftime(time, 100, "_boundary_%m%d%Y_%H%M%S", localtm) ;
		strcpy(boundary_file, name) ;
		strcat(boundary_file, time) ;
		cout << "boundary_file " << boundary_file << endl ;

		file_interior.open(interior_file) ;
		file_boundary.open(boundary_file) ;
		file_interior << "Problem " << name << endl ;
		file_boundary << "Problem " << name << endl ;
		file_interior << "h " << T ;
		file_boundary << "h " << T ;
#endif
	}


	~auto_tune()
	{
		//delete all zoids and clear the projections
		destroy_auto_tune_dag() ;
#ifdef WRITE_DAG
		file_interior.close() ;
		file_boundary.close() ;
#endif
	}

	template <typename F, typename BF>
    inline void do_power_of_two_time_cut(int t0, int t1,
        grid_info<N_RANK> const & grid, F const & f, BF const & bf,
		void * array, void * m_array)
	{
		//Pochoir_Array<double, N_RANK> * array = 
		//				(Pochoir_Array<double, N_RANK> *) arr ;
		int T = t1 - t0 ;
		int W = 0 ;  //max_width among all dimensions
		int slope ;
		unsigned long volume = 1 ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			volume *= m_algo.phys_length_ [i] ;
			if (m_algo.phys_length_ [i] > W)
			{
				W = m_algo.phys_length_ [i] ;
				slope = m_algo.slope_ [i] ;
			}
		}
		struct timespec start, end;
		double dag_time = 0. ;
		double expected_run_time = 0 ;
		//find index of most significant bit that is set
		int Wn = W / (slope << 1) ;
		int index_msb = (sizeof(int) << 3) - __builtin_clz(Wn) - 1 ;
		//h1 = 2^floor(lg(Wn)). The zoid with height h1 undergoes a space cut.
		int h1 = 1 << index_msb, h2, index ;
		if (T < h1)
		{
			index_msb = (sizeof(int) << 3) - __builtin_clz(T) - 1 ;
			h1 = 1 << index_msb ;
		}
		bool read_dag = false ;
		clock_gettime(CLOCK_MONOTONIC, &start) ;
		read_dag = read_dag_from_file(grid, T, h1, expected_run_time) ;
		clock_gettime(CLOCK_MONOTONIC, &end) ;
		if (read_dag)
		{
			dag_time = tdiff2(&end, &start) ;
			cout << "read dag from file " << endl ;
			cout << "# vertices " << m_num_vertices << endl ;
			cout << "DAG capacity " << m_zoids.capacity() << endl ;
			std::cout << "DAG : consumed time :" << 1.0e3 * dag_time
					<< "ms" << std::endl;
			cout << "Predicted run time " << expected_run_time * 1e3 << "ms" << endl;
		}
		else
		{
			initialize(grid, h1, h1, true) ;
			//back up data
			//copy_data(&(m_array[0]), array->data(), volume) ;
			//copy_data(m_array, array, volume) ;

			//do a dry run
    		//m_algo.power_of_two_time_cut(t0, t0 + h1, grid, f, bf) ;
			m_head.push_back (ULONG_MAX) ;
			clock_gettime(CLOCK_MONOTONIC, &start) ;
			build_auto_tune_dag_sawzoid(t0, t0 + h1, grid, f, bf, 0) ;

			int offset = t0 + T / h1 * h1 ;
			expected_run_time += m_zoids[m_head[0]].time * (int) (T / h1) ;
			h2 = t1 - offset ;
			index = 1 ;
			while (h2 >= 1)
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
				expected_run_time += m_zoids[m_head[index]].time ;
				offset += h ;
				h2 = t1 - offset ;
				index++ ;
			}
			clock_gettime(CLOCK_MONOTONIC, &end) ;
			//copy data back.
			//copy_data(array->data(), &(m_array[0]), volume) ;
			copy_data(array, m_array, volume) ;
			dag_time = tdiff2(&end, &start) ;
			cout << "# vertices " << m_num_vertices << endl ;
			cout << "DAG capacity " << m_zoids.capacity() << endl ;
			std::cout << "DAG : consumed time :" << 1.0e3 * dag_time
					<< "ms" << std::endl;
			cout << "Predicted run time " << expected_run_time * 1e3 << "ms" << endl;
			clear_projections() ;
#ifndef NDEBUG
			//print_dag() ;
#endif
			cout << "begin compress dag" << endl ;
			compress_dag () ;
			cout << "# vertices after compression " << m_num_vertices << endl ;
			cout << "DAG capacity after compression " << m_zoids.capacity() << endl ;
#ifndef NDEBUG
			//print_dag() ;
#endif
			write_dag_to_file(grid, T) ;
			create_simple_zoids() ;
		}
		double compute_time = 0. ;
		clock_gettime(CLOCK_MONOTONIC, &start) ;
		int m = T / h1 ;
		for (int i = 0 ; i < m ; i++)
		{
			cout << "t0 " << t0 << " t1 " << t1 << 
				" h1 " << h1 << " t0 + h1 " <<
				t0 + h1 << endl ;
			sawzoid_space_time_cut_boundary(t0, t0 + h1, grid, 
				&(m_simple_zoids [m_head [0]]), f, bf) ;
				//&(m_zoids [m_head [0]]), f, bf) ;
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
				&(m_simple_zoids [m_head [index]]), f, bf) ;
				//&(m_zoids [m_head [index]]), f, bf) ;
			t0 += h ;
			h2 = t1 - t0 ;
			index++ ;
		}
		clock_gettime(CLOCK_MONOTONIC, &end) ;
		compute_time = tdiff2(&end, &start) ;
		std::cout << "Compute time :" << 1.0e3 * compute_time
				<< "ms" << std::endl;
	}

	template <typename F, typename BF>
    inline void do_trap_space_time_cuts(int t0, int t1,
        grid_info<N_RANK> const & grid, F const & f, BF const & bf, 
		void * array, void * m_array)
		//Pochoir_Array<double, N_RANK> * array)
	{
		assert (t0 < t1) ;
		int T = t1 - t0 ;
		if (t0 >= t1)
		{
			return ;
		}
		cout << "t0 " << t0 << " t1 " << t1 << endl ;
		int W = 0 ;  //max_width among all dimensions
		int slope ;
		unsigned long volume = 1 ;
		for (int i = 0 ; i < N_RANK ; i++)
		{
			volume *= m_algo.phys_length_ [i] ;
			if (m_algo.phys_length_ [i] > W)
			{
				W = m_algo.phys_length_ [i] ;
				slope = m_algo.slope_ [i] ;
			}		
		}
		struct timespec start, end;
		double dag_time = 0. ;
		double expected_run_time = 0 ;
		if (W >= 2 * slope * T)
		{
			//max width is >= 2 * sigma * h implies the zoid is ready for 
			//space cuts
			bool read_dag = false ;
			clock_gettime(CLOCK_MONOTONIC, &start) ;
			read_dag = read_dag_from_file(grid, T, T, expected_run_time) ;
			clock_gettime(CLOCK_MONOTONIC, &end) ;
			if (read_dag)
			{
				dag_time = tdiff2(&end, &start) ;
				cout << "read dag from file " << endl ;
				cout << "# vertices " << m_num_vertices << endl ;
				cout << "DAG capacity " << m_zoids.capacity() << endl ;
				std::cout << "DAG : consumed time :" << 1.0e3 * dag_time
						<< "ms" << std::endl;
				cout << "Predicted run time " << expected_run_time * 1e3 << "ms" << endl;
			}
			else
			{
				initialize(grid, T, T, false) ;
				//back up data
				//copy_data(m_array, array, volume) ;
				//do a dry run
    			//m_algo.shorter_duo_sim_obase_bicut_p(t0, t1, grid, f, bf) ;
				m_head.push_back (ULONG_MAX) ;
				//cout << "m_head.size() " << m_head.size() << endl ;
				//cout << "mhead [0] " << m_head[0] << endl ;
				clock_gettime(CLOCK_MONOTONIC, &start) ;
				build_auto_tune_dag_trap(t0, t1, grid, f, bf, 0) ;
				clock_gettime(CLOCK_MONOTONIC, &end) ;
				dag_time = tdiff2(&end, &start) ;
				expected_run_time += m_zoids[m_head[0]].time ;
				cout << "# vertices " << m_num_vertices << endl ;
				cout << "DAG capacity " << m_zoids.capacity() << endl ;
				std::cout << "DAG : consumed time :" << 1.0e3 * dag_time
						<< "ms" << std::endl;
				cout << "Predicted run time " << expected_run_time * 1e3 << "ms" << endl;
				clear_projections() ;
#ifndef NDEBUG
				//print_dag() ;
#endif
				cout << "begin compress dag" << endl ;
				compress_dag () ;
				cout << "# vertices after compression" << m_num_vertices << endl ;
				cout << "DAG capacity after compression" << m_zoids.capacity() << endl ;
				//copy data back.
				//copy_data(array->data(), &(m_array[0]), volume) ;
				copy_data(array, m_array, volume) ;
				write_dag_to_file(grid, T) ;
				create_simple_zoids() ;
			}
			double compute_time = 0. ;
			clock_gettime(CLOCK_MONOTONIC, &start) ;
			trap_space_time_cut_boundary(t0, t1, grid, 
				&(m_simple_zoids [m_head [0]]), f, bf) ;
			clock_gettime(CLOCK_MONOTONIC, &end) ;
			compute_time = tdiff2(&end, &start) ;
			std::cout << "Compute time :" << 1.0e3 * compute_time
				<< "ms" << std::endl;
		}
		else
		{
			//the zoid will be cut in time until it becomes normal
			cout << "slope " << slope << endl ;
			//choose h1 to be the normalized width
			//if w < 2 * slope, no space cut is possible,
			//					choose height as 1.
			int h1 = W >= 2 * slope ? W / (2 * slope) : 1 ;
			int h2 = T - T / h1 * h1 ;
			
			bool read_dag = false ;
			clock_gettime(CLOCK_MONOTONIC, &start) ;
			read_dag = read_dag_from_file(grid, T, h1, expected_run_time) ;
			clock_gettime(CLOCK_MONOTONIC, &end) ;
			if (read_dag)
			{
				dag_time = tdiff2(&end, &start) ;
				cout << "read dag from file " << endl ;
				cout << "# vertices " << m_num_vertices << endl ;
				cout << "DAG capacity " << m_zoids.capacity() << endl ;
				std::cout << "DAG : consumed time :" << 1.0e3 * dag_time
						<< "ms" << std::endl;
				cout << "Predicted run time " << expected_run_time * 1e3 << "ms" << endl;
			}
			else
			{
				initialize(grid, h1, h2, false) ;
				//back up data
				//copy_data(m_array, array, volume) ;

				//do a dry run
    			//m_algo.shorter_duo_sim_obase_bicut_p(t0, t0 + h1, grid, f, bf) ;
				cout << "h1 " << h1 << " h2 " << h2 << endl ;
				clock_gettime(CLOCK_MONOTONIC, &start) ;
				m_head.push_back (ULONG_MAX) ;
				build_auto_tune_dag_trap(t0, t0 + h1, grid, f, bf, 0) ;
				expected_run_time += m_zoids[m_head[0]].time * (int) (T / h1) ;
				cout << " t0 + T / h1 * h1  " << t0 + T / h1 * h1 << endl ;
				if (h2 > 0)
				{
					m_head.push_back (ULONG_MAX) ;
					build_auto_tune_dag_trap(t0 + T / h1 * h1, t1, grid, f, bf, 1) ;
					expected_run_time += m_zoids[m_head[1]].time ;
				}
				clock_gettime(CLOCK_MONOTONIC, &end) ;
				dag_time = tdiff2(&end, &start) ;
				cout << "# vertices " << m_num_vertices << endl ;
				cout << "DAG capacity " << m_zoids.capacity() << endl ;
				cout << "DAG : consumed time :" << 1.0e3 * dag_time
						<< "ms" << std::endl;
				cout << "Predicted run time " << expected_run_time * 1e3 << "ms" << endl;
				clear_projections() ;
#ifndef NDEBUG
				//print_dag() ;
#endif
				cout << "begin compress dag " << endl ;
				compress_dag () ;
				cout << "# vertices after compression" << m_num_vertices << endl ;
				cout << "DAG capacity after compression" << m_zoids.capacity() << endl ;
				//copy data back.
				copy_data(array, m_array, volume) ;
				write_dag_to_file(grid, T) ;
				create_simple_zoids() ;
			}
		
			double compute_time = 0. ;
			clock_gettime(CLOCK_MONOTONIC, &start) ;
			int m = T / h1 ;
			for (int i = 0 ; i < m ; i++)
			{
				cout << "t0 " << t0 << " t1 " << t1 << 
					" h1 " << h1 << " t0 + h1 " <<
					t0 + h1 << endl ;
				trap_space_time_cut_boundary(t0, t0 + h1, grid, 
					&(m_simple_zoids [m_head [0]]), f, bf) ;
				t0 += h1 ;
			}
			if (h2 > 0)
			{
				cout << "t0 " << t0 << " t1 " << t1 << 
					" h2 " << h2 << " t0 + h2 " <<
					t0 + h2 << endl ;
				trap_space_time_cut_boundary(t0, t0 + h2, grid, 
					&(m_simple_zoids [m_head [1]]), f, bf) ;
			}
			clock_gettime(CLOCK_MONOTONIC, &end) ;
			compute_time = tdiff2(&end, &start) ;
			std::cout << "Compute time :" << 1.0e3 * compute_time
					<< "ms" << std::endl;
		}
	}
} ;

#endif
