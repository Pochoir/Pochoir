/*
 **********************************************************************************
 *  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 *  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 * 		                     Charles E. Leiserson <cel@mit.edu>
 * 	 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   Suggestsions:                  yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 ********************************************************************************/

#ifndef POCHOIR_WALK_H
#define POCHOIR_WALK_H

#include <cstdlib>
#include <cstdio>
#include <cstdio>
#include <cassert>
// #include <iostream>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h>
// #include "pochoir_types.hpp"
#include "pochoir_common.hpp"
#include <map>

using namespace std;

template <int N_RANK>
struct meta_grid_boundary {
    template <typename BF>
	static inline void single_step(int t, Grid_Info<N_RANK> const & grid, Grid_Info<N_RANK> const & initial_grid, BF const & bf); 
};

template <>
struct meta_grid_boundary <8>{
    template <typename BF>
	static inline void single_step(int t, Grid_Info<8> const & grid, Grid_Info<8> const & initial_grid, BF const & bf) {
        /* add cilk_for here will only lower the performance */
		for (int i = grid.x0[7]; i < grid.x1[7]; ++i) {
            int new_i = pmod_lu(i, initial_grid.x0[7], initial_grid.x1[7]);
			for (int j = grid.x0[6]; j < grid.x1[6]; ++j) {
                int new_j = pmod_lu(j, initial_grid.x0[6], initial_grid.x1[6]);
        for (int k = grid.x0[5]; k < grid.x1[5]; ++k) {
            int new_k = pmod_lu(k, initial_grid.x0[5], initial_grid.x1[5]);
            for (int l = grid.x0[4]; l < grid.x1[4]; ++l) {
                int new_l = pmod_lu(l, initial_grid.x0[4], initial_grid.x1[4]);
        for (int m = grid.x0[3]; m < grid.x1[3]; ++m) {
            int new_m = pmod_lu(m, initial_grid.x0[3], initial_grid.x1[3]);
            for (int n = grid.x0[2]; n < grid.x1[2]; ++n) {
                int new_n = pmod_lu(n, initial_grid.x0[2], initial_grid.x1[2]);
        for (int o = grid.x0[1]; o < grid.x1[1]; ++o) {
            int new_o = pmod_lu(o, initial_grid.x0[1], initial_grid.x1[1]);
            for (int p = grid.x0[0]; p < grid.x1[0]; ++p) {
                int new_p = pmod_lu(p, initial_grid.x0[0], initial_grid.x1[0]);
#ifdef CHECK_SHAPE
                if (inRun) {
                    home_cell_[8] = new_p; 
                    home_cell_[7] = new_o; home_cell_[6] = new_n; 
                    home_cell_[5] = new_m; home_cell_[4] = new_l; 
                    home_cell_[3] = new_k; home_cell_[2] = new_j; 
                    home_cell_[1] = new_i; home_cell_[0] = t;
                }
#endif
                bf(t, new_i, new_j, new_k, new_l, new_m, new_n, new_o, new_p);
            } } } } } } } }
    }
};

template <>
struct meta_grid_boundary <7>{
    template <typename BF>
	static inline void single_step(int t, Grid_Info<7> const & grid, Grid_Info<7> const & initial_grid, BF const & bf) {
        /* add cilk_for here will only lower the performance */
		for (int i = grid.x0[6]; i < grid.x1[6]; ++i) {
            int new_i = pmod_lu(i, initial_grid.x0[6], initial_grid.x1[6]);
			for (int j = grid.x0[5]; j < grid.x1[5]; ++j) {
                int new_j = pmod_lu(j, initial_grid.x0[5], initial_grid.x1[5]);
        for (int k = grid.x0[4]; k < grid.x1[4]; ++k) {
            int new_k = pmod_lu(k, initial_grid.x0[4], initial_grid.x1[4]);
            for (int l = grid.x0[3]; l < grid.x1[3]; ++l) {
                int new_l = pmod_lu(l, initial_grid.x0[3], initial_grid.x1[3]);
        for (int m = grid.x0[2]; m < grid.x1[2]; ++m) {
            int new_m = pmod_lu(m, initial_grid.x0[2], initial_grid.x1[2]);
            for (int n = grid.x0[1]; n < grid.x1[1]; ++n) {
                int new_n = pmod_lu(n, initial_grid.x0[1], initial_grid.x1[1]);
        for (int o = grid.x0[0]; o < grid.x1[0]; ++o) {
            int new_o = pmod_lu(o, initial_grid.x0[0], initial_grid.x1[0]);
#ifdef CHECK_SHAPE
                if (inRun) {
                    home_cell_[7] = new_o; home_cell_[6] = new_n; 
                    home_cell_[5] = new_m; home_cell_[4] = new_l; 
                    home_cell_[3] = new_k; home_cell_[2] = new_j; 
                    home_cell_[1] = new_i; home_cell_[0] = t;
                }
#endif
                bf(t, new_i, new_j, new_k, new_l, new_m, new_n, new_o);
            } } } } } } }
    }
};

template <>
struct meta_grid_boundary <6>{
    template <typename BF>
	static inline void single_step(int t, Grid_Info<6> const & grid, Grid_Info<6> const & initial_grid, BF const & bf) {
        /* add cilk_for here will only lower the performance */
		for (int i = grid.x0[5]; i < grid.x1[5]; ++i) {
            int new_i = pmod_lu(i, initial_grid.x0[5], initial_grid.x1[5]);
			for (int j = grid.x0[4]; j < grid.x1[4]; ++j) {
                int new_j = pmod_lu(j, initial_grid.x0[4], initial_grid.x1[4]);
        for (int k = grid.x0[3]; k < grid.x1[3]; ++k) {
            int new_k = pmod_lu(k, initial_grid.x0[3], initial_grid.x1[3]);
            for (int l = grid.x0[2]; l < grid.x1[2]; ++l) {
                int new_l = pmod_lu(l, initial_grid.x0[2], initial_grid.x1[2]);
        for (int m = grid.x0[1]; m < grid.x1[1]; ++m) {
            int new_m = pmod_lu(m, initial_grid.x0[1], initial_grid.x1[1]);
            for (int n = grid.x0[0]; n < grid.x1[0]; ++n) {
                int new_n = pmod_lu(n, initial_grid.x0[0], initial_grid.x1[0]);
#ifdef CHECK_SHAPE
                if (inRun) {
                    home_cell_[6] = new_n; 
                    home_cell_[5] = new_m; home_cell_[4] = new_l; 
                    home_cell_[3] = new_k; home_cell_[2] = new_j; 
                    home_cell_[1] = new_i; home_cell_[0] = t;
                }
#endif
                bf(t, new_i, new_j, new_k, new_l, new_m, new_n);
            } } } } } } 
    }
};

template <>
struct meta_grid_boundary <5>{
    template <typename BF>
	static inline void single_step(int t, Grid_Info<5> const & grid, Grid_Info<5> const & initial_grid, BF const & bf) {
        /* add cilk_for here will only lower the performance */
		for (int i = grid.x0[4]; i < grid.x1[4]; ++i) {
            int new_i = pmod_lu(i, initial_grid.x0[4], initial_grid.x1[4]);
			for (int j = grid.x0[3]; j < grid.x1[3]; ++j) {
                int new_j = pmod_lu(j, initial_grid.x0[3], initial_grid.x1[3]);
        for (int k = grid.x0[2]; k < grid.x1[2]; ++k) {
            int new_k = pmod_lu(k, initial_grid.x0[2], initial_grid.x1[2]);
            for (int l = grid.x0[1]; l < grid.x1[1]; ++l) {
                int new_l = pmod_lu(l, initial_grid.x0[1], initial_grid.x1[1]);
        for (int m = grid.x0[0]; m < grid.x1[0]; ++m) {
            int new_m = pmod_lu(m, initial_grid.x0[0], initial_grid.x1[0]);
#ifdef CHECK_SHAPE
                if (inRun) {
                    home_cell_[5] = new_m; home_cell_[4] = new_l; 
                    home_cell_[3] = new_k; home_cell_[2] = new_j; 
                    home_cell_[1] = new_i; home_cell_[0] = t;
                }
#endif
                bf(t, new_i, new_j, new_k, new_l, new_m);
            } } } } } 
    }
};

template <>
struct meta_grid_boundary <4>{
    template <typename BF>
	static inline void single_step(int t, Grid_Info<4> const & grid, Grid_Info<4> const & initial_grid, BF const & bf) {
        /* add cilk_for here will only lower the performance */
		for (int i = grid.x0[3]; i < grid.x1[3]; ++i) {
            int new_i = pmod_lu(i, initial_grid.x0[3], initial_grid.x1[3]);
			for (int j = grid.x0[2]; j < grid.x1[2]; ++j) {
                int new_j = pmod_lu(j, initial_grid.x0[2], initial_grid.x1[2]);
        for (int k = grid.x0[1]; k < grid.x1[1]; ++k) {
            int new_k = pmod_lu(k, initial_grid.x0[1], initial_grid.x1[1]);
            for (int l = grid.x0[0]; l < grid.x1[0]; ++l) {
                int new_l = pmod_lu(l, initial_grid.x0[0], initial_grid.x1[0]);
#ifdef CHECK_SHAPE
                if (inRun) {
                    home_cell_[4] = new_l; 
                    home_cell_[3] = new_k; home_cell_[2] = new_j; 
                    home_cell_[1] = new_i; home_cell_[0] = t;
                }
#endif
                bf(t, new_i, new_j, new_k, new_l);
            } } } } 
    } 
};

template <>
struct meta_grid_boundary <3>{
    template <typename BF>
	static inline void single_step(int t, Grid_Info<3> const & grid, Grid_Info<3> const & initial_grid, BF const & bf) {
        /* add cilk_for here will only lower the performance */
		for (int i = grid.x0[2]; i < grid.x1[2]; ++i) {
            int new_i = pmod_lu(i, initial_grid.x0[2], initial_grid.x1[2]);
			for (int j = grid.x0[1]; j < grid.x1[1]; ++j) {
                int new_j = pmod_lu(j, initial_grid.x0[1], initial_grid.x1[1]);
        for (int k = grid.x0[0]; k < grid.x1[0]; ++k) {
            int new_k = pmod_lu(k, initial_grid.x0[0], initial_grid.x1[0]);
#ifdef CHECK_SHAPE
                if (inRun) {
                    home_cell_[3] = new_k; home_cell_[2] = new_j; 
                    home_cell_[1] = new_i; home_cell_[0] = t;
                } 
#endif
                bf(t, new_i, new_j, new_k);
        } } }
	} 
};

template <>
struct meta_grid_boundary <2>{
    template <typename BF>
	static inline void single_step(int t, Grid_Info<2> const & grid, Grid_Info<2> const & initial_grid, BF const & bf) {
		for (int i = grid.x0[1]; i < grid.x1[1]; ++i) {
#if (KLEIN == 0)
            int new_i = pmod_lu(i, initial_grid.x0[1], initial_grid.x1[1]);
#endif
			for (int j = grid.x0[0]; j < grid.x1[0]; ++j) {
#if (KLEIN == 0)
                int new_j = pmod_lu(j, initial_grid.x0[0], initial_grid.x1[0]);
#else
                int new_i = i, new_j = j;
                klein(new_i, new_j, initial_grid);
#endif
#ifdef CHECK_SHAPE
                if (inRun) {
                    home_cell_[2] = new_j; 
                    home_cell_[1] = new_i; home_cell_[0] = t;
                } 
#endif
                bf(t, new_i, new_j);
			} }
	} 
};

template <>
struct meta_grid_boundary <1>{
    template <typename BF>
	static inline void single_step(int t, Grid_Info<1> const & grid, Grid_Info<1> const & initial_grid, BF const & bf) {
		for (int i = grid.x0[0]; i < grid.x1[0]; ++i) {
            int new_i = pmod_lu(i, initial_grid.x0[0], initial_grid.x1[0]);
#ifdef CHECK_SHAPE
            if (inRun) {
                home_cell_[1] = new_i; home_cell_[0] = t;
            }
#endif
		    bf(t, new_i);
        }
	} 
};

template <int N_RANK>
struct meta_grid_interior {
    template <typename F>
	static inline void single_step(int t, Grid_Info<N_RANK> const & grid, Grid_Info<N_RANK> const & initial_grid, F const & f); 
};

template <>
struct meta_grid_interior <8>{
    template <typename F>
	static inline void single_step(int t, Grid_Info<8> const & grid, Grid_Info<8> const & initial_grid, F const & f) {
        /* add cilk_for here will only lower the performance */
		for (int i = grid.x0[7]; i < grid.x1[7]; ++i) {
			for (int j = grid.x0[6]; j < grid.x1[6]; ++j) {
        for (int k = grid.x0[5]; k < grid.x1[5]; ++k) {
            for (int l = grid.x0[4]; l < grid.x1[4]; ++l) {
        for (int m = grid.x0[3]; m < grid.x1[3]; ++m) {
            for (int n = grid.x0[2]; n < grid.x1[2]; ++n) {
        for (int o = grid.x0[1]; o < grid.x1[1]; ++o) {
            for (int p = grid.x0[0]; p < grid.x1[0]; ++p) {
                f(t, i, j, k, l, m, n, o, p);
            } } } } } } } }
    }
};

template <>
struct meta_grid_interior <7>{
    template <typename F>
	static inline void single_step(int t, Grid_Info<7> const & grid, Grid_Info<7> const & initial_grid, F const & f) {
        /* add cilk_for here will only lower the performance */
		for (int i = grid.x0[6]; i < grid.x1[6]; ++i) {
			for (int j = grid.x0[5]; j < grid.x1[5]; ++j) {
        for (int k = grid.x0[4]; k < grid.x1[4]; ++k) {
            for (int l = grid.x0[3]; l < grid.x1[3]; ++l) {
        for (int m = grid.x0[2]; m < grid.x1[2]; ++m) {
            for (int n = grid.x0[1]; n < grid.x1[1]; ++n) {
        for (int o = grid.x0[0]; o < grid.x1[0]; ++o) {
                f(t, i, j, k, l, m, n, o);
            } } } } } } }
    }
};

template <>
struct meta_grid_interior <6>{
    template <typename F>
	static inline void single_step(int t, Grid_Info<6> const & grid, Grid_Info<6> const & initial_grid, F const & f) {
        /* add cilk_for here will only lower the performance */
		for (int i = grid.x0[5]; i < grid.x1[5]; ++i) {
			for (int j = grid.x0[4]; j < grid.x1[4]; ++j) {
        for (int k = grid.x0[3]; k < grid.x1[3]; ++k) {
            for (int l = grid.x0[2]; l < grid.x1[2]; ++l) {
        for (int m = grid.x0[1]; m < grid.x1[1]; ++m) {
            for (int n = grid.x0[0]; n < grid.x1[0]; ++n) {
                f(t, i, j, k, l, m, n);
            } } } } } }
    } 
};

template <>
struct meta_grid_interior <5>{
    template <typename F>
	static inline void single_step(int t, Grid_Info<5> const & grid, Grid_Info<5> const & initial_grid, F const & f) {
        /* add cilk_for here will only lower the performance */
		for (int i = grid.x0[4]; i < grid.x1[4]; ++i) {
			for (int j = grid.x0[3]; j < grid.x1[3]; ++j) {
        for (int k = grid.x0[2]; k < grid.x1[2]; ++k) {
            for (int l = grid.x0[1]; l < grid.x1[1]; ++l) {
        for (int m = grid.x0[0]; m < grid.x1[0]; ++m) {
                f(t, i, j, k, l, m);
            } } } } } 
    }
};

template <>
struct meta_grid_interior <4>{
    template <typename F>
	static inline void single_step(int t, Grid_Info<4> const & grid, Grid_Info<4> const & initial_grid, F const & f) {
        /* add cilk_for here will only lower the performance */
		for (int i = grid.x0[3]; i < grid.x1[3]; ++i) {
			for (int j = grid.x0[2]; j < grid.x1[2]; ++j) {
        for (int k = grid.x0[1]; k < grid.x1[1]; ++k) {
            for (int l = grid.x0[0]; l < grid.x1[0]; ++l) {
                f(t, i, j, k, l);
            } } } }
	} 
};

template <>
struct meta_grid_interior <3>{
    template <typename F>
	static inline void single_step(int t, Grid_Info<3> const & grid, Grid_Info<3> const & initial_grid, F const & f) {
        /* add cilk_for here will only lower the performance */
		for (int i = grid.x0[2]; i < grid.x1[2]; ++i) {
			for (int j = grid.x0[1]; j < grid.x1[1]; ++j) {
        for (int k = grid.x0[0]; k < grid.x1[0]; ++k) {
                f(t, i, j, k);
        } } }
	} 
};

template <>
struct meta_grid_interior <2>{
    template <typename F>
	static inline void single_step(int t, Grid_Info<2> const & grid, Grid_Info<2> const & initial_grid, F const & f) {
		for (int i = grid.x0[1]; i < grid.x1[1]; ++i) {
			for (int j = grid.x0[0]; j < grid.x1[0]; ++j) {
                f(t, i, j);
			} }
	} 
};

template <>
struct meta_grid_interior <1>{
    template <typename F>
	static inline void single_step(int t, Grid_Info<1> const & grid, Grid_Info<1> const & initial_grid, F const & f) {
		for (int i = grid.x0[0]; i < grid.x1[0]; ++i) {
		    f(t, i);
        }
	} 
};

static inline void set_worker_count(const char * nstr) 
{
#if 1
    if (0 != __cilkrts_set_param("nworkers", nstr)) {
        printf("Failed to set worker count\n");
    } else {
        printf("Successfully set worker count to %s\n", nstr);
    }
#endif
}

template <int BASE, int EXP>
struct power {
    enum { value = BASE * power<BASE, EXP-1>::value };
};

template <int BASE>
struct power<BASE, 1> {
    enum { value = BASE };
}; 

template <int N_RANK>
struct Algorithm {
	private:
        /* different stencils will have different slopes */
        /* We cut coarser in internal region, finer at boundary
         * to maximize the performance and reduce the region that
         * needs special treatment
        */
        int dx_recursive_[N_RANK];
        int dx_recursive_boundary_[N_RANK];
        int dt_recursive_;
        int dx_homo_[N_RANK]; /* threshold for largest inhomogeneous region */
        int dt_homo_;
        const int dt_recursive_boundary_;
        int Z;
        const int r_t; /* # of pieces cut in time dimension */
        int N_CORES;
        typedef int index_info[N_RANK];
        typedef struct {
            int level; /* level is how many dimensions we have cut so far */
            int t0, t1;
            Grid_Info<N_RANK> grid;
        } queue_info;

        /* we can use toggled circular queue! */
        Grid_Info<N_RANK> phys_grid_;
        int phys_length_[N_RANK];
        int slope_[N_RANK+1];
        int ulb_boundary[N_RANK], uub_boundary[N_RANK], lub_boundary[N_RANK];
        bool boundarySet, physGridSet, slopeSet, opksSet, ptsSet;
        int lcm_unroll_, time_shift_;
        Pochoir_Combined_Obase_Kernel<N_RANK> ** opks_;
		int num_kernels ; //number of kernels generated
        int sz_base_data_, sz_sync_data_;
        Spawn_Tree<N_RANK> * tree_;
        Color_Region<N_RANK> * color_region_;
        Vector_Info< Homogeneity > * homogeneity_vector_;
		int num_time_steps ; // # of time steps
#ifdef COUNT_PROJECTION
		multimap <unsigned long, unsigned long> map_1d ;
		typedef struct
		{
			int x [4] ;     //the 4 x co-ordinates
			int y [4] ;     //the 4 y co-ordinates
			unsigned long area ;        //area of the bigger rectangle
			char type ;     //'r' for rectangle, 'o' for octagon
			int larger_base ; //0 for bottom, 1 for top
		} projection_2d ;
		multimap <unsigned long, projection_2d *> map_2d ;
#endif
	public:
#if STAT
    /* sim_count_cut will be accessed outside Algorithm object */
    cilk::reducer_opadd<int> sim_count_cut[SUPPORT_RANK];
    cilk::reducer_opadd<int> interior_region_count, boundary_region_count;
    cilk::reducer_opadd<long long> interior_points_count, boundary_points_count;
#endif
    int num_kernel_, num_cond_kernel_, num_bkernel_, num_cond_bkernel_;
    typedef enum {TILE_NCORES, TILE_BOUNDARY, TILE_MP} algor_type;
    
    /* constructor */
    Algorithm (int const _slope[]) : dt_recursive_boundary_(1), r_t(1), lcm_unroll_(1), time_shift_(0) {
        for (int i = 0; i < N_RANK; ++i) {
            slope_[i] = _slope[i];
            dx_recursive_boundary_[i] = 1;
            dx_homo_[i] = 8;
            // dx_recursive_boundary_[i] = _slope[i];
//            dx_recursive_boundary_[i] = tune_dx_boundary;
            ulb_boundary[i] = uub_boundary[i] = lub_boundary[i] = 0;
            // dx_recursive_boundary_[i] = 10;
        }
        dt_homo_ = (int) dx_homo_[0]/2;
        Z = 10000;
        boundarySet = false;
        physGridSet = false;
        opksSet = ptsSet = false;
        slopeSet = true;
        sz_base_data_ = sz_sync_data_ = 0;
        num_kernel_ = num_cond_kernel_ = num_bkernel_ = num_cond_bkernel_ = 0;
        // pks_ = NULL; pts_ = NULL;
        opks_ = NULL; 
        // pure_region_ = NULL; 
        color_region_ = NULL; homogeneity_vector_ = NULL;
#define ALGOR_QUEUE_SIZE (power<5, N_RANK>::value)
#if STAT
//        for (int i = 0; i < SUPPORT_RANK; ++i) {
//            sim_count_cut[i] = 0;
//        }
#else
        // N_CORES = __cilkrts_get_nworkers();
#endif
//        cout << " N_CORES = " << N_CORES << endl;

    }

    ~Algorithm() {
        del_ele(color_region_);
        del_ele(homogeneity_vector_);
#ifdef COUNT_PROJECTION
		if (N_RANK == 1 && map_1d.size() > 0)
		{
			cout << "(N, T, # of projections) " 
			     << " (" << phys_length_ [0] << "," << num_time_steps << "," <<
				 map_1d.size() << ") " << endl ;
			multimap<unsigned long, unsigned long>::iterator pos = 
						map_1d.begin();
			cout << "slope " << slope_ [0] << endl ;
			cout << "smallest length " << pos->first << endl ;
			/*for (pos = map_1d.begin() ; pos != map_1d.end() ; pos++)
			{
				cout << pos->first << "," << pos->second << endl ;
			}*/
		}
		else if (N_RANK == 2 && map_2d.size() > 0)
		{
			cout << "(N1, N2, T, # of projections) " 
			     << " (" << phys_length_ [0] << "," <<  phys_length_ [1] << "," 
				<< num_time_steps << "," <<
				 map_2d.size() << ") " << endl ;
			cout << "slope [0, 1] " << slope_ [0] << "," << slope_ [1] 
				<< endl ;
			//cout << "# of projections " << map_2d.size() << endl ;

			typename multimap <unsigned long, projection_2d *>::iterator pos =  
			//multimap <unsigned long, unsigned long>::iterator pos =  
															map_2d.begin() ;
			cout << "smallest area " << pos->first << endl ; 
			for (pos = map_2d.begin() ; pos != map_2d.end() ; pos++)
			{
				//projection_2d & p = pos->second ;
				projection_2d * p = pos->second ;
				/*cout << "Area " << p->area << " type " << p->type << 
				" (x0, x1) " << "(" << p->x [0] << "," << p->x [1] << ") " <<
				" (x2, x3) " << "(" << p->x [2] << "," << p->x [3] << ") " <<
				" (y0, y1) " << "(" << p->y [0] << "," << p->y [1] << ") " <<
				" (y2, y3) " << "(" << p->y [2] << "," << p->y [3] << ") " <<
				 endl ;*/
				delete p ;
				pos->second = 0 ;
			}
		}
#endif
    }
    /* README!!!: set_phys_grid()/set_stride() must be called before call to 
     * - walk_adaptive 
     * - walk_ncores_hybrid
     * - walk_ncores_boundary
     */
    inline void set_tree(Spawn_Tree<N_RANK> * _tree) { tree_ = _tree; }
    inline int get_sz_base_data(void) { return sz_base_data_; }
    inline int get_sz_sync_data(void) { return sz_sync_data_; }
    inline void read_stat_kernel(int & _num_kernel_, int & _num_cond_kernel_, int & _num_bkernel_, int & _num_cond_bkernel_) {
        _num_kernel_ = num_kernel_; _num_cond_kernel_ = num_cond_kernel_;
        _num_bkernel_ = num_bkernel_; _num_cond_bkernel_ = num_cond_bkernel_;
    }
    inline void set_thres(int arr_type_size) {
#if 0
        /* following threshold for debugging only! */
        dt_recursive_ = 1;
        dx_recursive_[0] = 1;
        dx_homo_[0] = 8;
        for (int i = N_RANK-1; i >= 1; --i) {
            dx_recursive_[i] = 1;
            dx_homo_[i] = 8;
        }
        if (slopeSet) 
            dt_homo_ = (int)dx_homo_[0]/(2 * slope_[0]);
        else 
            dt_homo_ = (int)dx_homo_[0]/2;

#else
#if 1
        /* following threshold for performance run! */
        dx_recursive_[0] = (N_RANK == 2) ? (int)ceil(float((80 * sizeof(double))/arr_type_size)) : (int)floor(float((600 * sizeof(double))/arr_type_size));
        dx_homo_[0] = 8;
        for (int i = N_RANK-1; i >= 1; --i) {
            dx_recursive_[i] = (N_RANK == 2) ? (int)ceil(float(80 * sizeof(double))/arr_type_size): 10;
            dx_homo_[i] = 8;
        }

        assert(slope_[0] != 0);
        dt_recursive_ = (N_RANK == 1) ? floor(dx_recursive_[0]/(2 * slope_[0]) - 100) : ((N_RANK == 2) ? floor(dx_recursive_[0]/(2 * slope_[0])-10) : 5);
#ifdef COUNT_PROJECTION
		//avoid any optimizations for space/time cuts
        dx_recursive_[0] = 3 ;
        dx_homo_[0] = 1;
        for (int i = N_RANK-1; i >= 1; --i) {
            dx_recursive_[i] = 3 ;
            dx_homo_[i] = 1;
        }
        dt_recursive_ = 1 ;
#endif
#else
        dx_recursive_[0] = (N_RANK == 1) ? 10000 : (N_RANK == 2 ? 100 : (N_RANK == 3 ? 20 : 10));
        for (int i = N_RANK-1; i >= 1; --i)
            dx_recursive_[i] = (N_RANK == 1) ? 10000 : (N_RANK == 2 ? 100 : (N_RANK == 3 ? 20 : 10));
        dt_recursive_ = (N_RANK == 1) ? 4500 : (N_RANK == 2 ? 45 : (N_RANK == 3 ? 4 : 2));
#endif
#endif
#if DEBUG 
        printf("arr_type_size = %d\n", arr_type_size);
        printf("dt_thres = %d, ", dt_recursive_);
        for (int i = N_RANK-1; i >=1; --i)
            printf("dx_thres[%d] = %d, ", i, dx_recursive_[i]);
        printf("dx_thres[%d] = %d\n", 0, dx_recursive_[0]);
#endif
    }
    inline void push_queue(int dep, int level, int t0, int t1, Grid_Info<N_RANK> const & grid);
    inline queue_info & top_queue(int dep);
    inline void pop_queue(int dep);

    void set_phys_grid(Grid_Info<N_RANK> const & grid);
    // void set_stride(int const stride[]);
    void set_slope(int const slope[]);
    void set_pts(Vector_Info<Pochoir_Guard<N_RANK> *> & _pgs);
    void set_opks(Pochoir_Combined_Obase_Kernel<N_RANK> ** _opks);
    void set_unroll(int _lcm_unroll) { 
        lcm_unroll_ = _lcm_unroll; 
        // dt_recursive_boundary_ = _lcm_unroll; 
    }
    void set_time_shift(int _time_shift) { time_shift_ = _time_shift; }
	void set_time_step(int num_time_steps) 
	{
		this->num_time_steps = num_time_steps ;
	}

#ifdef COUNT_PROJECTION
	//print projection
	inline void print_projection(int t0, int t1, 
								Grid_Info<N_RANK> const & grid);
#endif

    Vector_Info< Homogeneity > & get_color_vector(void) { return (*homogeneity_vector_); }
    inline bool touch_boundary(int i, int lt, Grid_Info<N_RANK> & grid);
    inline bool within_boundary(int t0, int t1, Grid_Info<N_RANK> & grid);
    
    /*******************************************************************************/
    /* adaptive version for irregular stencil computation */
    inline void adaptive_space_bicut(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void adaptive_space_bicut_p(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void adaptive_bicut(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void adaptive_bicut_p(int t0, int t1, Grid_Info<N_RANK> const grid);
    /* adaptive version for irregular stencil computation */
    /*******************************************************************************/

    /*******************************************************************************/
    /* meta functions to generate the execute plan */
    inline void gen_plan_space_bicut_p(Node_Info<N_RANK> * parent, int t0, int t1, Grid_Info<N_RANK> const grid, int rec_level);
    inline void gen_plan_space_cut_p(Node_Info<N_RANK> * parent, int t0, int t1, Grid_Info<N_RANK> const grid, int rec_level);
    inline void gen_plan_bicut_p(Node_Info<N_RANK> * parent, int t0, int t1, Grid_Info<N_RANK> const grid, int rec_level);
    inline void gen_plan_cut_p(Node_Info<N_RANK> * parent, int t0, int t1, Grid_Info<N_RANK> const grid, int rec_level);
    /* meta functions to generate the execute plan */
    /*******************************************************************************/
    /* meta functions to run the plan 
     * -- This is the version for unroll-t-tile-* modes
     */
    inline void plan_space_bicut(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void space_bicut_interior(int t0, int t1, 
									Grid_Info<N_RANK> const grid);
    inline void plan_space_bicut_p(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void space_bicut_boundary(int t0, int t1, 
									Grid_Info<N_RANK> const grid);
    inline void plan_space_cut(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void space_cut_interior(int t0, int t1, 
									Grid_Info<N_RANK> const grid);
    inline void plan_space_cut_p(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void space_cut_boundary(int t0, int t1, 
									Grid_Info<N_RANK> const grid);
    inline void plan_bicut(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void bicut_interior(int t0, int t1, Grid_Info<N_RANK> const grid);
    inline void plan_bicut_p(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void bicut_boundary(int t0, int t1, Grid_Info<N_RANK> const grid);
    inline void plan_cut(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void space_time_cut_interior(int t0, int t1, 
										Grid_Info<N_RANK> const grid);
    inline void plan_cut_p(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void space_time_cut_boundary(int t0, int t1, 
										Grid_Info<N_RANK> const grid);

    /* meta functions to run the plan */
    /*******************************************************************************/
    /* meta functions to run the plan 
     * -- 'm' is the version for all-cond-tile-* modes
     */
    template <typename F>
    inline void plan_space_cut_m(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n, F const & f);
    template <typename F, typename BF>
    inline void plan_space_cut_mp(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n, F const & f, BF const & bf);
    template <typename F>
    inline void plan_bicut_m(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n, F const & f);
    template <typename F, typename BF>
    inline void plan_bicut_mp(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n, F const & f, BF const & bf);
    /* meta functions to run the plan */
    /*******************************************************************************/
    /* followings are the sim cut of both top and bottom bar */
    template <typename F>
    inline void shorter_duo_sim_obase_space_cut(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    template <typename F>
    inline void shorter_duo_sim_obase_bicut(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);

    template <typename F, typename BF>
    inline void shorter_duo_sim_obase_space_cut_p(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f, BF const & bf);
    template <typename F, typename BF>
    inline void shorter_duo_sim_obase_bicut_p(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f, BF const & bf);

    /* followings are the sim cut of both top and bottom bar */
    template <typename F>
    inline void duo_sim_obase_space_cut(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    template <typename F>
    inline void duo_sim_obase_bicut(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);

    template <typename F, typename BF>
    inline void duo_sim_obase_space_cut_p(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f, BF const & bf);
    template <typename F, typename BF>
    inline void duo_sim_obase_bicut_p(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f, BF const & bf);

    /* followings are sim cut only on bottom bar */
    template <typename F>
    inline void sim_obase_space_cut(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    template <typename F>
    inline void sim_obase_bicut(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);

    template <typename F, typename BF>
    inline void sim_obase_space_cut_p(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f, BF const & bf);
    template <typename F, typename BF>
    inline void sim_obase_bicut_p(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f, BF const & bf);

    template <typename F> 
	inline void walk_serial(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);

    /* all recursion-based algorithm */
    template <typename F> 
    inline void walk_adaptive(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    template <typename F> 
    inline void walk_bicut(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    /* recursive algorithm for obase */
    template <typename F> 
    inline void obase_m(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    template <typename F> 
    inline void obase_adaptive(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    template <typename F> 
    inline void obase_bicut(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    template <typename F, typename BF> 
    inline void walk_ncores_boundary_p(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f, BF const & bf);
    template <typename F, typename BF> 
    inline void walk_bicut_boundary_p(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f, BF const & bf);
    template <typename BF> 
    inline void obase_boundary_p(int t0, int t1, Grid_Info<N_RANK> const grid, BF const & bf);
    template <typename BF> 
    inline void obase_bicut_boundary_p(int t0, int t1, Grid_Info<N_RANK> const grid, BF const & bf);
    template <typename F, typename BF> 
    inline void obase_boundary_p(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f, BF const & bf);
    template <typename F, typename BF> 
    inline void obase_bicut_boundary_p(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f, BF const & bf);

    /* all loop-based algorithm */
    template <typename F> 
    inline void cut_time(algor_type algor, int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    template <typename F> 
    inline void naive_cut_space_mp(int dim, int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    template <typename F> 
    inline void naive_cut_space_ncores(int dim, int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    template <typename F> 
    inline void cut_space_ncores_boundary(int dim, int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    inline void single_step(int t, int i);
    inline void single_step(int t, int i, int j);
#if DEBUG_FACILITY 
	void print_grid(FILE * fp, int t0, int t1, Grid_Info<N_RANK> const & grid);
	void print_sync(FILE * fp);
	void print_index(int t, int const idx[]);
	void print_region(int t, int const head[], int const tail[]);
#endif

	//Find the homogeneity of the grid
    inline void find_homogeneity(Grid_Info<N_RANK> const & grid);
	
	void set_num_kernels(int num_kernels) 
	{
		assert (num_kernels > 0) ;
		this->num_kernels = num_kernels ;
	}
};

template <int N_RANK>
void Algorithm<N_RANK>::set_phys_grid(Grid_Info<N_RANK> const & grid)
{
    phys_grid_ = grid;
    for (int i = 0; i < N_RANK; ++i)
        phys_length_[i] = grid.x1[i] - grid.x0[i];
    physGridSet = true;
    if (slopeSet) {
        /* set up the lb/ub_boundary */
        for (int i = 0; i < N_RANK; ++i) {
            ulb_boundary[i] = phys_grid_.x1[i] - slope_[i];
            uub_boundary[i] = phys_grid_.x1[i] + slope_[i];
            lub_boundary[i] = phys_grid_.x0[i] + slope_[i];
        }
    }
}

template <int N_RANK>
void Algorithm<N_RANK>::set_slope(int const slope[])
{
    for (int i = 0; i < N_RANK; ++i)
        slope_[i] = slope[i];
    slopeSet = true;
    if (physGridSet) {
        /* set up the lb/ub_boundary */
        for (int i = 0; i < N_RANK; ++i) {
            ulb_boundary[i] = phys_grid_.x1[i] - slope_[i];
            uub_boundary[i] = phys_grid_.x1[i] + slope_[i];
            lub_boundary[i] = phys_grid_.x0[i] + slope_[i];
        }
    }
}

template <int N_RANK> 
void Algorithm<N_RANK>::set_pts(Vector_Info< Pochoir_Guard<N_RANK> *> & _pgs) {
    /* for pgs_, we use the container Vector_Info since it's called in Gen_Plan,
     * the time of which is not counted into the total running time;
     * for pts_, we keep using raw pointer for performance
     */
    assert(physGridSet);
    if (color_region_ == NULL) {
        color_region_ = new Color_Region<N_RANK>(_pgs, phys_grid_);
        homogeneity_vector_ = new Vector_Info< Homogeneity >();
    }
    ptsSet = true;
    return;
}

template <int N_RANK>
void Algorithm<N_RANK>::set_opks(Pochoir_Combined_Obase_Kernel<N_RANK> ** _opks) {
    /* for pgs_, we use the container Vector_Info since it's called in Gen_Plan,
     * the time of which is not counted into the total running time;
     * for opks_, we keep using raw pointer for performance
     */
    assert(physGridSet);
    opks_ = _opks;
    opksSet = true;
    return;
}

#if DEBUG_FACILITY 
template <int N_RANK>
void Algorithm<N_RANK>::print_grid(FILE *fp, int t0, int t1, Grid_Info<N_RANK> const & grid)
{
    int i;
    fprintf(fp, "{ BASE, ");
    fprintf(fp, "t = {%d, %d}, {", t0, t1);

    fprintf(fp, "x0 = {");
    for (i = 0; i < N_RANK-1; ++i) {
        /* print x0[3] */
        fprintf(fp, "%lu, ", grid.x0[i]);
    }
    fprintf(fp, "%lu}, ", grid.x0[i]);

    fprintf(fp, "x1 = {");
    for (i = 0; i < N_RANK-1; ++i) {
        /* print x1[3] */
        fprintf(fp, "%lu, ", grid.x1[i]);
    }
    fprintf(fp, "%lu}, ", grid.x1[i]);

    fprintf(fp, "dx0 = {");
    for (i = 0; i < N_RANK-1; ++i) {
        /* print dx0[3] */
        fprintf(fp, "%d, ", grid.dx0[i]);
    }
    fprintf(fp, "%d}, ", grid.dx0[i]);

    fprintf(fp, "dx1 = {");
    for (i = 0; i < N_RANK-1; ++i) {
        /* print dx1[3] */
        fprintf(fp, "%d, ", grid.dx1[i]);
    }
    fprintf(fp, "%d}}}, \n", grid.dx1[i]);
    fflush(fp);
    return;
}

template <int N_RANK>
void Algorithm<N_RANK>::print_sync(FILE * fp)
{
    int i;
    fprintf(fp, "{ SYNC, ");
    fprintf(fp, "t = {0, 0}, {");

    fprintf(fp, "x0 = {");
    for (i = 0; i < N_RANK-1; ++i) {
        /* print x0[3] */
        fprintf(fp, "0, ");
    }
    fprintf(fp, "0}, ");

    fprintf(fp, "x1 = {");
    for (i = 0; i < N_RANK-1; ++i) {
        /* print x1[3] */
        fprintf(fp, "0, ");
    }
    fprintf(fp, "0}, ");

    fprintf(fp, "dx0 = {");
    for (i = 0; i < N_RANK-1; ++i) {
        /* print dx0[3] */
        fprintf(fp, "0, ");
    }
    fprintf(fp, "0}, ");

    fprintf(fp, "dx1 = {");
    for (i = 0; i < N_RANK-1; ++i) {
        /* print dx1[3] */
        fprintf(fp, "0, ");
    }
    fprintf(fp, "0}}}, \n");
    fflush(fp);
    return;
}

template <int N_RANK>
void Algorithm<N_RANK>::print_index(int t, int const idx[])
{
    printf("U(t=%lu, {", t);
    for (int i = 0; i < N_RANK; ++i) {
        printf("%lu ", idx[i]);
    }
    printf("}) ");
    fflush(stdout);
}

template <int N_RANK>
void Algorithm<N_RANK>::print_region(int t, int const head[], int const tail[])
{
    printf("%s:%lu t=%lu, {", __FUNCTION__, __LINE__, t);
    for (int i = 0; i < N_RANK; ++i) {
        printf("{%lu, %lu} ", head[i], tail[i]);
    }
    printf("}\n");
    fflush(stdout);

}

#endif /* end if DEBUG */
#endif /* POCHOIR_WALK_H */
