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
                int new_l = pmod_lu(l, initial_gird.x0[4], initial_grid.x1[4]);
        for (int m = grid.x0[3]; m < grid.x1[3]; ++m) {
            int new_m = pmod_lu(m, initial_grid.x0[3], initial_grid.x1[3]);
            for (int n = grid.x0[2]; n < grid.x1[2]; ++n) {
                int new_n = pmod_lu(n, initial_grid.x0[2], initial_grid.x1[2]);
        for (int o = grid.x0[1]; o < grid.x1[1]; ++o) {
            int new_o = pmod_lu(o, initial_grid.x0[1], initial_grid.x1[1]);
            for (int p = grid.x0[0]; p < grid.x1[0]; ++p) {
                int new_p = pmod_lu(p, initial_grid.x0[0], initial_grid.x1[0]);
                if (inRun) {
                    home_cell_[8] = new_p; home_cell_[7] = new_o;
                    home_cell_[6] = new_n; home_cell_[5] = new_m;
                    home_cell_[4] = new_l; home_cell_[3] = new_k;
                    home_cell_[2] = new_j; home_cell_[1] = new_i;
                }
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
                int new_l = pmod_lu(l, initial_gird.x0[3], initial_grid.x1[3]);
        for (int m = grid.x0[2]; m < grid.x1[2]; ++m) {
            int new_m = pmod_lu(m, initial_grid.x0[2], initial_grid.x1[2]);
            for (int n = grid.x0[1]; n < grid.x1[1]; ++n) {
                int new_n = pmod_lu(n, initial_grid.x0[1], initial_grid.x1[1]);
        for (int o = grid.x0[0]; o < grid.x1[0]; ++o) {
            int new_o = pmod_lu(o, initial_grid.x0[0], initial_grid.x1[0]);
                if (inRun) {
                    home_cell_[7] = new_o;
                    home_cell_[6] = new_n; home_cell_[5] = new_m;
                    home_cell_[4] = new_l; home_cell_[3] = new_k;
                    home_cell_[2] = new_j; home_cell_[1] = new_i;
                }
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
                int new_l = pmod_lu(l, initial_gird.x0[2], initial_grid.x1[2]);
        for (int m = grid.x0[1]; m < grid.x1[1]; ++m) {
            int new_m = pmod_lu(m, initial_grid.x0[1], initial_grid.x1[1]);
            for (int n = grid.x0[0]; n < grid.x1[0]; ++n) {
                int new_n = pmod_lu(n, initial_grid.x0[0], initial_grid.x1[0]);
                if (inRun) {
                    home_cell_[6] = new_n; home_cell_[5] = new_m;
                    home_cell_[4] = new_l; home_cell_[3] = new_k;
                    home_cell_[2] = new_j; home_cell_[1] = new_i;
                }
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
                int new_l = pmod_lu(l, initial_gird.x0[1], initial_grid.x1[1]);
        for (int m = grid.x0[0]; m < grid.x1[0]; ++m) {
            int new_m = pmod_lu(m, initial_grid.x0[0], initial_grid.x1[0]);
                if (inRun) {
                    home_cell_[5] = new_m;
                    home_cell_[4] = new_l; home_cell_[3] = new_k;
                    home_cell_[2] = new_j; home_cell_[1] = new_i;
                }
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
                int new_l = pmod_lu(l, initial_gird.x0[0], initial_grid.x1[0]);
                if (inRun) {
                    home_cell_[4] = new_l; home_cell_[3] = new_k;
                    home_cell_[2] = new_j; home_cell_[1] = new_i;
                }
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
                do {
                    home_cell_[3] = inRun ? new_k : 0;
                    home_cell_[2] = inRun ? new_j : 0;
                    home_cell_[1] = inRun ? new_i : 0;
                } while (0);
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
                do {
                    home_cell_[2] = inRun ? new_j : 0;
                    home_cell_[1] = inRun ? new_i : 0;
                } while (0);
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
            do {
                home_cell_[1] = inRun ? new_i : 0;
            } while (0);
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

template <int N_RANK>
struct power {
    enum { value = 5 * power<N_RANK-1>::value };
};

template <>
struct power<1> {
    enum {value = 5};
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
        const int dt_recursive_boundary_;
        int Z;
        const int r_t; /* # of pieces cut in time dimension */
        const int pad_point_;
        int N_CORES;
        typedef int index_info[N_RANK];
        typedef struct {
            int level; /* level is how many dimensions we have cut so far */
            int t0, t1;
            Grid_Info<N_RANK> grid;
        } queue_info;

        int ALGOR_QUEUE_SIZE;

        /* we can use toggled circular queue! */
        Grid_Info<N_RANK> phys_grid_;
        int phys_length_[N_RANK];
        int slope_[N_RANK];
        int ulb_boundary[N_RANK], uub_boundary[N_RANK], lub_boundary[N_RANK];
        bool boundarySet, physGridSet, slopeSet, pgkSet, opgkSet;
        int sz_pgk_;
        int lcm_unroll_, time_shift_;
        Pochoir_Guard_Kernel<N_RANK> * pgk_;
        Pochoir_Obase_Guard_Kernel<N_RANK> * opgk_;
#if PURE_REGION_ALL
        Pure_Region_All<N_RANK> * pure_region_;
#else
        Pure_Region_Corners<N_RANK> * pure_region_;
#endif
        int sz_base_data_, sz_sync_data_;
        Spawn_Tree<N_RANK> * tree_;
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
    Algorithm (int const _slope[]) : dt_recursive_boundary_(1), r_t(1), lcm_unroll_(1), time_shift_(0), pad_point_(2) {
        for (int i = 0; i < N_RANK; ++i) {
            slope_[i] = _slope[i];
            dx_recursive_boundary_[i] = 1;
            // dx_recursive_boundary_[i] = _slope[i];
//            dx_recursive_boundary_[i] = tune_dx_boundary;
            ulb_boundary[i] = uub_boundary[i] = lub_boundary[i] = 0;
            // dx_recursive_boundary_[i] = 10;
        }
        Z = 10000;
        boundarySet = false;
        physGridSet = false;
        pgkSet = opgkSet = false;
        slopeSet = true;
        sz_base_data_ = sz_sync_data_ = 0;
        num_kernel_ = num_cond_kernel_ = num_bkernel_ = num_cond_bkernel_ = 0;
        /* ALGOR_QUEUE_SIZE = 3^N_RANK */
        // ALGOR_QUEUE_SIZE = power<N_RANK>::value;
#define ALGOR_QUEUE_SIZE (power<N_RANK>::value)
#if STAT
//        for (int i = 0; i < SUPPORT_RANK; ++i) {
//            sim_count_cut[i] = 0;
//        }
#else
        N_CORES = __cilkrts_get_nworkers();
#endif
//        cout << " N_CORES = " << N_CORES << endl;

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
        dt_recursive_ = 1;
        dx_recursive_[0] = 1;
        for (int i = N_RANK-1; i >= 1; --i)
            dx_recursive_[i] = 1;
#else
#if 1
        dx_recursive_[0] = (N_RANK == 2) ? (int)ceil(float((80 * sizeof(double))/arr_type_size)) : (int)floor(float((600 * sizeof(double))/arr_type_size));
//        dx_recursive_[0] = 30;
        for (int i = N_RANK-1; i >= 1; --i)
            dx_recursive_[i] = (N_RANK == 2) ? (int)ceil(float(80 * sizeof(double))/arr_type_size): 10;
        assert(slope_[0] != 0);
        dt_recursive_ = (N_RANK == 1) ? floor(dx_recursive_[0]/(2 * slope_[0]) - 100) : ((N_RANK == 2) ? floor(dx_recursive_[0]/(2 * slope_[0])-10) : 5);
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
    void set_pgk(int _sz_pgk, Pochoir_Guard_Kernel<N_RANK> * _pgk);
    void set_obase_pgk(int _sz_pgk, Pochoir_Obase_Guard_Kernel<N_RANK> * _opgk);
    void set_unroll(int _lcm_unroll) { lcm_unroll_ = _lcm_unroll; }
    void set_time_shift(int _time_shift) { time_shift_ = _time_shift; }

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
    inline void gen_plan_space_bicut_p(Node_Info<N_RANK> * parent, int t0, int t1, Grid_Info<N_RANK> const grid);
    inline void gen_plan_bicut_p(Node_Info<N_RANK> * parent, int t0, int t1, Grid_Info<N_RANK> const grid);
    /* meta functions to generate the execute plan */
    /*******************************************************************************/
    /* meta functions to run the plan */
    inline void plan_space_bicut(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void plan_space_bicut_p(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void plan_bicut(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
    inline void plan_bicut_p(int t0, int t1, Grid_Info<N_RANK> const grid, int region_n);
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
	inline void base_case_kernel_interior(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f);
    template <typename BF> 
	inline void base_case_kernel_boundary(int t0, int t1, Grid_Info<N_RANK> const grid, BF const & bf);
    template <typename G1, typename F1, typename G2, typename F2> 
	inline void base_case_kernel_stagger(int t0, int t1, Grid_Info<N_RANK> const grid, G1 const & g1, F1 const & f1, G2 const & g2, F2 const & f2);
    inline void base_case_kernel_guard(int t0, int t1, Grid_Info<N_RANK> const grid);
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
    if (pgkSet || opgkSet) {
        pure_region_->set_phys_grid(phys_grid_);
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
void Algorithm<N_RANK>::set_pgk(int _sz_pgk, Pochoir_Guard_Kernel<N_RANK> * _pgk) {
    sz_pgk_ = _sz_pgk; pgk_ = _pgk; 
    pgkSet = true;
#if PURE_REGION_ALL
    pure_region_ = new Pure_Region_All<N_RANK>(_sz_pgk, _pgk);
#else
    pure_region_ = new Pure_Region_Corners<N_RANK>(_sz_pgk, _pgk);
#endif
    if (physGridSet) {
        pure_region_->set_phys_grid(phys_grid_);
    }
    return;
}

template <int N_RANK> 
void Algorithm<N_RANK>::set_obase_pgk(int _sz_pgk, Pochoir_Obase_Guard_Kernel<N_RANK> * _opgk) {
    sz_pgk_ = _sz_pgk; opgk_ = _opgk; 
    opgkSet = true;
#if PURE_REGION_ALL
    pure_region_ = new Pure_Region_All<N_RANK>(_sz_pgk, _opgk);
#else
    pure_region_ = new Pure_Region_Corners<N_RANK>(_sz_pgk, _opgk);
#endif
    if (physGridSet) {
        pure_region_->set_phys_grid(phys_grid_);
    }
    return;
}

template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::base_case_kernel_interior(int t0, int t1, Grid_Info<N_RANK> const grid, F const & f) {
	Grid_Info<N_RANK> l_grid = grid;
	for (int t = t0; t < t1; ++t) {
        /* execute one single time step */
        meta_grid_interior<N_RANK>::single_step(t, l_grid, phys_grid_, f);

        /* because the shape is trapezoid! */
        for (int i = 0; i < N_RANK; ++i) {
            l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
        }
	}
}

template <int N_RANK> template <typename BF>
inline void Algorithm<N_RANK>::base_case_kernel_boundary(int t0, int t1, Grid_Info<N_RANK> const grid, BF const & bf) {
	Grid_Info<N_RANK> l_grid = grid;
	for (int t = t0; t < t1; ++t) {
        home_cell_[0] = t;
        /* execute one single time step */
        meta_grid_boundary<N_RANK>::single_step(t, l_grid, phys_grid_, bf);

        /* because the shape is trapezoid! */
        for (int i = 0; i < N_RANK; ++i) {
            l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
        }
	}
}

template <int N_RANK> 
inline void Algorithm<N_RANK>::base_case_kernel_guard(int t0, int t1, Grid_Info<N_RANK> const grid) {
    /* Each kernel update the entire region independently and entirely!!! */
	Grid_Info<N_RANK> l_grid = grid;
    Pochoir_Generic_Kernel<N_RANK> l_kernel(sz_pgk_, pgk_);
    for (int t = t0; t < t1; ) {
        home_cell_[0] = t;
        meta_grid_boundary<N_RANK>::single_step(t, l_grid, phys_grid_, l_kernel);
        /* adjust the trapezoids */
        for (int i = 0; i < N_RANK; ++i) {
            l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
        } 
        ++t;
        l_kernel.shift_pointer();
    } /* end for 't' */
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
