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
 *********************************************************************************
 */

#ifndef EXPR_WALK_LOOPS_H
#define EXPR_WALK_LOOPS_H

#include "pochoir_walk.hpp"

#define MAX(a, b) ((a) >= (b) ? (a) : (b))

template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::naive_cut_space_mp(int dim, int t0, int t1, Grid_Info<N_RANK> const grid, F const & f)
{
	/* This is the version that cut into as many pieces as we can */
	/* cut into Space dimension one after another */
	int i;
	int lt = t1 - t0;
	int bl = MAX(2*slope_[dim]*lt, dx_recursive_);
	int lx = (dim < N_RANK) ? (grid.x1[dim] - grid.x0[dim]) : 0;
	bool can_cut = (dim < N_RANK) ? (lx/bl >= 2) : false;

#if DEBUG
//	printf("dim = %d :", dim);
//	print_grid(stdout, t0, t1, grid);
//	fflush(stdout);
#endif
	if (!can_cut || dim == N_RANK) {
		if (dim < N_RANK)
			naive_cut_space_mp(dim+1, t0, t1, grid, f);
		else {
			assert(dim == N_RANK);
#if DEBUG
	//		printf("%s:%d base_case_kernel\n", __FUNCTION__, __LINE__);
	//		fflush(stdout);
#endif
//			base_case_kernel(t0, t1, grid);
			base_case_kernel(t0, t1, grid, f);
		}
		return;
	} else {
		assert(can_cut);
		Grid_Info<N_RANK> l_grid = grid;
		int r = lx / bl;
		int sep = bl;
		for (i = 0; i < r - 1; i++) {
			l_grid.x0[dim] = grid.x0[dim] + i * sep;
			l_grid.dx0[dim] = slope_[dim];
			l_grid.x1[dim] = grid.x0[dim] + (i + 1) * sep;
			l_grid.dx1[dim] = -slope_[dim];
			cilk_spawn naive_cut_space_mp(dim+1, t0, t1, l_grid, f);
		}
		l_grid.x0[dim] = grid.x0[dim] + i * sep;
		l_grid.dx0[dim] = slope_[dim];
		l_grid.x1[dim] = grid.x1[dim];
		l_grid.dx1[dim] = -slope_[dim];
		naive_cut_space_mp(dim+1, t0, t1, l_grid, f);
		cilk_sync;

		if (grid.dx0[dim] != slope_[dim]) {
			l_grid.x0[dim] = grid.x0[dim];
			l_grid.dx0[dim] = grid.dx0[dim];
			l_grid.x1[dim] = grid.x0[dim];
			l_grid.dx1[dim] = slope_[dim];
			cilk_spawn naive_cut_space_mp(dim+1, t0, t1, l_grid, f);
		}
		for (i = 1; i < r; i++) {
			l_grid.x0[dim] = grid.x0[dim] + i * sep;
			l_grid.dx0[dim] = -slope_[dim];
			l_grid.x1[dim] = grid.x0[dim] + i * sep;
			l_grid.dx1[dim] = slope_[dim];
			cilk_spawn naive_cut_space_mp(dim+1, t0, t1, l_grid, f);
		}
		if (grid.dx1[dim] != -slope_[dim]) {
			l_grid.x0[dim] = grid.x1[dim];
			l_grid.dx0[dim] = -slope_[dim];
			l_grid.x1[dim] = grid.x1[dim];
			l_grid.dx1[dim] = grid.dx1[dim];
			cilk_spawn naive_cut_space_mp(dim+1, t0, t1, l_grid, f);
		}
		return;
	}
}

template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::naive_cut_space_ncores(int dim, int t0, int t1, Grid_Info<N_RANK> const grid, F const & f)
{
	/* This version cut into exactly N_CORES pieces */
	/* cut into Space dimension one after another */
	int i;
	int lt = t1 - t0;
	int lx = (dim < N_RANK) ? (grid.x1[dim] - grid.x0[dim]) : 0;
	bool can_cut = (dim < N_RANK) ? (N_CORES * 2 * slope_[dim] * lt <= lx) : false;
	//printf("TILE_NCORES\n");
#if DEBUG
//	printf("dim = %d :", dim);
//	print_grid(stdout, t0, t1, grid);
//	fflush(stdout);
#endif
	if (!can_cut || dim == N_RANK) {
		if (dim < N_RANK)
			naive_cut_space_ncores(dim+1, t0, t1, grid, f);
		else {
			assert(dim == N_RANK);
#if DEBUG
//			print_grid(stdout, t0, t1, grid);
#endif
//			base_case_kernel(t0, t1, grid);
			base_case_kernel(t0, t1, grid, f);

		}
		return;
	} else {
		assert(can_cut);
		Grid_Info<N_RANK> l_grid = grid;
		int sep = lx / N_CORES;
		for (i = 0; i < N_CORES - 1; i++) {
			l_grid.x0[dim] = grid.x0[dim] + i * sep;
			l_grid.dx0[dim] = slope_[dim];
			l_grid.x1[dim] = grid.x0[dim] + (i + 1) * sep;
			l_grid.dx1[dim] = -slope_[dim];
			cilk_spawn naive_cut_space_ncores(dim+1, t0, t1, l_grid, f);
		}
		l_grid.x0[dim] = grid.x0[dim] + i * sep;
		l_grid.dx0[dim] = slope_[dim];
		l_grid.x1[dim] = grid.x1[dim];
		l_grid.dx1[dim] = -slope_[dim];
		naive_cut_space_ncores(dim+1, t0, t1, l_grid, f);
#if DEBUG
//		fprintf(stdout, "cilk_sync\n");
//		fflush(stdout);
#endif
		cilk_sync;

		if (grid.dx0[dim] != slope_[dim]) {
			l_grid.x0[dim] = grid.x0[dim];
			l_grid.dx0[dim] = grid.dx0[dim];
			l_grid.x1[dim] = grid.x0[dim];
			l_grid.dx1[dim] = slope_[dim];
			cilk_spawn naive_cut_space_ncores(dim+1, t0, t1, l_grid, f);
		}
		for (i = 1; i < N_CORES; i++) {
			l_grid.x0[dim] = grid.x0[dim] + i * sep;
			l_grid.dx0[dim] = -slope_[dim];
			l_grid.x1[dim] = grid.x0[dim] + i * sep;
			l_grid.dx1[dim] = slope_[dim];
			cilk_spawn naive_cut_space_ncores(dim+1, t0, t1, l_grid, f);
		}
		if (grid.dx1[dim] != -slope_[dim]) {
			l_grid.x0[dim] = grid.x1[dim];
			l_grid.dx0[dim] = -slope_[dim];
			l_grid.x1[dim] = grid.x1[dim];
			l_grid.dx1[dim] = grid.dx1[dim];
			cilk_spawn naive_cut_space_ncores(dim+1, t0, t1, l_grid, f);
		}
		return;
	}
}

template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::cut_space_ncores_boundary(int dim, int t0, int t1, Grid_Info<N_RANK> const grid, F const & f)
{
	/* This version cut into exactly NCORES pieces */
	/* cut into Space dimension one after another */
	int i;
	int lt = t1 - t0;
	int lx;  
	bool can_cut, call_boundary = false;

	if (dim < N_RANK) {
		lx = grid.x1[dim] - grid.x0[dim];
	} else 
		lx = 0;
	can_cut = (dim < N_RANK) ? (N_CORES * 2 * slope_[dim] * lt <= lx) : false;
	//printf("TILE_NCORES\n");
#if DEBUG
//	printf("dim = %d :", dim);
//	print_grid(stdout, t0, t1, grid);
//	fflush(stdout);
#endif
	if (!can_cut || dim == N_RANK) {
		if (dim < N_RANK) {
			cut_space_ncores_boundary(dim+1, t0, t1, grid, f);
		} else {
			assert(dim == N_RANK);
#if DEBUG
//			print_grid(stdout, t0, t1, grid);
#endif
			call_boundary = false;
			for (int i = 0; i < N_RANK; i++) {
				call_boundary |= (grid.x0[i] == initial_grid_.x0[i] || grid.x1[i] == initial_grid_.x1[i]);
			}
			if (call_boundary) 
                //we will defer the processing of boundary condition later
				//base_case_kernel_boundary(t0, t1, grid, f);
				base_case_kernel(t0, t1, grid, f);
			else
				base_case_kernel(t0, t1, grid, f);
		}
		return;
	} else {
		assert(can_cut);
		Grid_Info<N_RANK> l_grid = grid;
		int sep = lx / N_CORES;
		int l_start = (grid.x0[dim]);
		int l_end = (grid.x1[dim]);
		for (i = 0; i < N_CORES - 1; i++) {
			l_grid.x0[dim] = l_start + i * sep;
			l_grid.dx0[dim] = slope_[dim];
			l_grid.x1[dim] = l_start + (i + 1) * sep;
			l_grid.dx1[dim] = -slope_[dim];
			cilk_spawn cut_space_ncores_boundary(dim+1, t0, t1, l_grid, f);
		}
		l_grid.x0[dim] = l_start + i * sep;
		l_grid.dx0[dim] = slope_[dim];
		l_grid.x1[dim] = l_end;
		l_grid.dx1[dim] = -slope_[dim];
		cut_space_ncores_boundary(dim+1, t0, t1, l_grid, f);
#if DEBUG
//		fprintf(stdout, "cilk_sync\n");
//		fflush(stdout);
#endif
		cilk_sync;

		if (grid.dx0[dim] != slope_[dim]) {
			l_grid.x0[dim] = grid.x0[dim];
			l_grid.dx0[dim] = grid.dx0[dim];
			l_grid.x1[dim] = grid.x0[dim];
			l_grid.dx1[dim] = slope_[dim];
			cilk_spawn cut_space_ncores_boundary(dim+1, t0, t1, l_grid, f);
		}
		for (i = 1; i < N_CORES; i++) {
			l_grid.x0[dim] = grid.x0[dim] + i * sep;
			l_grid.dx0[dim] = -slope_[dim];
			l_grid.x1[dim] = grid.x0[dim] + i * sep;
			l_grid.dx1[dim] = slope_[dim];
			cilk_spawn cut_space_ncores_boundary(dim+1, t0, t1, l_grid, f);
		}
		if (grid.dx1[dim] != -slope_[dim]) {
			l_grid.x0[dim] = grid.x1[dim];
			l_grid.dx0[dim] = -slope_[dim];
			l_grid.x1[dim] = grid.x1[dim];
			l_grid.dx1[dim] = grid.dx1[dim];
			cilk_spawn cut_space_ncores_boundary(dim+1, t0, t1, l_grid, f);
		}

		return;
	}
}

template <int N_RANK> template <typename F>
inline void Algorithm<N_RANK>::cut_time(algor_type algor, int t0, int t1, Grid_Info<N_RANK> const grid, F const & f)
{
	/* cut into Time dimension */
	int i;
	int r_t = (t1 - t0)/dt_recursive_;
	if (r_t < 2) {
		switch(algor) {
		case TILE_NCORES: 
			naive_cut_space_ncores(0, t0, t1, grid, f);
			break;
		case TILE_BOUNDARY:
			cut_space_ncores_boundary(0, t0, t1, grid, f);
			break;
		case TILE_MP:
			naive_cut_space_mp(0, t0, t1, grid, f);
			break;
		default:
			break;
		}
		return;
	} else {
		assert(r_t >= 2);
		for (i = 0; i < r_t; i++) {
			switch(algor) {
			case TILE_NCORES: 
				naive_cut_space_ncores(0, t0+i*dt_recursive_, t0+(i+1)*dt_recursive_, grid, f);
#if DEBUG
//				fprintf(stdout, "cilk_sync\n");
//				fflush(stdout);
#endif
				break;
			case TILE_BOUNDARY:
				cut_space_ncores_boundary(0, t0+i*dt_recursive_, t0+(i+1)*dt_recursive_, grid, f);
				break;
			case TILE_MP:
				naive_cut_space_mp(0, t0+i*dt_recursive_, t0+(i+1)*dt_recursive_, grid, f);
				break;
			default:
				break;
			}
		}
		if (t1 > t0+i*dt_recursive_) {
			switch(algor) {
			case TILE_NCORES:
				naive_cut_space_ncores(0, t0+i*dt_recursive_, t1, grid, f);
				break;
			case TILE_BOUNDARY:
				cut_space_ncores_boundary(0, t0+i*dt_recursive_, t1, grid, f);
				break;
			case TILE_MP:
				naive_cut_space_mp(0, t0+i*dt_recursive_, t1, grid, f);
				break;
			default:
				break;
			}
		}
		return;
	}
}

#endif /* EXPR_WALK_LOOPS_H */
