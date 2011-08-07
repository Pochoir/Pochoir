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

#ifndef POCHOIR_WRAPPER_H
#define POCHOIR_WRAPPER_H

#if DEBUG
#include <cstdio>
#endif

#include "pochoir_common.hpp"
#include "pochoir_range.hpp"
#include "pochoir_walk.hpp"
#include "pochoir_walk_recursive.hpp"
#include "pochoir_walk_loops.hpp"

/* serial_loops() is not necessary because we can call base_case_kernel() to 
 * mimic the same behavior of serial_loops()
 */
template <typename F>
void serial_loops(Pochoir_Domain _tR, Pochoir_Domain _iR, Pochoir_Domain _jR, F const & f) { 
    size_t t_first = _tR.first(), t_last = _tR.last(), t_stride = _tR.stride();
    size_t i_first = _iR.first(), i_last = _iR.last(), i_stride = _iR.stride();
    size_t j_first = _jR.first(), j_last = _jR.last(), j_stride = _jR.stride(); 
    for (size_t t = t_first; t <= t_last ; t += t_stride) {
    for (size_t i = i_first; i <= i_last ; i += i_stride) {
    for (size_t j = j_first; j <= j_last ; j += j_stride) {
        f(t, i, j);
	} } }
} 

template <typename F>
void serial_loops(Pochoir_Domain _tR, Pochoir_Domain _iR, F const & f) { 
    size_t t_first = _tR.first(), t_last = _tR.last(), t_stride = _tR.stride();
    size_t i_first = _iR.first(), i_last = _iR.last(), i_stride = _iR.stride();
    for (size_t t = t_first; t <= t_last ; t += t_stride) {
    for (size_t i = i_first; i <= i_last ; i += i_stride) {
        f(t, i);
	} } 
} 

/* these are for those fall in full effective region and dont have any boundary conditions */
template <typename F>
void pochoir(Pochoir_Domain const & _tR, Pochoir_Domain const & _iR, Pochoir_Domain const & _jR, const size_t _slope[], F const f) {
	size_t l_t0 = _tR.first(), l_t1 = _tR.last();
    Grid_Info_2 l_grid;
    Algorithm<3, Grid_Info_2> algor(_slope);
    size_t l_stride[2];

    l_grid.x0[0] = _iR.first(); l_grid.x1[0] = _iR.first() + _iR.stride() * _iR.size();
    l_grid.x0[1] = _jR.first(); l_grid.x1[1] = _jR.first() + _jR.stride() * _jR.size();
    l_grid.dx0[0] = 0; l_grid.dx1[0] = 0;
    l_grid.dx0[1] = 0; l_grid.dx1[1] = 0;
    l_stride[0] = _iR.stride(); l_stride[1] = _jR.stride();
	//algor.walk_serial(l_t0, l_t1, l_grid, f);
	//algor.base_case_kernel(l_t0, l_t1, l_grid, f);
    algor.set_initial_grid(l_grid);
    algor.set_stride(l_stride);
    algor.walk_ncores_adaptive(l_t0, l_t1, l_grid, f);
}

template <typename F>
void pochoir(Pochoir_Domain const & _tR, Pochoir_Domain const & _iR, const size_t _slope[], F const f) {
	size_t l_t0 = _tR.first(), l_t1 = _tR.last();
    Algorithm<2, Grid_Info_1> algor(_slope);
    Grid_Info_1 l_grid;
    size_t l_stride[1];

    l_grid.x0[0] = _iR.first(); l_grid.x1[0] = _iR.first() + _iR.stride() * _iR.size();
    l_grid.dx0[0] = 0; l_grid.dx1[0] = 0;
    l_stride[0] = _iR.stride();
	//algor.walk_serial(l_t0, l_t1, l_grid, f);
	//algor.base_case_kernel(l_t0, l_t1, l_grid, f);
    algor.set_initial_grid(l_grid);
    algor.set_stride(l_stride);
    algor.walk_ncores_adaptive(l_t0, l_t1, l_grid, f);
}

/* Non-periodic: F is for internal region, and BF is for boundary condition processing */
template <typename F, typename BF>
void pochoir(Pochoir_Domain const & _tR, Pochoir_Domain const & _iR, Pochoir_Domain const & _jR, const size_t _slope[], F const & f, BF const & bf) {
	size_t l_t0 = _tR.first(), l_t1 = _tR.last();
    Grid_Info_2 l_grid;
    Algorithm<3, Grid_Info_2> algor(_slope);
    size_t l_stride[2];

    l_grid.x0[0] = _iR.first(); l_grid.x1[0] = _iR.first() + _iR.stride() * _iR.size();
    l_grid.x0[1] = _jR.first(); l_grid.x1[1] = _jR.first() + _jR.stride() * _jR.size();
    l_grid.dx0[0] = 0; l_grid.dx1[0] = 0;
    l_grid.dx0[1] = 0; l_grid.dx1[1] = 0;
    l_stride[0] = _iR.stride(); l_stride[1] = _jR.stride();
	//algor.walk_serial(l_t0, l_t1, l_grid, f);
	//algor.base_case_kernel(l_t0, l_t1, l_grid, f);
    algor.set_initial_grid(l_grid);
    algor.set_stride(l_stride);
    algor.walk_ncores_boundary(l_t0, l_t1, l_grid, f, bf);
}

template <typename F, typename BF>
void pochoir(Pochoir_Domain const & _tR, Pochoir_Domain const & _iR, const size_t _slope[], F const & f, BF const & bf) {
	size_t l_t0 = _tR.first(), l_t1 = _tR.last();
    Algorithm<2, Grid_Info_1> algor(_slope);
    Grid_Info_1 l_grid;
    size_t l_stride[1];

    l_grid.x0[0] = _iR.first(); l_grid.x1[0] = _iR.first() + _iR.stride() * _iR.size();
    l_grid.dx0[0] = 0; l_grid.dx1[0] = 0;
    l_stride[0] = _iR.stride();
	//algor.walk_serial(l_t0, l_t1, l_grid, f);
	//algor.base_case_kernel(l_t0, l_t1, l_grid, f);
    algor.set_initial_grid(l_grid);
    algor.set_stride(l_stride);
    algor.walk_ncores_boundary(l_t0, l_t1, l_grid, f, bf);
}

/* pochoir_p() will only call the Periodic stencil function,
 * while pochoir() may NOT -- Shall we merge the function with those of non-periodic,
 * and just use one parameter to distinguish them???
 */
template <typename F, typename BF>
void pochoir_p(Pochoir_Domain const & _tR, Pochoir_Domain const & _iR, Pochoir_Domain const & _jR, const size_t _slope[], F const & f, BF const & bf) {
	size_t l_t0 = _tR.first(), l_t1 = _tR.last();
    Grid_Info_2 l_grid;
    Algorithm<3, Grid_Info_2> algor(_slope);
    size_t l_stride[2];

    l_grid.x0[0] = _iR.first(); l_grid.x1[0] = _iR.first() + _iR.stride() * _iR.size();
    l_grid.x0[1] = _jR.first(); l_grid.x1[1] = _jR.first() + _jR.stride() * _jR.size();
    l_grid.dx0[0] = 0; l_grid.dx1[0] = 0;
    l_grid.dx0[1] = 0; l_grid.dx1[1] = 0;
    l_stride[0] = _iR.stride(); l_stride[1] = _jR.stride();
	//algor.walk_serial(l_t0, l_t1, l_grid, f);
	//algor.base_case_kernel(l_t0, l_t1, l_grid, f);
    algor.set_initial_grid(l_grid);
    algor.set_stride(l_stride);
    algor.walk_ncores_boundary_p(l_t0, l_t1, l_grid, f, bf);
}

template <typename F, typename BF>
void pochoir_p(Pochoir_Domain const & _tR, Pochoir_Domain const & _iR, const size_t _slope[], F const & f, BF const & bf) {
	size_t l_t0 = _tR.first(), l_t1 = _tR.last();
    Algorithm<2, Grid_Info_1> algor(_slope);
    Grid_Info_1 l_grid;
    Algorithm<2, Grid_Info_1>::index_info l_stride;

    l_grid.x0[0] = _iR.first(); l_grid.x1[0] = _iR.first() + _iR.stride() * _iR.size();
    l_grid.dx0[0] = 0; l_grid.dx1[0] = 0;
    l_stride[0] = _iR.stride();
	//algor.walk_serial(l_t0, l_t1, l_grid, f);
	//algor.base_case_kernel(l_t0, l_t1, l_grid, f);
    algor.set_initial_grid(l_grid);
    algor.set_stride(l_stride);
    algor.walk_ncores_boundary_p(l_t0, l_t1, l_grid, f, bf);
}

/* these are for those fall in full effective region and dont have any boundary conditions */
template <typename F>
void obase(Pochoir_Domain const & _tR, Pochoir_Domain const & _iR, Pochoir_Domain const & _jR, const size_t _slope[], F const f) {
	size_t l_t0 = _tR.first(), l_t1 = _tR.last();
    Grid_Info_2 l_grid;
    Algorithm<3, Grid_Info_2> algor(_slope);
    Algorithm<3, Grid_Info_2>::index_info l_stride;

    l_grid.x0[0] = _iR.first(); l_grid.x1[0] = _iR.first() + _iR.stride() * _iR.size();
    l_grid.x0[1] = _jR.first(); l_grid.x1[1] = _jR.first() + _jR.stride() * _jR.size();
    l_grid.dx0[0] = 0; l_grid.dx1[0] = 0;
    l_grid.dx0[1] = 0; l_grid.dx1[1] = 0;
    l_stride[0] = _iR.stride(); l_stride[1] = _jR.stride();
	//algor.walk_serial(l_t0, l_t1, l_grid, f);
	//algor.base_case_kernel(l_t0, l_t1, l_grid, f);
    algor.set_initial_grid(l_grid);
    algor.set_stride(l_stride);
    algor.obase_adaptive(l_t0, l_t1, l_grid, f);
}

template <typename F>
void obase(Pochoir_Domain const & _tR, Pochoir_Domain const & _iR, const size_t _slope[], F const f) {
	size_t l_t0 = _tR.first(), l_t1 = _tR.last();
    Algorithm<2, Grid_Info_1> algor(_slope);
    Grid_Info_1 l_grid;
    Algorithm<2, Grid_Info_1>::index_info l_stride;

    l_grid.x0[0] = _iR.first(); l_grid.x1[0] = _iR.first() + _iR.stride() * _iR.size();
    l_grid.dx0[0] = 0; l_grid.dx1[0] = 0;
    l_stride[0] = _iR.stride();
	//algor.walk_serial(l_t0, l_t1, l_grid, f);
	//algor.base_case_kernel(l_t0, l_t1, l_grid, f);
    algor.set_initial_grid(l_grid);
    algor.set_stride(l_stride);
    algor.obase_adaptive(l_t0, l_t1, l_grid, f);
}

/* Non-periodic: F is for internal region, and BF is for boundary condition processing */
template <typename F, typename BF>
void obase(Pochoir_Domain const & _tR, Pochoir_Domain const & _iR, Pochoir_Domain const & _jR, const size_t _slope[], F const & f, BF const & bf) {
	size_t l_t0 = _tR.first(), l_t1 = _tR.last();
    Grid_Info_2 l_grid;
    Algorithm<3, Grid_Info_2> algor(_slope);
    size_t l_stride[2];

    l_grid.x0[0] = _iR.first(); l_grid.x1[0] = _iR.first() + _iR.stride() * _iR.size();
    l_grid.x0[1] = _jR.first(); l_grid.x1[1] = _jR.first() + _jR.stride() * _jR.size();
    l_grid.dx0[0] = 0; l_grid.dx1[0] = 0;
    l_grid.dx0[1] = 0; l_grid.dx1[1] = 0;
    l_stride[0] = _iR.stride(); l_stride[1] = _jR.stride();
	//algor.walk_serial(l_t0, l_t1, l_grid, f);
	//algor.base_case_kernel(l_t0, l_t1, l_grid, f);
    algor.set_initial_grid(l_grid);
    algor.set_stride(l_stride);
    algor.obase_boundary(l_t0, l_t1, l_grid, f, bf);
}

template <typename F, typename BF>
void obase(Pochoir_Domain const & _tR, Pochoir_Domain const & _iR, const size_t _slope[], F const & f, BF const & bf) {
	size_t l_t0 = _tR.first(), l_t1 = _tR.last();
    Algorithm<2, Grid_Info_1> algor(_slope);
    Grid_Info_1 l_grid;
    size_t l_stride[1];

    l_grid.x0[0] = _iR.first(); l_grid.x1[0] = _iR.first() + _iR.stride() * _iR.size();
    l_grid.dx0[0] = 0; l_grid.dx1[0] = 0;
    l_stride[0] = _iR.stride();
	//algor.walk_serial(l_t0, l_t1, l_grid, f);
	//algor.base_case_kernel(l_t0, l_t1, l_grid, f);
    algor.set_initial_grid(l_grid);
    algor.set_stride(l_stride);
    algor.obase_boundary(l_t0, l_t1, l_grid, f, bf);
}

/* obase_p() will only call the Periodic stencil function,
 * while obase() may NOT -- Shall we merge the function with those of non-periodic,
 * and just use one parameter to distinguish them???
 */
template <typename F, typename BF>
void obase_p(Pochoir_Domain const & _tR, Pochoir_Domain const & _iR, Pochoir_Domain const & _jR, const size_t _slope[], F const & f, BF const & bf) {
	size_t l_t0 = _tR.first(), l_t1 = _tR.last();
    Grid_Info_2 l_grid;
    Algorithm<3, Grid_Info_2> algor(_slope);
    size_t l_stride[2];

    l_grid.x0[0] = _iR.first(); l_grid.x1[0] = _iR.first() + _iR.stride() * _iR.size();
    l_grid.x0[1] = _jR.first(); l_grid.x1[1] = _jR.first() + _jR.stride() * _jR.size();
    l_grid.dx0[0] = 0; l_grid.dx1[0] = 0;
    l_grid.dx0[1] = 0; l_grid.dx1[1] = 0;
    l_stride[0] = _iR.stride(); l_stride[1] = _jR.stride();
	//algor.walk_serial(l_t0, l_t1, l_grid, f);
	//algor.base_case_kernel(l_t0, l_t1, l_grid, f);
    algor.set_initial_grid(l_grid);
    algor.set_stride(l_stride);
    algor.obase_boundary_p(l_t0, l_t1, l_grid, f, bf);
}

template <typename F, typename BF>
void obase_p(Pochoir_Domain const & _tR, Pochoir_Domain const & _iR, const size_t _slope[], F const & f, BF const & bf) {
	size_t l_t0 = _tR.first(), l_t1 = _tR.last();
    Algorithm<2, Grid_Info_1> algor(_slope);
    Grid_Info_1 l_grid;
    size_t l_stride[1];

    l_grid.x0[0] = _iR.first(); l_grid.x1[0] = _iR.first() + _iR.stride() * _iR.size();
    l_grid.dx0[0] = 0; l_grid.dx1[0] = 0;
    l_stride[0] = _iR.stride();
	//algor.walk_serial(l_t0, l_t1, l_grid, f);
	//algor.base_case_kernel(l_t0, l_t1, l_grid, f);
    algor.set_initial_grid(l_grid);
    algor.set_stride(l_stride);
    algor.obase_boundary_p(l_t0, l_t1, l_grid, f, bf);
}

#endif /* POCHOIR_WRAPPER_H */
