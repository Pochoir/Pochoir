/*
 **********************************************************************************
 *  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 *  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 *                           Charles E. Leiserson <cel@mit.edu>
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

#ifndef POCHOIR_KERNEL_HPP
#define POCHOIR_KERNEL_HPP

#include "pochoir_common.hpp"

using namespace std;

/* Pochoir_Kernel for Phase I */
template <int N_RANK>
struct Pochoir_Kernel {
    typedef typename Pochoir_Types<N_RANK>::T_Kernel T_Kernel;
    T_Kernel kernel_;
    /* we don't make a separate copy of Pochoir_Shape for each Pochoir_Kernel,
     * since it won't be accumulated for each single kernel,
     * so we just point to the shape, instead of copying it.
     */
    Pochoir_Shape<N_RANK> * shape_;
    int shape_size_, time_shift_, toggle_, slope_[N_RANK];
    Pochoir_Kernel(void) { shape_ = NULL; }
    template <int N_SIZE>
    Pochoir_Kernel(T_Kernel _kernel, Pochoir_Shape<N_RANK> (& _shape)[N_SIZE]) : kernel_(_kernel) {
        int l_min_time_shift=0, l_max_time_shift=0, depth=0;
        for (int r = 0; r < N_RANK; ++r) {
            slope_[r] = 0;
        }
        for (int i = 0; i < N_SIZE; ++i) {
            if (_shape[i].shift[0] < l_min_time_shift)
                l_min_time_shift = _shape[i].shift[0];
            if (_shape[i].shift[0] > l_max_time_shift)
                l_max_time_shift = _shape[i].shift[0];
        }
        depth = l_max_time_shift - l_min_time_shift;
        time_shift_ = 0 - l_min_time_shift;
        toggle_ = depth + 1;
        shape_ = _shape;
        for (int i = 0; i < N_SIZE; ++i) {
            for (int r = 0; r < N_RANK; ++r) {
                slope_[r] = max(slope_[r], abs((int)ceil((float)_shape[i].shift[r+1]/(l_max_time_shift - _shape[i].shift[0]))));
            }
        }
        shape_size_ = N_SIZE;
    }
    Pochoir_Shape<N_RANK> * Get_Shape() { return shape_; }
    int Get_Shape_Size() { return shape_size_; }
    T_Kernel & Get_Kernel(void) { return (kernel_); }
    template <typename ... IS>
    inline void operator() (int t, IS ... is) const { kernel_(t, is ...); }
    ~Pochoir_Kernel() { }
};

/* Pochoir_Obase_Kernel for Phase II */
/* define an abstract base class to act as a generic class pointer
 * to Pochoir_Obase_Kernel<F, N_RANK>
 */
template <int N_RANK>
struct Pochoir_Base_Kernel {
    // typedef std::function<void (int, int, const Grid_Info<N_RANK> &)> T_Kernel;
    // virtual T_Kernel & Get_Kernel (void) = 0;
    virtual Pochoir_Shape<N_RANK> * Get_Shape() = 0;
    virtual void operator()(int t0, int t1, Grid_Info<N_RANK> const & grid) const = 0;
    virtual int Get_Shape_Size() = 0;
};

/* latest version */
template <typename F, int N_RANK>
struct Pochoir_Obase_Kernel : public Pochoir_Base_Kernel<N_RANK> {
    F & kernel_;
    Pochoir_Shape<N_RANK> * shape_;
    int shape_size_, time_shift_, toggle_, slope_[N_RANK];
    Pochoir_Obase_Kernel(void) { shape_ = NULL; }
    template <int N_SIZE>
    Pochoir_Obase_Kernel(F & _kernel, Pochoir_Shape<N_RANK> (& _shape)[N_SIZE]) : kernel_(_kernel) {
        int l_min_time_shift=0, l_max_time_shift=0, depth=0;
        for (int r = 0; r < N_RANK; ++r) {
            slope_[r] = 0;
        }
        for (int i = 0; i < N_SIZE; ++i) {
            if (_shape[i].shift[0] < l_min_time_shift)
                l_min_time_shift = _shape[i].shift[0];
            if (_shape[i].shift[0] > l_max_time_shift)
                l_max_time_shift = _shape[i].shift[0];
        }
        depth = l_max_time_shift - l_min_time_shift;
        time_shift_ = 0 - l_min_time_shift;
        toggle_ = depth + 1;
        for (int i = 0; i < N_SIZE; ++i) {
            for (int r = 0; r < N_RANK; ++r) {
                slope_[r] = max(slope_[r], abs((int)ceil((float)_shape[i].shift[r+1]/(l_max_time_shift - _shape[i].shift[0]))));
            }
        }
        shape_size_ = N_SIZE;
        shape_ = _shape;
#ifdef DEBUG
        printf("time_shift_ = %d, toggle = %d\n", time_shift_, toggle_);
        for (int r = 0; r < N_RANK; ++r) {
            printf("slope[%d] = %d, ", r, slope_[r]);
        }
        printf("\n");
#endif
    }
    Pochoir_Shape<N_RANK> * Get_Shape() { return shape_; }
    int Get_Shape_Size() { return shape_size_; }
    // F & Get_Kernel (void) { return kernel_; }
    void operator() (int t0, int t1, Grid_Info<N_RANK> const & grid) const {
        kernel_(t0, t1, grid);
    }
    ~Pochoir_Obase_Kernel() { }
};

template <int N_RANK>
struct Pochoir_Guard {
    typedef typename Pochoir_Types<N_RANK>::T_Guard T_Guard;
    T_Guard guard_;
    Pochoir_Guard(void) { }
    Pochoir_Guard(T_Guard _guard) : guard_(_guard) { }
    template <typename ... IS>
    inline bool operator() (int t, IS ... is) const { return (guard_)(t, is ...); }
    T_Guard & Get_Guard (void) { return (guard_); }
};

#define Pochoir_Guard_3D_Begin(name, t, i, j, k) \
    auto __##name##__ = [&](int t, int i, int j, int k) -> bool {

#define Pochoir_Guard_2D_Begin(name, t, i, j) \
    auto __##name##__ = [&](int t, int i, int j) -> bool {

#define Pochoir_Guard_1D_Begin(name, t, i) \
    auto __##name##__ = [&](int t, int i) -> bool {

#define Pochoir_Guard_3D_End(name) \
    }; \
    Pochoir_Guard<3> name(__##name##__);

#define Pochoir_Guard_2D_End(name) \
    }; \
    Pochoir_Guard<2> name(__##name##__);

#define Pochoir_Guard_1D_End(name) \
    }; \
    Pochoir_Guard<1> name(__##name##__);

/* Version for Kernel in Phase I */
#define Pochoir_Kernel_1D_Begin(name, t, i) \
    auto __##name##__ = [&](int t, int i) { 

#define Pochoir_Kernel_2D_Begin(name, t, i, j) \
    auto __##name##__ = [&](int t, int i, int j) {

#define Pochoir_Kernel_3D_Begin(name, t, i, j, k) \
    auto __##name##__ = [&](int t, int i, int j, int k) {

#define Pochoir_Kernel_1D_End(name, shape) \
    }; \
    Pochoir_Kernel<1> name(__##name##__, shape);

#define Pochoir_Kernel_2D_End(name, shape) \
    }; \
    Pochoir_Kernel<2> name(__##name##__, shape);

#define Pochoir_Kernel_3D_End(name, shape) \
    }; \
    Pochoir_Kernel<3> name(__##name##__, shape);

/* Default Guard: Always return true */
Pochoir_Guard_3D_Begin(Default_Guard_3D, t, i, j, k)
    return true;
Pochoir_Guard_3D_End(Default_Guard_3D)

Pochoir_Guard_2D_Begin(Default_Guard_2D, t, i, j)
    return true;
Pochoir_Guard_2D_End(Default_Guard_2D)

Pochoir_Guard_1D_Begin(Default_Guard_1D, t, i)
    return true;
Pochoir_Guard_1D_End(Default_Guard_1D)

template <int N_RANK>
struct Pochoir_Tile_Kernel {
    /* for a N-dimensional stencil, it can have at most N+1-dimensional tile
     * because we have to count the time dimension
     */
    int size_[N_RANK+1], stride_[N_RANK+1], pointer_[N_RANK+1];
    int total_size_;
    /* Each Pochoir_Guard_Kernel can guard a tile of different kernel functions */
    Pochoir_Kernel<N_RANK> * kernel_;
    /* Each Pochoir_Guard_Kernel can hold one and only one guard function */
    Pochoir_Tile_Kernel() {
        for (int i = 0; i < N_RANK+1; ++i) {
            size_[i] = stride_[i] = pointer_[i] = 0;
        }
        total_size_ = 0;
        kernel_ = NULL;
    }
    /* if we use the default copy constructor, the kernel_ will only get a copy of
     * address instead of the real kernels
     */
    explicit Pochoir_Tile_Kernel(Pochoir_Tile_Kernel<N_RANK> const & t) {
        for (int i = 0; i < N_RANK + 1; ++i) {
            size_[i] = t.size_[i];
            stride_[i] = t.stride_[i];
            pointer_[i] = t.pointer_[i];
        }
        del_arr(kernel_);
        total_size_ = t.total_size_;
        kernel_ = new Pochoir_Kernel<N_RANK>[total_size_];
        for (int i = 0; i < total_size_; ++i) {
            kernel_[i] = t.kernel_[i];
        }
    }
    Pochoir_Tile_Kernel<N_RANK> & operator= (Pochoir_Tile_Kernel<N_RANK> const & t) {
        for (int i = 0; i < N_RANK + 1; ++i) {
            size_[i] = t.size_[i];
            stride_[i] = t.stride_[i];
            pointer_[i] = t.pointer_[i];
        }
        del_arr(kernel_);
        total_size_ = t.total_size_;
        kernel_ = new Pochoir_Kernel<N_RANK>[total_size_];
        for (int i = 0; i < total_size_; ++i) {
            kernel_[i] = t.kernel_[i];
        }
        return (*this);
    }
    ~Pochoir_Tile_Kernel() {
        del_arr(kernel_);
    }
};

template <int N_RANK>
struct Pochoir_Combined_Obase_Kernel {
    /* For obased kernel, all staggered kernels should already be unrolled
     * by Pochoir compiler
     */
    int unroll_;
    Pochoir_Base_Kernel<N_RANK> * kernel_;
    Pochoir_Base_Kernel<N_RANK> * cond_kernel_;
    Pochoir_Base_Kernel<N_RANK> * bkernel_;
    Pochoir_Base_Kernel<N_RANK> * cond_bkernel_;
};      

template <int N_RANK>
struct Pochoir_Run_Regional_Guard_Tile_Kernel {
    int const time_shift_;
    Pochoir_Guard<N_RANK> * pg_;
    Pochoir_Tile_Kernel<N_RANK> * pt_;
    Pochoir_Run_Regional_Guard_Tile_Kernel(int _time_shift, Pochoir_Guard<N_RANK> * _pg, Pochoir_Tile_Kernel<N_RANK> * _pt) : time_shift_(_time_shift), pg_(_pg), pt_(_pt) { }
    ~Pochoir_Run_Regional_Guard_Tile_Kernel() { }
    template <typename I>
    int set_pointer(int dim, I i) const { 
        assert(dim == 0);
        int l_idx = pt_->size_[0] == 0 ? 0 : (i % pt_->size_[0]) * pt_->stride_[0];
        return (l_idx);
    }
    template <typename I, typename ... IS>
    int set_pointer(int dim, I i, IS ... is) const {
        if (dim == 0) {
            int l_idx = pt_->size_[0] == 0 ? 0 : (i % pt_->size_[0]) * pt_->stride_[0];
            return (l_idx);
        } else {
            if (pt_->size_[dim] == 0) {
                return set_pointer(dim-1, i, is ...);
            } else {
                return ((i % pt_->size_[dim]) * pt_->stride_[dim] + set_pointer(dim-1, is ...));
            }
        }
    }
    template <typename ... IS>
    inline void operator() (int t, IS ... is) const {
        if ((*pg_)(t, is ...)) {
            int l_kernel_pointer = set_pointer(N_RANK, t, is ...);
            (pt_->kernel_[l_kernel_pointer])(t, is ...);
        }
    }
};

/* So far, we assume that there is no inhomogeneity in time,
 * the inhomogeneity only occurs in all spatial dimensions.
 */
template <int N_RANK>
struct Color_Region {
    /* sizeof(int) = 4 bytes = 32 bits -- let's assume that user won't
     * register more than 32 different kernels for now
     */
    int sz_pgk_;
    Vector_Info< Pochoir_Guard<N_RANK> * > & pgs_;
    Grid_Info<N_RANK> & phys_grid_;
    Color_Region(Vector_Info< Pochoir_Guard<N_RANK> * > & _pgs, Grid_Info<N_RANK> & _grid) : sz_pgk_(_pgs.size()), pgs_(_pgs), phys_grid_(_grid) { }
    T_color get_color() { return 0; }
    T_color operator() (int t0, int t1, Grid_Info<N_RANK> const & grid) { 
        /* return 0 to make g++ happy! */
        return 0; 
    }
};

template <>
struct Color_Region<1> {
    int sz_pgk_;
    Vector_Info< Pochoir_Guard<1> * > & pgs_;
    Grid_Info<1> & phys_grid_;
    Color_Region(Vector_Info< Pochoir_Guard<1> * > & _pgs, Grid_Info<1> & _grid) : sz_pgk_(_pgs.size()), pgs_(_pgs), phys_grid_(_grid) { }

    T_color get_color(int t, int i) {
        T_color l_color = 0;
        for (int pt = 0; pt < sz_pgk_; ++pt) {
            l_color <<= 1;
            bool is_color = (*pgs_[pt])(t, i);
            if (is_color)
                l_color |= 1;
        }
        return l_color;
    }

    Homogeneity operator() (int t0, int t1, Grid_Info<1> const & grid) {
        Grid_Info<1> l_grid = grid;
        int start_i = pmod_lu(l_grid.x0[0], phys_grid_.x0[0], phys_grid_.x1[0]);
        T_color start_color = get_color(t0, start_i);
        T_color l_o = start_color, l_a = start_color;
        for (int t = t0; t < t1; ++t) {
            for (int i = l_grid.x0[0]; i < l_grid.x1[0]; ++i) {
                const int new_i = pmod_lu(i, phys_grid_.x0[0], phys_grid_.x1[0]);
                T_color l_color = get_color(t, new_i);
                l_o |= l_color; l_a &= l_color;
            }
            /* Adjust trapezoid */
            l_grid.x0[0] += l_grid.dx0[0]; l_grid.x1[0] += l_grid.dx1[0];
        }
        return Homogeneity(l_o, l_a, sz_pgk_);
    }
};

template <>
struct Color_Region<2> {
    int sz_pgk_;
    Vector_Info< Pochoir_Guard<2> * > & pgs_;
    Grid_Info<2> & phys_grid_;
    Color_Region(Vector_Info< Pochoir_Guard<2> * > & _pgs, Grid_Info<2> & _grid) : sz_pgk_(_pgs.size()), pgs_(_pgs), phys_grid_(_grid) { }

    T_color get_color(int t, int i, int j) {
        T_color l_color = 0;
        for (int pt = 0; pt < sz_pgk_; ++pt) {
            l_color <<= 1;
            bool is_color = (*pgs_[pt])(t, i, j);
            if (is_color)
                l_color |= 1;
        }
        return l_color;
    }

    Homogeneity operator() (int t0, int t1, Grid_Info<2> const & grid) {
        Grid_Info<2> l_grid = grid;
        int start_i = pmod_lu(l_grid.x0[1], phys_grid_.x0[1], phys_grid_.x1[1]);
        int start_j = pmod_lu(l_grid.x0[0], phys_grid_.x0[0], phys_grid_.x1[0]);
        T_color start_color = get_color(t0, start_i, start_j);
        T_color l_o = start_color, l_a = start_color;
        for (int t = t0; t < t1; ++t) {
            for (int i = l_grid.x0[1]; i < l_grid.x1[1]; ++i) {
                const int new_i = pmod_lu(i, phys_grid_.x0[1], phys_grid_.x1[1]);
        for (int j = l_grid.x0[0]; j < l_grid.x1[0]; ++j) {
            const int new_j = pmod_lu(j, phys_grid_.x0[0], phys_grid_.x1[0]);
                T_color l_color = get_color(t, new_i, new_j);
                l_o |= l_color; l_a &= l_color;
        }
            }
            /* Adjust trapezoid */
            l_grid.x0[0] += l_grid.dx0[0]; l_grid.x1[0] += l_grid.dx1[0];
            l_grid.x0[1] += l_grid.dx0[1]; l_grid.x1[1] += l_grid.dx1[1];
        }
        return Homogeneity(l_o, l_a, sz_pgk_);
    }
};

template <>
struct Color_Region<3> {
    int sz_pgk_;
    Vector_Info< Pochoir_Guard<3> * > & pgs_;
    Grid_Info<3> & phys_grid_;
    Color_Region(Vector_Info< Pochoir_Guard<3> * > & _pgs, Grid_Info<3> & _grid) : sz_pgk_(_pgs.size()), pgs_(_pgs), phys_grid_(_grid) { }

    T_color get_color(int t, int i, int j, int k) {
        T_color l_color = 0;
        for (int pt = 0; pt < sz_pgk_; ++pt) {
            l_color <<= 1;
            bool is_color = (*pgs_[pt])(t, i, j, k);
            if (is_color)
                l_color |= 1;
        }
        return l_color;
    }

    Homogeneity operator() (int t0, int t1, Grid_Info<3> const & grid) {
        Grid_Info<3> l_grid = grid;
        int start_i = pmod_lu(l_grid.x0[2], phys_grid_.x0[2], phys_grid_.x1[2]);
        int start_j = pmod_lu(l_grid.x0[1], phys_grid_.x0[1], phys_grid_.x1[1]);
        int start_k = pmod_lu(l_grid.x0[0], phys_grid_.x0[0], phys_grid_.x1[0]);
        T_color start_color = get_color(t0, start_i, start_j, start_k);
        T_color l_o = start_color, l_a = start_color;
        for (int t = t0; t < t1; ++t) {
            for (int i = l_grid.x0[2]; i < l_grid.x1[2]; ++i) {
                const int new_i = pmod_lu(i, phys_grid_.x0[2], phys_grid_.x1[2]);
        for (int j = l_grid.x0[1]; j < l_grid.x1[1]; ++j) {
            const int new_j = pmod_lu(j, phys_grid_.x0[1], phys_grid_.x1[1]);
            for (int k = l_grid.x0[0]; k < l_grid.x1[0]; ++k) {
                const int new_k = pmod_lu(k, phys_grid_.x0[0], phys_grid_.x1[0]);
                T_color l_color = get_color(t, new_i, new_j, new_k);
                l_o |= l_color; l_a &= l_color;
            }
        }
            }
            /* Adjust trapezoid */
            l_grid.x0[0] += l_grid.dx0[0]; l_grid.x1[0] += l_grid.dx1[0];
            l_grid.x0[1] += l_grid.dx0[1]; l_grid.x1[1] += l_grid.dx1[1];
            l_grid.x0[2] += l_grid.dx0[2]; l_grid.x1[2] += l_grid.dx1[2];
        }
        return Homogeneity(l_o, l_a, sz_pgk_);
    }
};

#endif /* POCHOIR_KERNEL_HPP */
