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

#ifndef EXPR_STENCIL_HPP
#define EXPR_STENCIL_HPP

#include "pochoir_common.hpp"
#include "pochoir_array.hpp"
/* assuming there won't be more than 10 Pochoir_Array in one Pochoir object! */
#define ARRAY_SIZE 10

template <int N_RANK>
class Pochoir {
    private:
        int slope_[N_RANK];
        grid_info<N_RANK> logic_grid_;
        grid_info<N_RANK> phys_grid_;
        int time_shift_;
        int toggle_;
        int unroll_;
        int timestep_;
        bool regArrayFlag, regLogicDomainFlag, regPhysDomainFlag, regShapeFlag;
        void checkFlag(bool flag, char const * str);
        void checkFlags(void);
        template <typename T_Array>
        void getPhysDomainFromArray(T_Array & arr);
        template <typename T_Array>
        void cmpPhysDomainFromArray(T_Array & arr);
        template <size_t N_SIZE>
        void Register_Shape(Pochoir_Shape<N_RANK> (& shape)[N_SIZE]);
        template <size_t N_SIZE1, size_t N_SIZE2>
        void Register_Shape(Pochoir_Shape<N_RANK> (& shape1)[N_SIZE1], Pochoir_Shape<N_RANK> (& shape2)[N_SIZE2]);
        Pochoir_Shape<N_RANK> * shape_;
        int shape_size_;
        int num_arr_;
        int arr_type_size_;
        int size_pochoir_guard_kernel_;
        /* assuming that the number of distinct sub-regions is less than 10 */
        Pochoir_Guard_Kernel<N_RANK> pochoir_guard_kernel_[10];

        /* Private Register Kernel Function */
        template <typename K>
        void reg_kernel(int pt, K k);
        template <typename K, typename ... KS>
        void reg_kernel(int pt, K k, KS ... ks);
        void reg_guard(typename Pochoir_Types<N_RANK>::T_Guard g);

    public:
    template <typename ... KS>
    void Register_Default_Kernel(KS ... ks);
    template <typename ... KS>
    void Register_Kernel(typename Pochoir_Types<N_RANK>::T_Guard g, KS ... ks);
    // get slope(s)
    int slope(int const _idx) { return slope_[_idx]; }
    template <size_t N_SIZE>
    Pochoir(Pochoir_Shape<N_RANK> (& shape)[N_SIZE]) {
        for (int i = 0; i < N_RANK; ++i) {
            slope_[i] = 0;
            logic_grid_.x0[i] = logic_grid_.x1[i] = logic_grid_.dx0[i] = logic_grid_.dx1[i] = 0;
            phys_grid_.x0[i] = phys_grid_.x1[i] = phys_grid_.dx0[i] = phys_grid_.dx1[i] = 0;
        }
        timestep_ = 0;
        regArrayFlag = regLogicDomainFlag = regPhysDomainFlag = regShapeFlag = false;
        Register_Shape(shape);
        regShapeFlag = true;
        num_arr_ = 0;
        arr_type_size_ = 0;
        size_pochoir_guard_kernel_ = 0;
    }

    template <size_t N_SIZE1, size_t N_SIZE2>
    Pochoir(Pochoir_Shape<N_RANK> (& shape1)[N_SIZE1], Pochoir_Shape<N_RANK> (& shape2)[N_SIZE2]) {
        for (int i = 0; i < N_RANK; ++i) {
            slope_[i] = 0;
            logic_grid_.x0[i] = logic_grid_.x1[i] = logic_grid_.dx0[i] = logic_grid_.dx1[i] = 0;
            phys_grid_.x0[i] = phys_grid_.x1[i] = phys_grid_.dx0[i] = phys_grid_.dx1[i] = 0;
        }
        timestep_ = 0;
        regArrayFlag = regLogicDomainFlag = regPhysDomainFlag = regShapeFlag = false;
        Register_Shape(shape1, shape2);
        regShapeFlag = true;
        num_arr_ = 0;
        arr_type_size_ = 0;
        size_pochoir_guard_kernel_ = 0;
    }
    /* currently, we just compute the slope[] out of the shape[] */
    /* We get the grid_info out of arrayInUse */
    template <typename T>
    void Register_Array(Pochoir_Array<T, N_RANK> & arr);

    /* We should still keep the Register_Domain for zero-padding!!! */
    template <typename D>
    void Register_Domain(D const & d);
    template <typename D, typename ... DS>
    void Register_Domain(D const & d, DS ... ds);
    /* register boundary value function with corresponding Pochoir_Array object directly */
    template <typename T_Array, typename RET>
    void registerBoundaryFn(T_Array & arr, RET (*_bv)(T_Array &, int, int, int)) {
        arr.Register_Boundary(_bv);
        Register_Array(arr);
    } 
    grid_info<N_RANK> get_phys_grid(void);

    /* Executable Spec */
    template <typename BF>
    void Run(int timestep, BF const & bf);
    void Run(int timestep);
    /* obase for zero-padded region */
    template <typename F>
    void Run_Obase(int timestep, F const & f);
    /* obase for interior and ExecSpec for boundary */
    template <typename F, typename BF>
    void Run_Obase(int timestep, F const & f, BF const & bf);
};

template <int N_RANK> template <typename K>
void Pochoir<N_RANK>::reg_kernel(int pt, K k) {
    pochoir_guard_kernel_[size_pochoir_guard_kernel_].kernel_[pt] = k; 
}

template <int N_RANK> template <typename K, typename ... KS>
void Pochoir<N_RANK>::reg_kernel(int pt, K k, KS ... ks) {
    pochoir_guard_kernel_[size_pochoir_guard_kernel_].kernel_[pt] = k;
    reg_kernel(pt+1, ks ...);
}

template <int N_RANK>
void Pochoir<N_RANK>::reg_guard(typename Pochoir_Types<N_RANK>::T_Guard g) {
    pochoir_guard_kernel_[size_pochoir_guard_kernel_].guard_ = g;
    return;
}
    
template <int N_RANK> template <typename ... KS>
void Pochoir<N_RANK>::Register_Kernel(typename Pochoir_Types<N_RANK>::T_Guard g, KS ... ks) {
    int l_size = sizeof...(KS);
    typedef typename Pochoir_Types<N_RANK>::T_Kernel T_Kernel;
    assert(size_pochoir_guard_kernel_ < 10);
    if (size_pochoir_guard_kernel_ >= 10) {
        printf("Pochoir Error: Register_Kernel > %d\n", size_pochoir_guard_kernel_);
        exit(1);
    }
    pochoir_guard_kernel_[size_pochoir_guard_kernel_].size_ = l_size;
    pochoir_guard_kernel_[size_pochoir_guard_kernel_].pointer_ = 0;
    pochoir_guard_kernel_[size_pochoir_guard_kernel_].kernel_ = (T_Kernel *) calloc(l_size, sizeof(T_Kernel));
    reg_guard(g);
    reg_kernel(0, ks ...);
    ++size_pochoir_guard_kernel_;
}

template <int N_RANK>
void Pochoir<N_RANK>::checkFlag(bool flag, char const * str) {
    if (!flag) {
        printf("\nPochoir registration error:\n");
        printf("You forgot to register %s.\n", str);
        exit(1);
    }
}

template <int N_RANK>
void Pochoir<N_RANK>::checkFlags(void) {
    checkFlag(regArrayFlag, "Pochoir array");
    checkFlag(regLogicDomainFlag, "Logic Domain");
    checkFlag(regPhysDomainFlag, "Physical Domain");
    checkFlag(regShapeFlag, "Shape");
    return;
}

template <int N_RANK> template <typename T_Array> 
void Pochoir<N_RANK>::getPhysDomainFromArray(T_Array & arr) {
    /* get the physical grid */
    for (int i = 0; i < N_RANK; ++i) {
        phys_grid_.x0[i] = 0; phys_grid_.x1[i] = arr.size(i);
        /* if logic domain is not set, let's set it the same as physical grid */
        if (!regLogicDomainFlag) {
            logic_grid_.x0[i] = 0; logic_grid_.x1[i] = arr.size(i);
        }
    }

    regPhysDomainFlag = true;
    regLogicDomainFlag = true;
}

template <int N_RANK> template <typename T_Array> 
void Pochoir<N_RANK>::cmpPhysDomainFromArray(T_Array & arr) {
    /* check the consistency of all engaged Pochoir_Array */
    for (int j = 0; j < N_RANK; ++j) {
        if (arr.size(j) != phys_grid_.x1[j]) {
            printf("Pochoir array size mismatch error:\n");
            printf("Registered Pochoir arrays have different sizes!\n");
            exit(1);
        }
    }
}

template <int N_RANK> template <typename T>
void Pochoir<N_RANK>::Register_Array(Pochoir_Array<T, N_RANK> & arr) {
    if (!regShapeFlag) {
        printf("Please register Shape before register Array!\n");
        exit(1);
    }

    if (num_arr_ == 0) {
        arr_type_size_ = sizeof(T);
        ++num_arr_;
#if DEBUG
        printf("arr_type_size = %d\n", arr_type_size_);
#endif
    } 
    if (!regPhysDomainFlag) {
        getPhysDomainFromArray(arr);
    } else {
        cmpPhysDomainFromArray(arr);
    }
    arr.Register_Shape(shape_, shape_size_, unroll_);
#if 0
    arr.set_slope(slope_);
    arr.set_toggle(toggle_);
    arr.alloc_mem();
#endif
    regArrayFlag = true;
}

template <int N_RANK> template <size_t N_SIZE>
void Pochoir<N_RANK>::Register_Shape(Pochoir_Shape<N_RANK> (& shape)[N_SIZE]) {
    /* currently we just get the slope_[] and toggle_ out of the shape[] */
    shape_ = new Pochoir_Shape<N_RANK>[N_SIZE];
    shape_size_ = N_SIZE;
    int l_min_time_shift=0, l_max_time_shift=0, depth=0;
    for (int i = 0; i < N_SIZE; ++i) {
        if (shape[i].shift[0] < l_min_time_shift)
            l_min_time_shift = shape[i].shift[0];
        if (shape[i].shift[0] > l_max_time_shift)
            l_max_time_shift = shape[i].shift[0];
        for (int r = 0; r < N_RANK+1; ++r) {
            shape_[i].shift[r] = shape[i].shift[r];
        }
    }
    depth = l_max_time_shift - l_min_time_shift;
    time_shift_ = 0 - l_min_time_shift;
    toggle_ = depth + 1;
    unroll_ = 1;
    for (int i = 0; i < N_SIZE; ++i) {
        for (int r = 1; r < N_RANK+1; ++r) {
            slope_[N_RANK-r] = max(slope_[N_RANK-r], abs((int)ceil((float)shape_[i].shift[r]/(l_max_time_shift - shape_[i].shift[0]))));
        }
    }
#if DEBUG 
    printf("time_shift_ = %d, toggle = %d\n", time_shift_, toggle_);
    for (int r = 0; r < N_RANK; ++r) {
        printf("slope[%d] = %d, ", r, slope_[r]);
    }
    printf("\n");
#endif
    regShapeFlag = true;
}

template <int N_RANK> template <size_t N_SIZE1, size_t N_SIZE2>
void Pochoir<N_RANK>::Register_Shape(Pochoir_Shape<N_RANK> (& shape1)[N_SIZE1], Pochoir_Shape<N_RANK> (& shape2)[N_SIZE2]) {
    /* currently we just get the slope_[] and toggle_ out of the shape[] */
    shape_ = new Pochoir_Shape<N_RANK>[N_SIZE1+N_SIZE2];
    shape_size_ = N_SIZE1+N_SIZE2;
    int l_min_time_shift=0, l_max_time_shift=0, depth=0;
    int i;
    for (i = 0; i < N_SIZE1; ++i) {
        if (shape1[i].shift[0] < l_min_time_shift)
            l_min_time_shift = shape1[i].shift[0];
        if (shape1[i].shift[0] > l_max_time_shift)
            l_max_time_shift = shape1[i].shift[0];
        for (int r = 0; r < N_RANK+1; ++r) {
            shape_[i].shift[r] = shape1[i].shift[r];
        }
    }
    for (i = 0; i < N_SIZE2; ++i) {
        if (shape2[i].shift[0] < l_min_time_shift)
            l_min_time_shift = shape2[i].shift[0];
        if (shape2[i].shift[0] > l_max_time_shift)
            l_max_time_shift = shape2[i].shift[0];
        for (int r = 0; r < N_RANK+1; ++r) {
            shape_[i+N_SIZE1].shift[r] = shape2[i].shift[r];
        }
    }
    depth = l_max_time_shift - l_min_time_shift;
    time_shift_ = 0 - l_min_time_shift;
    toggle_ = depth + 1;
    unroll_ = 2;
    for (i = 0; i < N_SIZE1+N_SIZE2; ++i) {
        for (int r = 1; r < N_RANK+1; ++r) {
            slope_[N_RANK-r] = max(slope_[N_RANK-r], abs((int)ceil((float)shape_[i].shift[r]/(l_max_time_shift - shape_[i].shift[0]))));
        }
    }
#if DEBUG 
    printf("time_shift_ = %d, toggle = %d\n", time_shift_, toggle_);
    for (int r = 0; r < N_RANK; ++r) {
        printf("slope[%d] = %d, ", r, slope_[r]);
    }
    printf("\n");
#endif
    regShapeFlag = true;
}

template <int N_RANK> template <typename D>
void Pochoir<N_RANK>::Register_Domain(D const & d) {
    logic_grid_.x0[0] = d.first();
    logic_grid_.x1[0] = d.first() + d.size();
    regLogicDomainFlag = true;
}

template <int N_RANK> template <typename D, typename ... DS>
void Pochoir<N_RANK>::Register_Domain(D const & d, DS ... ds) {
    int l_pointer = sizeof...(DS);
    logic_grid_.x0[l_pointer] = d.first();
    logic_grid_.x1[l_pointer] = d.first() + d.size();
}

template <int N_RANK> 
grid_info<N_RANK> Pochoir<N_RANK>::get_phys_grid(void) {
    return phys_grid_;
}

/* Run the kernel functions stored in array of function pointers */
template <int N_RANK>
void Pochoir<N_RANK>::Run(int timestep) {
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    algor.set_unroll(unroll_);
    timestep_ = timestep;
    /* base_case_kernel() will mimic exact the behavior of serial nested loop!
    */
    checkFlags();
    inRun = true;
    algor.base_case_kernel_guard(0 + time_shift_, timestep + time_shift_, logic_grid_, size_pochoir_guard_kernel_, pochoir_guard_kernel_);
    inRun = false;

}

/* obase for zero-padded area! */
template <int N_RANK> template <typename F>
void Pochoir<N_RANK>::Run_Obase(int timestep, F const & f) {
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    algor.set_unroll(unroll_);
    timestep_ = timestep;
    checkFlags();
#if BICUT
#if 0
    fprintf(stderr, "Call obase_bicut\n");
    algor.obase_bicut(0+time_shift_, timestep+time_shift_, logic_grid_, f);
#else
//     fprintf(stderr, "Call shorter_duo_sim_obase_bicut\n");
   // algor.sim_obase_bicut(0+time_shift_, timestep+time_shift_, logic_grid_, f);
    algor.shorter_duo_sim_obase_bicut(0+time_shift_, timestep+time_shift_, logic_grid_, f);
    // algor.duo_sim_obase_bicut(0+time_shift_, timestep+time_shift_, logic_grid_, f);
#if STAT
    for (int i = 1; i < SUPPORT_RANK; ++i) {
        fprintf(stderr, "sim_count_cut[%d] = %ld\n", i, algor.sim_count_cut[i].get_value());
    }
#endif
#endif
#else
    algor.obase_m(0+time_shift_, timestep+time_shift_, logic_grid_, f);
#endif
}

/* obase for interior and ExecSpec for boundary */
template <int N_RANK> template <typename F, typename BF>
void Pochoir<N_RANK>::Run_Obase(int timestep, F const & f, BF const & bf) {
    int l_total_points = 1;
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    algor.set_unroll(unroll_);
    /* this version uses 'f' to compute interior region, 
     * and 'bf' to compute boundary region
     */
    timestep_ = timestep;
    checkFlags();
#if BICUT
#if 0
    fprintf(stderr, "Call obase_bicut_boundary_P\n");
    algor.obase_bicut_boundary_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#else
//    fprintf(stderr, "Call sim_obase_bicut_P\n");
//    hyper-space cut
    // algor.sim_obase_bicut_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
    // cutting based on shorter bar
    algor.shorter_duo_sim_obase_bicut_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
    // cutting based on longer bar
    // algor.duo_sim_obase_bicut_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
    // serial space cut
    // algor.obase_bicut_boundary_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#if STAT
    for (int i = 1; i < SUPPORT_RANK; ++i) {
        fprintf(stderr, "sim_count_cut[%d] = %ld\n", i, algor.sim_count_cut[i].get_value());
    }
#endif
#endif
#else
    algor.obase_boundary_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#endif
}

#endif
