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
//#include "pochoir_modified_cuts.hpp"
//#include "sawzoid2.hpp"
#include "sawzoid.hpp"
#include "pochoir_walk_recursive.hpp"
#ifdef COUNT_PROJECTIONS
//#include "projections.hpp"
#endif
#include "pochoir_array.hpp"
//#include "clones_2d_heat.hpp"
//#include "clones_2d_diffusion.hpp"
//#include "clones_2d_output.hpp"
//#include "clones_2dwave.hpp"
//#include "clones.hpp"

#ifdef KERNEL_SELECTION
#include "kernel_selection_trap.hpp"
//#include "kernel_selection_sawzoid.hpp"
#include "kernel_selection_sawzoid_middle.hpp"
#elif defined GENEITY_TEST
//#include "symbolic_walk_better_memory.hpp"
//#include "geneity_problem_trap.hpp"
//#include "pochoir_modified_cuts_heterogeneity.hpp"
#include "sawzoid_middle_heterogeneity.hpp"
//#include "symbolic_walk.hpp"
#elif defined AUTO_TUNE
//#include "auto_tuning_trap.hpp"
//#include "auto_tuning_sawzoid.hpp"
//#include "auto_tuning_sawzoid_middle.hpp"
#include "auto_tuning_arbitrary_cuts_sawzoid.hpp"
#include "auto_tuning_arbitrary_cuts_trap.hpp"
//#include "auto_tuning_arbitrary_cuts_sawzoid_middle.hpp"
//#include "dag_sawzoid.hpp"
#endif

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
        Pochoir_Shape<N_RANK> * shape_;
        int shape_size_;
        int num_arr_;
        int arr_type_size_;
		//eka - adding a pointer to pochoir array
		//Pochoir_Array<double, N_RANK> * arr_ ;
		void * arr_ ;
		int resolution_ ;
		ofstream * outputFile_ ;
		char * problem_name ; //name of the problem we are solving
    public:

	void set_resolution(int r)
	{
		resolution_ = r ;
	}

	void set_output_file(ofstream * file)
	{
		outputFile_ = file ;
	}

	void set_problem_name(char * name)
	{
		problem_name = name ;
	}

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
    }
    /* currently, we just compute the slope[] out of the shape[] */
    /* We get the grid_info out of arrayInUse */
    template <typename T>
    void Register_Array(Pochoir_Array<T, N_RANK> & arr);

    /* We should still keep the Register_Domain for zero-padding!!! */
    template <typename Domain>
    void Register_Domain(Domain const & i);
    template <typename Domain>
    void Register_Domain(Domain const & i, Domain const & j);
    template <typename Domain>
    void Register_Domain(Domain const & i, Domain const & j, Domain const & k);
    template <typename Domain>
    void Register_Domain(Domain const & i, Domain const & j, Domain const & k, Domain const & l);
    template <typename Domain>
    void Register_Domain(Domain const & i, Domain const & j, Domain const & k, Domain const & l, Domain const & m);
    template <typename Domain>
    void Register_Domain(Domain const & i, Domain const & j, Domain const & k, Domain const & l, Domain const & m, Domain const & n);
    template <typename Domain>
    void Register_Domain(Domain const & i, Domain const & j, Domain const & k, Domain const & l, Domain const & m, Domain const & n, Domain const & o);
    template <typename Domain>
    void Register_Domain(Domain const & i, Domain const & j, Domain const & k, Domain const & l, Domain const & m, Domain const & n, Domain const & o, Domain const & p);

    /* register boundary value function with corresponding Pochoir_Array object directly */
    template <typename T_Array, typename RET>
    void registerBoundaryFn(T_Array & arr, RET (*_bv)(T_Array &, int, int, int)) {
        arr.Register_Boundary(_bv);
        Register_Array(arr);
    } 
    /* Executable Spec */
    template <typename BF>
    void Run(int timestep, BF const & bf);
    /* safe/unsafe Executable Spec */
    template <typename F, typename BF>
    void Run(int timestep, F const & f, BF const & bf);
    /* obase for zero-padded region */
    template <typename F>
    void Run_Obase(int timestep, F const & f);
    /* obase for interior and ExecSpec for boundary */
    template <typename F, typename BF>
    void Run_Obase(int timestep, F const & f, BF const & bf);
};

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
        cout << "Please register Shape before register Array!" << endl;
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
    arr.Register_Shape(shape_, shape_size_);
#if 0
    arr.set_slope(slope_);
    arr.set_toggle(toggle_);
    arr.alloc_mem();
#endif
    regArrayFlag = true;
	//eka - storing a pointer to pochoir array
	//arr_ = &(arr) ;
	arr_ = arr.data() ;
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
    for (int i = 0; i < N_SIZE; ++i) {
        for (int r = 0; r < N_RANK; ++r) {
            slope_[r] = max(slope_[r], abs((int)ceil((float)shape[i].shift[r+1]/(l_max_time_shift - shape[i].shift[0]))));
        }
    }
#if DEBUG 
    cout << "time_shift_ = " << time_shift_ << ", toggle = " << toggle_ << endl;
    for (int r = 0; r < N_RANK; ++r) {
        printf("slope[%d] = %d, ", r, slope_[r]);
    }
    printf("\n");
#endif
    regShapeFlag = true;
}

template <int N_RANK> template <typename Domain>
void Pochoir<N_RANK>::Register_Domain(Domain const & r_i, Domain const & r_j, Domain const & r_k, Domain const & r_l, Domain const & r_m, Domain const & r_n, Domain const & r_o, Domain const & r_p) {
    logic_grid_.x0[7] = r_i.first();
    logic_grid_.x1[7] = r_i.first() + r_i.size();
    logic_grid_.x0[6] = r_j.first();
    logic_grid_.x1[6] = r_j.first() + r_j.size();
    logic_grid_.x0[5] = r_k.first();
    logic_grid_.x1[5] = r_k.first() + r_k.size();
    logic_grid_.x0[4] = r_l.first();
    logic_grid_.x1[4] = r_l.first() + r_l.size();
    logic_grid_.x0[3] = r_m.first();
    logic_grid_.x1[3] = r_m.first() + r_m.size();
    logic_grid_.x0[2] = r_n.first();
    logic_grid_.x1[2] = r_n.first() + r_n.size();
    logic_grid_.x0[1] = r_o.first();
    logic_grid_.x1[1] = r_o.first() + r_o.size();
    logic_grid_.x0[0] = r_p.first();
    logic_grid_.x1[0] = r_p.first() + r_p.size();
    regLogicDomainFlag = true;
}

template <int N_RANK> template <typename Domain>
void Pochoir<N_RANK>::Register_Domain(Domain const & r_i, Domain const & r_j, Domain const & r_k, Domain const & r_l, Domain const & r_m, Domain const & r_n, Domain const & r_o) {
    logic_grid_.x0[6] = r_i.first();
    logic_grid_.x1[6] = r_i.first() + r_i.size();
    logic_grid_.x0[5] = r_j.first();
    logic_grid_.x1[5] = r_j.first() + r_j.size();
    logic_grid_.x0[4] = r_k.first();
    logic_grid_.x1[4] = r_k.first() + r_k.size();
    logic_grid_.x0[3] = r_l.first();
    logic_grid_.x1[3] = r_l.first() + r_l.size();
    logic_grid_.x0[2] = r_m.first();
    logic_grid_.x1[2] = r_m.first() + r_m.size();
    logic_grid_.x0[1] = r_n.first();
    logic_grid_.x1[1] = r_n.first() + r_n.size();
    logic_grid_.x0[0] = r_o.first();
    logic_grid_.x1[0] = r_o.first() + r_o.size();
    regLogicDomainFlag = true;
}

template <int N_RANK> template <typename Domain>
void Pochoir<N_RANK>::Register_Domain(Domain const & r_i, Domain const & r_j, Domain const & r_k, Domain const & r_l, Domain const & r_m, Domain const & r_n) {
    logic_grid_.x0[5] = r_i.first();
    logic_grid_.x1[5] = r_i.first() + r_i.size();
    logic_grid_.x0[4] = r_j.first();
    logic_grid_.x1[4] = r_j.first() + r_j.size();
    logic_grid_.x0[3] = r_k.first();
    logic_grid_.x1[3] = r_k.first() + r_k.size();
    logic_grid_.x0[2] = r_l.first();
    logic_grid_.x1[2] = r_l.first() + r_l.size();
    logic_grid_.x0[1] = r_m.first();
    logic_grid_.x1[1] = r_m.first() + r_m.size();
    logic_grid_.x0[0] = r_n.first();
    logic_grid_.x1[0] = r_n.first() + r_n.size();
    regLogicDomainFlag = true;
}

template <int N_RANK> template <typename Domain>
void Pochoir<N_RANK>::Register_Domain(Domain const & r_i, Domain const & r_j, Domain const & r_k, Domain const & r_l, Domain const & r_m) {
    logic_grid_.x0[4] = r_i.first();
    logic_grid_.x1[4] = r_i.first() + r_i.size();
    logic_grid_.x0[3] = r_j.first();
    logic_grid_.x1[3] = r_j.first() + r_j.size();
    logic_grid_.x0[2] = r_k.first();
    logic_grid_.x1[2] = r_k.first() + r_k.size();
    logic_grid_.x0[1] = r_l.first();
    logic_grid_.x1[1] = r_l.first() + r_l.size();
    logic_grid_.x0[0] = r_m.first();
    logic_grid_.x1[0] = r_m.first() + r_m.size();
    regLogicDomainFlag = true;
}

template <int N_RANK> template <typename Domain>
void Pochoir<N_RANK>::Register_Domain(Domain const & r_i, Domain const & r_j, Domain const & r_k, Domain const & r_l) {
    logic_grid_.x0[3] = r_i.first();
    logic_grid_.x1[3] = r_i.first() + r_i.size();
    logic_grid_.x0[2] = r_j.first();
    logic_grid_.x1[2] = r_j.first() + r_j.size();
    logic_grid_.x0[1] = r_k.first();
    logic_grid_.x1[1] = r_k.first() + r_k.size();
    logic_grid_.x0[0] = r_l.first();
    logic_grid_.x1[0] = r_l.first() + r_l.size();
    regLogicDomainFlag = true;
}

template <int N_RANK> template <typename Domain>
void Pochoir<N_RANK>::Register_Domain(Domain const & r_i, Domain const & r_j, Domain const & r_k) {
    logic_grid_.x0[2] = r_i.first();
    logic_grid_.x1[2] = r_i.first() + r_i.size();
    logic_grid_.x0[1] = r_j.first();
    logic_grid_.x1[1] = r_j.first() + r_j.size();
    logic_grid_.x0[0] = r_k.first();
    logic_grid_.x1[0] = r_k.first() + r_k.size();
    regLogicDomainFlag = true;
}

template <int N_RANK> template <typename Domain>
void Pochoir<N_RANK>::Register_Domain(Domain const & r_i, Domain const & r_j) {
    logic_grid_.x0[1] = r_i.first();
    logic_grid_.x1[1] = r_i.first() + r_i.size();
    logic_grid_.x0[0] = r_j.first();
    logic_grid_.x1[0] = r_j.first() + r_j.size();
    regLogicDomainFlag = true;
}

template <int N_RANK> template <typename Domain>
void Pochoir<N_RANK>::Register_Domain(Domain const & r_i) {
    logic_grid_.x0[0] = r_i.first();
    logic_grid_.x1[0] = r_i.first() + r_i.size();
    regLogicDomainFlag = true;
}

/* Executable Spec */
template <int N_RANK> template <typename BF>
void Pochoir<N_RANK>::Run(int timestep, BF const & bf) {
    /* this version uses 'f' to compute interior region, 
     * and 'bf' to compute boundary region
     */
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    timestep_ = timestep;
    /* base_case_kernel() will mimic exact the behavior of serial nested loop!
    */
    checkFlags();
    inRun = true;
    algor.base_case_kernel_boundary(0 + time_shift_, timestep + time_shift_, logic_grid_, bf);
    inRun = false;
    // algor.sim_bicut_zero(0 + time_shift_, timestep + time_shift_, logic_grid_, bf);
    /* obase_boundary_p() is a parallel divide-and-conquer algorithm, which checks
     * boundary for every point
     */
    // algor.obase_boundary_p(0, timestep, logic_grid_, bf);
}

/* safe/non-safe ExecSpec */
template <int N_RANK> template <typename F, typename BF>
void Pochoir<N_RANK>::Run(int timestep, F const & f, BF const & bf) {
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    /* this version uses 'f' to compute interior region, 
     * and 'bf' to compute boundary region
     */
    timestep_ = timestep;
    checkFlags();
#pragma isat marker M2_begin
#if BICUT
#if 1
    algor.walk_bicut_boundary_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#else
    algor.sim_bicut_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#endif
#else
    algor.walk_ncores_boundary_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#endif
#pragma isat marker M2_end
}

/* obase for zero-padded area! */
template <int N_RANK> template <typename F>
void Pochoir<N_RANK>::Run_Obase(int timestep, F const & f) {
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    timestep_ = timestep;
    checkFlags();
#if BICUT
#if 0
    fprintf(stderr, "Call obase_bicut\n");
#pragma isat marker M2_begin
    algor.obase_bicut(0+time_shift_, timestep+time_shift_, logic_grid_, f);
#pragma isat marker M2_end
#else
//     fprintf(stderr, "Call shorter_duo_sim_obase_bicut\n");
#pragma isat marker M2_begin
   // algor.sim_obase_bicut(0+time_shift_, timestep+time_shift_, logic_grid_, f);
	cout << "no boundary kernel version " << endl ;
#if 1
    printf("shorter_duo_sim_obase_bicut!\n");
    algor.shorter_duo_sim_obase_bicut(0+time_shift_, timestep+time_shift_, logic_grid_, f);
#else
    printf("stevenj!\n");
    algor.stevenj(0+time_shift_, timestep+time_shift_, logic_grid_, f);
#endif
    // algor.duo_sim_obase_bicut(0+time_shift_, timestep+time_shift_, logic_grid_, f);
#pragma isat marker M2_end
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
    /* this version uses 'f' to compute interior region, 
     * and 'bf' to compute boundary region
     */
    timestep_ = timestep;
    checkFlags();
#if BICUT
#if 0
    fprintf(stderr, "Call obase_bicut_boundary_P\n");
#pragma isat marker M2_begin
    algor.obase_bicut_boundary_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#pragma isat marker M2_end
#else
//    fprintf(stderr, "Call sim_obase_bicut_P\n");
#pragma isat marker M2_begin
#if 1
    algor.set_time_step(timestep_);
	algor.set_time_shift(time_shift_) ;
	cout << "time_shift_ " << time_shift_ << " timestep " << timestep << endl ;
	grid_info <N_RANK> grid = logic_grid_ ;
	for (int i = N_RANK - 1 ; i >= 0 ; i--)
	{
		cout << "grid.dx0[i] " << grid.dx0[i] << endl ;
		cout << "grid.dx1[i] " << grid.dx1[i] << endl ;
		cout << " x0 [" << i << "] " << grid.x0 [i] 
			<< " x1 [" << i << "] " << grid.x1 [i] 
			<< " x2 [" << i << "] " << grid.x0[i] + grid.dx0[i] * timestep_
			<< " x3 [" << i << "] " << grid.x1[i] + grid.dx1[i] * timestep_
			<< endl ; 
		cout << "slope_ [i] " << slope_ [i] << endl ;
	}
	
#ifdef COUNT_PROJECTIONS
//    algor.compute_projections(0+time_shift_, timestep+time_shift_, logic_grid_) ;
	algor.set_thres_auto_tuning() ;
#endif
#ifdef COARSEN_BASE_CASE_WRT_BOTTOM_SIDE
	cout << "coarsen base case wrt bottom side " << endl ;
#else
	cout << "coarsen base case wrt shorter side " << endl ;
#endif
#ifdef KERNEL_SELECTION
	cout << "kernel selection" << endl ;
#endif
#ifdef GENEITY_TEST
	cout << "geneity testing" << endl ;
#endif
#ifdef DEFAULT_SPACE_CUT
	cout << "default space cut " << endl ;
#else
	cout << "modified space cut " << endl ;
#endif
#ifdef AUTO_TUNE
	algor.set_thres_auto_tuning() ;
	cout << "auto tune" << endl ;
#ifdef TIME_INVARIANCE_INTERIOR
	cout << "time invariance interior " << endl ;
#else
	cout << "space-time invariance interior " << endl ;
#endif

#ifdef TIME_INVARIANCE_BOUNDARY
	cout << "time invariance boundary " << endl ;
#else
	cout << "space-time invariance boundary " << endl ;
#endif

#ifdef FIXED_TIME_CUT
	cout << "fixed time cut" << endl ;
#else
	cout << "arbitrary time cut" << endl ;
#endif

#ifdef FIXED_SPACE_CUT
	cout << "fixed space cut" << endl ;
#else
	cout << "arbitrary space cut " << endl ;
#endif
#endif

#ifdef TRAP
	cout << "default time cut " << endl ;
	//algor.set_thres_auto_tuning() ;
#ifndef USE_PROJECTION
    algor.shorter_duo_sim_obase_bicut_p(time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#else
#ifdef KERNEL_SELECTION
	predicate <N_RANK> p (phys_grid_) ; 
	p.set_resolution(resolution_) ;
	p.set_output_file(outputFile_) ;
	pochoir_clone_array <N_RANK> c1(*arr_, p) ; //may not compile
	kernel_selection_trap<N_RANK> ks(algor, phys_grid_, 0) ;
	ks.set_clone_array(&c1) ;
	ks.do_default_space_time_cuts(time_shift_, timestep+time_shift_,
							logic_grid_, f, bf, p) ;
#elif defined GENEITY_TEST
	predicate <N_RANK> p (phys_grid_) ; 
	p.set_resolution(resolution_) ;
	//p.set_output_file(outputFile_) ;
	pochoir_clone_array <N_RANK> c1(*arr_, p) ; //may not compile
	heterogeneity<N_RANK> hg(algor, phys_grid_, 0) ;
	//geneity_problem<N_RANK> hg(algor, phys_grid_, 0) ;
	hg.set_clone_array(&c1) ;
	//hg.build_heterogeneity_dag(0, timestep, logic_grid_, p, 0) ;
	hg.do_default_space_time_cuts(time_shift_, timestep+time_shift_,
							logic_grid_, f, bf, p) ;
#ifndef NDEBUG
	cout << "calling print dag " << endl ;
	hg.print_dag() ;
	hg.print_heterogeneity() ;
#endif
#elif defined AUTO_TUNE
	//cout << "address of home cell " << &home_cell_ << endl ;
	auto_tune<N_RANK> at(algor, phys_grid_, 1, problem_name, timestep_,
						 arr_type_size_) ;
	struct timeval start, end;
	double compute_time = 0. ;
	gettimeofday(&start, 0);
	at.do_trap_space_time_cuts(time_shift_, timestep+time_shift_,
								logic_grid_, f, bf, arr_) ;
	gettimeofday(&end, 0);
	compute_time = tdiff(&end, &start) ;
	//std::cout << "compute time :" << 1.0e3 * compute_time << "ms" << std::endl;
	//at.print_dag() ;
#endif

#endif
#else
	cout << "pow2 time cut " << endl ;
	//algor.set_thres_auto_tuning() ;
#ifndef USE_PROJECTION
	struct timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);
    algor.power_of_two_time_cut(time_shift_, timestep+time_shift_, logic_grid_, f, bf);
	clock_gettime(CLOCK_MONOTONIC, &end) ;
	cout << "compute time " << tdiff2(&end, &start) * 1e3 << "ms" << endl ;
#else
#ifdef KERNEL_SELECTION
	predicate <N_RANK> p (phys_grid_) ; 
	p.set_resolution(resolution_) ;
	p.set_output_file(outputFile_) ;
	pochoir_clone_array <N_RANK> c1(*arr_, p) ; //may not compile
    kernel_selection_sawzoid<N_RANK> ks(algor, phys_grid_, 0) ;
    ks.set_clone_array(&c1) ;
    ks.do_power_of_two_time_cut(time_shift_, timestep+time_shift_,
                            logic_grid_, f, bf, p) ;
#elif defined GENEITY_TEST
	algor.set_thres_auto_tuning() ;
	predicate <N_RANK> p (phys_grid_) ; 
	p.set_resolution(resolution_) ;
	p.set_output_file(outputFile_) ;
	pochoir_clone_array <N_RANK> c1(*arr_, p) ; //may not compile
	heterogeneity<N_RANK> hg(algor, phys_grid_, 1) ;
	hg.set_clone_array(&c1) ;
	struct timeval start, end;
	double compute_time = 0. ;
	//hg.build_heterogeneity_dag_modified(0, timestep, logic_grid_, p, 0) ;
	gettimeofday(&start, 0);
	hg.do_power_of_two_time_cut(time_shift_, timestep+time_shift_,
								logic_grid_, f, bf, p) ;
	gettimeofday(&end, 0);
	compute_time = tdiff(&end, &start) ;
	//std::cout << "compute time :" << 1.0e3 * compute_time << "ms" << std::endl;
#ifndef NDEBUG
	//cout << "calling print dag " << endl ;
	//hg.print_dag() ;
	//hg.print_heterogeneity() ;
#endif
#elif defined AUTO_TUNE
	//cout << "address of home cell " << &home_cell_ << endl ;
	auto_tune<N_RANK> at(algor, phys_grid_, 1, problem_name, timestep_,
						arr_type_size_) ;
	struct timeval start, end;
	double compute_time = 0. ;
	gettimeofday(&start, 0);
	at.do_power_of_two_time_cut(time_shift_, timestep+time_shift_,
								logic_grid_, f, bf, arr_) ;
	gettimeofday(&end, 0);
	compute_time = tdiff(&end, &start) ;
	//std::cout << "compute time :" << 1.0e3 * compute_time << "ms" << std::endl;

	//at.print_dag() ;
#endif

#endif
#endif
#else
    printf("stevenj_p!\n");
    algor.stevenj_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#endif
#pragma isat marker M2_end
#if STAT
    for (int i = 1; i < SUPPORT_RANK; ++i) {
        fprintf(stderr, "sim_count_cut[%d] = %ld\n", i, algor.sim_count_cut[i].get_value());
    }
#endif
#endif
#else
#pragma isat marker M2_begin
    algor.obase_boundary_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#pragma isat marker M2_end
#endif
}

#endif
