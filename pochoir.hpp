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
#include "pochoir_kernel.hpp"
#include "pochoir_dloader.hpp"
#include "pochoir_walk_recursive.hpp"
#include "pochoir_walk_loops.hpp"
#include "rbq/rbq_2d_and.hpp"
#include "rbq/rbq_2d_or.hpp"

/* assuming there won't be more than 10 Pochoir_Array in one Pochoir object! */
// #define ARRAY_SIZE 10

template <int N_RANK>
class Pochoir {
    private:
        int slope_[N_RANK];
        Grid_Info<N_RANK> logic_grid_;
        Grid_Info<N_RANK> phys_grid_;
        int time_shift_;  //depth of the stencil
        int toggle_;
        int timestep_;	//# of time steps
        bool regArrayFlag_, regLogicDomainFlag_, regPhysDomainFlag_, regShapeFlag_;
        void checkFlag(bool flag, char const * str);
        void checkFlags(void);
        template <typename T_Array>
        void getPhysDomainFromArray(T_Array & arr);
        template <typename T_Array>
        void cmpPhysDomainFromArray(T_Array & arr);
        void Register_Shape(Pochoir_Shape<N_RANK> * shape, int N_SIZE);
        Vector_Info< Pochoir_Shape<N_RANK> > shape_;
        int num_arr_;
        int arr_type_size_;
        int lcm_unroll_;
        int order_num_;
        double pochoir_time_;
        Pochoir_Mode pmode_;
		
		//pointer to rbq
		rbq_2d<POCHOIR_WORD_TYPE> * rbq ;
        Vector_Info< Pochoir_Tile_Kernel<N_RANK> * > reg_tile_kernels_;
        Vector_Info< Pochoir_Combined_Obase_Kernel<N_RANK> * > reg_obase_kernels_;
        Vector_Info< Pochoir_Guard<N_RANK> * > reg_guards_;
        Vector_Info< Pochoir_Guard<N_RANK> * > reg_obase_guards_;

        /* Private Register Kernel Function */
        template <typename I>
        void reg_tile_dim(Pochoir_Tile_Kernel<N_RANK> & pt, int dim, I i); 
        template <typename I, typename ... IS>
        void reg_tile_dim(Pochoir_Tile_Kernel<N_RANK> & pt, int dim, I i, IS ... is);

    public:
     /* Register_Tile_Kernels, by default, are inclusive */
    void Register_Tile_Kernels(Pochoir_Guard<N_RANK> & g, Pochoir_Kernel<N_RANK> & k);
    template <int N_SIZE1>
    void Register_Tile_Kernels(Pochoir_Guard<N_RANK> & g, Pochoir_Kernel<N_RANK> (& tile)[N_SIZE1]);
    template <int N_SIZE1, int N_SIZE2>
    void Register_Tile_Kernels(Pochoir_Guard<N_RANK> & g, Pochoir_Kernel<N_RANK> (& tile)[N_SIZE1][N_SIZE2]);
    template <int N_SIZE1, int N_SIZE2, int N_SIZE3>
    void Register_Tile_Kernels(Pochoir_Guard<N_RANK> & g, Pochoir_Kernel<N_RANK> (& tile)[N_SIZE1][N_SIZE2][N_SIZE3]);

    /* Register_Tile_Obase_Kernels are exclusive */
    void Register_Tile_Obase_Kernels(Pochoir_Guard<N_RANK> & g, int unroll, Pochoir_Base_Kernel<N_RANK> & k, Pochoir_Base_Kernel<N_RANK> & bk);
    void Register_Tile_Obase_Kernels(Pochoir_Guard<N_RANK> & g, int unroll, Pochoir_Base_Kernel<N_RANK> & k, Pochoir_Base_Kernel<N_RANK> & cond_k, Pochoir_Base_Kernel<N_RANK> & bk, Pochoir_Base_Kernel<N_RANK> & cond_bk); 
    // get slope(s)
    int slope(int const _idx) { return slope_[_idx]; }
    Pochoir() {
        for (int i = 0; i < N_RANK; ++i) {
            slope_[i] = 0;
            logic_grid_.x0[i] = logic_grid_.x1[i] = logic_grid_.dx0[i] = logic_grid_.dx1[i] = 0;
            phys_grid_.x0[i] = phys_grid_.x1[i] = phys_grid_.dx0[i] = phys_grid_.dx1[i] = 0;
        }
        time_shift_ = 0; toggle_ = 0;
        timestep_ = 0;
        regArrayFlag_ = regLogicDomainFlag_ = regPhysDomainFlag_ = regShapeFlag_ = false;
        num_arr_ = 0;
        arr_type_size_ = 0;
        lcm_unroll_ = 1;
        pmode_ = Pochoir_Null;
        order_num_ = 0;
        pochoir_time_ = INF;
    }
    ~Pochoir() {
        regArrayFlag_ = regLogicDomainFlag_ = regPhysDomainFlag_ = regShapeFlag_ = false;
        size_t l_size = reg_tile_kernels_.size();
        for (int i = 0; i < l_size; ++i) {
            auto k = reg_tile_kernels_[i];
            del_ele (k);
        }
        l_size = reg_obase_kernels_.size();
        for (int i = 0; i < l_size; ++i) {
            auto k = reg_obase_kernels_[i];
            del_ele (k);
        }
        l_size = reg_guards_.size();
        for (int i = 0; i < l_size; ++i) {
            auto g = reg_guards_[i];
            del_ele (g);
        }
        l_size = reg_obase_guards_.size();
        for (int i = 0; i < l_size; ++i) {
            auto g = reg_obase_guards_[i];
            del_ele (g);
        }
    }

    /* currently, we just compute the slope[] out of the shape[] */
    /* We get the Grid_Info out of arrayInUse */
    template <typename T>
    void Register_Array(Pochoir_Array<T, N_RANK> & a);
    template <typename T, typename ... TS>
    void Register_Array(Pochoir_Array<T, N_RANK> & a, Pochoir_Array<TS, N_RANK> ... as);

    /* We should still keep the Register_Domain for zero-padding!!! */
    template <typename D>
    void Register_Domain(D const & d);
    template <typename D, typename ... DS>
    void Register_Domain(D const & d, DS ... ds);
    /* register boundary value function with corresponding Pochoir_Array object directly */
    template <typename T_Array, typename T_RET>
    void registerBoundaryFn(T_Array & arr, T_RET (*_bv)(T_Array &, int, int, int)) {
        arr.Register_Boundary(_bv);
        Register_Array(arr);
    } 
    Grid_Info<N_RANK> & get_phys_grid(void);

    /* Executable Spec */
    /* to remove: Run(timestep) */
    //void Run(int timestep);
    //void Run_Obase(int timestep);
    Pochoir_Plan<N_RANK> & Gen_Plan(int timepstep);
    template <typename T_Array>
    Pochoir_Plan<N_RANK> & Gen_Plan_Obase(int timepstep, const char * pochoir_mode, const char * fname, T_Array & a);
    template <typename T_Array>
    Pochoir_Plan<N_RANK> & Gen_Plan_Obase_internal(int timepstep, const char * pochoir_mode, const char * fname, T_Array & a);
    template <typename T_Array>
    Pochoir_Plan<N_RANK> & generate_clones(int timepstep, const char * pochoir_mode, const char * fname, T_Array & a);
    template <typename T_Array, typename ... T_ArrayS>
    Pochoir_Plan<N_RANK> & Gen_Plan_Obase(int timepstep, const char * pochoir_mode, const char * fname, T_Array & a, T_ArrayS ... as);
    Pochoir_Plan<N_RANK> & Load_Plan(const char * file_name);
    template <typename T_Array>
    Pochoir_Plan<N_RANK> & Load_Plan_Obase(const char * file_name, T_Array & a);
    template <typename T_Array, typename ... T_ArrayS>
    Pochoir_Plan<N_RANK> & Load_Plan_Obase(const char * file_name, T_Array & a, T_ArrayS ... as);
    void Store_Plan(const char * file_name, Pochoir_Plan<N_RANK> & _plan);
    void Destroy_Plan(Pochoir_Plan<N_RANK> & _plan); 
	
	//new run functions that use a plan
    void Run(Pochoir_Plan<N_RANK> & _plan);
    void Run_Obase_All_Cond_with_plan(Pochoir_Plan<N_RANK> & _plan);
    void Run_Obase_Unroll_T_with_plan(Pochoir_Plan<N_RANK> & _plan);

	//wrapper functions to switch between the new functions that use plan
	//and the old functions that do not use a plan
    void Run_Obase_All_Cond(Pochoir_Plan<N_RANK> & _plan);
    void Run_Obase_Unroll_T(Pochoir_Plan<N_RANK> & _plan);

	//run functions that do not use a plan
	void Run_Obase_Unroll_T_without_plan() ;
};

template <int N_RANK> template <typename I>
void Pochoir<N_RANK>::reg_tile_dim(Pochoir_Tile_Kernel<N_RANK> & pt, int dim, I i) {
    pt.size_[0] = i; 
    pt.pointer_[0] = 0;
    pt.stride_[0] = 1;
    for (int i = 1; i < dim; ++i) {
        pt.stride_[i] = pt.stride_[i - 1] * pt.size_[i - 1];
    }
    int l_total_size_ = 1;
    for (int i = 0; i < dim; ++i) {
        l_total_size_ *= pt.size_[i];
    }
    pt.total_size_ = l_total_size_;
    return;
}

template <int N_RANK> template <typename I, typename ... IS>
void Pochoir<N_RANK>::reg_tile_dim(Pochoir_Tile_Kernel<N_RANK> & pt, int dim, I i, IS ... is) {
    int l_dim = sizeof ... (IS);
    pt.size_[l_dim] = i;
    pt.pointer_[l_dim] = 0;
    reg_tile_dim(pt, dim, is ...);
    return;
}

template <int N_RANK>
void Pochoir<N_RANK>::Register_Tile_Kernels(Pochoir_Guard<N_RANK> & g, Pochoir_Kernel<N_RANK> & k) {
    typedef Pochoir_Kernel<N_RANK> T_Kernel;
    pmode_ = Pochoir_Tile;
    Pochoir_Guard<N_RANK> * l_pig = new Pochoir_Guard<N_RANK>();
    *l_pig = g;
    Pochoir_Tile_Kernel<N_RANK> * l_pit = new Pochoir_Tile_Kernel<N_RANK>();
    l_pit->kernel_ = new T_Kernel[1];
    l_pit->kernel_[0] = k;
    Register_Shape(k.Get_Shape(), k.Get_Shape_Size());
    reg_tile_dim(*l_pit, 1, 1);

    lcm_unroll_ = lcm(lcm_unroll_, 1);
    reg_guards_.push_back(l_pig);
    reg_tile_kernels_.push_back(l_pit);
    return;
}

template <int N_RANK> template <int N_SIZE1>
void Pochoir<N_RANK>::Register_Tile_Kernels(Pochoir_Guard<N_RANK> & g, Pochoir_Kernel<N_RANK> (& tile)[N_SIZE1]) {
    typedef Pochoir_Kernel<N_RANK> T_Kernel;
    pmode_ = Pochoir_Tile;
    Pochoir_Guard<N_RANK> * l_pig = new Pochoir_Guard<N_RANK>();
    *l_pig = g;
    Pochoir_Tile_Kernel<N_RANK> * l_pit = new Pochoir_Tile_Kernel<N_RANK>();
    l_pit->kernel_ = new T_Kernel[N_SIZE1];
    for (int i = 0; i < N_SIZE1; ++i) {
        l_pit->kernel_[i] = tile[i];
        Register_Shape(tile[i].Get_Shape(), tile[i].Get_Shape_Size());
    }
    reg_tile_dim(*l_pit, 1, N_SIZE1);

    lcm_unroll_ = lcm(lcm_unroll_, N_SIZE1);
    reg_guards_.push_back(l_pig);
    reg_tile_kernels_.push_back(l_pit);
    return;
}

template <int N_RANK> template <int N_SIZE1, int N_SIZE2>
void Pochoir<N_RANK>::Register_Tile_Kernels(Pochoir_Guard<N_RANK> & g, Pochoir_Kernel<N_RANK> (& tile)[N_SIZE1][N_SIZE2]) {
    typedef Pochoir_Kernel<N_RANK> T_Kernel;
    pmode_ = Pochoir_Tile;
    Pochoir_Guard<N_RANK> * l_pig = new Pochoir_Guard<N_RANK>();
    *l_pig = g;
    Pochoir_Tile_Kernel<N_RANK> * l_pit = new Pochoir_Tile_Kernel<N_RANK>();
    l_pit->kernel_ = new T_Kernel[N_SIZE1 * N_SIZE2];
    for (int i = 0; i < N_SIZE1; ++i) {
        for (int j = 0; j < N_SIZE2; ++j) {
            l_pit->kernel_[i * N_SIZE2 + j] = tile[i][j];
            Register_Shape(tile[i][j].Get_Shape(), tile[i][j].Get_Shape_Size());
        }
    }
    reg_tile_dim(*l_pit, 2, N_SIZE1, N_SIZE2);

    lcm_unroll_ = lcm(lcm_unroll_, N_SIZE1);
    reg_guards_.push_back(l_pig);
    reg_tile_kernels_.push_back(l_pit);
    return;
}

template <int N_RANK> template <int N_SIZE1, int N_SIZE2, int N_SIZE3>
void Pochoir<N_RANK>::Register_Tile_Kernels(Pochoir_Guard<N_RANK> & g, Pochoir_Kernel<N_RANK> (& tile)[N_SIZE1][N_SIZE2][N_SIZE3]) {
    typedef Pochoir_Kernel<N_RANK> T_Kernel;
    pmode_ = Pochoir_Tile;
    Pochoir_Guard<N_RANK> * l_pig = new Pochoir_Guard<N_RANK>();
    *l_pig = g;
    Pochoir_Tile_Kernel<N_RANK> * l_pit = new Pochoir_Tile_Kernel<N_RANK>();
    l_pit->kernel_ = new T_Kernel[N_SIZE1 * N_SIZE2 * N_SIZE3];
    for (int i = 0; i < N_SIZE1; ++i) {
        for (int j = 0; j < N_SIZE2; ++j) {
    for (int k = 0; k < N_SIZE3; ++k) {
        l_pit->kernel_[i * N_SIZE2 * N_SIZE3 + j * N_SIZE3 + k] = tile[i][j][k];
        Register_Shape(tile[i][j][k].Get_Shape(), tile[i][j][k].Get_Shape_Size());
    }
        }
    }
    reg_tile_dim(*l_pit, 3, N_SIZE1, N_SIZE2, N_SIZE3);

    lcm_unroll_ = lcm(lcm_unroll_, N_SIZE1);
    reg_guards_.push_back(l_pig);
    reg_tile_kernels_.push_back(l_pit);
    return;
}

template <int N_RANK> 
void Pochoir<N_RANK>::Register_Tile_Obase_Kernels(Pochoir_Guard<N_RANK> & g, int unroll, Pochoir_Base_Kernel<N_RANK> & k, Pochoir_Base_Kernel<N_RANK> & bk) {
    // typedef Pochoir_Base_Kernel<N_RANK> T_Kernel;
    pmode_ = Pochoir_Obase_Tile;
    Pochoir_Guard<N_RANK> * l_opg = new Pochoir_Guard<N_RANK>();
    *l_opg = g;
    Pochoir_Combined_Obase_Kernel<N_RANK> * l_opk = new Pochoir_Combined_Obase_Kernel<N_RANK>();
    l_opk->unroll_ = unroll;
    l_opk->kernel_ = &k;
    Register_Shape(k.Get_Shape(), k.Get_Shape_Size());
    l_opk->bkernel_ = &bk;
    Register_Shape(bk.Get_Shape(), bk.Get_Shape_Size());
    l_opk->cond_kernel_ = NULL;
    l_opk->cond_bkernel_ = NULL;
    lcm_unroll_ = lcm(lcm_unroll_, unroll);

    reg_obase_guards_.push_back(l_opg);
    reg_obase_kernels_.push_back(l_opk);
    return;
}

template <int N_RANK> 
void Pochoir<N_RANK>::Register_Tile_Obase_Kernels(Pochoir_Guard<N_RANK> & g, int unroll, Pochoir_Base_Kernel<N_RANK> & k, Pochoir_Base_Kernel<N_RANK> & cond_k, Pochoir_Base_Kernel<N_RANK> & bk, Pochoir_Base_Kernel<N_RANK> & cond_bk) {
    // typedef Pochoir_Base_Kernel<N_RANK> T_Kernel;
    pmode_ = Pochoir_Obase_Tile;
    Pochoir_Guard<N_RANK> * l_opg = new Pochoir_Guard<N_RANK>();
    *l_opg = g;
    Pochoir_Combined_Obase_Kernel<N_RANK> * l_opk = new Pochoir_Combined_Obase_Kernel<N_RANK>();
    l_opk->unroll_ = unroll;
    l_opk->kernel_ = &k;
    Register_Shape(k.Get_Shape(), k.Get_Shape_Size());
    l_opk->cond_kernel_ = &cond_k;
    Register_Shape(cond_k.Get_Shape(), cond_k.Get_Shape_Size());
    l_opk->bkernel_ = &bk;
    Register_Shape(bk.Get_Shape(), bk.Get_Shape_Size());
    l_opk->cond_bkernel_ = &cond_bk;
    Register_Shape(cond_bk.Get_Shape(), cond_bk.Get_Shape_Size());
    lcm_unroll_ = lcm(lcm_unroll_, unroll);

    reg_obase_guards_.push_back(l_opg);
    reg_obase_kernels_.push_back(l_opk);
    return;
}

template <int N_RANK>
void Pochoir<N_RANK>::checkFlag(bool flag, char const * str) {
    if (!flag) {
        ERROR_ARGS("\nPochoir registration error:\n" "You forgot to register %s.\n", str);
    }
}

template <int N_RANK>
void Pochoir<N_RANK>::checkFlags(void) {
    checkFlag(regArrayFlag_, "Pochoir array");
    checkFlag(regLogicDomainFlag_, "Logic Domain");
    checkFlag(regPhysDomainFlag_, "Physical Domain");
    checkFlag(regShapeFlag_, "Shape");
    return;
}

template <int N_RANK> template <typename T_Array> 
void Pochoir<N_RANK>::getPhysDomainFromArray(T_Array & arr) {
    /* get the physical grid */
    for (int i = 0; i < N_RANK; ++i) {
        phys_grid_.x0[i] = 0; phys_grid_.x1[i] = arr.size(i);
        /* if logic domain is not set, let's set it the same as physical grid */
        if (!regLogicDomainFlag_) {
            logic_grid_.x0[i] = 0; logic_grid_.x1[i] = arr.size(i);
        }
    }

    regPhysDomainFlag_ = true;
    regLogicDomainFlag_ = true;
}

template <int N_RANK> template <typename T_Array> 
void Pochoir<N_RANK>::cmpPhysDomainFromArray(T_Array & arr) {
    /* check the consistency of all engaged Pochoir_Array */
    for (int j = 0; j < N_RANK; ++j) {
        if (arr.size(j) != phys_grid_.x1[j]) {
            ERROR("Pochoir array size mismatch error:\n" 
                    "Registered Pochoir arrays have different sizes!\n");
        }
    }
}

template <int N_RANK> template <typename T>
void Pochoir<N_RANK>::Register_Array(Pochoir_Array<T, N_RANK> & a) {
    if (!regShapeFlag_) {
        ERROR("Please register Shape before register Array!\n");
    }

    if (num_arr_ == 0) {
        arr_type_size_ = sizeof(T);
        // arr_type_size_ = sizeof(double);
        LOG_ARGS(0, "arr_type_size = %d\n", arr_type_size_);
        ++num_arr_;
    } 
    if (!regPhysDomainFlag_) {
        getPhysDomainFromArray(a);
    } else {
        cmpPhysDomainFromArray(a);
    }
    a.Register_Shape(shape_);
    regArrayFlag_ = true;
}

template <int N_RANK> template <typename T, typename ... TS>
void Pochoir<N_RANK>::Register_Array(Pochoir_Array<T, N_RANK> & a, Pochoir_Array<TS, N_RANK> ... as) {
    if (!regShapeFlag_) {
        ERROR("Please register Shape before register Array!\n");
    }

    if (num_arr_ == 0) {
        arr_type_size_ = sizeof(T);
        LOG_ARGS(0, "arr_type_size = %d\n", arr_type_size_);
        ++num_arr_;
    } 
    if (!regPhysDomainFlag_) {
        getPhysDomainFromArray(a);
    } else {
        cmpPhysDomainFromArray(a);
    }
    a.Register_Shape(shape_);
    Register_Array(as ...);
    regArrayFlag_ = true;
}

template <int N_RANK> 
void Pochoir<N_RANK>::Register_Shape(Pochoir_Shape<N_RANK> * shape, int N_SIZE) {
    /* we extract time_shift_, toggle_, slope_ out of Shape */
    /* we make a copy of the Pochoir_Shape since it can be accumulated
     * across multiple kernels
     */
    regShapeFlag_ = true;
    /* currently we just get the slope_[]  and toggle_ out of the shape[]  */
    int l_min_time_shift=0, l_max_time_shift=0, depth=0;
    for (int i = 0; i < N_SIZE; ++i)  {
        if (shape[i] .shift[0] < l_min_time_shift)
            l_min_time_shift = shape[i].shift[0];
        if (shape[i].shift[0] > l_max_time_shift)
            l_max_time_shift = shape[i].shift[0];
    }
    depth = l_max_time_shift - l_min_time_shift;
    time_shift_ = max(time_shift_, 0 - l_min_time_shift);
    toggle_ = max(toggle_, depth + 1);
    for (int i = 0; i < N_SIZE; ++i) {
        shape_.push_back_unique(shape[i]);
    }
    for (int i = 0; i < N_SIZE; ++i) {
        for (int r = 0; r < N_RANK; ++r) {
            slope_[r] = max(slope_[r], abs((int)ceil((float)shape[i].shift[r+1]/(l_max_time_shift - shape[i].shift[0]))));
        }
    }
    LOG_ARGS(0, "toggle = %d\n", toggle_);
    for (int r = 0; r < N_RANK; ++r) {
        LOG_ARGS(0, "slope[%d] = %d, ", r, slope_[r]);
    }
    LOG(0, "\n");
}

template <int N_RANK> template <typename D>
void Pochoir<N_RANK>::Register_Domain(D const & d) {
    logic_grid_.x0[0] = d.first();
    logic_grid_.x1[0] = d.first() + d.size();
    regLogicDomainFlag_ = true;
}

template <int N_RANK> template <typename D, typename ... DS>
void Pochoir<N_RANK>::Register_Domain(D const & d, DS ... ds) {
    int l_pointer = sizeof...(DS);
    logic_grid_.x0[l_pointer] = d.first();
    logic_grid_.x1[l_pointer] = d.first() + d.size();
}

template <int N_RANK> 
Grid_Info<N_RANK> & Pochoir<N_RANK>::get_phys_grid(void) {
    return phys_grid_;
}

template <int N_RANK> 
Pochoir_Plan<N_RANK> & Pochoir<N_RANK>::Gen_Plan(int timestep) {
    Pochoir_Plan<N_RANK> * l_plan = new Pochoir_Plan<N_RANK>();
    int l_sz_base_data, l_sz_sync_data;
    Spawn_Tree<N_RANK> * l_tree;
    Node_Info<N_RANK> * l_root, * l_internal;

    timestep_ = timestep;
    checkFlags();
    l_tree = new Spawn_Tree<N_RANK>();
    l_tree->set_add_empty_region(true);
    l_root = l_tree->get_root();
    Homogeneity white_clone(0);
    l_internal = new Node_Info<N_RANK>(0 + time_shift_, timestep + time_shift_, logic_grid_, white_clone);
    l_tree->add_node(l_root, l_internal, IS_SPAWN);
    l_sz_base_data = 2;
    l_sz_sync_data = 0;
    l_plan->alloc_base_data(l_sz_base_data);
    l_plan->alloc_sync_data(l_sz_sync_data);
    int l_tree_size_begin = l_tree->size();
    LOG_ARGS(0, "sz_base_data = %d, sz_sync_data = %d\n", l_sz_base_data, l_sz_sync_data);
    LOG_ARGS(0, "tree size = %d\n", l_tree_size_begin);
    /* after remove all nodes, the only remaining node will be the 'root' */
    while (l_tree_size_begin > 1) {
        l_tree->dfs_until_sync(l_root->left, (*(l_plan->base_data_)));
        int l_tree_size_end = l_tree->size();
        if (l_tree_size_begin - l_tree_size_end > 0) {
            l_plan->sync_data_->push_back(l_tree_size_begin - l_tree_size_end);
        }
        l_tree->dfs_rm_sync(l_root->left);
        l_tree_size_begin = l_tree->size();
        LOG_ARGS(0, "tree size = %d\n", l_tree_size_begin);
    }
    l_plan->sync_data_->scan();
    l_plan->sync_data_->push_back(END_SYNC);
    l_plan->set_order_num(order_num_);
    ++order_num_;
    del_ele(l_tree);
    return (*l_plan);
}

template <int N_RANK> template <typename T_Array> 
Pochoir_Plan<N_RANK> & Pochoir<N_RANK>::generate_clones(int timestep, 
				const char * pochoir_mode, const char * fname, T_Array & a) 
{
    struct timeval l_start, l_end;
    double l_min_tdiff = INF;

    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    algor.set_time_shift(time_shift_);
    /* we unroll at least twice in out-most loop to promote performance */
    lcm_unroll_ = (lcm_unroll_ > 1) ? lcm_unroll_ : 2;
    algor.set_unroll(lcm_unroll_);
#ifdef DEBUG
    printf("lcm_unroll_ = %d\n", lcm_unroll_);
#endif
    /* set individual unroll factor from opgk_ */
    if (pmode_ == Pochoir_Tile) {
        algor.set_pts(reg_guards_);
    } else {
        ERROR("Something is wrong in generate_clones(Timestep)!\n");
    }
    timestep_ = timestep;
    checkFlags();
	
	//initialize the rbq pointer
	//rbq = new rbq_2d<POCHOIR_WORD_TYPE> ()
	algor.find_homogeneity(logic_grid_) ;

    char color_vector_fname[FNAME_LENGTH], kernel_info_fname[FNAME_LENGTH];
    sprintf(color_vector_fname, "%s_%d_color.dat", fname, order_num_);
    sprintf(kernel_info_fname, "%s_kernel_info.cpp", fname);

    Vector_Info< Homogeneity > & l_color_vector = algor.get_color_vector();
    Homogeneity * white_clone = NULL;

    if (l_color_vector.size() > 0) {
        white_clone = new Homogeneity(l_color_vector[0].size());
    } else {
        white_clone = new Homogeneity(0);
    }
    /* sort the color vector according to the member 'measure_' */
    l_color_vector.sort();
    l_color_vector.set_size(55);
    /* make the 'rec_level' of white_clone to be 0, so it won't be eliminated out */
    l_color_vector.push_back_unique(*white_clone, 0);
    std::ofstream os_color_vector;
    os_color_vector.open(color_vector_fname, ofstream::out | ofstream::trunc);
    if (os_color_vector.is_open()) {
        os_color_vector << l_color_vector << std::endl;
    } else {
        ERROR("os_color_vector is NOT open! exit!\n");
    }
    os_color_vector.close();

    l_min_tdiff = INF;
    gettimeofday(&l_start, 0);
    char cmd[200];
    sprintf(cmd, "./genstencils -order %d %s %s %s", order_num_, pochoir_mode, color_vector_fname, kernel_info_fname);
    LOG_ARGS(INF, "%s\n", cmd);
    int ret = system(cmd);
    if (ret == -1) {
        ERROR("system() call to genstencils failed!\n");
    }
    LOG(INF, "./genstencils exits!\n");
    gettimeofday(&l_end, 0);
    l_min_tdiff = min (l_min_tdiff, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(INF, "Kernel Generation (genstencils) time : %.6f milliseconds\n", l_min_tdiff);

    l_min_tdiff = INF;
    gettimeofday(&l_start, 0);
    char gen_kernel_fname [strlen(fname) + 20];
    sprintf(gen_kernel_fname, "./%s_%d_gen_kernel", fname, order_num_);
    LOG_ARGS(0, "gen_kernel_fname = %s\n", gen_kernel_fname);
    char cpp_filename[strlen(gen_kernel_fname) + 10], so_filename[strlen(gen_kernel_fname) + 10];
    sprintf(cpp_filename, "%s.cpp", gen_kernel_fname);
    sprintf(so_filename, "%s.so", gen_kernel_fname);
#if 0 
    /* This branch is for debugging purpose */
    sprintf(cmd, "icpc -o %s -shared -nostartfiles -fPIC -O0 -g -std=c++0x -I${POCHOIR_LIB_PATH} %s", so_filename, cpp_filename);
#else
    /* This branch is for best performance */
    // sprintf(cmd, "icpc -o %s -shared -nostartfiles -fPIC -O2 -Wall -Werror -unroll-agressive -funroll-loops -xHost -fno-alias -fno-fnalias -fp-model precise -std=c++0x -I${POCHOIR_LIB_PATH} %s", so_filename, cpp_filename);
    sprintf(cmd, "icpc -o %s -shared -nostartfiles -fPIC -O2 -Wall -Werror -unroll-agressive -funroll-loops -xHost -fno-alias -fno-fnalias -fp-model precise -std=c++0x -I${POCHOIR_LIB_PATH} %s", so_filename, cpp_filename);
    // sprintf(cmd, "g++-cilk -o %s -shared -nostartfiles -fPIC -O2 -std=c++0x -I${POCHOIR_LIB_PATH} %s", so_filename, cpp_filename);
#endif

    LOG_ARGS(INF, "%s\n", cmd);
    ret = system(cmd);
    if (ret == -1) {
        ERROR("system() call failed!");
    }
    gettimeofday(&l_end, 0);
    l_min_tdiff = min (l_min_tdiff, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(INF, "Kernel Compilation (icpc) time : %.6f milliseconds\n", l_min_tdiff);

    l_min_tdiff = INF;
    gettimeofday(&l_start, 0);
    Pochoir_Plan<N_RANK> * l_plan = new Pochoir_Plan<N_RANK>();
    l_plan->set_order_num(order_num_);
    l_plan->set_fname(fname);
    l_plan->load_kernels(*this, a); 
    gettimeofday(&l_end, 0);
    l_min_tdiff = min (l_min_tdiff, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(INF, "Dynamic loading time : %.6f milliseconds\n", l_min_tdiff);

    ++order_num_;
    /* clean up the memory to avoid memory leak */
    del_ele(white_clone);
    return (*l_plan);
}

template <int N_RANK> template <typename T_Array> 
Pochoir_Plan<N_RANK> & Pochoir<N_RANK>::Gen_Plan_Obase(int timestep, const char * pochoir_mode, const char * fname, T_Array & a) 
{
	//return Gen_Plan_Obase_internal(timestep, pochoir_mode, fname, a) ;
	return generate_clones(timestep, pochoir_mode, fname, a) ;
}



template <int N_RANK> template <typename T_Array> 
Pochoir_Plan<N_RANK> & Pochoir<N_RANK>::Gen_Plan_Obase_internal(int timestep, const char * pochoir_mode, const char * fname, T_Array & a) {
    Pochoir_Plan<N_RANK> * l_plan = new Pochoir_Plan<N_RANK>();
    int l_sz_base_data, l_sz_sync_data;
    Spawn_Tree<N_RANK> * l_tree;
    Node_Info<N_RANK> * l_root;
    struct timeval l_start, l_end;
    double l_min_tdiff = INF;

    gettimeofday(&l_start, 0);
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    algor.set_time_shift(time_shift_);
    /* we unroll at least twice in out-most loop to promote performance */
    lcm_unroll_ = (lcm_unroll_ > 1) ? lcm_unroll_ : 2;
    algor.set_unroll(lcm_unroll_);
#ifdef DEBUG
    printf("lcm_unroll_ = %d\n", lcm_unroll_);
#endif
    /* set individual unroll factor from opgk_ */
    if (pmode_ == Pochoir_Tile) {
        algor.set_pts(reg_guards_);
    } else {
        ERROR("Something is wrong in Gen_Plan_Obase(Timestep)!\n");
    }
    timestep_ = timestep;
    checkFlags();
    l_tree = new Spawn_Tree<N_RANK>();
    l_tree->set_add_empty_region(false);
    l_root = l_tree->get_root();
    algor.set_tree(l_tree);
	
    algor.gen_plan_cut_p(l_root, 0 + time_shift_, timestep + time_shift_, logic_grid_, 0);
	l_sz_base_data = algor.get_sz_base_data();
    l_sz_sync_data = max(1, algor.get_sz_sync_data());
    l_plan->alloc_base_data(l_sz_base_data);
    l_plan->alloc_sync_data(l_sz_sync_data);
    int l_tree_size_begin = l_tree->size();
    LOG_ARGS(0, "sz_base_data = %d, sz_sync_data = %d\n", l_sz_base_data, l_sz_sync_data);
    LOG_ARGS(0, "tree size = %d\n", l_tree_size_begin);
    /* after remove all nodes, the only remaining node will be the 'root' */
    while (l_tree_size_begin > 1) {
        l_tree->dfs_until_sync(l_root->left, (*(l_plan->base_data_)));
        int l_tree_size_end = l_tree->size();
        if (l_tree_size_begin - l_tree_size_end > 0) {
            l_plan->sync_data_->push_back(l_tree_size_begin - l_tree_size_end);
        }
        l_tree->dfs_rm_sync(l_root->left);
        l_tree_size_begin = l_tree->size();
        LOG_ARGS(0, "tree size = %d\n", l_tree_size_begin);
    }
    l_plan->sync_data_->scan();
    l_plan->sync_data_->push_back(END_SYNC);
    l_plan->set_order_num(order_num_);
    l_plan->set_fname(fname);

    char color_vector_fname[FNAME_LENGTH], kernel_info_fname[FNAME_LENGTH];
    sprintf(color_vector_fname, "%s_%d_color.dat", fname, order_num_);
    sprintf(kernel_info_fname, "%s_kernel_info.cpp", fname);

    Vector_Info< Homogeneity > & l_color_vector = algor.get_color_vector();
    Homogeneity * white_clone = NULL;

    if (l_color_vector.size() > 0) {
        white_clone = new Homogeneity(l_color_vector[0].size());
        // white_clone = (Homogeneity *)calloc(1, sizeof(Homogeneity)); 
        // if (white_clone != NULL) 
        //     white_clone->set_size(l_color_vector[0].size());
    } else {
        white_clone = new Homogeneity(0);
    }
    /* sort the color vector according to the member 'measure_' */
    l_color_vector.sort();
    l_color_vector.set_size(55);
    /* make the 'rec_level' of white_clone to be 0, so it won't be eliminated out */
    l_color_vector.push_back_unique(*white_clone, 0);
    l_plan->change_region_n(l_color_vector);
    std::ofstream os_color_vector;
    os_color_vector.open(color_vector_fname, ofstream::out | ofstream::trunc);
    if (os_color_vector.is_open()) {
        os_color_vector << l_color_vector << std::endl;
    } else {
        ERROR("os_color_vector is NOT open! exit!\n");
    }
    os_color_vector.close();
    gettimeofday(&l_end, 0);
    l_min_tdiff = min (l_min_tdiff, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(INF, "Pochoir_Plan generation time : %.6f milliseconds\n", l_min_tdiff);

    l_min_tdiff = INF;
    gettimeofday(&l_start, 0);
    char cmd[200];
    sprintf(cmd, "./genstencils -order %d %s %s %s", order_num_, pochoir_mode, color_vector_fname, kernel_info_fname);
    LOG_ARGS(INF, "%s\n", cmd);
    int ret = system(cmd);
    if (ret == -1) {
        ERROR("system() call to genstencils failed!\n");
    }
    LOG(INF, "./genstencils exits!\n");
    gettimeofday(&l_end, 0);
    l_min_tdiff = min (l_min_tdiff, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(INF, "Kernel Generation (genstencils) time : %.6f milliseconds\n", l_min_tdiff);

    l_min_tdiff = INF;
    gettimeofday(&l_start, 0);
    char gen_kernel_fname [strlen(fname) + 20];
    sprintf(gen_kernel_fname, "./%s_%d_gen_kernel", fname, order_num_);
    LOG_ARGS(0, "gen_kernel_fname = %s\n", gen_kernel_fname);
    char cpp_filename[strlen(gen_kernel_fname) + 10], so_filename[strlen(gen_kernel_fname) + 10];
    sprintf(cpp_filename, "%s.cpp", gen_kernel_fname);
    sprintf(so_filename, "%s.so", gen_kernel_fname);
#if 0 
    /* This branch is for debugging purpose */
    sprintf(cmd, "icpc -o %s -shared -nostartfiles -fPIC -O0 -g -std=c++0x -I${POCHOIR_LIB_PATH} %s", so_filename, cpp_filename);
#else
    /* This branch is for best performance */
    // sprintf(cmd, "icpc -o %s -shared -nostartfiles -fPIC -O2 -Wall -Werror -unroll-agressive -funroll-loops -xHost -fno-alias -fno-fnalias -fp-model precise -std=c++0x -I${POCHOIR_LIB_PATH} %s", so_filename, cpp_filename);
    sprintf(cmd, "icpc -o %s -shared -nostartfiles -fPIC -O2 -Wall -Werror -unroll-agressive -funroll-loops -xHost -fno-alias -fno-fnalias -fp-model precise -std=c++0x -I${POCHOIR_LIB_PATH} %s", so_filename, cpp_filename);
    // sprintf(cmd, "g++-cilk -o %s -shared -nostartfiles -fPIC -O2 -std=c++0x -I${POCHOIR_LIB_PATH} %s", so_filename, cpp_filename);
#endif

    LOG_ARGS(INF, "%s\n", cmd);
    ret = system(cmd);
    if (ret == -1) {
        ERROR("system() call failed!");
    }
    gettimeofday(&l_end, 0);
    l_min_tdiff = min (l_min_tdiff, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(INF, "Kernel Compilation (icpc) time : %.6f milliseconds\n", l_min_tdiff);

    l_min_tdiff = INF;
    gettimeofday(&l_start, 0);
    l_plan->load_kernels(*this, a); 
    gettimeofday(&l_end, 0);
    l_min_tdiff = min (l_min_tdiff, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(INF, "Dynamic loading time : %.6f milliseconds\n", l_min_tdiff);

    ++order_num_;
    /* clean up the memory to avoid memory leak */
    del_ele(white_clone);
    del_ele(l_tree); 
    return (*l_plan);
}

template <int N_RANK> template <typename T_Array, typename ... T_ArrayS> 
Pochoir_Plan<N_RANK> & Pochoir<N_RANK>::Gen_Plan_Obase(int timestep, const char * pochoir_mode, const char * fname, T_Array & a, T_ArrayS ... as) {
    Pochoir_Plan<N_RANK> * l_plan = new Pochoir_Plan<N_RANK>();
    int l_sz_base_data, l_sz_sync_data;
    Spawn_Tree<N_RANK> * l_tree;
    Node_Info<N_RANK> * l_root;
    struct timeval l_start, l_end;
    double l_min_tdiff = INF;
    gettimeofday(&l_start, 0);
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    algor.set_time_shift(time_shift_);
    /* we unroll at least twice in out-most loop to promote performance */
    lcm_unroll_ = (lcm_unroll_ > 1) ? lcm_unroll_ : 2;
    algor.set_unroll(lcm_unroll_);
#ifdef DEBUG
    printf("lcm_unroll_ = %d\n", lcm_unroll_);
#endif
    /* set individual unroll factor from opgk_ */
    if (pmode_ == Pochoir_Tile) {
        // set color_region
        algor.set_pts(reg_guards_);
    } else {
        ERROR("Something is wrong in Gen_Plan_Obase(Timestep)!\n");
    }
    timestep_ = timestep;
    checkFlags();
    l_tree = new Spawn_Tree<N_RANK>();
    l_tree->set_add_empty_region(false);
    l_root = l_tree->get_root();
    algor.set_tree(l_tree);
    algor.gen_plan_cut_p(l_root, 0 + time_shift_, timestep + time_shift_, logic_grid_, 0);
    l_sz_base_data = algor.get_sz_base_data();
    l_sz_sync_data = max(1, algor.get_sz_sync_data());
    l_plan->alloc_base_data(l_sz_base_data);
    l_plan->alloc_sync_data(l_sz_sync_data);
    int l_tree_size_begin = l_tree->size();
    LOG_ARGS(0, "sz_base_data = %d, sz_sync_data = %d\n", l_sz_base_data, l_sz_sync_data);
    LOG_ARGS(0, "tree size = %d\n", l_tree_size_begin);
    /* after remove all nodes, the only remaining node will be the 'root' */
    while (l_tree_size_begin > 1) {
        l_tree->dfs_until_sync(l_root->left, (*(l_plan->base_data_)));
        int l_tree_size_end = l_tree->size();
        if (l_tree_size_begin - l_tree_size_end > 0) {
            l_plan->sync_data_->push_back(l_tree_size_begin - l_tree_size_end);
        }
        l_tree->dfs_rm_sync(l_root->left);
        l_tree_size_begin = l_tree->size();
        LOG_ARGS(0, "tree size = %d\n", l_tree_size_begin);
    }
    l_plan->sync_data_->scan();
    l_plan->sync_data_->push_back(END_SYNC);
    l_plan->set_order_num(order_num_);
    l_plan->set_fname(fname);

    char color_vector_fname[100], kernel_info_fname[100];
    sprintf(color_vector_fname, "%s_%d_color.dat\0", fname, order_num_);
    sprintf(kernel_info_fname, "%s_kernel_info.cpp\0", fname);

    Vector_Info< Homogeneity > & l_color_vector = algor.get_color_vector();
    Homogeneity * white_clone = NULL;

    if (l_color_vector.size() > 0)
        white_clone = new Homogeneity(l_color_vector[0].size());
    else
        white_clone = new Homogeneity(0);
    /* sort the color vector according to the member 'measure_' */
    l_color_vector.sort();
    l_color_vector.set_size(55);
    /* make the 'rec_level' of white_clone to be 0, so it won't be eliminated out */
    l_color_vector.push_back_unique(*white_clone, 0);
    l_plan->change_region_n(l_color_vector);
    std::ofstream os_color_vector;
    os_color_vector.open(color_vector_fname, ofstream::out | ofstream::trunc);
    if (os_color_vector.is_open()) {
        os_color_vector << l_color_vector << std::endl;
    } else {
        ERROR("os_color_vector is NOT open! exit!\n");
    }
    os_color_vector.close();

    gettimeofday(&l_end, 0);
    l_min_tdiff = min (l_min_tdiff, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(INF, "Pochoir_Plan generation time : %.6f milliseconds\n", l_min_tdiff);

    l_min_tdiff = INF;
    gettimeofday(&l_start, 0);
    char cmd[200];
    sprintf(cmd, "./genstencils -order %d %s %s %s\0", order_num_, pochoir_mode, color_vector_fname, kernel_info_fname);
    LOG_ARGS(INF, "%s\n", cmd);
    int ret = system(cmd);
    if (ret == -1) {
        ERROR("system() call to genstencils failed!\n");
    }
    LOG(INF, "./genstencils exits!\n");
    gettimeofday(&l_end, 0);
    l_min_tdiff = min (l_min_tdiff, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(INF, "Kernel Generation (genstencils) time : %.6f milliseconds\n", l_min_tdiff);

    l_min_tdiff = INF;
    gettimeofday(&l_start, 0);
    char gen_kernel_fname [strlen(fname) + 20];
    sprintf(gen_kernel_fname, "./%s_%d_gen_kernel", fname, order_num_);
    LOG_ARGS(0, "gen_kernel_fname = %s\n", gen_kernel_fname);
    char cpp_filename[strlen(gen_kernel_fname) + 10], so_filename[strlen(gen_kernel_fname) + 10];
    sprintf(cpp_filename, "%s.cpp", gen_kernel_fname);
    sprintf(so_filename, "%s.so", gen_kernel_fname);
#if 0
    sprintf(cmd, "icpc -o %s -shared -nostartfiles -fPIC -O0 -g -std=c++0x -I${POCHOIR_LIB_PATH} %s", so_filename, cpp_filename);
#else
    sprintf(cmd, "icpc -o %s -shared -nostartfiles -fPIC -O2 -Wall -Werror -unroll-agressive -funroll-loops -xHost -fno-alias -fno-fnalias -fp-model precise -std=c++0x -I${POCHOIR_LIB_PATH} %s\0", so_filename, cpp_filename);
    // sprintf(cmd, "g++-cilk -o %s -shared -nostartfiles -fPIC -O2 -std=c++0x -I${POCHOIR_LIB_PATH} %s", so_filename, cpp_filename);
#endif

    LOG_ARGS(INF, "%s\n", cmd);
    ret = system(cmd);
    if (ret == -1) {
        ERROR("system() call failed!");
    }
    gettimeofday(&l_end, 0);
    l_min_tdiff = min (l_min_tdiff, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(INF, "Kernel Compilation (icpc) time : %.6f milliseconds\n", l_min_tdiff);

    l_min_tdiff = INF;
    gettimeofday(&l_start, 0);
    l_plan->load_kernels(*this, a, as ...); 
    gettimeofday(&l_end, 0);
    l_min_tdiff = min (l_min_tdiff, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(INF, "Dynamic loading time : %.6f milliseconds\n", l_min_tdiff);

    ++order_num_;
    /* clean up the temporary local memory */
    del_ele(white_clone);
    del_ele(l_tree);
    return (*l_plan);
}

template <int N_RANK>
void Pochoir<N_RANK>::Store_Plan(const char * file_name, Pochoir_Plan<N_RANK> & _plan) {
    char * l_base_file_name = new char[100];
    sprintf(l_base_file_name, "base_%s", file_name);
    char * l_sync_file_name = new char[100];
    sprintf(l_sync_file_name, "sync_%s", file_name);
	if (_plan.sz_base_data_ > 0)
	{
    	_plan.store_plan(l_base_file_name, l_sync_file_name);
	}
    del_arr(l_base_file_name);
    del_arr(l_sync_file_name);
    return;
}

template <int N_RANK> 
Pochoir_Plan<N_RANK> & Pochoir<N_RANK>::Load_Plan(const char * file_name) {
    char * l_base_file_name = new char[100];
    sprintf(l_base_file_name, "base_%s", file_name);
    char * l_sync_file_name = new char[100];
    sprintf(l_sync_file_name, "sync_%s", file_name);
    Pochoir_Plan<N_RANK> * l_plan = new Pochoir_Plan<N_RANK>();
    l_plan->load_plan(l_base_file_name, l_sync_file_name);
    return (*l_plan);
}

template <int N_RANK> template <typename T_Array>
Pochoir_Plan<N_RANK> & Pochoir<N_RANK>::Load_Plan_Obase(const char * file_name, T_Array & a) {
    char * l_base_file_name = new char[100];
    sprintf(l_base_file_name, "base_%s", file_name);
    char * l_sync_file_name = new char[100];
    sprintf(l_sync_file_name, "sync_%s", file_name);
    Pochoir_Plan<N_RANK> * l_plan = new Pochoir_Plan<N_RANK>();
    l_plan->load_plan(l_base_file_name, l_sync_file_name);

    /* assuming the *.so is already compiled by a previous run of Gen_Plan_Obase() */
    l_plan->load_kernels(*this, a); 
    return (*l_plan);
}

template <int N_RANK> template <typename T_Array, typename ... T_ArrayS>
Pochoir_Plan<N_RANK> & Pochoir<N_RANK>::Load_Plan_Obase(const char * file_name, T_Array & a, T_ArrayS ... as) {
    char * l_base_file_name = new char[100];
    sprintf(l_base_file_name, "base_%s", file_name);
    char * l_sync_file_name = new char[100];
    sprintf(l_sync_file_name, "sync_%s", file_name);
    Pochoir_Plan<N_RANK> * l_plan = new Pochoir_Plan<N_RANK>();
    l_plan->load_plan(l_base_file_name, l_sync_file_name);

    /* assuming the *.so is already compiled by a previous run of Gen_Plan_Obase() */
    l_plan->load_kernels(*this, a, as ...);
    return (*l_plan);
}

template <int N_RANK>
void Pochoir<N_RANK>::Destroy_Plan(Pochoir_Plan<N_RANK> & _plan) {
    struct timeval l_start, l_end;
    double l_min_tdiff = INF;
    gettimeofday(&l_start, 0);
    _plan.unload_kernels();
    Pochoir_Plan<N_RANK> * l_plan = &_plan;
    del_ele(l_plan);
    gettimeofday(&l_end, 0);
    l_min_tdiff = min (l_min_tdiff, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(0, "Dynamic Unloading time : %.6f milliseconds\n", l_min_tdiff);
    return;
}

template <int N_RANK>
void Pochoir<N_RANK>::Run(Pochoir_Plan<N_RANK> & _plan) {
#ifdef CHECK_SHAPE
    inRun = true;
#endif
    int offset = 0, l_sz_pigk = reg_tile_kernels_.size();
    for (int j = 0; _plan.sync_data_->region_[j] != END_SYNC; ++j) {
        for (int i = offset; i < _plan.sync_data_->region_[j]; ++i) {
            //int l_region_n = _plan.base_data_->region_[i].region_n_;
            int l_t0 = _plan.base_data_->region_[i].t0_;
            int l_t1 = _plan.base_data_->region_[i].t1_;
            Grid_Info<N_RANK> l_grid = _plan.base_data_->region_[i].grid_;
            for (int t = l_t0; t < l_t1; ++t) {
                for (int i = 0; i < l_sz_pigk; ++i) {
                    Pochoir_Run_Regional_Guard_Tile_Kernel<N_RANK> l_kernel(time_shift_, reg_guards_[i], reg_tile_kernels_[i]);
                    meta_grid_boundary<N_RANK>::single_step(t, l_grid, phys_grid_, l_kernel);
                }
                for (int i = 0; i < N_RANK; ++i) {
                    l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
                }
            }
        }
        offset = _plan.sync_data_->region_[j];
    }
#ifdef CHECK_SHAPE
    inRun = false;
#endif
    return;
}

template <int N_RANK>
void Pochoir<N_RANK>::Run_Obase_Unroll_T_with_plan(Pochoir_Plan<N_RANK> & _plan){
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    algor.set_unroll(lcm_unroll_);
    algor.set_time_shift(time_shift_);
    if (pmode_ == Pochoir_Obase_Tile) {
        algor.set_opks(reg_obase_kernels_.get_root());
    } else {
        ERROR("Something is wrong in Run_Obase(Plan)!\n");
    }
    checkFlags();
    struct timeval l_start, l_end;
    gettimeofday(&l_start, 0);
#if USE_CILK_FOR
    int offset = 0;
    for (int j = 0; _plan.sync_data_->region_[j] != END_SYNC; ++j) {
        cilk_for (int i = offset; i < _plan.sync_data_->region_[j]; ++i) {
            int l_region_n = _plan.base_data_->region_[i].region_n_;
            assert(l_region_n >= 0);
            int l_t0 = _plan.base_data_->region_[i].t0_;
            int l_t1 = _plan.base_data_->region_[i].t1_;
            Grid_Info<N_RANK> l_grid = _plan.base_data_->region_[i].grid_;
            algor.plan_cut_p(l_t0, l_t1, l_grid, l_region_n);
        }
        offset = _plan.sync_data_->region_[j];
    }
#else
    int offset = 0, i;
    for (int j = 0; _plan.sync_data_->region_[j] != END_SYNC; ++j) {
        for (i = offset; i < _plan.sync_data_->region_[j]-1; ++i) {
            int l_region_n = _plan.base_data_->region_[i].region_n_;
            assert(l_region_n >= 0);
            int l_t0 = _plan.base_data_->region_[i].t0_;
            int l_t1 = _plan.base_data_->region_[i].t1_;
            Grid_Info<N_RANK> l_grid = _plan.base_data_->region_[i].grid_;
            cilk_spawn algor.plan_cut_p(l_t0, l_t1, l_grid, l_region_n);
        }
        int l_region_n = _plan.base_data_->region_[i].region_n_;
        assert(l_region_n >= 0);
        int l_t0 = _plan.base_data_->region_[i].t0_;
        int l_t1 = _plan.base_data_->region_[i].t1_;
        Grid_Info<N_RANK> l_grid = _plan.base_data_->region_[i].grid_;
        algor.plan_cut_p(l_t0, l_t1, l_grid, l_region_n);
        cilk_sync;
        offset = _plan.sync_data_->region_[j];
    }
#endif
    gettimeofday(&l_end, 0);
    pochoir_time_ = min(pochoir_time_, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(0, "Pochoir time = %.6f milliseconds\n", pochoir_time_);
#if DEBUG
    int l_num_kernel = 0, l_num_cond_kernel = 0, l_num_bkernel = 0, l_num_cond_bkernel = 0;
    algor.read_stat_kernel(l_num_kernel, l_num_cond_kernel, l_num_bkernel, l_num_cond_bkernel);
    LOG_ARGS(0, "kernel = %d, cond_kernel = %d, bkernel = %d, cond_bkernel = %d\n",
            l_num_kernel, l_num_cond_kernel, l_num_bkernel, l_num_cond_bkernel);
#endif
    return;
}

template <int N_RANK> 
void Pochoir<N_RANK>::Run_Obase_All_Cond_with_plan(Pochoir_Plan<N_RANK> & _plan) {
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    algor.set_unroll(lcm_unroll_);
    algor.set_time_shift(time_shift_);
    size_t l_size = reg_obase_kernels_.size();
    if (pmode_ == Pochoir_Obase_Tile) {
        algor.set_opks(reg_obase_kernels_.get_root());
    } else {
        ERROR("Something is wrong in Run_Obase(Plan)!\n");
    }
    checkFlags();
    struct timeval l_start, l_end;
    gettimeofday(&l_start, 0);
#if USE_CILK_FOR
    int offset = 0;
    for (int j = 0; _plan.sync_data_->region_[j] != END_SYNC; ++j) {
        cilk_for (int i = offset; i < _plan.sync_data_->region_[j]; ++i) {
            int l_region_n = _plan.base_data_->region_[i].region_n_;
            assert(l_region_n >= 0);
            int l_t0 = _plan.base_data_->region_[i].t0_;
            int l_t1 = _plan.base_data_->region_[i].t1_;
            Grid_Info<N_RANK> l_grid = _plan.base_data_->region_[i].grid_;
            auto f = (*reg_obase_kernels_[l_region_n]).kernel_;
            auto bf = (*reg_obase_kernels_[l_region_n]).bkernel_;
            algor.plan_bicut_mp(l_t0, l_t1, l_grid, l_region_n, (*f), (*bf));
        }
        offset = _plan.sync_data_->region_[j];
    }
#else
    int offset = 0, i;
    for (int j = 0; _plan.sync_data_->region_[j] != END_SYNC; ++j) {
        for (i = offset; i < _plan.sync_data_->region_[j]-1; ++i) {
            int l_region_n = _plan.base_data_->region_[i].region_n_;
            assert(l_region_n >= 0);
            int l_t0 = _plan.base_data_->region_[i].t0_;
            int l_t1 = _plan.base_data_->region_[i].t1_;
            Grid_Info<N_RANK> l_grid = _plan.base_data_->region_[i].grid_;
            auto f = (*reg_obase_kernels_[l_region_n]).kernel_;
            auto bf = (*reg_obase_kernels_[l_region_n]).bkernel_;
            cilk_spawn algor.plan_bicut_mp(l_t0, l_t1, l_grid, l_region_n, (*f), (*bf));
        }
        int l_region_n = _plan.base_data_->region_[i].region_n_;
        assert(l_region_n >= 0);
        int l_t0 = _plan.base_data_->region_[i].t0_;
        int l_t1 = _plan.base_data_->region_[i].t1_;
        Grid_Info<N_RANK> l_grid = _plan.base_data_->region_[i].grid_;
        auto f = (*reg_obase_kernels_[l_region_n]).kernel_;
        auto bf = (*reg_obase_kernels_[l_region_n]).bkernel_;
        algor.plan_bicut_mp(l_t0, l_t1, l_grid, l_region_n, (*f), (*bf));
        cilk_sync;
        offset = _plan.sync_data_->region_[j];
    }
#endif
    gettimeofday(&l_end, 0);
    pochoir_time_ = min (pochoir_time_, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(0, "Pochoir time = %.6f milliseconds\n", pochoir_time_);
#if DEBUG
    int l_num_kernel = 0, l_num_cond_kernel = 0, l_num_bkernel = 0, l_num_cond_bkernel = 0;
    algor.read_stat_kernel(l_num_kernel, l_num_cond_kernel, l_num_bkernel, l_num_cond_bkernel);
    LOG_ARGS(0, "kernel = %d, cond_kernel = %d, bkernel = %d, cond_bkernel = %d\n",
            l_num_kernel, l_num_cond_kernel, l_num_bkernel, l_num_cond_bkernel);
#endif
    return;
}


template <int N_RANK>
void Pochoir<N_RANK>::Run_Obase_Unroll_T(Pochoir_Plan<N_RANK> & _plan) 
{
	//Run_Obase_Unroll_T_with_plan(_plan) ;
	Run_Obase_Unroll_T_without_plan() ;
}
template <int N_RANK>
void Pochoir<N_RANK>::Run_Obase_All_Cond(Pochoir_Plan<N_RANK> & _plan)
{
	Run_Obase_All_Cond_with_plan(_plan) ;
}

template <int N_RANK>
void Pochoir<N_RANK>::Run_Obase_Unroll_T_without_plan()
{
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    algor.set_unroll(lcm_unroll_);
    algor.set_time_shift(time_shift_);
    algor.set_time_step(timestep_);
    if (pmode_ == Pochoir_Obase_Tile) {
		cout << "Not using plans at runtime" << endl ;
        algor.set_opks(reg_obase_kernels_.get_root());
		algor.set_num_kernels(reg_obase_kernels_.pointer_) ;
    } else {
        ERROR("Something is wrong in Run_Obase(Plan)!\n");
    }
    checkFlags();
    struct timeval l_start, l_end;
    gettimeofday(&l_start, 0);
	algor.space_time_cut_boundary(0 + time_shift_, timestep_ + time_shift_, 
								logic_grid_);
    gettimeofday(&l_end, 0);
    pochoir_time_ = min(pochoir_time_, (1.0e3 * tdiff(&l_end, &l_start)));
    LOG_ARGS(0, "Pochoir time = %.6f milliseconds\n", pochoir_time_);
}
#endif
