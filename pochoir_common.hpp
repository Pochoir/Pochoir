/*
 **********************************************************************************
 *  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 *  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 * 		                     Charles E. Leiserson <cel@mit.edu>
 * 	 
 *   This program is delete software: you can redistribute it and/or modify
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

#ifndef POCHOIR_COMMON_H
#define POCHOIR_COMMON_H

#include <sys/time.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>

static inline double tdiff (struct timeval *a, struct timeval *b)
{
	    return a->tv_sec - b->tv_sec + 1e-6 * (a->tv_usec - b->tv_usec);
}

static inline int StrToInt(const char * s)
{
  return atoi(s);
}

/* greatest common divisor */
static inline int gcd(int a, int b) {
    if (b == 0)
        return a;
    else
        return gcd(b, a % b);
}

/* lowest common multiple of 'a' and 'b' */
static inline int lcm(int a, int b) {
    return ((a * b) / gcd(a, b));
}

#define ARRAY_LENGTH(x) (int)(sizeof(x)/sizeof(x[0]))

/* due to the fact that bit trick is much slower than conditional instruction,
 * let's disable it for now!!!
 */
#define BIT_TRICK 0
#define INF 100000000
#define SUPPORT_RANK 9
#define DEBUG_FACILITY 1
// #define DEBUG 0
#define END_SYNC -1

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
/* a bit tricky version of modulo operation, assuming a < 2 * b */
#define pmod(a, b) ((a) - ((b) & -((a)>=(b))))
#define pmod_lu(a, lb, ub) ((a) - (((ub)-(lb)) & -((a)>=(ub))))

#define pCond(b, x, y) (x&(-b)) | (y&-(!b))

static inline bool select(bool b, bool x, bool y) {
    return (x&(-b)) | (y&-(!b));
}
static inline int select(bool b, int x, int y) {
    return (x&(-b)) | (y&-(!b));
}
static inline float select(bool b, float x, float y) {
    int __ir__ = ((*(int*)&x) & (-b)) | ((*(int*)&y) & -(!b)); 
    return *(float*)&__ir__; 
}   
static inline double select(bool b, double x, double y) {
    long __ir__ = ((*(long*)&x) & (-b)) | ((*(long*)&y) & -(!b));
    return *(double*)&__ir__;
}

typedef int T_dim;
typedef int T_index;

template <int N_RANK>
struct Grid_Info {
    int x0[N_RANK], x1[N_RANK];
    int dx0[N_RANK], dx1[N_RANK];
};

template <int N_RANK>
struct Pochoir_Shape {
    /* N_RANK + 1 because we probably have to include the time dimension
     * to correctly calculate the slope[]
     */
    int shift[N_RANK+1];
};
 
enum Meta_Op { IS_ROOT, IS_SPAWN, IS_SYNC, IS_INTERNAL };

template <int N_RANK>
struct Region_Info {
    int t0, t1;
    Grid_Info<N_RANK> grid;
    int region_n;
    Region_Info() {
        t0 = 0; t1 = 0; region_n = -1;
    }

    Region_Info(int _t0, int _t1, Grid_Info<N_RANK> _grid, int _region_n) {
        t0 = _t0; t1 = _t1; grid = _grid;
        region_n = _region_n;
    }

    friend std::ofstream & operator<<(std::ofstream & fs, Region_Info<N_RANK> const & r) {
        int i;
        fs << "{ BASE, ";
        fs << "t = {" << r.t0 << ", " << r.t1 << "}, {";

        fs << "x0 = {";
        for (i = 0; i < N_RANK-1; ++i) {
            fs << r.grid.x0[i] << ", ";
        }
        fs << r.grid.x0[i] << "}, ";

        fs << "x1 = {";
        for (i = 0; i < N_RANK-1; ++i) {
            fs << r.grid.x1[i] << ", ";
        }
        fs << r.grid.x1[i] << "}, ";

        fs << "dx0 = {";
        for (i = 0; i < N_RANK-1; ++i) {
            fs << r.grid.dx0[i] << ", ";
        }
        fs << r.grid.dx0[i] << "}, ";

        fs << "dx1 = {";
        for (i = 0; i < N_RANK-1; ++i) {
            fs << r.grid.dx1[i] << ", ";
        }
        fs << r.grid.dx1[i] << "}}, " << r.region_n << "}";
        return fs;
    }

    void pscanf(FILE * fs) {
        int i;
        fscanf(fs, "{ BASE, ");
        fscanf(fs, "t = {%d, %d}, {",  &(t0), &(t1));

        fscanf(fs, "x0 = {");
        for (i = 0; i < N_RANK-1; ++i) {
            fscanf(fs, "%d, ", &(grid.x0[i]));
        }
        fscanf(fs, "%d}, ", &(grid.x0[i]));

        fscanf(fs, "x1 = {");
        for (i = 0; i < N_RANK-1; ++i) {
            fscanf(fs, "%d, ", &(grid.x1[i]));
        }
        fscanf(fs, "%d}, ", &(grid.x1[i]));

        fscanf(fs, "dx0 = {");
        for (i = 0; i < N_RANK-1; ++i) {
            fscanf(fs, "%d, ", &(grid.dx0[i]));
        }
        fscanf(fs, "%d}, ", &(grid.dx0[i]));

        fscanf(fs, "dx1 = {");
        for (i = 0; i < N_RANK-1; ++i) {
            fscanf(fs, "%d, ", &(grid.dx1[i]));
        }
        fscanf(fs, "%d}}, %d}\n", &(grid.dx1[i]), &(region_n));
        return;
    }

    friend std::ostream & operator<<(std::ostream & fs, Region_Info<N_RANK> const & r) {
        int i;
        fs << "{ BASE, ";
        fs << "t = {" << r.t0 << ", " << r.t1 << "}, {";

        fs << "x0 = {";
        for (i = 0; i < N_RANK-1; ++i) {
            fs << r.grid.x0[i] << ", ";
        }
        fs << r.grid.x0[i] << "}, ";

        fs << "x1 = {";
        for (i = 0; i < N_RANK-1; ++i) {
            fs << r.grid.x1[i] << ", ";
        }
        fs << r.grid.x1[i] << "}, ";

        fs << "dx0 = {";
        for (i = 0; i < N_RANK-1; ++i) {
            fs << r.grid.dx0[i] << ", ";
        }
        fs << r.grid.dx0[i] << "}, ";

        fs << "dx1 = {";
        for (i = 0; i < N_RANK-1; ++i) {
            fs << r.grid.dx1[i] << ", ";
        }
        fs << r.grid.dx1[i] << "}}, " << r.region_n << "}";
        return fs;
    }
};

template <typename T>
struct Vector_Info {
    T * region_;
    int pointer_, size_;
    Vector_Info(int size) {
        region_ = (T *) calloc(size, sizeof(T));
        pointer_ = 0; size_ = size;
#if DEBUG
        printf("init size = %d\n", size_);
#endif
    }
    ~Vector_Info() {
        pointer_ = size_ = 0;
        delete region_;
    } 
    void add_element(T ele) {
#if DEBUG
//        std::cerr << "add_element " << ele << std::endl;
#endif
        if (pointer_ < size_) {
            region_[pointer_] = ele;
            ++pointer_;
        } else {
#if DEBUG
            printf("realloc memory size = %d -> %d!\n", size_, 2*size_);
#endif
            T * l_region = (T *) calloc(2*size_, sizeof(T));
            if (l_region != NULL) {
                for (int i = 0; i < size_; ++i) {
                    l_region[i] = region_[i];
                }
                delete(region_);
                region_ = l_region;
                region_[pointer_] = ele;
                ++pointer_;
                size_ = 2 * size_;
            } else {
                printf("realloc wrong!\n");
                exit(1);
            }
        }
        return;
    }
    void scan () {
        for (int i = 1; i < pointer_; ++i) {
            region_[i] = region_[i] + region_[i-1];
        }
    }
    int size() { return pointer_; }

    friend std::ofstream & operator<<(std::ofstream & fs, Vector_Info<T> const & v) {
        for (int i = 0; i < v.pointer_; ++i) {
            fs << v.region_[i] << std::endl;
            // std::cerr << region_[i] << "\n";
        }
        return fs;
    }
    friend std::ostream & operator<<(std::ostream & fs, Vector_Info<T> const & v) {
        for (int i = 0; i < v.pointer_; ++i) {
            fs << v.region_[i] << std::endl;
            // std::cerr << region_[i] << "\n";
        }
        return fs;
    }
};

template <int N_RANK>
struct Pochoir_Plan {
    Vector_Info< Region_Info<N_RANK> > * base_data_;
    Vector_Info<int> * sync_data_;
    int sz_base_data_, sz_sync_data_;
    Pochoir_Plan (int _sz_base_data, int _sz_sync_data) : sz_base_data_(_sz_base_data), sz_sync_data_(_sz_sync_data) {
        base_data_ = new Vector_Info< Region_Info<N_RANK> >(_sz_base_data);
        sync_data_ = new Vector_Info<int>(_sz_sync_data);
    }
    Pochoir_Plan () {
        sz_base_data_ = sz_sync_data_ = 0;
        base_data_ = NULL;
        sync_data_ = NULL;
    }
    void alloc_base_data(int _sz_base_data) {
        sz_base_data_ = _sz_base_data;
        base_data_ = new Vector_Info< Region_Info<N_RANK> >(_sz_base_data);
    }
    void alloc_sync_data(int _sz_sync_data) {
        sz_sync_data_ = _sz_sync_data;
        sync_data_ = new Vector_Info<int>(_sz_sync_data);
    }
    ~Pochoir_Plan() {
        sz_base_data_ = sz_sync_data_ = 0;
        delete base_data_;
        delete sync_data_;
    }
    void store_plan(const char * base_file_name, const char * sync_file_name) {
        /* extract data and store plan */
        std::ofstream os_base_data(base_file_name);
        std::ofstream os_sync_data(sync_file_name);
        if (os_base_data.is_open()) {
#if DEBUG
            printf("os_base_data is open!\n");
#endif
            os_base_data << (*base_data_);
        } else {
            printf("os_base_data is NOT open! exit!\n");
            exit(1);
        }
        if (os_sync_data.is_open()) {
#if DEBUG
            printf("os_sync_data is open!\n");
#endif
            os_sync_data << (*sync_data_);
        } else {
            printf("os_sync_data is NOT open! exit!\n");
            exit(1);
        }
        os_base_data.close();
        os_sync_data.close();
        return;
    }
    Pochoir_Plan<N_RANK> & load_plan(const char * base_file_name, const char * sync_file_name) {
        if (sz_base_data_ != 0) {
            delete base_data_;
            sz_base_data_ = 0;
        }
        if (sz_sync_data_ != 0) {
            delete sync_data_;
            sz_sync_data_ = 0;
        }
        alloc_base_data(10);
        alloc_sync_data(10);

        FILE * is_base_data = fopen(base_file_name, "r");

        if (is_base_data != NULL) {
#if DEBUG
            printf("is_base_data is open!\n");
#endif
            while (!feof(is_base_data)) {
                Region_Info<N_RANK> l_region;
                l_region.pscanf(is_base_data);
                base_data_->add_element(l_region);
            }
        } else {
            printf("is_base_data is NOT open! exit!\n");
            exit(1);
        }
        fclose(is_base_data);
#if DEBUG
        std::cerr << "base_data : \n" << (*base_data_) << std::endl;
#endif

        FILE * is_sync_data = fopen(sync_file_name, "r");
        if (is_sync_data != NULL) {
#if DEBUG
            printf("is_sync_data is open!\n");
#endif
            while (!feof(is_sync_data)) {
                int l_sync;
                fscanf(is_sync_data, "%d\n", &l_sync);
                sync_data_->add_element(l_sync);
            }
        } else {
            printf("is_sync_data is NOT open! exit!\n");
            exit(1);
        }
        fclose(is_sync_data);
#if DEBUG
        std::cerr << "sync_data : \n" << (*sync_data_) << std::endl;
#endif

        return (*this);
    }
};

template <int N_RANK>
struct Node_Info {
    Region_Info<N_RANK> region_;
    Node_Info<N_RANK> *parent, *left, *right;
    enum Meta_Op op;

    Node_Info() { parent = left = right = NULL; }

    Node_Info(enum Meta_Op _op) {
        /* constructor */
        op = _op;
        parent = left = right = NULL;
    }

    Node_Info(int _t0, int _t1, Grid_Info<N_RANK> & _grid) : region_(_t0, _t1, _grid, -1) {
        /* constructor */
        op = IS_INTERNAL; 
        parent = left = right = NULL;
    }

    ~Node_Info() {
        /* destructor */
        parent = left = right = NULL;
    }   
};

template <int N_RANK>
struct Spawn_Tree {
    Node_Info<N_RANK> * root_;
    int size_;
    Spawn_Tree() { 
        /* constructor */
        root_ = new Node_Info<N_RANK>(IS_ROOT); size_ = 1; 
    }

    Spawn_Tree(int _t0, int _t1, Grid_Info<N_RANK> & _grid) {
        /* constructor */
        root_ = new Node_Info<N_RANK>(_t0, _t1, _grid);
        size_ = 1;
    }
    ~Spawn_Tree() {
        /* free the entire tree */
        if (root_->left == NULL) {
            delete root_;
            return;
        } else {
            dfs_rm_tree(root_);
            return;
        }
    }
    void dfs_rm_tree(Node_Info<N_RANK> * parent) {
        Node_Info<N_RANK> * l_node;

        if (parent->left == NULL) {
            rm_node(parent);
            return;
        } else {
            l_node = parent->right;
            dfs_rm_tree(parent->left);
            assert(parent->left == NULL);
            rm_node(parent);
            dfs_rm_tree(l_node);
            return;
        }
    }
    Node_Info<N_RANK> * get_root() {  return root_; }
    int size() { return size_; }
    void add_node(Node_Info<N_RANK> * parent, Node_Info<N_RANK> * child, enum Meta_Op _op) {
        if (parent->left == NULL) {
            parent->left = child;
            child->parent = parent;
            child->op = _op;
            ++size_;
            return;
        }
        Node_Info<N_RANK> * youngest_child = parent->left;
        while (youngest_child->right != NULL)
            youngest_child = youngest_child->right;
        youngest_child->right = child;
        child->parent = parent;
        child->op = _op;
        ++size_;
        return;
    }

    void add_node(Node_Info<N_RANK> * parent, Node_Info<N_RANK> * child, enum Meta_Op _op, int _region_n) {
        if (parent->left == NULL) {
            parent->left = child;
            child->parent = parent;
            child->op = _op; (child->region_).region_n = _region_n;
            ++size_;
            return;
        }
        Node_Info<N_RANK> * youngest_child = parent->left;
        while (youngest_child->right != NULL)
            youngest_child = youngest_child->right;
        youngest_child->right = child;
        child->parent = parent;
        child->op = _op; (child->region_).region_n = _region_n;
        ++size_;
        return;
    }

    void rm_node(Node_Info<N_RANK> * node) {
        if (node == node->parent->left) {
            /* I am the biggest brother */
            node->parent->left = node->right;
        } else {
            /* I am NOT the biggest brother */
            Node_Info<N_RANK> * l_node = node->parent->left;
            while (l_node->right != node) {
                l_node = l_node->right;
            }
            assert(l_node->right == node);
            l_node->right = node->right;
        }
        delete(node);
        if (node != NULL)
            node = NULL;
        --size_;
    }
    void dfs_until_sync(Node_Info<N_RANK> * node, Vector_Info< Region_Info<N_RANK> > & base_data) {
        if (node == NULL) {
            return;
        }
        /* visit node */
        if (node->op == IS_SYNC) {
            return;
        }
        Node_Info<N_RANK> * l_node = NULL;
        if (node->op == IS_SPAWN) {
            assert(node->left == NULL);
            l_node = node->right;
            base_data.add_element(node->region_);
            rm_node(node);
        } else if (node->op == IS_INTERNAL) {
            l_node = node->right;
            dfs_until_sync(node->left, base_data);
        }
        /* visit its brothers */
        dfs_until_sync(l_node, base_data);
        return;
    }
    void dfs_rm_sync(Node_Info<N_RANK> * node) {
        if (node == NULL) {
            return;
        }
        if (node->op == IS_SPAWN) {
            return;
        }
        if (node->op == IS_SYNC && node == node->parent->left) {
            /* remove the biggest sync */
            rm_node(node);
            return;
        }
        Node_Info<N_RANK> * l_node = NULL;
        if (node->op == IS_INTERNAL && node->left == NULL) {
            /* remove the biggest empty internal node */
            l_node = node->right;
            rm_node(node);
        } else if (node->op == IS_INTERNAL && node->left != NULL) {
            l_node = node->right;
            dfs_rm_sync(node->left);
            if (node->left == NULL) {
                rm_node(node);
            }
        }
        dfs_rm_sync(l_node);
        return;
    }
};

template <int N_RANK, size_t N>
size_t ArraySize (Pochoir_Shape<N_RANK> (& arr)[N]) { return N; }

#define KLEIN 0
#define USE_CILK_FOR 0
#define BICUT 1
#define STAT 0
static bool inRun = false;
static int home_cell_[9];

static inline void klein(int & new_i, int & new_j, Grid_Info<2> const & grid) {
    int l_arr_size_1 = grid.x1[1] - grid.x0[1];
    int l_arr_size_0 = grid.x1[0] - grid.x0[0];

    if (new_i < grid.x0[1])
        new_i += l_arr_size_1;
    else if (new_i >= grid.x1[1])
        new_i -= l_arr_size_1;
    if (new_j < grid.x0[0]) {
        new_j += l_arr_size_0;
        new_i  = grid.x0[1] + (grid.x1[1] - 1 - new_i);
    } else if (new_j >= grid.x1[0]) {
        new_j -= l_arr_size_0;
        new_i  = grid.x0[1] + (grid.x1[1] - 1 - new_i);
    }
    return;
}

static inline void klein_region(Grid_Info<2> & grid, Grid_Info<2> const & initial_grid) {
    Grid_Info<2> orig_grid;
    const int l_arr_size_1 = initial_grid.x1[1] - initial_grid.x0[1];
    const int l_arr_size_0 = initial_grid.x1[0] - initial_grid.x0[0];

    if (grid.x0[1] >= initial_grid.x1[1]) {
        grid.x0[1] -= l_arr_size_1;
        grid.x1[1] -= l_arr_size_1;
    } else if (grid.x1[1] < initial_grid.x0[1]) {
        grid.x0[1] += l_arr_size_1;
        grid.x1[1] += l_arr_size_1;
    } 
    orig_grid = grid;
    if (grid.x0[0] >= initial_grid.x1[0]) {
        grid.x0[0] -= l_arr_size_0;
        grid.x1[0] -= l_arr_size_0;
        grid.x0[1] = initial_grid.x0[1] + (initial_grid.x1[1] - orig_grid.x1[1]);
        grid.x1[1] = initial_grid.x0[1] + (initial_grid.x1[1] - orig_grid.x0[1]);
        grid.dx0[1] = -orig_grid.dx1[1];
        grid.dx1[1] = -orig_grid.dx0[1];
    } else if (grid.x1[0] < initial_grid.x0[0]) {
        grid.x0[0] += l_arr_size_0;
        grid.x1[0] += l_arr_size_0;
        grid.x0[1] = initial_grid.x0[1] + (initial_grid.x1[1] - orig_grid.x1[1]);
        grid.x1[1] = initial_grid.x0[1] + (initial_grid.x1[1] - orig_grid.x0[1]);
        grid.dx0[1] = -orig_grid.dx1[1];
        grid.dx1[1] = -orig_grid.dx0[1];
    }
    return;
}

#define Pochoir_1D Pochoir<1>
#define Pochoir_2D Pochoir<2>
#define Pochoir_3D Pochoir<3>
#define Pochoir_4D Pochoir<4>
#define Pochoir_5D Pochoir<5>
#define Pochoir_6D Pochoir<6>
#define Pochoir_7D Pochoir<7>
#define Pochoir_8D Pochoir<8>

#define Pochoir_Array_1D(type) Pochoir_Array<type, 1>
#define Pochoir_Array_2D(type) Pochoir_Array<type, 2>
#define Pochoir_Array_3D(type) Pochoir_Array<type, 3>
#define Pochoir_Array_4D(type) Pochoir_Array<type, 4>
#define Pochoir_Array_5D(type) Pochoir_Array<type, 5>
#define Pochoir_Array_6D(type) Pochoir_Array<type, 6>
#define Pochoir_Array_7D(type) Pochoir_Array<type, 7>
#define Pochoir_Array_8D(type) Pochoir_Array<type, 8>

#define Pochoir_Shape_1D Pochoir_Shape<1>
#define Pochoir_Shape_2D Pochoir_Shape<2>
#define Pochoir_Shape_3D Pochoir_Shape<3>
#define Pochoir_Shape_4D Pochoir_Shape<4>
#define Pochoir_Shape_5D Pochoir_Shape<5>
#define Pochoir_Shape_6D Pochoir_Shape<6>
#define Pochoir_Shape_7D Pochoir_Shape<7>
#define Pochoir_Shape_8D Pochoir_Shape<8>

/* these lambda functions are for computing internal/boundary region,
 * the original 'f'/'bf'
 */

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

/* - these function templates are for computing boundary values, currently
 *   icc doesn't support capturing the lambda function by function objects,
 *   so, we have to utilize the function pointers!
 * - because these functions will be called inside T & operator() functions,
 *   so we have to return a value of T&
 */
#define Pochoir_Boundary_1D(name, arr, t, i) \
    template <typename T> \
    T name (Pochoir_Array<T, 1> & arr, int t, int i) { 

#define Pochoir_Boundary_2D(name, arr, t, i, j) \
    template <typename T> \
    T name (Pochoir_Array<T, 2> & arr, int t, int i, int j) { 

#define Pochoir_Boundary_3D(name, arr, t, i, j, k) \
    template <typename T> \
    T name (Pochoir_Array<T, 3> & arr, int t, int i, int j, int k) { 

#define Pochoir_Boundary_4D(name, arr, t, i, j, k, l) \
    template <typename T> \
    T name (Pochoir_Array<T, 4> & arr, int t, int i, int j, int k, int l) { 

#define Pochoir_Boundary_5D(name, arr, t, i, j, k, l, m) \
    template <typename T> \
    T name (Pochoir_Array<T, 5> & arr, int t, int i, int j, int k, int l, int m) { 

#define Pochoir_Boundary_6D(name, arr, t, i, j, k, l, m, n) \
    template <typename T> \
    T name (Pochoir_Array<T, 6> & arr, int t, int i, int j, int k, int l, int m, int n) { 

#define Pochoir_Boundary_7D(name, arr, t, i, j, k, l, m, n, o) \
    template <typename T> \
    T name (Pochoir_Array<T, 7> & arr, int t, int i, int j, int k, int l, int m, int n, int o) { 

#define Pochoir_Boundary_8D(name, arr, t, i, j, k, l, m, n, o, p) \
    template <typename T> \
    T name (Pochoir_Array<T, 8> & arr, int t, int i, int j, int k, int l, int m, int n, int o, int p) { 

#define Pochoir_Boundary_End }

#endif /* POCHOIR_COMMON_H */
