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

#ifndef POCHOIR_COMMON_H
#define POCHOIR_COMMON_H

#include <sys/time.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <functional>
#include <vector>
#include <boost/type_traits.hpp> /* to decide if a typename T is a pointer? */
#include "pochoir_errmsg.hpp"
#include "pochoir_dloader.hpp"

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

#if 0
#define cilk_for for
#define cilk_spawn 
#define cilk_sync
#endif

/* due to the fact that bit trick is much slower than conditional instruction,
 * let's disable it for now!!!
 */
#define BIT_TRICK 0
// #define INF 100000000
#define SUPPORT_RANK 9
#define DEBUG_FACILITY 1
// #define DEBUG 0
#define END_SYNC -1
#define KLEIN 0
#define USE_CILK_FOR 0
#define NONE_EXCLUSIVE_IFS -1
#define CROSS_REGION -2
#define BICUT 1
#define STAT 0
#define VECTOR_SIZE 10
#define FNAME_LENGTH 200

enum Pochoir_Mode {
    Pochoir_Null,
    Pochoir_Stagger,
    Pochoir_Tile,
    Pochoir_Obase_Tile
};

// define an alias to the array of array of Pochoir_Kernel
#define POCHOIR_TILE Pochoir_Kernel
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define pabs(a, b) ((a) > (b) ? ((a) - (b)) : ((b) - (a)))
/* a bit tricky version of modulo operation, assuming a < 2 * b */
#define pmod(a, b) ((a) - ((b) & -((a)>=(b))))
#define pmod_lu(a, lb, ub) ((a) - (((ub)-(lb)) & -((a)>=(ub))))

#define pCond(b, x, y) (x&(-b)) | (y&-(!b))
#define del_ele(obj) \
    do { \
    if ((obj) != NULL) {\
        delete (obj); \
        (obj) = NULL; \
    } } while (0) 

#define del_arr(obj) \
    do { \
    if ((obj) != NULL) {\
        delete [] (obj); \
        (obj) = NULL; \
    } } while (0)

static inline bool select(bool b, bool x, bool y) {
    return (x&(-b)) | (y&-(!b));
}
static inline int select(bool b, int x, int y) {
    return (x&(-b)) | (y&-(!b));
}
static inline unsigned int log2_floor(unsigned long a) {
    return 63-__builtin_clzl(a);
}
static inline unsigned int log2_floor(unsigned int a) {
    return 31-__builtin_clz(a);
}
static inline unsigned int pow2(unsigned int e) {
    assert(e < 32);
    return (0x1 << e);
}
static inline unsigned long pow2l(unsigned int e) {
    assert(e < 64);
    return (0x1 << e);
}

typedef int T_dim;
typedef int T_index;
/* T_color could be of type int, long, ..., which could be able to
 * perform bit-wise operations 
 */
typedef int T_color;

struct Homogeneity {
    T_color o_, a_;
    int size_; /* # registered kernel in total */
    Homogeneity() { o_ = 0; a_ = 0; size_ = 0; }
    Homogeneity(T_color _o, T_color _a) : o_(_o), a_(_a), size_(0) { }
    Homogeneity(T_color _o, T_color _a, int _size) : o_(_o), a_(_a), size_(_size) { }
    Homogeneity(Homogeneity const & _h) : o_(_h.o_), a_(_h.a_), size_(_h.size_) { }
    Homogeneity(int _size) : size_(_size) {
        /* set up a white color with size '_size' */
        o_ = 0; a_ = 0;
        for (int i = 0; i < _size; ++i) {
            o_ <<= 1; o_ |= 1;
        }
    }
    inline void set_size(int _size) {
        /* set up a white color with size '_size' */
        size_ = _size;
        o_ = 0; a_ = 0;
        for (int i = 0; i < _size; ++i) {
            o_ <<= 1; o_ |= 1;
        }
        return;
    }
    inline int size(void) { return size_; }
    inline int size(void) const { return size_; }
    inline bool operator<= (Homogeneity const & h) {
        /* define the partial order between Homogeneities */
        return ((a_ & h.a_) == a_ && (~o_ & ~h.o_) == ~o_);
    }
    inline Homogeneity operator+ (Homogeneity const & h) {
        /* to combine regions, or to get the maximum-norm kernel
         * <= both this and h
         */
        T_color l_o = o_ | h.o_; T_color l_a = a_ & h.a_;

        return Homogeneity(l_o, l_a, size_);
    }
    inline bool operator== (Homogeneity const & h) {
        return (size_ == h.size_ && o_ == h.o_ && a_ == h.a_);
    }
    inline Homogeneity & operator= (Homogeneity const & h) {
        size_ = h.size_; o_ = h.o_; a_ = h.a_;
        return (*this);
    }
    inline bool is_homogeneous(void) { return (o_ == a_); }
    T_color norm(void) {
        /* norm |h| = |(o, a)| = # 0 bits in o-a = Hamming weight |~(o-a)| 
         * = # eliminated if's in an (o, a) clone
         */
        T_color v = ~(o_ ^ a_);
        T_color c = 0;
        int l_size = size_;
        for (c = 0; l_size > 0 && v; v >>= 1, --l_size)
            c += v & 1;
        return c;
    }
    T_color reverse_bits(T_color v){
#if 0
        T_color r = v;
        int s = size_ - 1; // extra shift needed at end

        for (v >>= 1; v; v >>= 1) {
            r <<= 1;
            r |= v & 1;
            s--;
        }
        r <<= s; // shift when v's highest bits are zero
        return r;
#else
        T_color r = 0;
        int s = size_;
        for ( ; v; v >>= 1, r <<= 1, --s) {
            r |= (v & 1);
        }
        r <<= s;
        return r;
#endif
    }
    
    friend std::ofstream & operator<<(std::ofstream & fs, Homogeneity const & h) {
        fs << "(" ;
        T_color l_mask = 0x1 << (h.size_ - 1);
        T_color l_o = h.o_, l_a = h.a_;
        for (int i = 0; i < h.size_; ++i, l_mask >>= 1) {
            if (l_o & l_mask)
                fs << "1";
            else
                fs << "0";
        }
        fs << ", ";
        l_mask = 0x1 << (h.size_ - 1);
        for (int i = 0; i < h.size_; ++i, l_mask >>= 1) {
            if (l_a & l_mask)
                fs << "1";
            else
                fs << "0";
        }
        fs << ")";
        return fs;
    }

    friend std::ostream & operator<<(std::ostream & fs, Homogeneity const & h) {
        fs << "(" ;
        T_color l_mask = 0x1 << (h.size_ - 1);
        T_color l_o = h.o_, l_a = h.a_;
        for (int i = 0; i < h.size_; ++i, l_mask >>= 1) {
            if (l_o & l_mask)
                fs << "1";
            else
                fs << "0";
        }
        fs << ", ";
        l_mask = 0x1 << (h.size_ - 1);
        for (int i = 0; i < h.size_; ++i, l_mask >>= 1) {
            if (l_a & l_mask)
                fs << "1";
            else
                fs << "0";
        }
        fs << ")";
        return fs;
    }

    void pscanf (FILE * _fs) {
        char ch;
        T_color l_o = 0, l_a = 0;
        int l_o_sz = 0, l_a_sz = 0;
        bool l_valid_digit = false;

        fscanf(_fs, "(");
        do {
            l_valid_digit = false;
            fscanf(_fs, "%c", &ch);
            if (ch == '0' || ch == '1') {
                l_o <<= 1;
                l_o |= ch - '0';
                l_valid_digit = true;
                ++l_o_sz;
            }
        } while (l_valid_digit);
        // fscanf(_fs, ", ");
        assert(ch == ',');
        if (ch != ',') {
            ERROR("delimitation ',' is missing!");
        }
        fscanf(_fs, " ");
        do {
            l_valid_digit = false;
            fscanf(_fs, "%c", &ch);
            if (ch == '0' || ch == '1') {
                l_a <<= 1;
                l_a |= ch - '0';
                l_valid_digit = true;
                ++l_a_sz;
            }
        } while (l_valid_digit);
        assert(ch == ')');
        if (l_o_sz == l_a_sz) {
            o_ = l_o; a_ = l_a;
            size_ = l_o_sz;
        } else {
            ERROR_ARGS("o_sz (%d) != a_sz (%d)", l_o_sz, l_a_sz);
        }
    }
};

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
    inline bool operator== (Pochoir_Shape<N_RANK> const & _shape) {
        for (int i = 0; i < N_RANK+1; ++i) {
            if (shift[i] != _shape.shift[i])
                return false;
        }
        return true;
    }
};
 
template <int N_RANK>
struct Pochoir_Types {
    typedef std::function<void (int, int)> T_Kernel;
    typedef std::function<void (int, int, Grid_Info<N_RANK> const &)> T_Obase_Kernel;
    typedef std::function<bool (int, int)> T_Guard;
};
template <>
struct Pochoir_Types<1> {
    typedef std::function<void (int, int)> T_Kernel;
    typedef std::function<void (int, int, Grid_Info<1> const &)> T_Obase_Kernel;
    typedef std::function<bool (int, int)> T_Guard;
};
template <>
struct Pochoir_Types<2> {
#if 1
    typedef std::function<void (int, int, int)> T_Kernel;
    typedef std::function<void (int, int, Grid_Info<2> const &)> T_Obase_Kernel;
    typedef std::function<bool (int, int, int)> T_Guard;
#else
    typedef void (*T_Kernel)(int, int, int);
    typedef void (*T_Obase_Kernel)(int, int, Grid_Info<2> const &);
    typedef bool (*T_Guard)(int, int, int);
#endif
};
template <>
struct Pochoir_Types<3> {
    typedef std::function<void (int, int, int, int)> T_Kernel;
    typedef std::function<void (int, int, Grid_Info<3> const &)> T_Obase_Kernel;
    typedef std::function<bool (int, int, int, int)> T_Guard;
};

enum Meta_Op { IS_ROOT, IS_SPAWN, IS_SYNC, IS_INTERNAL };

template <int N_RANK>
struct Region_Info {
    int t0_, t1_;
    Grid_Info<N_RANK> grid_;
    Homogeneity color_;
    int region_n_;
    Region_Info() : t0_(0), t1_(0), region_n_(NONE_EXCLUSIVE_IFS) { }

    Region_Info(int _t0, int _t1, Grid_Info<N_RANK> const & _grid, Homogeneity _color, int _region_n) : t0_(_t0), t1_(_t1), grid_(_grid), color_(_color), region_n_(_region_n) { }

    inline Region_Info<N_RANK> & operator= (Region_Info<N_RANK> const & r) {
        t0_ = r.t0_; t1_ = r.t1_;
        grid_ = r.grid_; color_ = r.color_; region_n_ = r.region_n_;
        return (*this);
    }
    friend std::ofstream & operator<<(std::ofstream & _fs, Region_Info<N_RANK> const & _r) {
        int i;
        _fs << "{ BASE, ";
        _fs << "t = {" << _r.t0_ << ", " << _r.t1_ << "}, {";

        _fs << "x0 = {";
        for (i = N_RANK-1; i > 0; --i) {
            _fs << _r.grid_.x0[i] << ", ";
        }
        _fs << _r.grid_.x0[i] << "}, ";

        _fs << "x1 = {";
        for (i = N_RANK-1; i > 0; --i) {
            _fs << _r.grid_.x1[i] << ", ";
        }
        _fs << _r.grid_.x1[i] << "}, ";

        _fs << "dx0 = {";
        for (i = N_RANK-1; i > 0; --i) {
            _fs << _r.grid_.dx0[i] << ", ";
        }
        _fs << _r.grid_.dx0[i] << "}, ";

        _fs << "dx1 = {";
        for (i = N_RANK-1; i > 0; --i) {
            _fs << _r.grid_.dx1[i] << ", ";
        }
        _fs << _r.grid_.dx1[i] << "}}, ";
        _fs << _r.color_ << ", "; 
        _fs << _r.region_n_ << "}";
        return _fs;
    }

    void pscanf(FILE * _fs) {
        int i;
        fscanf(_fs, "{ BASE, ");
        fscanf(_fs, "t = {%d, %d}, {",  &(t0_), &(t1_));

        fscanf(_fs, "x0 = {");
        for (i = N_RANK-1; i > 0; --i) {
            fscanf(_fs, "%d, ", &(grid_.x0[i]));
        }
        fscanf(_fs, "%d}, ", &(grid_.x0[i]));

        fscanf(_fs, "x1 = {");
        for (i = N_RANK-1; i > 0; --i) {
            fscanf(_fs, "%d, ", &(grid_.x1[i]));
        }
        fscanf(_fs, "%d}, ", &(grid_.x1[i]));

        fscanf(_fs, "dx0 = {");
        for (i = N_RANK-1; i > 0; --i) {
            fscanf(_fs, "%d, ", &(grid_.dx0[i]));
        }
        fscanf(_fs, "%d}, ", &(grid_.dx0[i]));

        fscanf(_fs, "dx1 = {");
        for (i = N_RANK-1; i > 0; --i) {
            fscanf(_fs, "%d, ", &(grid_.dx1[i]));
        }
        fscanf(_fs, "%d}}, ", &(grid_.dx1[i]));
        color_.pscanf(_fs);
        fscanf(_fs, ", %d}\n", &(region_n_));
        return;
    }

    friend std::ostream & operator<<(std::ostream & _fs, Region_Info<N_RANK> const & _r) {
        int i;
        _fs << "{ BASE, ";
        _fs << "t = {" << _r.t0_ << ", " << _r.t1_ << "}, {";

        _fs << "x0 = {";
        for (i = N_RANK-1; i > 0; --i) {
            _fs << _r.grid_.x0[i] << ", ";
        }
        _fs << _r.grid_.x0[i] << "}, ";

        _fs << "x1 = {";
        for (i = N_RANK-1; i > 0; --i) {
            _fs << _r.grid_.x1[i] << ", ";
        }
        _fs << _r.grid_.x1[i] << "}, ";

        _fs << "dx0 = {";
        for (i = N_RANK-1; i > 0; --i) {
            _fs << _r.grid_.dx0[i] << ", ";
        }
        _fs << _r.grid_.dx0[i] << "}, ";

        _fs << "dx1 = {";
        for (i = N_RANK-1; i > 0; --i) {
            _fs << _r.grid_.dx1[i] << ", ";
        }
        _fs << _r.grid_.dx1[i] << "}}, ";
        _fs << _r.color_ << ", "; 
        _fs << _r.region_n_ << "}";
        return _fs;
    }
};

template <typename T>
bool is_basic_data_type() { return false; };
template <> bool is_basic_data_type<char>() { return true; };
template <> bool is_basic_data_type<int>() { return true; };
template <> bool is_basic_data_type<long>() { return true; };
template <> bool is_basic_data_type<float>() { return true; };
template <> bool is_basic_data_type<double>() { return true; };
template <> bool is_basic_data_type< Homogeneity >() { return true; }
template <> bool is_basic_data_type<char*>() { return true; };
template <> bool is_basic_data_type<int*>() { return true; };
template <> bool is_basic_data_type<long*>() { return true; };
template <> bool is_basic_data_type<float*>() { return true; };
template <> bool is_basic_data_type<double*>() { return true; };
template <> bool is_basic_data_type<Homogeneity*>() { return true; }

typedef double T_measure;

template <typename T>
struct Vector_Info {
    T * region_;
    /* the measure_ is the unit for sorting the elements stored in the Vector */
    T_measure * measure_;
    int pointer_, size_;

    typedef boost::integral_constant<bool, ::boost::is_pointer<T>::value> region_type_;
    typedef boost::integral_constant<bool, ::boost::is_pointer<T_measure>::value> measure_type_;

    template <typename TT>
    inline TT * setup_region(int _size, const boost::false_type&) {
        /* set up the region which is a value, rather than a pointer */
        TT * l_tt = new TT[_size];
        return l_tt;
    }
    template <typename TT>
    inline TT * setup_region(int _size, const boost::true_type&) {
        /* set up the region which is a pointer 
         * - so l_tt is a pointer to pointer
         * - an array of pointer
         */
        TT * l_tt = new TT[_size];
        for (int i = 0; i < _size; ++i) {
            l_tt[i] = NULL;
        }
        return l_tt;
    }
    template <typename TT>
    inline int release_region(TT * _region, int _size, bool _is_basic_type, const boost::false_type&) {
        /* release region containing only values */
        del_arr(_region);
        return 0;
    }
    template <typename TT>
    inline int release_region(TT * _region, int _size, bool _is_basic_type, const boost::true_type&) {
        /* release region containing pointers: 
         * this method also release the memory region each entry pointer
         * points to
         */
        del_arr(_region);
        return 0;
    }
    Vector_Info() {
        region_ = setup_region<T>(VECTOR_SIZE, region_type_());
        measure_ = setup_region<T_measure>(VECTOR_SIZE, measure_type_());
        for (int i = 0; i < VECTOR_SIZE; ++i)
            measure_[i] = 0;
        pointer_ = 0; size_ = VECTOR_SIZE;
#if DEBUG
        printf("init size = %d\n", size_);
#endif
    }
    Vector_Info(int size) {
        region_ = setup_region<T>(size, region_type_());
        measure_ = setup_region<T_measure>(size, measure_type_());
        for (int i = 0; i < size; ++i)
            measure_[i] = 0;
        pointer_ = 0; size_ = size;
#if DEBUG
        printf("init size = %d\n", size_);
#endif
    }
    Vector_Info(Vector_Info<T> const & rhs) {
        int l_rhs_size = rhs.size();
        region_ = setup_region<T>(l_rhs_size, region_type_());
        measure_ = setup_region<T_measure>(l_rhs_size, measure_type_());
        size_ = l_rhs_size;
        for (int i = 0; i < l_rhs_size; ++i) {
            region_[i] = rhs[i];
            measure_[i] = rhs.measure_[i];
        }
        pointer_ = l_rhs_size;
    }

    ~Vector_Info() {
        /* how to write a destructor if the region_[] is an array of function objects? 
         */
        LOG(0, "Call destructor of Vector_Info"); 
        release_region<T>(region_, pointer_, is_basic_data_type<T>(), region_type_());
        release_region<T_measure>(measure_, pointer_, is_basic_data_type<T_measure>(), measure_type_());
        pointer_ = size_ = 0;
    } 
    inline void sort(void) {
        /* bubble sort */
        for (int i = 0; i < pointer_ - 1; ++i) {
            T_measure l_max = measure_[i];
            int l_max_idx = i;
            for (int j = i; j < pointer_ - 1; ++j) {
                if (measure_[j] > l_max) {
                    l_max_idx = j;
                    l_max = measure_[j];
                }
            }
            if (l_max_idx != i) {
                /* exchange the elements 'i' and 'l_max_idx' */
                T_measure tmp = measure_[i];
                measure_[i] = measure_[l_max_idx];
                measure_[l_max_idx] = tmp;
                T r_tmp = region_[i];
                region_[i] = region_[l_max_idx];
                region_[l_max_idx] = r_tmp;
            }
        }
    }
    /* push_back() without measurement */
    inline void push_back(T const & ele) {
        /* by push_back, the element added into the vector will always
         * occupy a seperate slot in the vector regardless if it's a duplicate
         * of existing elements
         */
#if DEBUG
        std::cerr << "push_back " << ele << std::endl;
#endif
        if (pointer_ < size_) {
            region_[pointer_] = ele;
            ++pointer_;
        } else {
#if DEBUG
            printf("realloc memory size = %d -> %d!\n", size_, 2*size_);
#endif
            T * l_region = setup_region<T>(2 * size_, region_type_());
            T_measure * l_measure = setup_region<T_measure>(2 * size_, measure_type_());
            if (l_region != NULL && l_measure != NULL) {
                for (int i = 0; i < size_; ++i) {
                    l_region[i] = region_[i];
                }
                release_region<T>(region_, pointer_, is_basic_data_type<T>(), region_type_());
                release_region<T_measure>(measure_, pointer_, is_basic_data_type<T_measure>(), measure_type_());
                region_ = l_region;
                measure_ = l_measure;
                region_[pointer_] = ele;
                ++pointer_;
                size_ = 2 * size_;
            } else {
                ERROR("realloc wrong!");
            }
        }
        return;
    }
    inline void push_back_unique(T const & ele) {
        /* by push_back_unique, it's guaranteed that the element added
         * into the vector will be union'ed against the existing elements
         */
#if DEBUG
        std::cerr << "push_back_unique " << ele << std::endl;
#endif
        /* to make sure every element in this vector is unique */
        for (int i = 0; i < pointer_; ++i) {
            if (region_[i] == ele)
                return;
        }
        push_back(ele); 
        return;
    }
    /* push_back() with measurement */
    inline void push_back(T const & ele, int rec_level) {
        /* by push_back, the element added into the vector will always
         * occupy a seperate slot in the vector regardless if it's a duplicate
         * of existing elements
         */
        T_measure l_inc = (T_measure)1/(0x1 << rec_level);
#if DEBUG
        std::cerr << "push_back " << ele << std::endl;
#endif
        if (pointer_ < size_) {
            region_[pointer_] = ele;
            measure_[pointer_] = l_inc;
            ++pointer_;
        } else {
#if DEBUG
            LOG_ARGS("realloc memory size = %d -> %d!\n", size_, 2 * size_);
#endif
            T * l_region = setup_region<T>(2 * size_, region_type_());
            T_measure * l_measure = setup_region<T_measure>(2 * size_, measure_type_());
            if (l_region != NULL && l_measure != NULL) {
                for (int i = 0; i < size_; ++i) {
                    l_region[i] = region_[i];
                    l_measure[i] = measure_[i];
                }
                release_region(region_, pointer_, is_basic_data_type<T>(), region_type_());
                release_region(measure_, pointer_, is_basic_data_type<T_measure>(), measure_type_());
                region_ = l_region;
                measure_ = l_measure;
                region_[pointer_] = ele;
                measure_[pointer_] = l_inc;
                ++pointer_;
                size_ = 2 * size_;
            } else {
                ERROR("realloc wrong!");
            }
        }
        return;
    }
    inline void push_back_unique(T const & ele, int rec_level) {
        /* by push_back_unique, it's guaranteed that the element added
         * into the vector will be union'ed against the existing elements
         */
        double l_inc = (double)1/(0x1 << rec_level);
#if DEBUG
        std::cerr << "push_back_unique " << ele << std::endl;
#endif
        /* to make sure every element in this vector is unique */
        for (int i = 0; i < pointer_; ++i) {
            if (region_[i] == ele) {
                measure_[i] += l_inc;
                return;
            }
        }
        push_back(ele, rec_level);
        return;
    }
    inline int get_index(T const & ele) {
        for (int i = 0; i < pointer_; ++i) {
            if (region_[i] == ele)
                return i;
        }
        return -1;
    }
    inline int get_largest_lb (T const & ele) {
        T const & l_orig = ele;
        /* l_subsume is initialized to white_clone */
        T l_subsume(ele.size());
        /* l_idx is initialized to the idx of white_clone */
        int l_idx = get_index(l_subsume);
        for (int i = 0; i < pointer_; ++i) {
            T & l_region = region_[i];
            if (l_region <= l_orig && l_subsume <= l_region) {
                l_idx = i;
                l_subsume = l_region;
            }
        }
        return l_idx;
    }
    inline void scan () {
        for (int i = 1; i < pointer_; ++i) {
            region_[i] = region_[i] + region_[i-1];
        }
    }
    T * get_root() { return region_; }
    inline T & operator[] (int _idx) { return region_[_idx]; }
    inline T_measure measure(int _idx) { return measure_[_idx]; }
    int size() { return pointer_; }
    int set_size(int _size) { 
        /* So, don't Over-size the vector */
        if (_size < pointer_) {
            pointer_ = _size;
        }
        return pointer_;
    }
    T & operator= (T const & rhs) {
        const int l_rhs_size = rhs.size();
        if (l_rhs_size <= size_) {
            for (int i = 0; i < l_rhs_size; ++i) {
                region_[i] = rhs[i];
                measure_[i] = rhs.measure(i);
            }
            pointer_ = l_rhs_size;
            return (*this);
        }
        /* l_rhs_size > size_ */
        T * l_region = setup_region(l_rhs_size, region_type_());
        T_measure * l_measure = setup_region(l_rhs_size, measure_type_());
        if (l_region == NULL || l_measure == NULL) {
            ERROR("Run out of memory!");
        }
        for (int i = 0; i < l_rhs_size; ++i) {
            l_region[i] = rhs[i];
            l_measure[i] = rhs.measure(i);
        }
        release_region(region_, pointer_, is_basic_data_type<T>(), region_type_());
        release_region(measure_, pointer_, is_basic_data_type<T_measure>(), measure_type_());
        region_ = l_region;
        measure_ = l_measure;
        size_ = l_rhs_size;
        pointer_ = l_rhs_size;
        return (*this);
    }

    friend std::ofstream & operator<<(std::ofstream & fs, Vector_Info<T> const & v) {
        for (int i = 0; i < v.pointer_; ++i) {
            fs << v.region_[i] << std::endl;
        }
        return fs;
    }
    friend std::ostream & operator<<(std::ostream & fs, Vector_Info<T> const & v) {
        for (int i = 0; i < v.pointer_; ++i) {
            fs << v.region_[i] << std::endl;
        }
        return fs;
    }
};

template <int N_RANK>
struct Pochoir_Plan {
    Vector_Info< Region_Info<N_RANK> > * base_data_;
    Vector_Info<int> * sync_data_;
    int sz_base_data_, sz_sync_data_;
    int order_num_;
    char * fname_;
    DynamicLoader * dloader_;
    Pochoir_Plan (int _sz_base_data, int _sz_sync_data) : sz_base_data_(_sz_base_data), sz_sync_data_(_sz_sync_data) {
        base_data_ = new Vector_Info< Region_Info<N_RANK> >(_sz_base_data);
        sync_data_ = new Vector_Info<int>(_sz_sync_data);
        order_num_ = -1;
        fname_ = NULL; 
        dloader_ = NULL;
    }
    Pochoir_Plan () {
        sz_base_data_ = sz_sync_data_ = 0;
        order_num_ = -1;
        fname_ = NULL; 
        dloader_ = NULL;
        base_data_ = NULL;
        sync_data_ = NULL;
    }
    ~Pochoir_Plan() {
        sz_base_data_ = sz_sync_data_ = 0;
        order_num_ = -1;
        del_ele(base_data_);
        del_ele(sync_data_);
        del_ele(dloader_);
        del_arr(fname_);
    }
    void change_region_n(Vector_Info< Homogeneity > & color_vectors) {
        int const l_size = base_data_->size();
        for (int i = 0; i < l_size; ++i) {
            Region_Info<N_RANK> & l_region = (*base_data_)[i];
            Homogeneity & l_homo = l_region.color_;
            int l_idx = color_vectors.get_largest_lb(l_homo);
            l_region.region_n_ = l_idx;
            l_region.color_ = color_vectors[l_idx];
        }
    }
    void alloc_base_data(int _sz_base_data) {
        sz_base_data_ = _sz_base_data;
        base_data_ = new Vector_Info< Region_Info<N_RANK> >(_sz_base_data);
    }
    void alloc_sync_data(int _sz_sync_data) {
        sz_sync_data_ = _sz_sync_data;
        // sync_data_ = new Vector_Info<int>(_sz_sync_data);
        sync_data_ = new Vector_Info<int>;
    }
    void set_order_num(int _order_num) { order_num_ = _order_num; }
    int get_order_num(void) { return order_num_; }
    void set_fname(const char * _fname) {
        // fname_ = new char[strlen(_fname)+2];
        fname_ = new char[120];
        strcpy(fname_, _fname);
    }
    char * get_fname(void) {
        return fname_;
    }
    template <typename T_Pochoir, typename T_Array>
    int load_kernels(T_Pochoir & _pochoir, T_Array & _a) {
        if (dloader_ != NULL) {
            WARNING("dloader != NULL!");
            return 0;
        }
        LOG(0, "<DLoader> starts loading!\n");
        /***************************************************************************************/
        char gen_kernel_fname [120];
        sprintf(gen_kernel_fname, "./%s_%d_gen_kernel", fname_, order_num_);
        LOG_ARGS(0, "gen_kernel_fname = %s\n", gen_kernel_fname);
        dloader_ = new DynamicLoader(gen_kernel_fname);
#if 1
        typedef int (*T_Create_Lambda)(T_Pochoir &, T_Array &);
        typedef int (*T_Register_Lambda)(T_Pochoir &);
        T_Create_Lambda create_lambdas = dloader_->load < T_Create_Lambda > ("Create_Lambdas");
        T_Register_Lambda register_lambdas = dloader_->load < T_Register_Lambda > ("Register_Lambdas");
#else
        // std::function < int (T_Pochoir &, T_Array &) > create_lambdas = dloader_->load < int (T_Pochoir &, T_Array &) > ("Create_Lambdas");
        // std::function < int (T_Pochoir &) > register_lambdas = dloader_->load < int (T_Pochoir &) > ("Register_Lambdas");
#endif

        create_lambdas(_pochoir, _a);
        register_lambdas(_pochoir);
        /***************************************************************************************/
        LOG(0, "<DLoader> ends loading!\n");
        return 0;
    }
    template <typename T_Pochoir, typename T_Array, typename ... T_ArrayS>
    int load_kernels(T_Pochoir & _pochoir, T_Array & _a, T_ArrayS ... _as) {
        if (dloader_ != NULL) {
            WARNING("dloader != NULL!");
            return 0;
        }
        LOG(0, "<DLoader> starts loading!\n");
        /***************************************************************************************/
        char gen_kernel_fname [strlen(fname_) + 20];
        sprintf(gen_kernel_fname, "./%s_%d_gen_kernel", fname_, order_num_);
        LOG_ARGS(0, "gen_kernel_fname = %s\n", gen_kernel_fname);
        dloader_ = new DynamicLoader(gen_kernel_fname);
#if 1
        typedef int (*T_Create_Lambda)(T_Pochoir &, T_Array &, T_ArrayS ...);
        typedef int (*T_Register_Lambda)(T_Pochoir &);
        T_Create_Lambda create_lambdas = dloader_->load < T_Create_Lambda > ("Create_Lambdas");
        T_Register_Lambda register_lambdas = dloader_->load < T_Register_Lambda > ("Register_Lambdas");
#else
        std::function < int (T_Pochoir &, T_Array &, T_ArrayS ...) > create_lambdas = dloader_->load < int (T_Pochoir &, T_Array &, T_ArrayS ...) > ("Create_Lambdas");
        std::function < int (T_Pochoir &) > register_lambdas = dloader_->load < int (T_Pochoir &) > ("Register_Lambdas");
#endif

        create_lambdas(_pochoir, _a, _as ...);
        register_lambdas(_pochoir);
        /***************************************************************************************/
        LOG(0, "<DLoader> ends loading!\n");
        return 0;
    }
    int unload_kernels(void) {
        LOG(0, "<DLoader> starts Deloading!\n");
        /***************************************************************************************/
        if (dloader_ == NULL) {
            WARNING("dloader == NULL");
            return 0;
        }
#if 1
        typedef int (*T_Destroy_Lambda)(void);
        T_Destroy_Lambda destroy_lambdas = dloader_->load < T_Destroy_Lambda > ("Destroy_Lambdas");
#else
        std::function < int (void) > destroy_lambdas = dloader_->load < int (void) > ("Destroy_Lambdas");
#endif
        destroy_lambdas();
        dloader_->close();
        del_ele(dloader_);
        /***************************************************************************************/
        LOG(0, "<DLoader> ends Deloading!\n");
        return 0;
    }
    void store_plan(const char * base_file_name, const char * sync_file_name) {
        /* extract data and store plan */
        std::ofstream os_base_data(base_file_name);
        std::ofstream os_sync_data(sync_file_name);
        if (os_base_data.is_open()) {
#if DEBUG
            printf("os_base_data is open!\n");
#endif
            os_base_data << "order_num = " << order_num_ << std::endl;
            os_base_data << "fname = " << fname_ << std::endl;
            os_base_data << (*base_data_);
        } else {
            printf("os_base_data is NOT open! exit!\n");
            exit(1);
        }
        if (os_sync_data.is_open()) {
#if DEBUG
            printf("os_sync_data is open!\n");
#endif
            os_sync_data << "order_num = " << order_num_ << std::endl;
            os_sync_data << "fname = " << fname_ << std::endl;
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
            del_ele(base_data_);
            sz_base_data_ = 0;
        }
        if (sz_sync_data_ != 0) {
            del_ele(sync_data_);
            sz_sync_data_ = 0;
        }
        alloc_base_data(10);
        alloc_sync_data(10);

        FILE * is_base_data = fopen(base_file_name, "r");
        int l_order_from_base_file = -1, l_order_from_sync_file = -1;
        char l_base_fname[FNAME_LENGTH], l_sync_fname[FNAME_LENGTH];

        if (is_base_data != NULL) {
#if DEBUG
            printf("is_base_data is open!\n");
#endif
            if (!feof(is_base_data)) {
                fscanf(is_base_data, "order_num = %d\n", &l_order_from_base_file);
            }
            if (!feof(is_base_data)) {
                fscanf(is_base_data, "fname = %s\n", l_base_fname);
            }
            while (!feof(is_base_data)) {
                Region_Info<N_RANK> l_region;
                l_region.pscanf(is_base_data);
                base_data_->push_back(l_region);
            }
            sz_base_data_ = base_data_->size();
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
            if (!feof(is_sync_data)) {
                fscanf(is_sync_data, "order_num = %d\n", &l_order_from_sync_file);
            }
            if (!feof(is_sync_data)) {
                fscanf(is_sync_data, "fname = %s\n", l_sync_fname);
            }
            while (!feof(is_sync_data)) {
                int l_sync;
                fscanf(is_sync_data, "%d\n", &l_sync);
                sync_data_->push_back(l_sync);
            }
            sz_sync_data_ = sync_data_->size();
        } else {
            printf("is_sync_data is NOT open! exit!\n");
            exit(1);
        }
        fclose(is_sync_data);
#if DEBUG
        std::cerr << "sync_data : \n" << (*sync_data_) << std::endl;
#endif

        if (l_order_from_base_file != l_order_from_sync_file) {
            ERROR_ARGS("order from base <%d>, order from sync <%d>", l_order_from_base_file, l_order_from_sync_file);
        } else {
            fprintf(stderr, "order from base = order from sync = %d\n", 
                    l_order_from_base_file);
            order_num_ = l_order_from_base_file;
        }
        if (strcmp(l_base_fname, l_sync_fname) != 0) {
            ERROR_ARGS("base_fname <%s>, sync_fname <%s>", l_base_fname, l_sync_fname);
        } else {
            fprintf(stderr, "base_fname = sync_fname = %s\n", l_base_fname);
            fname_ = new char[strlen(l_base_fname)];
            strcpy(fname_, l_base_fname);
        }
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

#if 0
    Node_Info(int _t0, int _t1, Grid_Info<N_RANK> & _grid) : region_(_t0, _t1, _grid, -1) {
        /* constructor */
        op = IS_INTERNAL; 
        parent = left = right = NULL;
    }
#endif

    Node_Info(int _t0, int _t1, Grid_Info<N_RANK> const & _grid, Homogeneity _color) : region_(_t0, _t1, _grid, _color, -1) {
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
    bool add_empty_region_;
    Spawn_Tree() { 
        /* constructor */
        root_ = new Node_Info<N_RANK>(IS_ROOT); size_ = 1;
        add_empty_region_ = false;
    }

    Spawn_Tree(int _t0, int _t1, Grid_Info<N_RANK> const & _grid, Homogeneity const & _color) {
        /* constructor */
        root_ = new Node_Info<N_RANK>(_t0, _t1, _grid, _color);
        size_ = 1;
        add_empty_region_ = false;
    }
    ~Spawn_Tree() {
        /* free the entire tree */
        add_empty_region_ = false;
        if (root_->left == NULL) {
            del_ele(root_);
            return;
        } else {
            dfs_rm_tree(root_);
            return;
        }
    }
    void set_add_empty_region(bool _add_empty_region) { 
        add_empty_region_ = _add_empty_region;
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

    void add_node(Node_Info<N_RANK> * parent, Node_Info<N_RANK> * child, enum Meta_Op _op, Homogeneity _color, int _region_n) {
        if (_op == IS_SPAWN && _region_n == NONE_EXCLUSIVE_IFS && !add_empty_region_) {
            /* we don't add the empty region into the tree */
            return;
        }
        if (parent->left == NULL) {
            parent->left = child;
            child->parent = parent;
            child->op = _op; 
            (child->region_).region_n_ = _region_n;
            (child->region_).color_ = _color;
            ++size_;
            return;
        }
        Node_Info<N_RANK> * youngest_child = parent->left;
        while (youngest_child->right != NULL)
            youngest_child = youngest_child->right;
        youngest_child->right = child;
        child->parent = parent;
        child->op = _op; 
        (child->region_).region_n_ = _region_n;
        (child->region_).color_ = _color;
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
        del_ele(node);
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
            assert(node->region_.region_n >= 0);
            l_node = node->right;
            base_data.push_back(node->region_);
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
