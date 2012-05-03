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

#ifndef POCHOIR_ARRAY_H
#define POCHOIR_ARRAY_H

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>

#include "pochoir_range.hpp"
#include "pochoir_common.hpp"
#include "pochoir_proxy.hpp"

using namespace std;

template <typename T, int N_RANK>
class Pochoir_Array;

template <typename T, int N_RANK>
struct Boundary {
    typedef T (*Types) (Pochoir_Array<T, N_RANK> &, int, int, int);
};

template <typename T>
struct Boundary<T, 1> {
    typedef T (*Types) (Pochoir_Array<T, 1> &, int, int);
};

template <typename T>
struct Boundary<T, 2> {
    typedef T (*Types) (Pochoir_Array<T, 2> &, int, int, int);
};

template <typename T>
struct Boundary<T, 3> {
    typedef T (*Types) (Pochoir_Array<T, 3> &, int, int, int, int);
};

template <typename T>
class Storage {
	private:
		T * storage_;
		int ref_;
	public:
		inline Storage(int _sz) {
#if 1
			storage_ = new T[_sz];
			ref_ = 1;
			for (int i = 0; i < _sz; ++i)
				storage_[i] = T();
#else
            storage_ = (T *)calloc(_sz, sizeof(T));
#endif
		}

		inline ~Storage() {
            del_arr(storage_);
		}

		inline void inc_ref() { 
			++ref_; 
		}

		inline void dec_ref() { 
			--ref_; 
		}

		inline int ref() { 
			return ref_; 
		}

		inline T & operator[] (int _idx) {
			return storage_[_idx];
		}

		inline T const & operator[] (int _idx) const {
			return storage_[_idx];
		}

		T * data() { return storage_; }
};

template <typename T, int N_RANK>
class Pochoir_Array {
	private:
		Storage<T> * view_; // real storage of elements
        T * data_; /* begining data pointer of view_, reserved for iterator! */
		typedef int size_info[N_RANK];
		size_info logic_size_; // logical of elements in each dimension
		size_info logic_start_, logic_end_; 
		size_info phys_size_; // physical of elements in each dimension
		size_info stride_; // stride of each dimension
        bool allocMemFlag_;
		int total_size_;
        int slope_[N_RANK], toggle_;
        Vector_Info< Pochoir_Shape<N_RANK> > shape_;
        int time_shift_;

        typedef typename Boundary<T, N_RANK>::Types BValue;
        BValue bv_;

        /* disable the copy constructor! */
		Pochoir_Array (Pochoir_Array<T, N_RANK> const & orig) { }
        /* disable the assignment operator! */
		Pochoir_Array<T, N_RANK> & operator= (Pochoir_Array<T, N_RANK> const & orig) { }

	public:
		/* create array with initial size 
         * - Following dimensions for constructors are spatial dimension
         * - all spatial dimensions are row-majored
         */

        template <typename I>
        inline void init (I sz) {
            logic_size_[0] = phys_size_[0] = sz;
            logic_start_[0] = 0; logic_end_[0] = sz;
            stride_[0] = 1; slope_[0] = 0; 
            toggle_ = 1; time_shift_ = 0;
            view_ = NULL; data_ = NULL;
            bv_ = NULL;
            total_size_ = 1; stride_[0] = 1;
            for (int i = 0; i < N_RANK; ++i) {
                total_size_ *= phys_size_[i];
                if (i < N_RANK - 1) {
                    stride_[i+1] = stride_[i] * phys_size_[i];
                }
            }
            allocMemFlag_ = false;
        }

        template <typename I, typename ... IS>
		inline void init (I sz, IS ... szs) {
            int l_dim = sizeof...(IS);
            logic_size_[l_dim] = sz; logic_start_[l_dim] = 0; logic_end_[l_dim] = sz;
            phys_size_[l_dim] = sz; 
            slope_[l_dim] = 0;
            init(szs...);
		}

        template <typename ... IS>
        explicit Pochoir_Array(IS ... szs) {
            init(szs ... );
        }

		/* Copy constructor -- create another view of the
		 * same array
         * -- let's try using the default copy constructor??
		 */
		/* destructor : free memory */
		~Pochoir_Array() {
            if (view_ != NULL) {
                view_->dec_ref();
                del_ele(view_);
            }
        
            allocMemFlag_ = false;
		}

		inline Storage<T> * view() { return view_; }

        inline T * data() { return data_; }
        /* return the function pointer which generates the boundary value! */
        BValue bv(void) { return bv_; }

        /* guarantee that only one version of boundary function is registered ! */
        void Register_Boundary(BValue _bv) { bv_ = _bv; }
        void unRegister_Boundary(void) { bv_ = NULL; }

        void Register_Domain(Grid_Info<N_RANK> initial_grid) {
            for (int i = 0; i < N_RANK; ++i) {
                logic_start_[i] = initial_grid.x0[i];
                logic_end_[i] = initial_grid.x1[i];
                logic_size_[i] = initial_grid.x1[i] - initial_grid.x0[i];
            }
        }

        /* This function will be called by Pochoir::Register_Array
         * from pochoir.hpp
         */
        void Register_Shape(Vector_Info< Pochoir_Shape<N_RANK> > & shape) {
            /* currently we just get the slope_[] and toggle_ out of the shape[] */
            int l_min_time_shift=0, l_max_time_shift=0, depth=0;
            for (int r = 0; r < N_RANK; ++r) {
                slope_[r] = 0;
            }
            int l_shape_size = shape.size();
            for (int i = 0; i < l_shape_size; ++i) {
                if (shape[i].shift[0] < l_min_time_shift)
                    l_min_time_shift = shape[i].shift[0];
                if (shape[i].shift[0] > l_max_time_shift)
                    l_max_time_shift = shape[i].shift[0];
            }
            depth = l_max_time_shift - l_min_time_shift;
            time_shift_ = max(time_shift_, 0 - l_min_time_shift);
            toggle_ = depth + 1;
            for (int i = 0; i < l_shape_size; ++i) {
                shape_.push_back_unique(shape[i]);
                for (int r = 0; r < N_RANK; ++r) {
                    slope_[r] = max(slope_[r], abs((int)ceil((float)shape[i].shift[r+1]/(l_max_time_shift - shape[i].shift[0]))));
                }
            }
#if DEBUG 
            printf("<%s> toggle = %d\n", __FUNCTION__, toggle_);
            for (int r = 0; r < N_RANK; ++r) {
                printf("<%s> slope[%d] = %d, ", __FUNCTION__, r, slope_[r]);
            }
            printf("\n");
#endif
            if (!allocMemFlag_) {
                alloc_mem();
            }
        }

        /* This function could be called directly from user's app to 
         * register a shape with Pochoir_Array
         */
        template <size_t N_SIZE>
        void Register_Shape(Pochoir_Shape<N_RANK> (& shape)[N_SIZE]) {
            /* currently we just get the slope_[] and toggle_ out of the shape[] */
            int l_min_time_shift=0, l_max_time_shift=0, depth=0;
            for (int r = 0; r < N_RANK; ++r) {
                slope_[r] = 0;
            }
            for (int i = 0; i < N_SIZE; ++i) {
                if (shape[i].shift[0] < l_min_time_shift)
                    l_min_time_shift = shape[i].shift[0];
                if (shape[i].shift[0] > l_max_time_shift)
                    l_max_time_shift = shape[i].shift[0];
            }
            depth = l_max_time_shift - l_min_time_shift;
            time_shift_ = max(time_shift_, 0 - l_min_time_shift);
            toggle_ = depth + 1;
            for (int i = 0; i < N_SIZE; ++i) {
                shape_.push_back_unique(shape[i]);
                for (int r = 0; r < N_RANK; ++r) {
                    slope_[r] = max(slope_[r], abs((int)ceil((float)shape[i].shift[r+1]/(l_max_time_shift - shape[i].shift[0]))));
                }
            }
#if DEBUG 
            printf("toggle = %d\n", toggle_);
            for (int r = 0; r < N_RANK; ++r) {
                printf("slope[%d] = %d, ", r, slope_[r]);
            }
            printf("\n");
#endif
            if (!allocMemFlag_) {
                alloc_mem();
            }
        }

        inline void print_shape(void) {
            printf("\nInput Pochoir_Shape<%d> = \n{", N_RANK);
            int l_shape_size = shape_.size();
            for (int i = 0; i < l_shape_size-1; ++i) {
                printf("{");
                for (int r = 0; r < N_RANK; ++r) {
                    printf("%d, ", shape_[i].shift[r]);
                }
                for (int r = N_RANK; r < N_RANK+1; ++r) {
                    printf("%d", shape_[i].shift[r]);
                }
                printf("}, ");
            }

            for (int i = l_shape_size-1; i < l_shape_size; ++i) {
                printf("{");
                for (int r = 0; r < N_RANK; ++r) {
                    printf("%d, ", shape_[i].shift[r]);
                }
                for (int r = N_RANK; r < N_RANK+1; ++r) {
                    printf("%d", shape_[i].shift[r]);
                }
                printf("}");
            }
            printf("}\n");
        }

        inline bool check_shape_shift(int const (& _shift) [N_RANK+1]) {
            bool l_match;
            int l_shape_size = shape_.size();
            for (int i = 0; i < l_shape_size; ++i) {
                l_match = true;
                for (int r = 0; l_match && r < N_RANK+1; ++r) {
                    if (r == 0) {
                        if (shape_[i].shift[0] != _shift[r]) {
                            l_match = false;
                            break;
                        }
                    } else if (shape_[i].shift[r] != _shift[r]) {
                        l_match = false;
                        break;
                    }
                }
                if (l_match)
                    return true;
            }
            return false;
        }

        void set_slope(int _slope[N_RANK]) { 
            for (int i = 0; i < N_RANK; ++i) 
                slope_[i] = _slope[i]; 
        }
        void set_toggle(int _toggle) { toggle_ = _toggle; }
        void alloc_mem(void) {
            if (!allocMemFlag_) {
                view_ = new Storage<T>(toggle_*total_size_) ;
                data_ = view_->data();
                allocMemFlag_ = true;
            }
        }
		/* return size */
		int phys_size(int _dim) const { return phys_size_[_dim]; }
		int logic_size(int _dim) const { return logic_size_[_dim]; }
        /* the size() function is for user's convenience! */
		int size(int _dim) const { return phys_size_[_dim]; }
		int slope(int _dim) const { return slope_[_dim]; }

		/* return total_size_ */
		int total_size() const { return total_size_; }

		/* return stride */
		int stride (int _dim) const { return stride_[_dim]; }

        template <typename I>
		inline int cal_addr (I _idx) {
			int l_idx = _idx * stride_[0];
			return l_idx;
		}

        template <typename I, typename ... IS>
        inline int cal_addr (I _idx, IS ... _idxs) {
            int l_dim = sizeof...(IS);
            return (_idx * stride_[l_dim] + cal_addr(_idxs ...));
        }

        template <typename I>
		inline bool check_bound (I _idx) {
            return (_idx < logic_start_[0] || _idx >= logic_end_[0]);
		}

        template <typename I, typename ... IS>
        inline bool check_bound (I _idx, IS ... _idxs) {
            int l_dim = sizeof...(IS);
            return (_idx < logic_start_[l_dim] || _idx >= logic_end_[l_dim] || check_bound(_idxs ...));
        }

        template <typename I>
		inline void print_idx (I _idx) {
            printf("%d", _idx);
		}

        template <typename I, typename ... IS>
        inline void print_idx (I _idx, IS ... _idxs) {
            printf("%d, ", _idx);
            print_idx(_idxs ...);
        }

        inline void print_index(int _sz, int * _idx) {
            for (int i = 0; i < _sz-1; ++i) {
                printf("%d, ", _idx[i]);
            }
            printf("%d", _idx[_sz-1]);
        }
        /* 
         * orig_value() is reserved for "ostream" : cout << Pochoir_Array
         */
        inline T orig_value (int _idx_t, size_info & _idx_s) {
            int l_idx = (_idx_t % toggle_) * total_size_;
            for (int i = 0; i < N_RANK; ++i) {
                l_idx += _idx_s[i] * stride_[i];
            }
            return (*view_)[l_idx];
        }

		/* index operator() for the format of a(i, j, k) 
         * - The highest dimension is always time dimension
         * - this is the uninterior version
         */
        template <typename I>
        void check_shape(int (&_shift)[N_RANK+1], int (&_index)[N_RANK+1], I _idx) {
            _shift[N_RANK] = _idx - home_cell_[N_RANK];
            _index[N_RANK] = _idx;
            bool l_within_shape = check_shape_shift(_shift);
            if (!l_within_shape) {
               printf("Pochoir off-shape access error:\n");
               printf("Pochoir array index (");
               print_index(N_RANK+1, _index);
#if DEBUG
               printf("), home_cell (");
               print_index(N_RANK+1, home_cell_);
#endif
               printf("), Shape shift index {");
               print_index(N_RANK+1, _shift);
               printf("}\n");
               print_shape();
               exit(1);
            }
        }
        template <typename I, typename ... IS>
        void check_shape(int (&_shift)[N_RANK+1], int (&_index)[N_RANK+1], I _idx, IS ... _idxs) {
            int l_dim = N_RANK - sizeof...(IS);
            _shift[l_dim] = _idx - home_cell_[l_dim];
            _index[l_dim] = _idx;
            check_shape(_shift, _index, _idxs ...);
        }

        template <typename ... IS>
		inline Pochoir_Proxy<T> operator() (int _idx_t, IS ... _idxs) {
            if (!allocMemFlag_) {
                printf("Pochoir array access error:\n");
                printf("A Pochoir array is accessed without being registered with a Pochoir object.\n");
                exit(1);
            }
#ifdef CHECK_SHAPE
            if (inRun) {
                int l_shift[N_RANK+1];
                int l_index[N_RANK+1];
                check_shape(l_shift, l_index, _idx_t, _idxs ...);
            }
#endif
            bool l_boundary = check_bound(_idxs ...);
            bool set_boundary = (l_boundary && bv_ != NULL);
            if (set_boundary) 
                return Pochoir_Proxy<T>(bv_(*this, _idx_t, _idxs...));
			int l_idx = cal_addr(_idxs ...) + (_idx_t % toggle_) * total_size_;
            return Pochoir_Proxy<T>(&((*view_)[l_idx]));
		}

        /* set()/get() pair to set/get boundary value in user supplied bvalue function */
        template <typename ... IS>
        inline T & set (int _idx_t, IS ... _idxs) {
            int l_idx = (_idx_t % toggle_) * total_size_ + cal_addr(_idx_t, _idxs ...);
            return (*view_)[l_idx];
        }

        template <typename ... IS>
        inline T get (int _idx_t, IS ... _idxs) {
            bool l_out_bound = check_bound(_idxs ...);
            if (l_out_bound) {
                printf("Pochoir illegal access by boundary function error:\n");
                printf("Out-of-range access by boundary function at index (");
                print_idx(_idx_t, _idxs ...);
                printf(")\n");
                exit(1);
            }
            int l_idx = (_idx_t % toggle_) * total_size_ + cal_addr(_idxs ...);
            return (*view_)[l_idx];
        }

		/* index operator() for the format of a.interior(i, j, k) 
         * - The highest dimension is always time dimension
         * - this is the interior (non-checking) version
         */
        template <typename ... IS>
        inline T & interior (int _idx_t, IS ... _idxs) {
            int l_idx = (_idx_t % toggle_) * total_size_ + cal_addr(_idxs ...);
            return (*view_)[l_idx];
        }

        template <typename ... IS>
		inline Pochoir_Proxy<T> boundary (int _idx_t, IS ... _idxs) {
            bool l_boundary = check_bound(_idxs ...);
            bool set_boundary = (l_boundary && bv_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv_(*this, _idx_t, _idxs ...));
			int l_idx = cal_addr(_idxs ...) + (_idx_t % toggle_) * total_size_;
            return Pochoir_Proxy<T>(&((*view_)[l_idx]));
		}

		/* size_info is of type int[] */
		static inline bool update_index(int * index, bool & line_break, int const * head_index, int const * tail_index)
		{
			int i = 0;
			bool done = false, whole_done = false;
			while (!done && i < N_RANK) {
				if (index[i] == (tail_index[i] - 1)) {
					index[i] = head_index[i];
					line_break = true;
					if (i == N_RANK-1)
						whole_done = true;
					i++;
				} else {
					index[i]++;
					done = true;
				}
			}
			return whole_done;
		}

#if 1
		template <typename T2, int N2>
		friend ostream& operator<<(ostream& os, Pochoir_Array<T2, N2> const & x); 
#endif
};

#if 1
template<typename T2, int N2>
ostream& operator<<(ostream& os, Pochoir_Array<T2, N2> const & x) { 
	typedef int size_info[N2];
	size_info l_index, l_head_index, l_tail_index;
	bool done = false, line_break = false;
	int i = 0;

	os << " Pochoir_Array : "; 
	for (int i = 0; i < N2; ++i) {
		l_index[i] = 0;
		l_head_index[i] = 0;
		l_tail_index[i] = x.phys_size(i);
		os << "Dim " << i << ", size<" << x.phys_size(i) << "> ; ";
	}
	os << endl;

	while (!done) {
		T2 x0, x1;
		x0 = const_cast<Pochoir_Array<T2, N2> &>(x).orig_value(0, l_index);
		x1 = const_cast<Pochoir_Array<T2, N2> &>(x).orig_value(1, l_index);
		os << setw(9) << x0 << " (" << x1 << ")" << " "; 
		done = const_cast<Pochoir_Array<T2, N2> &>(x).update_index(l_index, line_break, l_head_index, l_tail_index);
		if (line_break) {
			os << endl;
			line_break = false;
		}
	}
	return os; 
}
#endif
#endif // POCHOIR_ARRAY_H
