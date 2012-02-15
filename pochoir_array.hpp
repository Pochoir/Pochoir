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

template <int DIM>
inline int cal_index(int const * _idx, int const * _stride) {
	return (_idx[DIM] * _stride[DIM]) + cal_index<DIM-1>(_idx, _stride);
}

template <>
inline int cal_index<0>(int const * _idx, int const * _stride) {
	/* 0-dim is always the time dimension */
	return (_idx[0] * _stride[0]);
}

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

#ifdef CHECK_SHAPE
#undef CHECK_SHAPE
#endif

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
        Pochoir_Shape<N_RANK> * shape_;
        int shape_size_;
        typedef T (*BValue_1D)(Pochoir_Array<T, 1> &, int, int);
        typedef T (*BValue_2D)(Pochoir_Array<T, 2> &, int, int, int);
        typedef T (*BValue_3D)(Pochoir_Array<T, 3> &, int, int, int, int);
        typedef T (*BValue_4D)(Pochoir_Array<T, 4> &, int, int, int, int, int);
        typedef T (*BValue_5D)(Pochoir_Array<T, 5> &, int, int, int, int, int, int);
        typedef T (*BValue_6D)(Pochoir_Array<T, 6> &, int, int, int, int, int, int, int);
        typedef T (*BValue_7D)(Pochoir_Array<T, 7> &, int, int, int, int, int, int, int, int);
        typedef T (*BValue_8D)(Pochoir_Array<T, 8> &, int, int, int, int, int, int, int, int, int);
        BValue_1D bv1_;
        BValue_2D bv2_;
        BValue_3D bv3_;
        BValue_4D bv4_;
        BValue_5D bv5_;
        BValue_6D bv6_;
        BValue_7D bv7_;
        BValue_8D bv8_;
	public:
		/* create array with initial size 
         * - Following dimensions for constructors are spatial dimension
         * - all spatial dimensions are row-majored
         */
        explicit Pochoir_Array (int sz0) {
            logic_size_[0] = phys_size_[0] = sz0;
            logic_start_[0] = 0; logic_end_[0] = sz0;
            stride_[0] = 1; 
            total_size_ = sz0;
            shape_ = NULL;
            view_ = NULL;
            bv1_ = NULL; bv2_ = NULL; bv3_ = NULL;
            allocMemFlag_ = false;
//            view_ = new Storage<T>(TOGGLE * total_size_);
//            data_ = view_->data();
        }

		explicit Pochoir_Array (int sz1, int sz0) {
			logic_size_[1] = sz1; logic_size_[0] = sz0; 
			phys_size_[1] = sz1; phys_size_[0] = sz0; 
            logic_start_[0] = 0; logic_end_[0] = sz0;
            logic_start_[1] = 0; logic_end_[1] = sz1;
			stride_[1] = sz0; stride_[0] = 1; 
            shape_ = NULL;
			total_size_ = phys_size_[0] * phys_size_[1];
			view_ = NULL;
            bv1_ = NULL; bv2_ = NULL; bv3_ = NULL;
            allocMemFlag_ = false;
//			  view_ = new Storage<T>(TOGGLE * total_size_) ;
//            data_ = view_->data();
		}

		explicit Pochoir_Array (int sz2, int sz1, int sz0) {
			logic_size_[2] = sz2; logic_size_[1] = sz1; logic_size_[0] = sz0; 
			phys_size_[2] = sz2; phys_size_[1] = sz1; phys_size_[0] = sz0; 
            logic_start_[0] = 0; logic_end_[0] = sz0;
            logic_start_[1] = 0; logic_end_[1] = sz1;
            logic_start_[2] = 0; logic_end_[2] = sz2;
			stride_[0] = 1;  
			total_size_ = phys_size_[2];
            shape_ = NULL;
			for (int i = 0; i < 2; ++i) {
				total_size_ *= phys_size_[i];
				stride_[i+1] = stride_[i] * phys_size_[i];
			}
			view_ = NULL;
			/* double the total_size_ because we are using toggle array */
            bv1_ = NULL; bv2_ = NULL; bv3_ = NULL;
            allocMemFlag_ = false;
//  		  view_ = new Storage<T>(TOGGLE*total_size_) ;
//            data_ = view_->data();
		}

		explicit Pochoir_Array (int sz3, int sz2, int sz1, int sz0) {
			logic_size_[3] = sz3; logic_size_[2] = sz2; logic_size_[1] = sz1; logic_size_[0] = sz0; 
			phys_size_[3] = sz3; phys_size_[2] = sz2; phys_size_[1] = sz1; phys_size_[0] = sz0; 
            logic_start_[0] = 0; logic_end_[0] = sz0;
            logic_start_[1] = 0; logic_end_[1] = sz1;
            logic_start_[2] = 0; logic_end_[2] = sz2;
            logic_start_[3] = 0; logic_end_[3] = sz3;
			stride_[0] = 1;  
			total_size_ = phys_size_[3];
            shape_ = NULL;
			for (int i = 0; i < 3; ++i) {
				total_size_ *= phys_size_[i];
				stride_[i+1] = stride_[i] * phys_size_[i];
			}
			view_ = NULL;
			/* double the total_size_ because we are using toggle array */
            bv1_ = NULL; bv2_ = NULL; bv3_ = NULL;
            allocMemFlag_ = false;
//			  view_ = new Storage<T>(TOGGLE*total_size_) ;
//            data_ = view_->data();
		}

		explicit Pochoir_Array (int sz4, int sz3, int sz2, int sz1, int sz0) {
			logic_size_[4] = sz4; logic_size_[3] = sz3; logic_size_[2] = sz2; logic_size_[1] = sz1; logic_size_[0] = sz0; 
			phys_size_[4] = sz4; phys_size_[3] = sz3; phys_size_[2] = sz2; phys_size_[1] = sz1; phys_size_[0] = sz0; 
            logic_start_[0] = 0; logic_end_[0] = sz0;
            logic_start_[1] = 0; logic_end_[1] = sz1;
            logic_start_[2] = 0; logic_end_[2] = sz2;
            logic_start_[3] = 0; logic_end_[3] = sz3;
            logic_start_[4] = 0; logic_end_[4] = sz4;
			stride_[0] = 1;  
			total_size_ = phys_size_[4];
            shape_ = NULL;
			for (int i = 0; i < 4; ++i) {
				total_size_ *= phys_size_[i];
				stride_[i+1] = stride_[i] * phys_size_[i];
			}
			view_ = NULL;
			/* double the total_size_ because we are using toggle array */
            bv1_ = NULL; bv2_ = NULL; bv3_ = NULL; bv4_ = NULL;
            allocMemFlag_ = false;
//			  view_ = new Storage<T>(TOGGLE*total_size_) ;
//            data_ = view_->data();
		}

		explicit Pochoir_Array (int sz5, int sz4, int sz3, int sz2, int sz1, int sz0) {
			logic_size_[5] = sz5; logic_size_[4] = sz4; logic_size_[3] = sz3; logic_size_[2] = sz2; logic_size_[1] = sz1; logic_size_[0] = sz0; 
			phys_size_[5] = sz5; phys_size_[4] = sz4; phys_size_[3] = sz3; phys_size_[2] = sz2; phys_size_[1] = sz1; phys_size_[0] = sz0; 
            logic_start_[0] = 0; logic_end_[0] = sz0;
            logic_start_[1] = 0; logic_end_[1] = sz1;
            logic_start_[2] = 0; logic_end_[2] = sz2;
            logic_start_[3] = 0; logic_end_[3] = sz3;
            logic_start_[4] = 0; logic_end_[4] = sz4;
            logic_start_[5] = 0; logic_end_[5] = sz5;
			stride_[0] = 1;  
			total_size_ = phys_size_[5];
            shape_ = NULL;
			for (int i = 0; i < 5; ++i) {
				total_size_ *= phys_size_[i];
				stride_[i+1] = stride_[i] * phys_size_[i];
			}
			view_ = NULL;
			/* double the total_size_ because we are using toggle array */
            bv1_ = NULL; bv2_ = NULL; bv3_ = NULL; bv4_ = NULL; bv5_ = NULL;
            allocMemFlag_ = false;
//			  view_ = new Storage<T>(TOGGLE*total_size_) ;
//            data_ = view_->data();
		}

		explicit Pochoir_Array (int sz6, int sz5, int sz4, int sz3, int sz2, int sz1, int sz0) {
			logic_size_[6] = sz6; logic_size_[5] = sz5; logic_size_[4] = sz4; logic_size_[3] = sz3; logic_size_[2] = sz2; logic_size_[1] = sz1; logic_size_[0] = sz0; 
			phys_size_[6] = sz6; phys_size_[5] = sz5; phys_size_[4] = sz4; phys_size_[3] = sz3; phys_size_[2] = sz2; phys_size_[1] = sz1; phys_size_[0] = sz0; 
            logic_start_[0] = 0; logic_end_[0] = sz0;
            logic_start_[1] = 0; logic_end_[1] = sz1;
            logic_start_[2] = 0; logic_end_[2] = sz2;
            logic_start_[3] = 0; logic_end_[3] = sz3;
            logic_start_[4] = 0; logic_end_[4] = sz4;
            logic_start_[5] = 0; logic_end_[5] = sz5;
            logic_start_[6] = 0; logic_end_[6] = sz6;
			stride_[0] = 1;  
			total_size_ = phys_size_[6];
            shape_ = NULL;
			for (int i = 0; i < 6; ++i) {
				total_size_ *= phys_size_[i];
				stride_[i+1] = stride_[i] * phys_size_[i];
			}
			view_ = NULL;
			/* double the total_size_ because we are using toggle array */
            bv1_ = NULL; bv2_ = NULL; bv3_ = NULL; bv4_ = NULL; bv5_ = NULL; bv6_ = NULL;
            allocMemFlag_ = false;
//			  view_ = new Storage<T>(TOGGLE*total_size_) ;
//            data_ = view_->data();
		}

		explicit Pochoir_Array (int sz7, int sz6, int sz5, int sz4, int sz3, int sz2, int sz1, int sz0) {
			logic_size_[7] = sz7; logic_size_[6] = sz6; logic_size_[5] = sz5; logic_size_[4] = sz4; logic_size_[3] = sz3; logic_size_[2] = sz2; logic_size_[1] = sz1; logic_size_[0] = sz0; 
			phys_size_[7] = sz7; phys_size_[6] = sz6; phys_size_[5] = sz5; phys_size_[4] = sz4; phys_size_[3] = sz3; phys_size_[2] = sz2; phys_size_[1] = sz1; phys_size_[0] = sz0; 
            logic_start_[0] = 0; logic_end_[0] = sz0;
            logic_start_[1] = 0; logic_end_[1] = sz1;
            logic_start_[2] = 0; logic_end_[2] = sz2;
            logic_start_[3] = 0; logic_end_[3] = sz3;
            logic_start_[4] = 0; logic_end_[4] = sz4;
            logic_start_[5] = 0; logic_end_[5] = sz5;
            logic_start_[6] = 0; logic_end_[6] = sz6;
            logic_start_[7] = 0; logic_end_[7] = sz7;
			stride_[0] = 1;  
			total_size_ = phys_size_[7];
            shape_ = NULL;
			for (int i = 0; i < 7; ++i) {
				total_size_ *= phys_size_[i];
				stride_[i+1] = stride_[i] * phys_size_[i];
			}
			view_ = NULL;
			/* double the total_size_ because we are using toggle array */
            bv1_ = NULL; bv2_ = NULL; bv3_ = NULL; bv4_ = NULL; bv5_ = NULL; bv6_ = NULL; bv7_ = NULL;
            allocMemFlag_ = false;
//			  view_ = new Storage<T>(TOGGLE*total_size_) ;
//            data_ = view_->data();
		}

		/* Copy constructor -- create another view of the
		 * same array
		 */
		Pochoir_Array (Pochoir_Array<T, N_RANK> const & orig) {
			total_size_ = orig.total_size();
			for (int i = 0; i < N_RANK; ++i) {
				phys_size_[i] = orig.phys_size(i);
				logic_size_[i] = orig.logic_size(i);
				stride_[i] = orig.stride(i);
                logic_start_[i] = 0; logic_end_[i] = logic_size_[i];
			}
			view_ = NULL;
			view_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).view();
			view_->inc_ref();
            /* We also get the BValue function pointer from orig */
            bv1_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_1D(); 
            bv2_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_2D(); 
            bv3_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_3D(); 
            bv4_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_4D(); 
            bv5_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_5D(); 
            bv6_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_6D(); 
            bv7_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_7D(); 
            bv8_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_8D(); 
            data_ = view_->data();
            allocMemFlag_ = true;
            shape_ = NULL;
		}

        /* assignment operator for vector<> */
		Pochoir_Array<T, N_RANK> & operator= (Pochoir_Array<T, N_RANK> const & orig) {
			total_size_ = orig.total_size();
			for (int i = 0; i < N_RANK; ++i) {
				phys_size_[i] = orig.phys_size(i);
				logic_size_[i] = orig.logic_size(i);
				stride_[i] = orig.stride(i);
			}
			view_ = NULL;
			view_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).view();
			view_->inc_ref();
            /* We also get the BValue function pointer from orig */
            bv1_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_1D(); 
            bv2_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_2D(); 
            bv3_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_3D(); 
            bv4_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_4D(); 
            bv5_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_5D(); 
            bv6_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_6D(); 
            bv7_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_7D(); 
            bv8_ = const_cast<Pochoir_Array<T, N_RANK> &>(orig).bv_8D(); 
            data_ = view_->data();
            allocMemFlag_ = true;
            shape_ = NULL;
            return *this;
		}

		/* destructor : free memory */
		~Pochoir_Array() {
            if (view_ != NULL) {
                view_->dec_ref();
                del_ele(view_);
            }
        
            allocMemFlag_ = false;
#if 0
            del_arr(shape_);
#endif
		}

		inline Storage<T> * view() {
			return view_;
		}

        inline T * data() { return data_; }
        /* return the function pointer which generates the boundary value! */
        BValue_1D bv_1D(void) { return bv1_; }
        BValue_2D bv_2D(void) { return bv2_; }
        BValue_3D bv_3D(void) { return bv3_; }
        BValue_4D bv_4D(void) { return bv4_; }
        BValue_5D bv_5D(void) { return bv5_; }
        BValue_6D bv_6D(void) { return bv6_; }
        BValue_7D bv_7D(void) { return bv7_; }
        BValue_8D bv_8D(void) { return bv8_; }

        /* guarantee that only one version of boundary function is registered ! */
        void Register_Boundary(BValue_1D _bv1) { bv1_ = _bv1;  bv2_ = NULL; bv3_ = NULL; bv4_ = NULL; bv5_ = NULL; bv6_ = NULL; bv7_ = NULL; bv8_ = NULL;}
        void Register_Boundary(BValue_2D _bv2) { bv2_ = _bv2;  bv1_ = NULL; bv3_ = NULL; bv4_ = NULL; bv5_ = NULL; bv6_ = NULL; bv7_ = NULL; bv8_ = NULL;}
        void Register_Boundary(BValue_3D _bv3) { bv3_ = _bv3;  bv1_ = NULL; bv2_ = NULL; bv4_ = NULL; bv5_ = NULL; bv6_ = NULL; bv7_ = NULL; bv8_ = NULL;}
        void Register_Boundary(BValue_4D _bv4) { bv4_ = _bv4;  bv1_ = NULL; bv2_ = NULL; bv3_ = NULL; bv5_ = NULL; bv6_ = NULL; bv7_ = NULL; bv8_ = NULL;}
        void Register_Boundary(BValue_5D _bv5) { bv5_ = _bv5;  bv1_ = NULL; bv2_ = NULL; bv3_ = NULL; bv4_ = NULL; bv6_ = NULL; bv7_ = NULL; bv8_ = NULL;}
        void Register_Boundary(BValue_6D _bv6) { bv6_ = _bv6;  bv1_ = NULL; bv2_ = NULL; bv3_ = NULL; bv4_ = NULL; bv5_ = NULL; bv7_ = NULL; bv8_ = NULL;}
        void Register_Boundary(BValue_7D _bv7) { bv7_ = _bv7;  bv1_ = NULL; bv2_ = NULL; bv3_ = NULL; bv4_ = NULL; bv5_ = NULL; bv6_ = NULL; bv8_ = NULL;}
        void Register_Boundary(BValue_8D _bv8) { bv8_ = _bv8;  bv1_ = NULL; bv2_ = NULL; bv3_ = NULL; bv4_ = NULL; bv5_ = NULL; bv6_ = NULL; bv7_ = NULL;}

        void unRegister_Boundary(void) { bv1_ = NULL;  bv2_ = NULL; bv3_ = NULL; bv4_ = NULL; ; bv5_ = NULL; bv6_ = NULL; bv7_ = NULL; bv8_ = NULL; return; }

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
        void Register_Shape(Pochoir_Shape<N_RANK> * shape, int shape_size) {
            /* currently we just get the slope_[] and toggle_ out of the shape[] */
            int l_min_time_shift=0, l_max_time_shift=0, depth=0;
#if 0
            assert(shape_ == NULL);
            shape_ = new Pochoir_Shape<N_RANK>[shape_size];
#endif
            shape_size_ = shape_size;
            for (int r = 0; r < N_RANK; ++r) {
                slope_[r] = 0;
            }
            for (int i = 0; i < shape_size; ++i) {
                if (shape[i].shift[0] < l_min_time_shift)
                    l_min_time_shift = shape[i].shift[0];
                if (shape[i].shift[0] > l_max_time_shift)
                    l_max_time_shift = shape[i].shift[0];
#if 0
                for (int r = 0; r < N_RANK+1; ++r) {
                    shape_[i].shift[r] = shape[i].shift[r];
                }
#endif
            }
            depth = l_max_time_shift - l_min_time_shift;
            toggle_ = depth + 1;
            shape_ = shape;
            for (int i = 0; i < shape_size; ++i) {
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
#if 0
            assert(shape_ == NULL);
            shape_ = new Pochoir_Shape<N_RANK>[N_SIZE];
#endif
            shape_size_ = N_SIZE;
            for (int r = 0; r < N_RANK; ++r) {
                slope_[r] = 0;
            }
            for (int i = 0; i < N_SIZE; ++i) {
                if (shape[i].shift[0] < l_min_time_shift)
                    l_min_time_shift = shape[i].shift[0];
                if (shape[i].shift[0] > l_max_time_shift)
                    l_max_time_shift = shape[i].shift[0];
#if 0
                for (int r = 0; r < N_RANK+1; ++r) {
                    shape_[i].shift[r] = shape[i].shift[r];
                }
#endif
            }
            depth = l_max_time_shift - l_min_time_shift;
            toggle_ = depth + 1;
            shape_ = shape;
            for (int i = 0; i < N_SIZE; ++i) {
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

#if 0
        /* This function could be called directly from user's app to 
         * register a shape with Pochoir_Array
         */
        template <size_t N_SIZE1, size_t N_SIZE2>
        void Register_Shape(Pochoir_Shape<N_RANK> (& shape1)[N_SIZE1], Pochoir_Shape<N_RANK> (& shape2)[N_SIZE2]) {
            /* currently we just get the slope_[] and toggle_ out of the shape[] */
            int l_min_time_shift=0, l_max_time_shift=0, depth=0;
            assert(shape_ == NULL);
            shape_ = new Pochoir_Shape<N_RANK>[N_SIZE1+N_SIZE2];
            shape_size_ = N_SIZE1+N_SIZE2;
            int i;
            for (int r = 0; r < N_RANK; ++r) {
                slope_[r] = 0;
            }
            for (i = 0; i < N_SIZE1; ++i) {
                if (shape1[i].shift[0] < l_min_time_shift)
                    l_min_time_shift = shape1[i].shift[0];
                if (shape1[i].shift[0] > l_max_time_shift)
                    l_max_time_shift = shape1[i].shift[0];
                for (int r = 0; r < N_RANK+1; ++i) {
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
            toggle_ = depth + 1;
            for (i = 0; i < N_SIZE1+N_SIZE2; ++i) {
                for (int r = 0; r < N_RANK; ++r) {
                    slope_[r] = max(slope_[r], abs((int)ceil((float)shape_[i].shift[r+1]/(l_max_time_shift - shape_[i].shift[0]))));
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
#endif

        inline void print_shape(void) {
            printf("\nInput Pochoir_Shape<%d> = \n{", N_RANK);
            for (int i = 0; i < shape_size_-1; ++i) {
                printf("{");
                for (int r = 0; r < N_RANK; ++r) {
                    printf("%d, ", shape_[i].shift[r]);
                }
                for (int r = N_RANK; r < N_RANK+1; ++r) {
                    printf("%d", shape_[i].shift[r]);
                }
                printf("}, ");
            }

            for (int i = shape_size_-1; i < shape_size_; ++i) {
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

        inline bool check_shape(int const (& l_shift) [N_RANK+1]) {
            bool shape_match;
            int const l_home_time_cord = shape_[0].shift[0];
            for (int i = 0; i < shape_size_; ++i) {
                shape_match = true;
                for (int r = 0; shape_match && r < N_RANK+1; ++r) {
                    if (r == 0) {
                        if (shape_[i].shift[0] != l_shift[r]) {
                            shape_match = false;
                            break;
                        }
                    } else if (shape_[i].shift[r] != l_shift[r]) {
                        shape_match = false;
                        break;
                    }
                }
                if (shape_match)
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

        inline bool check_boundary(size_info const & _idx) const {
            bool touch_boundary = false;
            for (int i = 0; i < N_RANK; ++i) {
                touch_boundary |= (_idx[i] < logic_start_[i]
                                | _idx[i] >= logic_end_[i]);
            }
            return touch_boundary;
        }

#define check_boundary1(_idx1, _idx0) \
            (_idx0 < logic_start_[0] || _idx0 >= logic_end_[0]) \
        

#define check_boundary2(_idx2, _idx1, _idx0) \
            (_idx0 < logic_start_[0] || _idx0 >= logic_end_[0] \
                 || _idx1 < logic_start_[1] || _idx1 >= logic_end_[1])

#define check_boundary3(_idx3, _idx2, _idx1, _idx0) \
            (_idx0 < logic_start_[0] || _idx0 >= logic_end_[0] \
                 || _idx1 < logic_start_[1] || _idx1 >= logic_end_[1] \
                 || _idx2 < logic_start_[2] || _idx2 >= logic_end_[2])

#define check_boundary4(_idx4, _idx3, _idx2, _idx1, _idx0) \
            (_idx0 < logic_start_[0] || _idx0 >= logic_end_[0] \
                 || _idx1 < logic_start_[1] || _idx1 >= logic_end_[1] \
                 || _idx2 < logic_start_[2] || _idx2 >= logic_end_[2] \
                 || _idx3 < logic_start_[3] || _idx3 >= logic_end_[3])

#define check_boundary5(_idx5, _idx4, _idx3, _idx2, _idx1, _idx0) \
            (_idx0 < logic_start_[0] || _idx0 >= logic_end_[0] \
                 || _idx1 < logic_start_[1] || _idx1 >= logic_end_[1] \
                 || _idx2 < logic_start_[2] || _idx2 >= logic_end_[2] \
                 || _idx3 < logic_start_[3] || _idx3 >= logic_end_[3] \
                 || _idx4 < logic_start_[4] || _idx4 >= logic_end_[4])

#define check_boundary6(_idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0) \
            (_idx0 < logic_start_[0] || _idx0 >= logic_end_[0] \
                 || _idx1 < logic_start_[1] || _idx1 >= logic_end_[1] \
                 || _idx2 < logic_start_[2] || _idx2 >= logic_end_[2] \
                 || _idx3 < logic_start_[3] || _idx3 >= logic_end_[3] \
                 || _idx4 < logic_start_[4] || _idx4 >= logic_end_[4] \
                 || _idx5 < logic_start_[5] || _idx5 >= logic_end_[5])

#define check_boundary7(_idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0) \
            (_idx0 < logic_start_[0] || _idx0 >= logic_end_[0] \
                 || _idx1 < logic_start_[1] || _idx1 >= logic_end_[1] \
                 || _idx2 < logic_start_[2] || _idx2 >= logic_end_[2] \
                 || _idx3 < logic_start_[3] || _idx3 >= logic_end_[3] \
                 || _idx4 < logic_start_[4] || _idx4 >= logic_end_[4] \
                 || _idx5 < logic_start_[5] || _idx5 >= logic_end_[5] \
                 || _idx6 < logic_start_[6] || _idx5 >= logic_end_[6])

#define check_boundary8(_idx8, _idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0) \
            (_idx0 < logic_start_[0] || _idx0 >= logic_end_[0] \
                 || _idx1 < logic_start_[1] || _idx1 >= logic_end_[1] \
                 || _idx2 < logic_start_[2] || _idx2 >= logic_end_[2] \
                 || _idx3 < logic_start_[3] || _idx3 >= logic_end_[3] \
                 || _idx4 < logic_start_[4] || _idx4 >= logic_end_[4] \
                 || _idx5 < logic_start_[5] || _idx5 >= logic_end_[5] \
                 || _idx6 < logic_start_[6] || _idx5 >= logic_end_[6] \
                 || _idx7 < logic_start_[7] || _idx5 >= logic_end_[7])
        /* 
         * orig_value() is reserved for "ostream" : cout << Pochoir_Array
         */
        inline T orig_value (int _timestep, size_info & _idx) {
            bool l_boundary = check_boundary(_idx);
            bool set_boundary = false;
            T l_bvalue = 0;
            if (l_boundary && bv1_ != NULL) {
                l_bvalue = bv1_(*this, _timestep, _idx[0]);
                set_boundary = true;
            } else if (l_boundary && bv2_ != NULL) {
                l_bvalue = bv2_(*this, _timestep, _idx[1], _idx[0]);
                set_boundary = true;
            } else if (l_boundary && bv3_ != NULL) {
                l_bvalue = bv3_(*this, _timestep, _idx[2], _idx[1], _idx[0]);
                set_boundary = true;
            } else if (l_boundary && bv4_ != NULL) {
                l_bvalue = bv4_(*this, _timestep, _idx[3], _idx[2], _idx[1], _idx[0]);
                set_boundary = true;
            } else if (l_boundary && bv5_ != NULL) {
                l_bvalue = bv5_(*this, _timestep, _idx[4], _idx[3], _idx[2], _idx[1], _idx[0]);
                set_boundary = true;
            } else if (l_boundary && bv6_ != NULL) {
                l_bvalue = bv6_(*this, _timestep, _idx[5], _idx[4], _idx[3], _idx[2], _idx[1], _idx[0]);
                set_boundary = true;
            } else if (l_boundary && bv7_ != NULL) {
                l_bvalue = bv7_(*this, _timestep, _idx[6], _idx[5], _idx[4], _idx[3], _idx[2], _idx[1], _idx[0]);
                set_boundary = true;
            } else if (l_boundary && bv8_ != NULL) {
                l_bvalue = bv8_(*this, _timestep, _idx[7], _idx[6], _idx[5], _idx[4], _idx[3], _idx[2], _idx[1], _idx[0]);
                set_boundary = true;
            }

            /* the highest dimension is time dimension! */
            int l_idx = cal_index<N_RANK>(_idx, stride_) + (_timestep % toggle_) * total_size_;
            return (set_boundary) ? l_bvalue : (*view_)[l_idx];
        }

		/* index operator() for the format of a(i, j, k) 
         * - The highest dimension is always time dimension
         * - this is the uninterior version
         */

		inline Pochoir_Proxy<T> operator() (int _idx1, int _idx0) {
            if (!allocMemFlag_) {
                printf("Pochoir array access error:\n");
                printf("A Pochoir array is accessed without being registered with a Pochoir object.\n");
                exit(1);
            }
#ifdef CHECK_SHAPE
            if (inRun) {
                int l_shift[N_RANK+1];
                int const l_home_time_cord = shape_[0].shift[0];
                l_shift[0] = _idx1 - home_cell_[0];
                l_shift[1] = _idx0 - home_cell_[1];
                bool l_within_shape = check_shape(l_shift);
                if (!l_within_shape) {
                    printf("Pochoir off-shape access error:\n");
                    printf("Pochoir array index (%d, %d)\nShape index {%d, %d}\n",
                            _idx1, _idx0, l_shift[0], l_shift[1]);
                    print_shape();
                    exit(1);
                }
            }
#endif
            bool l_boundary = check_boundary1(_idx1, _idx0);
            bool set_boundary = (l_boundary && bv1_ != NULL);
            if (set_boundary) 
                return Pochoir_Proxy<T>(bv1_(*this, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + (_idx1 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> operator() (int _idx2, int _idx1, int _idx0) {
            if (!allocMemFlag_) {
                printf("Pochoir array access error:\n");
                printf("A Pochoir array is accessed without being registered with a Pochoir object.\n");
                exit(1);
            }
#ifdef CHECK_SHAPE
            if (inRun) {
                int l_shift[N_RANK+1];
                int const l_home_time_cord = shape_[0].shift[0];
                l_shift[0] = _idx2 - home_cell_[0];
                l_shift[1] = _idx1 - home_cell_[1]; l_shift[2] = _idx0 - home_cell_[2];
                bool l_within_shape = check_shape(l_shift);
                // printf("branch called!\n");
                /* very weird!!! This branch is never called in optimized run,
                 * but if comment in/out this branch will affect the correctness
                 * of periodic boundary condition!!!
                 */
                if (!l_within_shape) {
#if 1
                    printf("Pochoir off-shape access error:\n");
                    printf("Pochoir array index (%d, %d, %d)\nShape index {%d, %d, %d}\n", _idx2, _idx1, _idx0, l_shift[0], l_shift[1], l_shift[2]);
                    print_shape();
                    exit(1);
#endif
                }
            }
#endif
            bool l_boundary = check_boundary2(_idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv2_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv2_(*this, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + (_idx2 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> operator() (int _idx3, int _idx2, int _idx1, int _idx0) {
            if (!allocMemFlag_) {
                printf("Pochoir array access error:\n");
                printf("A Pochoir array is accessed without being registered with a Pochoir object.\n");
                exit(1);
            }
#ifdef CHECK_SHAPE
            if (inRun) {
                int l_shift[N_RANK+1];
                int const l_home_time_cord = shape_[0].shift[0];
                l_shift[0] = _idx3 - home_cell_[0];
                l_shift[1] = _idx2 - home_cell_[1]; l_shift[2] = _idx1 - home_cell_[2];
                l_shift[3] = _idx0 - home_cell_[3];
                bool l_within_shape = check_shape(l_shift);
                if (!l_within_shape) {
                    printf("Pochoir off-shape access error:\n");
                    printf("Pochoir array index (%d, %d, %d, %d)\nShape index {%d, %d, %d, %d}\n",
                            _idx3, _idx2, _idx1, _idx0,
                            l_shift[0], l_shift[1], l_shift[2], l_shift[3]);
                    print_shape();
                    exit(1);
                }
            }
#endif
            bool l_boundary = check_boundary3(_idx3, _idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv3_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv3_(*this, _idx3, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + (_idx3 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> operator() (int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            if (!allocMemFlag_) {
                printf("Pochoir array access error:\n");
                printf("A Pochoir array is accessed without being registered with a Pochoir object.\n");
                exit(1);
            }
#ifdef CHECK_SHAPE
            if (inRun) {
                int l_shift[N_RANK+1];
                int const l_home_time_cord = shape_[0].shift[0];
                l_shift[0] = _idx4 - home_cell_[0];
                l_shift[1] = _idx3 - home_cell_[1]; l_shift[2] = _idx2 - home_cell_[2];
                l_shift[3] = _idx1 - home_cell_[3]; l_shift[4] = _idx0 - home_cell_[4];
                bool l_within_shape = check_shape(l_shift);
                if (!l_within_shape) {
                    printf("Pochoir off-shape access error:\n");
                    printf("Pochoir array index (%d, %d, %d, %d, %d)\nShape index {%d, %d, %d, %d, %d}\n",
                            _idx4, _idx3, _idx2, _idx1, _idx0,
                            l_shift[0], l_shift[1], l_shift[2], l_shift[3], l_shift[4]);
                    print_shape();
                    exit(1);
                }
            }
#endif
            bool l_boundary = check_boundary4(_idx4, _idx3, _idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv4_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv4_(*this, _idx4, _idx3, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + (_idx4 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx); 
		}

		inline Pochoir_Proxy<T> operator() (int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            if (!allocMemFlag_) {
                printf("Pochoir array access error:\n");
                printf("A Pochoir array is accessed without being registered with a Pochoir object.\n");
                exit(1);
            }
#ifdef CHECK_SHAPE
            if (inRun) {
                int l_shift[N_RANK+1];
                int const l_home_time_cord = shape_[0].shift[0];
                l_shift[0] = _idx5 - home_cell_[0];
                l_shift[1] = _idx4 - home_cell_[1]; l_shift[2] = _idx3 - home_cell_[2];
                l_shift[3] = _idx2 - home_cell_[3]; l_shift[4] = _idx1 - home_cell_[4];
                l_shift[5] = _idx0 - home_cell_[5];
                bool l_within_shape = check_shape(l_shift);
                if (!l_within_shape) {
                    printf("Pochoir off-shape access error:\n");
                    printf("Pochoir array index (%d, %d, %d, %d, %d, %d)\nShape index {%d, %d, %d, %d, %d, %d}\n",
                            _idx5, _idx4, _idx3, _idx2, _idx1, _idx0,
                            l_shift[0], l_shift[1], l_shift[2], l_shift[3], l_shift[4], l_shift[5]);
                    print_shape();
                    exit(1);
                }
            }
#endif
            bool l_boundary = check_boundary5(_idx5, _idx4, _idx3, _idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv5_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv5_(*this, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + (_idx5 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> operator() (int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            if (!allocMemFlag_) {
                printf("Pochoir array access error:\n");
                printf("A Pochoir array is accessed without being registered with a Pochoir object.\n");
                exit(1);
            }
#ifdef CHECK_SHAPE
            if (inRun) {
                int l_shift[N_RANK+1];
                int const l_home_time_cord = shape_[0].shift[0];
                l_shift[0] = _idx6 - home_cell_[0];
                l_shift[1] = _idx5 - home_cell_[1]; l_shift[2] = _idx4 - home_cell_[2];
                l_shift[3] = _idx3 - home_cell_[3]; l_shift[4] = _idx2 - home_cell_[4];
                l_shift[5] = _idx1 - home_cell_[5]; l_shift[6] = _idx0 - home_cell_[6];
                bool l_within_shape = check_shape(l_shift);
                if (!l_within_shape) {
                    printf("Pochoir off-shape access error:\n");
                    printf("Pochoir array index (%d, %d, %d, %d, %d, %d, %d)\nShape index {%d, %d, %d, %d, %d, %d, %d}\n",
                            _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0,
                            l_shift[0], l_shift[1], l_shift[2], l_shift[3], l_shift[4], l_shift[5], l_shift[6]);
                    print_shape();
                    exit(1);
                }
            }
#endif
            bool l_boundary = check_boundary6(_idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv6_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv6_(*this, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + (_idx6 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> operator() (int _idx7, int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            if (!allocMemFlag_) {
                printf("Pochoir array access error:\n");
                printf("A Pochoir array is accessed without being registered with a Pochoir object.\n");
                exit(1);
            }
#ifdef CHECK_SHAPE
            if (inRun) {
                int l_shift[N_RANK+1];
                int const l_home_time_cord = shape_[0].shift[0];
                l_shift[0] = _idx7 - home_cell_[0];
                l_shift[1] = _idx6 - home_cell_[1]; l_shift[2] = _idx5 - home_cell_[2];
                l_shift[3] = _idx4 - home_cell_[3]; l_shift[4] = _idx3 - home_cell_[4];
                l_shift[5] = _idx2 - home_cell_[5]; l_shift[6] = _idx1 - home_cell_[6];
                l_shift[7] = _idx0 - home_cell_[7];
                bool l_within_shape = check_shape(l_shift);
                if (!l_within_shape) {
                    printf("Pochoir off-shape access error:\n");
                    printf("Pochoir array index (%d, %d, %d, %d, %d, %d, %d, %d)\nShape index {%d, %d, %d, %d, %d, %d, %d, %d}\n",
                            _idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0,
                            l_shift[0], l_shift[1], l_shift[2], l_shift[3], l_shift[4], l_shift[5], l_shift[6], l_shift[7]);
                    print_shape();
                    exit(1);
                }
            }
#endif
            bool l_boundary = check_boundary7(_idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv7_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv7_(*this, _idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + _idx6 * stride_[6] + (_idx7 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> operator() (int _idx8, int _idx7, int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            if (!allocMemFlag_) {
                printf("Pochoir array access error:\n");
                printf("A Pochoir array is accessed without being registered with a Pochoir object.\n");
                exit(1);
            }
#ifdef CHECK_SHAPE
            if (inRun) {
                int l_shift[N_RANK+1];
                int const l_home_time_cord = shape_[0].shift[0];
                l_shift[0] = _idx8 - home_cell_[0];
                l_shift[1] = _idx7 - home_cell_[1]; l_shift[2] = _idx6 - home_cell_[2];
                l_shift[3] = _idx5 - home_cell_[3]; l_shift[4] = _idx4 - home_cell_[4];
                l_shift[5] = _idx3 - home_cell_[5]; l_shift[6] = _idx2 - home_cell_[6];
                l_shift[7] = _idx1 - home_cell_[7]; l_shift[8] = _idx0 - home_cell_[8];
                bool l_within_shape = check_shape(l_shift);
                if (!l_within_shape) {
                    printf("Pochoir off-shape access error:\n");
                    printf("Pochoir array index (%d, %d, %d, %d, %d, %d, %d, %d, %d)\nShape{%d, %d, %d, %d, %d, %d, %d, %d, %d}\n",
                            _idx8, _idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0,
                            l_shift[0], l_shift[1], l_shift[2], l_shift[3], l_shift[4], l_shift[5], l_shift[6], l_shift[7], l_shift[8]);
                    print_shape();
                    exit(1);
                }
            }
#endif
            bool l_boundary = check_boundary8(_idx8, _idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv8_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv8_(*this, _idx8, _idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + _idx6 * stride_[6] + _idx7 * stride_[7] + (_idx8 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

        /* set()/get() pair to set/get boundary value in user supplied bvalue function */
		inline T & set (int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + (_idx1 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & set (int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + (_idx2 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & set (int _idx3, int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + (_idx3 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & set (int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + (_idx4 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & set (int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + (_idx5 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & set (int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + (_idx6 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & set (int _idx7, int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + _idx6 * stride_[6] + (_idx7 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & set (int _idx8, int _idx7, int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + _idx6 * stride_[6] + _idx7 * stride_[7] + (_idx8 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T get (int _idx1, int _idx0) {
            if (check_boundary1(_idx1, _idx0)) {
                printf("Pochoir illegal access by boundary function error:\n");
                printf("Out-of-range access by boundary function at index (%d, %d)\n", _idx1, _idx0);
                exit(1);
            }
			int l_idx = _idx0 * stride_[0] + (_idx1 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T get (int _idx2, int _idx1, int _idx0) {
            if (check_boundary2(_idx2, _idx1, _idx0)) {
                printf("Pochoir illegal access by boundary function error:\n");
                printf("Out-of-range access by boundary function at index (%d, %d, %d)\n", _idx2, _idx1, _idx0);
                exit(1);
            }
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + (_idx2 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T get (int _idx3, int _idx2, int _idx1, int _idx0) {
            if (check_boundary3(_idx3, _idx2, _idx1, _idx0)) {
                printf("Pochoir illegal access by boundary function error:\n");
                printf("Out-of-range access by boundary function at index (%d, %d, %d, %d)\n", _idx3, _idx2, _idx1, _idx0);
                exit(1);
            }
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + (_idx3 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T get (int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            if (check_boundary4(_idx4, _idx3, _idx2, _idx1, _idx0)) {
                printf("Pochoir illegal access by boundary function error:\n");
                printf("Out-of-range access by boundary function at index (%d, %d, %d, %d, %d)\n", _idx4, _idx3, _idx2, _idx1, _idx0);
                exit(1);
            }
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + (_idx4 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T get (int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            if (check_boundary5(_idx5, _idx4, _idx3, _idx2, _idx1, _idx0)) {
                printf("Pochoir illegal access by boundary function error:\n");
                printf("Out-of-range access by boundary function at index (%d, %d, %d, %d, %d, %d)\n", _idx5, _idx4, _idx3, _idx2, _idx1, _idx0);
                exit(1);
            }
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + (_idx5 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T get (int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            if (check_boundary6(_idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0)) {
                printf("Pochoir illegal access by boundary function error:\n");
                printf("Out-of-range accesss by boundary function at index (%d, %d, %d, %d, %d, %d, %d)\n", _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0);
                exit(1);
            }
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + (_idx6 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T get (int _idx7, int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            if (check_boundary7(_idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0)) {
                printf("Pochoir illegal access by boundary function error:\n");
                printf("Out-of-range access by boundary function at index (%d, %d, %d, %d, %d, %d, %d, %d)\n", _idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0);
                exit(1);
            }
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + _idx6 * stride_[6] + (_idx7 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T get (int _idx8, int _idx7, int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            if (check_boundary8(_idx8, _idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0)) {
                printf("Pochoir illegal access by boundary function error:\n");
                printf("Out-of-range access by boundary function at index (%d, %d, %d, %d, %d, %d, %d, %d, %d)\n", _idx8, _idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0);
                exit(1);
            }
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + _idx6 * stride_[6] + _idx7 * stride_[7] + (_idx8 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		/* index operator() for the format of a.interior(i, j, k) 
         * - The highest dimension is always time dimension
         * - this is the interior (non-checking) version
         */

		inline T & interior (int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + (_idx1 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & interior (int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + (_idx2 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & interior (int _idx3, int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + (_idx3 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & interior (int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + (_idx4 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & interior (int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + (_idx5 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & interior (int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + (_idx6 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & interior (int _idx7, int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + _idx6 * stride_[6] + (_idx7 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline T & interior (int _idx8, int _idx7, int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + _idx6 * stride_[6] + _idx7 * stride_[7] + (_idx8 % toggle_) * total_size_;
			return (*view_)[l_idx];
		}

		inline Pochoir_Proxy<T> boundary (int _idx1, int _idx0) {
            bool l_boundary = check_boundary1(_idx1, _idx0);
            bool set_boundary = (l_boundary && bv1_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv1_(*this, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + (_idx1 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> boundary (int _idx2, int _idx1, int _idx0) {
            bool l_boundary = check_boundary2(_idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv2_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv2_(*this, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + (_idx2 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> boundary (int _idx3, int _idx2, int _idx1, int _idx0) {
            bool l_boundary = check_boundary3(_idx3, _idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv3_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv3_(*this, _idx3, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + (_idx3 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> boundary (int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            bool l_boundary = check_boundary4(_idx4, _idx3, _idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv4_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv4_(*this, _idx4, _idx3, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + (_idx4 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> boundary (int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            bool l_boundary = check_boundary5(_idx5, _idx4, _idx3, _idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv5_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv5_(*this, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + (_idx5 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> boundary (int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            bool l_boundary = check_boundary6(_idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv6_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv6_(*this, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + (_idx6 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> boundary (int _idx7, int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            bool l_boundary = check_boundary7(_idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv7_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv7_(*this, _idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + _idx6 * stride_[6] + (_idx7 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
		}

		inline Pochoir_Proxy<T> boundary (int _idx8, int _idx7, int _idx6, int _idx5, int _idx4, int _idx3, int _idx2, int _idx1, int _idx0) {
            bool l_boundary = check_boundary8(_idx8, _idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0);
            bool set_boundary = (l_boundary && bv8_ != NULL);
            if (set_boundary)
                return Pochoir_Proxy<T>(bv8_(*this, _idx8, _idx7, _idx6, _idx5, _idx4, _idx3, _idx2, _idx1, _idx0));
			int l_idx = _idx0 * stride_[0] + _idx1 * stride_[1] + _idx2 * stride_[2] + _idx3 * stride_[3] + _idx4 * stride_[4] + _idx5 * stride_[5] + _idx6 * stride_[6] + _idx7 * stride_[7] + (_idx8 % toggle_) * total_size_;
            return Pochoir_Proxy<T>(data_ + l_idx);
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
