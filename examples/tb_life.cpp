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

/* Test bench - 2D Game of Life, Periodic version */
#include <cstdio>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>

#include <pochoir.hpp>

using namespace std;
#define TIMES 1
//#define CHECK_RESULT

void check_result(int t, int j, int i, bool a, bool b)
{
	if (a == b) {
//		printf("a(%d, %d, %d) == b(%d, %d, %d) == %s : passed!\n", t, j, i, t, j, i, a ? "True" : "False");
	} else {
		printf("a(%d, %d, %d) = %s, b(%d, %d, %d) = %s : FAILED!\n", t, j, i, a ? "True" : "False", t, j, i, b ? "True" : "False");
	}

}

Pochoir_Boundary_2D(life_bv_2D, arr, t, i, j)
    const int arr_size_1 = arr.size(1);
    const int arr_size_0 = arr.size(0);

    int new_i = (i >= arr_size_1) ? (i - arr_size_1) : (i < 0 ? i + arr_size_1 : i);
    int new_j = (j >= arr_size_0) ? (j - arr_size_0) : (j < 0 ? j + arr_size_0 : j);

    /* we use arr.get(...) instead of arr(...) to implement different boundary
     * checking strategy: In arr(...), if off-boundary access occurs, we call
     * boundary function to supply a value; in arr.get(...), if off-boundary 
     * access occurs, we will print the off-boundary access and quit!
     */
    return arr.get(t, new_i, new_j);
    // return arr.get(t, new_i, -1);
Pochoir_Boundary_End

int main(int argc, char * argv[])
{
	int t;
	struct timeval start, end;
    //int N_SIZE = 0, T_SIZE = 0;
    int N1 = 500, N2 = 100, T_SIZE = 731;

    if (argc < 4) {
        printf("argc < 4, quit! \n");
        exit(1);
    }
    N1 = StrToInt(argv[1]);
    N2 = StrToInt(argv[2]);
    T_SIZE = StrToInt(argv[3]);
	
    /*N_SIZE = StrToInt(argv[1]);
    T_SIZE = StrToInt(argv[2]);
    printf("N_SIZE = %d, T_SIZE = %d\n", N_SIZE, T_SIZE);*/
    Pochoir_Shape_2D life_shape_2D[] = {{0, 0, 0}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 1}, {-1, 0, -1}, {-1, 1, 1}, {-1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}, {-1, 0, 0}};
    Pochoir_2D life_2D(life_shape_2D) ;
	Pochoir_Array_2D(bool) a(N1, N2) ;

    a.Register_Boundary(life_bv_2D);
#ifdef CHECK_RESULT
	Pochoir_Array_2D(bool) b(N1, N2) ;
    b.Register_Shape(life_shape_2D);
#endif

    life_2D.Register_Array(a);

#ifdef LIFE_BIT_TRICK
    Pochoir_2D bt_life_2D(life_shape_2D);
	Pochoir_Array_2D(bool) c(N1, N2);
    c.Register_Boundary(life_bv_2D);
    bt_life_2D.Register_Array(c);
#endif
	for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
		a(0, i, j) = (rand() & 0x1) ? true : false;
		a(1, i, j) = 0; 
#ifdef CHECK_RESULT
        b(0, i, j) = a(0, i, j);
        b(1, i, j) = 0;
#endif
#ifdef LIFE_BIT_TRICK
        c(0, i, j) = a(0, i, j);
        c(1, i, j) = 0;
#endif
	} }


    Pochoir_Kernel_2D(life_2D_fn, t, i, j)
    int neighbors = a(t-1, i-1, j-1) + a(t-1, i-1, j) + a(t-1, i-1, j+1) +
                    a(t-1, i, j-1)                  + a(t-1, i, j+1) +
                    a(t-1, i+1, j-1) + a(t-1, i+1, j) + a(t-1, i+1, j+1);
    if (a(t-1, i, j) == true && neighbors < 2)
        a(t, i, j) = false;
    else if (a(t-1, i, j) == true && neighbors > 3) { 
        a(t, i, j) = false;
    } else if (a(t-1, i, j) == true && (neighbors == 2 || neighbors == 3)) {
        a(t, i, j) = a(t-1, i, j);
    } else if (a(t-1, i, j) == false && neighbors == 3) 
        a(t, i, j) = true;
    else
        a(t, i, j) = a(t-1, i, j);
    Pochoir_Kernel_End

	char name [100] ;
	sprintf(name, "life") ;
	life_2D.set_problem_name(name) ;

	gettimeofday(&start, 0);
    for (int times = 0; times < TIMES; ++times) {
        life_2D.Run(T_SIZE, life_2D_fn);
    }
	gettimeofday(&end, 0);
	std::cout << "Pochoir : consumed time :" << 1.0e3 * tdiff(&end, &start)/TIMES << "ms" << std::endl;
    printf("N1 = %d, N2 = %d, T_SIZE = %d\n", N1, N2, T_SIZE);
    printf("Game of Life : %d x %d, %d time steps\n", N1, N2, T_SIZE);

#ifdef LIFE_BIT_TRICK
    Pochoir_Kernel_2D(bt_life_2D_fn, t, i, j)
    int neighbors = c(t-1, i-1, j-1) + c(t-1, i-1, j) + c(t-1, i-1, j+1) +
                    c(t-1, i, j-1)                  + c(t-1, i, j+1) +
                    c(t-1, i+1, j-1) + c(t-1, i+1, j) + c(t-1, i+1, j+1);
     int x = ((neighbors-2)>>31) | ~((neighbors-3)>>31); // x==0 iff neighbors==2 (otherwise -1)
     int y = ((neighbors-3)>>31) | ~((neighbors-4)>>31); // y==0 iff neighbors==3 (otherwise -1)
     //c(t,i,j) = 1&((~x&c(t-1,i,j))|~y);
     c(t,i,j) = 1&((~x&c(t-1,i,j)) | ~y);
    bool set0 = (c(t-1, i, j) == false && neighbors == 3); /* set true */
    bool set1 = (c(t-1, i, j) == true && neighbors < 2) 
             || (c(t-1, i, j) == true && neighbors > 3); /* set false */
    c(t, i, j) = set0 ? true : (set1 ? false : c(t-1, i, j));
    Pochoir_Kernel_End

	sprintf(name, "life_bit_trick") ;
	bt_life_2D.set_problem_name(name) ;

	gettimeofday(&start, 0);
    for (int times = 0; times < TIMES; ++times) {
        bt_life_2D.Run(T_SIZE, bt_life_2D_fn);
    }
	gettimeofday(&end, 0);
	std::cout << "Pochoir (Bit Trick): consumed time :" << 1.0e3 * tdiff(&end, &start)/TIMES << "ms" << std::endl;
#endif
//#if 1
#ifdef CHECK_RESULT
    b.Register_Boundary(life_bv_2D);
	gettimeofday(&start, 0);
    for (int times = 0; times < TIMES; ++times) {
	for (int t = 1; t < T_SIZE+1; ++t) {
	cilk_for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
        int neighbors = b(t-1, i-1, j-1) + b(t-1, i-1, j) + b(t-1, i-1, j+1) +
                        b(t-1, i, j-1)                  + b(t-1, i, j+1) +
                        b(t-1, i+1, j-1) + b(t-1, i+1, j) + b(t-1, i+1, j+1);
        if (b(t-1, i, j) == true && neighbors < 2)
            b(t, i, j) = false;
        else if (b(t-1, i, j) == true && neighbors > 3) { 
            b(t, i, j) = false;
        } else if (b(t-1, i, j) == true && (neighbors == 2 || neighbors == 3)) {
            b(t, i, j) = b(t-1, i, j);
        } else if (b(t-1, i, j) == false && neighbors == 3) 
            b(t, i, j) = true;
        else
            b(t, i, j) = b(t-1, i, j);
   	} } }
    }

	gettimeofday(&end, 0);
	std::cout << "Naive Loop: consumed time :" << 1.0e3 * tdiff(&end, &start) / TIMES << "ms" << std::endl;
#endif
//#if 1
#ifdef CHECK_RESULT
    printf("compare a with b : ");
	for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
		check_result(t, i, j, a.interior(t, i, j), b.interior(t, i, j));
	} } 
    printf("passed!\n");

#ifdef LIFE_BIT_TRICK
	t = T_SIZE;
    printf("compare a with c : ");
	for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
		check_result(t, i, j, a.interior(t, i, j), c.interior(t, i, j));
	} } 
    printf("passed!\n");
    printf("compare c with b : ");
	for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
		check_result(t, i, j, c.interior(t, i, j), b.interior(t, i, j));
	} } 
    printf("passed!\n");
#endif
#endif

	return 0;
}
