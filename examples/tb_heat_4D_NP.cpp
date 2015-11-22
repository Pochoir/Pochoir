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

/* Test bench - 4D heat equation, Non-periodic version */
#include <cstdio>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>

#include <pochoir.hpp>

using namespace std;
#define TIMES 1
#define N_RANK 4
#define TOLERANCE (1e-6)
//#define CHECK_RESULT

void check_result(int t, int i, int j, int k, int l, double a, double b)
{
	if (abs(a - b) < TOLERANCE) {
//		printf("a(%d, %d, %d, %d) == b(%d, %d, %d, %d) == %f : passed!\n", t, i, j, k, t, i, j, k, a);
	} else {
		printf("a(%d, %d, %d, %d, %d) = %f, b(%d, %d, %d, %d, %d) = %f : FAILED!\n", t, i, j, k, l, a, t, i, j, k, l, b);
	}

}

Pochoir_Boundary_4D(heat_bv_4D, arr, t, i, j, k, l)
    return 0;
Pochoir_Boundary_End

int main(int argc, char * argv[])
{
	const int BASE = 1024;
	int t;
	struct timeval start, end;
    double min_tdiff = INF;
    //int N_SIZE = 0, T_SIZE = 0;
    int N1 = 0, N2 = 0, N3 = 0, N4 = 0, T_SIZE = 0;

    if (argc < 6) {
        printf("argc < 6, quit! \n");
        exit(1);
    }
    //N_SIZE = StrToInt(argv[1]);
    //T_SIZE = StrToInt(argv[2]);
    N1 = StrToInt(argv[1]);
    N2 = StrToInt(argv[2]);
    N3 = StrToInt(argv[3]);
    N4 = StrToInt(argv[4]);
    T_SIZE = StrToInt(argv[5]);
    //printf("N_SIZE = %d, T_SIZE = %d\n", N_SIZE, T_SIZE);
    //Pochoir_Shape_4D heat_shape_4D[] = {{0, 0, 0, 0}, {-1, 1, 0, 0}, {-1, -1, 0, 0}, {-1, 0, 0, 0}, {-1, 0, 0, -1}, {-1, 0, 0, 1}, {-1, 0, 1, 0}, {-1, 0, -1, 0}};
    Pochoir_Shape_4D heat_shape_4D[] = {{0, 0, 0, 0, 0}, {-1, 1, 0, 0, 0}, {-1, -1, 0, 0, 0}, {-1, 0, 0, 0, 0}, {-1, 0, 0, -1, 0}, {-1, 0, 0, 1, 0}, {-1, 0, 1, 0, 0}, {-1, 0, -1, 0, 0}, {-1, 0, 0, 0, 1}, {-1, 0, 0, 0, -1}};
    Pochoir_4D heat_4D(heat_shape_4D);
	//Pochoir_Array_4D(double) a(N_SIZE, N_SIZE, N_SIZE, N_SIZE) ;
	Pochoir_Array_4D(double) a(N1, N2, N3, N4) ;
    //Pochoir_Domain I(1, N_SIZE-1), J(1, N_SIZE-1), K(1, N_SIZE-1), L(1, N_SIZE - 1);
    Pochoir_Domain I(1, N1-1), J(1, N2-1), K(1, N3-1), L(1, N4 - 1);
    heat_4D.Register_Array(a);
    a.Register_Boundary(heat_bv_4D);
#ifdef CHECK_RESULT
	//Pochoir_Array_4D(double) b(N_SIZE, N_SIZE, N_SIZE, N_SIZE);
	Pochoir_Array_4D(double) b(N1, N2, N3, N4) ;
    b.Register_Shape(heat_shape_4D);
    b.Register_Boundary(heat_bv_4D);
#endif

	for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
    for (int k = 0; k < N3; ++k) {
    for (int l = 0; l < N4; ++l) {
        a(0, i, j, k, l) = 1.0 * (rand() % BASE); 
        a(1, i, j, k, l) = 0; 
#ifdef CHECK_RESULT
        b(0, i, j, k, l) = a(0, i, j, k, l);
        b(1, i, j, k, l) = 0;
#endif
	} } } }

    Pochoir_Kernel_4D(heat_4D_fn, t, i, j, k, l)
	   a(t, i, j, k, l) = 
           0.125 * (a(t-1, i+1, j, k, l) - 2.0 * a(t-1, i, j, k, l) + a(t-1, i-1, j, k, l)) 
         + 0.125 * (a(t-1, i, j+1, k, l) - 2.0 * a(t-1, i, j, k, l) + a(t-1, i, j-1, k, l)) 
         + 0.125 * (a(t-1, i, j, k+1, l) - 2.0 * a(t-1, i, j, k, l) + a(t-1, i, j, k-1, l))
         + 0.125 * (a(t-1, i, j, k, l + 1) - 2.0 * a(t-1, i, j, k, l) + a(t-1, i, j, k, l - 1))
         + a(t-1, i, j, k, l);
    Pochoir_Kernel_End

    /* we have to bind arrayInUse and Shape together 
     * => One arrayInUse, one shape[] => One slope[]
     * because each arrayInUse needs to know the slope to determine
     * the boundary region and when to call the user supplied boundary
     * value function
     */
    //heat_4D.Register_Domain(I, J, K, L);

	char name [100] ;
	sprintf(name, "heat_4D_NP") ;
	heat_4D.set_problem_name(name) ;

    for (int times = 0; times < TIMES; ++times) {
	    gettimeofday(&start, 0);
        heat_4D.Run(T_SIZE, heat_4D_fn);
	    gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
	std::cout << "Pochoir ET: consumed time :" << min_tdiff << "ms" << std::endl;
    printf("N1 = %d, N2 = %d, N3 = %d, N4 = %d, T_SIZE = %d\n", N1, N2, N3, N4,
			T_SIZE);

//#if 1
#ifdef CHECK_RESULT
    min_tdiff = INF;
    /* cilk_for + zero-padding */
    for (int times = 0; times < TIMES; ++times) {
	gettimeofday(&start, 0);
	for (int t = 0; t < T_SIZE; ++t) {
    cilk_for (int i = 0; i < N1 ; ++i) {
	for (int j = 0; j < N2 ; ++j) {
    for (int k = 0; k < N3 ; ++k) {
    for (int l = 0; l < N4 ; ++l) {
	   b(t + 1, i, j, k, l) = 
           0.125 * (b(t, i+1, j, k, l) - 2.0 * b(t, i, j, k, l) + b(t, i-1, j, k, l)) 
         + 0.125 * (b(t, i, j+1, k, l) - 2.0 * b(t, i, j, k, l) + b(t, i, j-1, k, l)) 
         + 0.125 * (b(t, i, j, k+1, l) - 2.0 * b(t, i, j, k, l) + b(t, i, j, k-1, l))
         + 0.125 * (b(t, i, j, k, l+1) - 2.0 * b(t, i, j, k, l) + b(t, i, j, k, l-1))
         + b(t, i, j, k, l);
    } } } } }
	gettimeofday(&end, 0);
    min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
	std::cout << "Naive Loop: consumed time :" << min_tdiff << "ms" << std::endl;

	t = T_SIZE;
	for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
    for (int k = 0; k < N3; ++k) {
    for (int l = 0; l < N4; ++l) {
		check_result(t, i, j, k, l, a.interior(t, i, j, k, l), b.interior(t, i, j, k, l));
	} } } }
#endif

	return 0;
}
