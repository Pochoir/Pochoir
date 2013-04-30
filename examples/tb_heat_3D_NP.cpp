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

/* Test bench - 3D heat equation, Non-periodic version */
#include <cstdio>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>

#include <pochoir.hpp>

using namespace std;
#define TIMES 1
#define N_RANK 3
#define TOLERANCE (1e-6)

void check_result(int t, int i, int j, int k, double a, double b)
{
	if (abs(a - b) < TOLERANCE) {
//		printf("a(%d, %d, %d, %d) == b(%d, %d, %d, %d) == %f : passed!\n", t, i, j, k, t, i, j, k, a);
	} else {
		printf("a(%d, %d, %d, %d) = %f, b(%d, %d, %d, %d) = %f : FAILED!\n", t, i, j, k, a, t, i, j, k, b);
	}

}

Pochoir_Boundary_3D(Pochoir_Aperiodic_3D, arr, t, i, j, k)
    return 0;
Pochoir_Boundary_End

int main(int argc, char * argv[])
{
	const int BASE = 1024;
	int t;
	struct timeval start, end;
    double min_tdiff = INF;
    int N_SIZE = 0, T_SIZE = 0;

    if (argc < 3) {
        printf("argc < 3, quit! \n");
        exit(1);
    }
    N_SIZE = StrToInt(argv[1]);
    T_SIZE = StrToInt(argv[2]);
    printf("N_SIZE = %d, T_SIZE = %d\n", N_SIZE, T_SIZE);
    Pochoir_Shape_3D heat_shape_3D[] = {{0, 0, 0, 0}, {-1, 1, 0, 0}, {-1, -1, 0, 0}, {-1, 0, 0, 0}, {-1, 0, 0, -1}, {-1, 0, 0, 1}, {-1, 0, 1, 0}, {-1, 0, -1, 0}};
    Pochoir_3D heat_3D(heat_shape_3D);
	Pochoir_Array_3D(double) a(N_SIZE, N_SIZE, N_SIZE), b(N_SIZE, N_SIZE, N_SIZE);
    Pochoir_Domain I(1, N_SIZE-1), J(1, N_SIZE-1), K(1, N_SIZE-1);
    heat_3D.Register_Array(a);
    b.Register_Shape(heat_shape_3D);

	for (int i = 0; i < N_SIZE; ++i) {
	for (int j = 0; j < N_SIZE; ++j) {
    for (int k = 0; k < N_SIZE; ++k) {
        if (i == 0 || i == N_SIZE-1
            || j == 0 || j == N_SIZE-1
            || k == 0 || k == N_SIZE-1) {
            a(0, i, j, k) = a(1, i, j, k) = 0;
        } else {
            a(0, i, j, k) = 1.0 * (rand() % BASE); 
            a(1, i, j, k) = 0; 
        }
        b(0, i, j, k) = a(0, i, j, k);
        b(1, i, j, k) = 0;
	} } }

    Pochoir_Kernel_3D(heat_3D_fn, t, i, j, k)
	   a(t, i, j, k) = 
           0.125 * (a(t-1, i+1, j, k) - 2.0 * a(t-1, i, j, k) + a(t-1, i-1, j, k)) 
         + 0.125 * (a(t-1, i, j+1, k) - 2.0 * a(t-1, i, j, k) + a(t-1, i, j-1, k)) 
         + 0.125 * (a(t-1, i, j, k+1) - 2.0 * a(t-1, i, j, k) + a(t-1, i, j, k-1))
         + a(t-1, i, j, k);
    Pochoir_Kernel_End

    /* we have to bind arrayInUse and Shape together 
     * => One arrayInUse, one shape[] => One slope[]
     * because each arrayInUse needs to know the slope to determine
     * the boundary region and when to call the user supplied boundary
     * value function
     */
    // a.Register_Boundary(Pochoir_Aperiodic_3D);
    heat_3D.Register_Domain(I, J, K);

#if 1
    for (int times = 0; times < TIMES; ++times) {
	    gettimeofday(&start, 0);
        heat_3D.Run(T_SIZE, heat_3D_fn);
	    gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
	std::cout << "Pochoir ET: consumed time :" << min_tdiff << "ms" << std::endl;

#endif
#if 1
    min_tdiff = INF;
    /* cilk_for + zero-padding */
    for (int times = 0; times < TIMES; ++times) {
	gettimeofday(&start, 0);
	for (int t = 1; t < T_SIZE+1; ++t) {
    cilk_for (int i = 1; i < N_SIZE-1; ++i) {
	for (int j = 1; j < N_SIZE-1; ++j) {
    for (int k = 1; k < N_SIZE-1; ++k) {
	   b.interior(t, i, j, k) = 
           0.125 * (b.interior(t-1, i+1, j, k) - 2.0 * b.interior(t-1, i, j, k) + b.interior(t-1, i-1, j, k)) 
         + 0.125 * (b.interior(t-1, i, j+1, k) - 2.0 * b.interior(t-1, i, j, k) + b.interior(t-1, i, j-1, k)) 
         + 0.125 * (b.interior(t-1, i, j, k+1) - 2.0 * b.interior(t-1, i, j, k) + b.interior(t-1, i, j, k-1))
         + b.interior(t-1, i, j, k);
    } } } }
	gettimeofday(&end, 0);
    min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
	std::cout << "Naive Loop: consumed time :" << min_tdiff << "ms" << std::endl;

	t = T_SIZE;
	for (int i = 1; i < N_SIZE-1; ++i) {
	for (int j = 1; j < N_SIZE-1; ++j) {
    for (int k = 1; k < N_SIZE-1; ++k) {
		check_result(t, i, j, k, a.interior(t, i, j, k), b.interior(t, i, j, k));
	} } }
#endif

	return 0;
}
