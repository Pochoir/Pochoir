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

/* Test bench - 1D heat equation, Non-periodic version */
#include <cstdio>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>

#include <pochoir.hpp>

using namespace std;
#define TIMES 1
#define N_RANK 2
#define TOLERANCE (1e-6)

bool check_result(int t, int i, double a, double b)
{
	if (abs(a - b) < TOLERANCE) {
        // printf("a(%d, %d) == b(%d, %d) == %f : PASSED!\n", t, i, t, i, a);
        return true;
	} else {
		printf("a(%d, %d) = %f, b(%d, %d) = %f : FAILED!\n", t, i, a, t, i, b);
        return false;
	}
}

Pochoir_Boundary_1D(heat_bv_1D, arr, t, i)
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
        printf("Usage: tb_heat_1D_NP_pochoir N_SIZE T_SIZE\n");
        exit(1);
    }

    N_SIZE = StrToInt(argv[1]);
    T_SIZE = StrToInt(argv[2]);
    printf("N_SIZE = %d, T_SIZE = %d\n", N_SIZE, T_SIZE);
	/* data structure of Pochoir - row major */
    Pochoir_Shape_1D heat_shape_1D[] = {{1, 0}, {0, 1}, {0, -1}, {0, 0}};
	Pochoir_Array_1D(double) a(N_SIZE), b(N_SIZE+2);
    Pochoir_1D heat_1D(heat_shape_1D);

	cout << "a(T+1, J, I) = 0.125 * (a(T, J+1, I) - 2.0 * a(T, J, I) + a(T, J-1, I)) + 0.125 * (a(T, J, I+1) - 2.0 * a(T, J, I) + a(T, J, I-1)) + a(T, J, I)" << endl;
    Pochoir_Kernel_1D(heat_1D_fn, t, i)
	   a(t+1, i) = 0.125 * (a(t, i+1) - 2.0 * a(t, i) + a(t, i-1));
    Pochoir_Kernel_End

    a.Register_Boundary(heat_bv_1D);
    heat_1D.Register_Array(a);
    b.Register_Shape(heat_shape_1D);

	for (int i = 0; i < N_SIZE; ++i) {
        a(0, i) = 1.0 * (rand() % BASE); 
        a(1, i) = 0; 
        b(0, i+1) = a(0, i);
        b(1, i+1) = 0;
	} 


#if 1
    for (int times = 0; times < TIMES; ++times) {
	    gettimeofday(&start, 0);
        heat_1D.Run(T_SIZE, heat_1D_fn);
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
    	for (int t = 0; t < T_SIZE; ++t) {
            cilk_for (int i = 1; i < N_SIZE+1; ++i) {
               b.interior(t+1, i) = 0.125 * (b.interior(t, i+1) - 2.0 * b.interior(t, i) + b.interior(t, i-1)); 
            } 
        }
    	gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
	std::cout << "Naive Loop: consumed time :" << min_tdiff << "ms" << std::endl;

	t = T_SIZE;
    bool passed = true;
	for (int i = 0; i < N_SIZE; ++i) {
		passed &= check_result(t, i, a.interior(t, i), b.interior(t, i+1));
	}

    if (passed) {
        printf("PASSED!\n");
    }
#endif

	return 0;
}
