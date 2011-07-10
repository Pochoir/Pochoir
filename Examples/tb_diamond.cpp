/*
 **********************************************************************************
 *  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 *  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 *                           Charles E. Leiserson <cel@mit.edu>
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

/* A diamond-shape irregular stencil computation in Pochoir
 */
#include <cstdio>
#include <cstddef>
// #include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>

#include <pochoir.hpp>

// using namespace std;
#define TIMES 1
#define TOLERANCE (1e-6)

void check_result(int t, int i, double a, double b)
{
    if (abs(a - b) < TOLERANCE) {
//      printf("a(%d, %d) == b(%d, %d) == %f : passed!\n", t, i, t, i, a);
    } else {
        printf("a(%d, %d) = %f, b(%d, %d) = %f : FAILED!\n", t, i, a, t, i, b);
    }
}

Pochoir_Boundary_1D(periodic_1D, arr, t, i)
    const int arr_size_0 = arr.size(0);

    int new_i = (i >= arr_size_0) ? (i - arr_size_0) : (i < 0 ? i + arr_size_0 : i);

    return arr.get(t, new_i);
Pochoir_Boundary_End

Pochoir_Boundary_1D(aperiodic_1D, arr, t, i)
    return 0;
Pochoir_Boundary_End

int main(int argc, char * argv[])
{
    const int BASE = 1024;
    int t;
    struct timeval start, end;
    double min_tdiff = INF;
    /* the 1D spatial dimension has 'N' points */
    int N = 0, T = 0;
    double umin, umax;

    if (argc < 3) {
        printf("argc < 3, quit! \n");
        exit(1);
    }
    N = StrToInt(argv[1]);
    T = StrToInt(argv[2]);
    printf("N = %d, T = %d\n", N, T);
    Pochoir_Shape_1D oned_3pt[] = {{0, 0}, {-1, 0}, {-1, -1}, {-1, 1}};
    Pochoir_Array_1D(double) a(N);
    Pochoir_Array_1D(double) b(N);
    Pochoir_1D leap_frog(oned_3pt);
    leap_frog.Register_Array(a);
    a.Register_Boundary(periodic_1D);
    b.Register_Shape(oned_3pt);
    b.Register_Boundary(periodic_1D);

    /* initialization */
    for (int i = 0; i < N; ++i) {
        a(0, i) = 1.0 * (rand() % BASE);
        b(0, i) = a(0, i);
    }

    Pochoir_Guard_1D(guard_interior, t, i)
        /* (T/2) / (N/2) = T / N */
        float _slope = float(T) / (N);
        if (float(t)/(i - N/2) < _slope ||
            float(T - t)/(i - N/2) < _slope ||
            float(T - t)/(N/2 - i) < _slope ||
            float(t)/(N/2 - i) < _slope)
            return false;
        else
            return true;
    Pochoir_Guard_End

    Pochoir_Guard_1D(guard_exterior, t, i)
        return (!guard_interior(t, i));
    Pochoir_Guard_End

    Pochoir_Kernel_1D(interior_0, t, i)
        a(t, i) = 0.1 * a(t-1, i-1) + 0.125 * a(t-1, i) + 0.189 * a(t-1, i+1);
    Pochoir_Kernel_End

    Pochoir_Kernel_1D(interior_1, t, i)
        a(t, i) = 0.02 * a(t-1, i-1) + 0.025 * a(t-1, i) + 0.0289 * a(t-1, i+1);
    Pochoir_Kernel_End

#if 0
    Pochoir_Kernel_1D(interior_2, t, i)
        a(t, i) = 0.03 * a(t-1, i-1) + 0.035 * a(t-1, i) + 0.0389 * a(t-1, i+1);
    Pochoir_Kernel_End
#endif

    Pochoir_Kernel_1D(exterior_0, t, i)
        a(t, i) = 0.1 * a(t-1, i-1) - 0.125 * a(t-1, i) - 0.189 * a(t-1, i+1);
    Pochoir_Kernel_End

    Pochoir_Kernel_1D(exterior_1, t, i)
        a(t, i) = 0.02 * a(t-1, i-1) - 0.025 * a(t-1, i) - 0.0289 * a(t-1, i+1);
    Pochoir_Kernel_End

    leap_frog.Register_Kernel(guard_interior, interior_0, interior_1);
    leap_frog.Register_Kernel(guard_exterior, exterior_0, exterior_1);

    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        leap_frog.Run(2*T);
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    printf("Pochoir time = %.6f ms\n", min_tdiff);
//    std::cout << "Pochoir time : " << min_tdiff << " ms" << std::endl;

    min_tdiff = INF;

    /* cilk_for */
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        for (int t = 1; t < 2 * T + 1; ++t) {
            cilk_for (int i = 0; i < N; ++i) {
            if (guard_interior(t, i)) {
                /* interior sub-region */
                if (t % 2 == 1) {
                    b(t, i) = 0.1 * b(t-1, i-1) + 0.125 * b(t-1, i) + 0.189 * b(t-1, i+1);
                }
                if (t % 2 == 0) {
                    b(t, i) = 0.02 * b(t-1, i-1) + 0.025 * b(t-1, i) + 0.0289 * b(t-1, i+1);
                }
            } else {
                /* exterior sub-region*/
                if (t % 2 == 1) {
                    b(t, i) = 0.1 * b(t-1, i-1) - 0.125 * b(t-1, i) - 0.189 * b(t-1, i+1);
                }
                if (t % 2 == 0) {
                    b(t, i) = 0.02 * b(t-1, i-1) - 0.025 * b(t-1, i) - 0.0289 * b(t-1, i+1);
                }
            }
            }
        }
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }

    printf("Parallel Loop time = %.6f ms\n", min_tdiff);
//    std::cout << "Parallel Loop time : " << min_tdiff << " ms" << std::endl;

    /* check results! */
    t = 2*T;
    for (int i = 0; i < N; ++i) {
        check_result(t, i, a(t, i), b(t, i));
    } 

    return 0;
}

