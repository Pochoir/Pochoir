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

using namespace std;

#define APP_DEBUG 0
#define TIMES 1
#define TOLERANCE (1e-6)

static double max_diff = 0;
static double max_a = 0, max_b = 0;

void check_result(int t, int i, double a, double b)
{
    double l_diff = abs(a - b);
    if (l_diff < TOLERANCE) {
    //    printf("a(%d, %d) == b(%d, %d) == %f : passed!\n", t, i, t, i, a);
    } else {
        if (l_diff > max_diff) {
            max_diff = l_diff;
            max_a = a; max_b = b;
        }
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

#define N 20000
#define T 1000

int main(int argc, char * argv[])
{
    const int BASE = 1024;
    int t;
    struct timeval start, end;
    double min_tdiff = INF;
    /* the 1D spatial dimension has 'N' points */
    // int N = 0, T = 0;
    double umin, umax;
    char pochoir_plan_file_name[100];

    if (argc < 3) {
        printf("argc < 3, quit! \n");
        exit(1);
    }
    // N = StrToInt(argv[1]);
    // T = StrToInt(argv[2]);
    printf("N = %d, T = %d\n", N, T);
    Pochoir_Shape_1D oned_3pt[] = {{0, 0}, {-1, 0}, {-1, -1}, {-1, 1}};
    Pochoir_Shape_1D shape_interior_0[] = {{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}};
    Pochoir_Shape_1D shape_interior_1[] = {{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}};
    Pochoir_Shape_1D shape_interior_2[] = {{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}};
    Pochoir_Shape_1D shape_exterior_0[] = {{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}};
    Pochoir_Shape_1D shape_exterior_1[] = {{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}};
    Pochoir_Array_1D(double) a(N);
    Pochoir_Array_1D(double) b(N);
    Pochoir_1D leap_frog;
    a.Register_Boundary(periodic_1D);
    b.Register_Shape(oned_3pt);
    b.Register_Boundary(periodic_1D);

    Pochoir_Guard_1D_Begin(guard_interior, t, i)
        /* (T/2) / (N/2) = T / N */
        float _slope = float(T) / (N);
        if ((t <= T/2 && i > N/2 && float(t)/(i - N/2) < _slope) ||
            (t > T/2 && i > N/2 && float(T - t)/(i - N/2) < _slope) ||
            (t > T/2 && i < N/2 && float(T - t)/(N/2 - i) < _slope) ||
            (t <= T/2 && i < N/2 && float(t)/(N/2 - i) < _slope))
            return false;
        else
            return true;
    Pochoir_Guard_1D_End(guard_interior)

    Pochoir_Guard_1D_Begin(guard_exterior, t, i)
        return (!__guard_interior__(t, i));
    Pochoir_Guard_1D_End(guard_exterior)

    Pochoir_Kernel_1D_Begin(interior_0, t, i)
#if APP_DEBUG
        fprintf(stderr, "<interior_0> a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.1 * a(t-1, i-1) + 0.15 * a(t-1, i) + 0.189 * a(t-1, i+1) + 0.8;
    Pochoir_Kernel_1D_End(interior_0, shape_interior_0)

    Pochoir_Kernel_1D_Begin(interior_1, t, i)
#if APP_DEBUG
        fprintf(stderr, "<interior_1> a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.2 * a(t-1, i-1) + 0.25 * a(t-1, i) + 0.289 * a(t-1, i+1) + 0.8;
    Pochoir_Kernel_1D_End(interior_1, shape_interior_1)

    Pochoir_Kernel_1D_Begin(interior_2, t, i)
#if APP_DEBUG
        fprintf(stderr, "<interior_2> a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.3 * a(t-1, i-1) + 0.35 * a(t-1, i) + 0.389 * a(t-1, i+1) + 0.8;
    Pochoir_Kernel_1D_End(interior_2, shape_interior_2)

    Pochoir_Kernel_1D_Begin(exterior_0, t, i)
#if APP_DEBUG
        fprintf(stderr, "<exterior_0> a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.1 * a(t-1, i-1) - 0.15 * a(t-1, i) - 0.189 * a(t-1, i+1) - 0.8;
    Pochoir_Kernel_1D_End(exterior_0, shape_exterior_0)

    Pochoir_Kernel_1D_Begin(exterior_1, t, i)
#if APP_DEBUG
        fprintf(stderr, "<exterior_1> a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.2 * a(t-1, i-1) - 0.25 * a(t-1, i) - 0.289 * a(t-1, i+1) - 0.8;
    Pochoir_Kernel_1D_End(exterior_1, shape_exterior_1)

    Pochoir_Kernel<1> tile_interior[3] = {interior_0, interior_1, interior_2};
    Pochoir_Kernel<1> tile_exterior[2] = {exterior_0, exterior_1};
    leap_frog.Register_Tile_Kernels(guard_interior, tile_interior);
    leap_frog.Register_Tile_Kernels(guard_exterior, tile_exterior);
    leap_frog.Register_Array(a);

    /* initialization */
    for (int i = 0; i < N; ++i) {
        a(0, i) = 1.0 * (rand() % BASE);
        b(0, i) = a(0, i);
    }

    Pochoir_Plan<1> & l_plan = leap_frog.Gen_Plan(T);
    sprintf(pochoir_plan_file_name, "pochoir_%d_%d.dat", N, T);
    leap_frog.Store_Plan(pochoir_plan_file_name, l_plan);
    // Pochoir_Plan<1> & ll_plan = leap_frog.Load_Plan(pochoir_plan_file_name);

    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        leap_frog.Run(l_plan);
        // leap_frog.Run(6*T);
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    printf("Pochoir time = %.6f ms\n", min_tdiff);
    leap_frog.Destroy_Plan(l_plan);
//    std::cout << "Pochoir time : " << min_tdiff << " ms" << std::endl;

    min_tdiff = INF;

    /* cilk_for */
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        for (int t = 1; t < T + 1; ++t) {
            cilk_for (int i = 0; i < N; ++i) {
            if (guard_interior(t, i)) {
                /* interior sub-region */
                if (t % 3 == 0) {
#if APP_DEBUG
                    fprintf(stderr, "<interior_0> b(%d, %d)\n", t, i);
#endif
                    b(t, i) = 0.1 * b(t-1, i-1) + 0.15 * b(t-1, i) + 0.189 * b(t-1, i+1) + 0.8;
                }
                if (t % 3 == 1) {
#if APP_DEBUG
                    fprintf(stderr, "<interior_1> b(%d, %d)\n", t, i);
#endif
                    b(t, i) = 0.2 * b(t-1, i-1) + 0.25 * b(t-1, i) + 0.289 * b(t-1, i+1) + 0.8;
                }
                if (t % 3 == 2) {
#if APP_DEBUG
                    fprintf(stderr, "<interior_2> b(%d, %d)\n", t, i);
#endif
                    b(t, i) = 0.3 * b(t-1, i-1) + 0.35 * b(t-1, i) + 0.389 * b(t-1, i+1) + 0.8;
                }
            } else {
                /* exterior sub-region*/
                if (t % 2 == 0) {
#if APP_DEBUG
                    fprintf(stderr, "<exterior_0> b(%d, %d)\n", t, i);
#endif
                    b(t, i) = 0.1 * b(t-1, i-1) - 0.15 * b(t-1, i) - 0.189 * b(t-1, i+1) - 0.8;
                }
                if (t % 2 == 1) {
#if APP_DEBUG
                    fprintf(stderr, "<exterior_1> b(%d, %d)\n", t, i);
#endif
                    b(t, i) = 0.2 * b(t-1, i-1) - 0.25 * b(t-1, i) - 0.289 * b(t-1, i+1) - 0.8;
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
    t = T;
    for (int i = 0; i < N; ++i) {
        check_result(t, i, a(t, i), b(t, i));
    } 
    printf("max_diff = %f, when a = %f, b = %f\n", max_diff, max_a, max_b);

    return 0;
}

