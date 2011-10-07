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
 *   Suggestions:                   yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 *********************************************************************************
 */

/* test bench for 2D checkerboard style stencil in Pochoir
 */
#include <cstdio>
#include <cstddef>
// #include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>

#include <pochoir.hpp>

// using namespace std;
#define APP_DEBUG 0
#define TIMES 1
#define TOLERANCE (1e-6)

void check_result(int t, int i, double a, double b)
{
    if (abs(a - b) < TOLERANCE) {
    //    printf("a(%d, %d) == b(%d, %d) == %f : passed!\n", t, i, t, i, a);
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
    char pochoir_plan_file_name[100];

    if (argc < 3) {
        printf("argc < 3, quit! \n");
        exit(1);
    }
    N = StrToInt(argv[1]);
    T = StrToInt(argv[2]);
    printf("N = %d, T = %d\n", N, T);
    Pochoir_Shape_1D oned_3pt[] = {{0, 0}, {-1, 0}, {-1, -1}, {-1, 1}};
    Pochoir_Shape_1D shape_exclusive_0[] = {{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}};
    Pochoir_Shape_1D shape_exclusive_1[] = {{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}};
    Pochoir_Shape_1D shape_inclusive_0[] = {{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}};
    Pochoir_Shape_1D shape_inclusive_1[] = {{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}};
    Pochoir_Shape_1D shape_tiny_inclusive_0[] = {{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}};
    Pochoir_Shape_1D shape_tiny_inclusive_1[] = {{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}};
    Pochoir_Array_1D(double) a(N);
    Pochoir_Array_1D(double) b(N);
    Pochoir_1D overlap;
    a.Register_Boundary(aperiodic_1D);
    b.Register_Shape(oned_3pt);
    b.Register_Boundary(aperiodic_1D);

    Pochoir_Guard_1D_Begin(g_exclusive_0, t, i)
        if (i < N/2)
            return true;
        else
            return false;
    Pochoir_Guard_1D_End(g_exclusive_0)

    Pochoir_Guard_1D_Begin(g_exclusive_1, t, i)
        if (i >= N/2) 
            return true;
        else 
            return false;
    Pochoir_Guard_1D_End(g_exclusive_1)

    Pochoir_Guard_1D_Begin(g_inclusive_0, t, i)
        if (i < N/2 && t < T/2)
            return true;
        else
            return false;
    Pochoir_Guard_1D_End(g_inclusive_0)

    Pochoir_Guard_1D_Begin(g_inclusive_1, t, i)
        if (i >= N/4 && i < N/3 && t < T*2/3)
            return true;
        else
            return false;
    Pochoir_Guard_1D_End(g_inclusive_1);

    Pochoir_Guard_1D_Begin(g_tiny_inclusive_0, t, i)
        if (i >= 1 && i < 5 && t > 1 && t <= 2)
            return true;
        else
            return false;
    Pochoir_Guard_1D_End(g_tiny_inclusive_0)

    Pochoir_Guard_1D_Begin(g_tiny_inclusive_1, t, i)
        if (i > 5 && i <= 10 && t > 2 && t < 4)
            return true;
        else
            return false;
    Pochoir_Guard_1D_End(g_tiny_inclusive_1)

    Pochoir_Kernel_1D_Begin(k_exclusive_0_0, t, i)
#if APP_DEBUG
        printf("<k_exclusive_0_0> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.11 * a(t-1, i-1) + 0.15 * a(t-1, i) + 0.185 * a(t-1, i+1) + 0.8;
    Pochoir_Kernel_1D_End(k_exclusive_0_0, shape_exclusive_0)

    Pochoir_Kernel_1D_Begin(k_exclusive_0_1, t, i)
#if APP_DEBUG
        printf("<k_exclusive_0_1> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.111 * a(t-1, i-1) + 0.151 * a(t-1, i) + 0.185 * a(t-1, i+1) + 0.81;
    Pochoir_Kernel_1D_End(k_exclusive_0_1, shape_exclusive_0)

    Pochoir_Kernel_1D_Begin(k_exclusive_0_2, t, i)
#if APP_DEBUG
        printf("<k_exclusive_0_2> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.112 * a(t-1, i-1) + 0.152 * a(t-1, i) + 0.185 * a(t-1, i+1) + 0.82;
    Pochoir_Kernel_1D_End(k_exclusive_0_2, shape_exclusive_0)

    Pochoir_Kernel_1D_Begin(k_exclusive_0_3, t, i)
#if APP_DEBUG
        printf("<k_exclusive_0_3> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.113 * a(t-1, i-1) + 0.153 * a(t-1, i) + 0.185 * a(t-1, i+1) + 0.83;
    Pochoir_Kernel_1D_End(k_exclusive_0_3, shape_exclusive_0)

    Pochoir_Kernel<1> tile_exclusive_0[2][2] = {{k_exclusive_0_0, k_exclusive_0_1}, {k_exclusive_0_2, k_exclusive_0_3}};

    Pochoir_Kernel_1D_Begin(k_exclusive_1_0, t, i)
#if APP_DEBUG
        printf("<k_exclusive_1_0> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.21 * a(t-1, i-1) + 0.25 * a(t-1, i) + 0.285 * a(t-1, i+1) + 0.8;
    Pochoir_Kernel_1D_End(k_exclusive_1_0, shape_exclusive_1)

    Pochoir_Kernel_1D_Begin(k_exclusive_1_1, t, i)
#if APP_DEBUG
        printf("<k_exclusive_1_1> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.212 * a(t-1, i-1) + 0.252 * a(t-1, i) + 0.285 * a(t-1, i+1) + 0.82;
    Pochoir_Kernel_1D_End(k_exclusive_1_1, shape_exclusive_1)

    Pochoir_Kernel<1> tile_exclusive_1[2] = {k_exclusive_1_0, k_exclusive_1_1};

    Pochoir_Kernel_1D_Begin(k_inclusive_0_0, t, i)
#if APP_DEBUG
        printf("<k_inclusive_0_0> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.31 * a(t-1, i-1) - 0.35 * a(t-1, i) + 0.385 * a(t-1, i+1) - 0.8;
    Pochoir_Kernel_1D_End(k_inclusive_0_0, shape_inclusive_0)

    Pochoir_Kernel_1D_Begin(k_inclusive_0_1, t, i)
#if APP_DEBUG
        printf("<k_inclusive_0_1> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.311 * a(t-1, i-1) - 0.351 * a(t-1, i) + 0.385 * a(t-1, i+1) - 0.81;
    Pochoir_Kernel_1D_End(k_inclusive_0_1, shape_inclusive_0)

    Pochoir_Kernel_1D_Begin(k_inclusive_0_2, t, i)
#if APP_DEBUG
        printf("<k_inclusive_0_2> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.312 * a(t-1, i-1) - 0.352 * a(t-1, i) + 0.385 * a(t-1, i+1) - 0.82;
    Pochoir_Kernel_1D_End(k_inclusive_0_2, shape_inclusive_0)

    Pochoir_Kernel_1D_Begin(k_inclusive_0_3, t, i)
#if APP_DEBUG
        printf("<k_inclusive_0_3> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.313 * a(t-1, i-1) - 0.353 * a(t-1, i) + 0.385 * a(t-1, i+1) - 0.83;
    Pochoir_Kernel_1D_End(k_inclusive_0_3, shape_inclusive_0)

    Pochoir_Kernel<1> tile_inclusive_0[2][2] = {{k_inclusive_0_0, k_inclusive_0_1}, {k_inclusive_0_2, k_inclusive_0_3}};

    Pochoir_Kernel_1D_Begin(k_inclusive_1_0, t, i)
#if APP_DEBUG
        printf("<k_inclusive_1_0> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.41 * a(t-1, i-1) - 0.45 * a(t-1, i) + 0.485 * a(t-1, i+1) - 0.8;
    Pochoir_Kernel_1D_End(k_inclusive_1_0, shape_inclusive_1)

    Pochoir_Kernel_1D_Begin(k_inclusive_1_1, t, i)
#if APP_DEBUG
        printf("<k_inclusive_1_1> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.414 * a(t-1, i-1) - 0.454 * a(t-1, i) + 0.485 * a(t-1, i+1) - 0.84;
    Pochoir_Kernel_1D_End(k_inclusive_1_1, shape_inclusive_1)

    Pochoir_Kernel<1> tile_inclusive_1[2] = {k_inclusive_1_0, k_inclusive_1_1};

    Pochoir_Kernel_1D_Begin(k_tiny_inclusive_0_0, t, i)
#if APP_DEBUG
        printf("<k_tiny_inclusive_0_0> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.51 * a(t-1, i-1) + 0.55 * a(t-1, i) - 0.585 * a(t-1, i+1) + 0.8;
    Pochoir_Kernel_1D_End(k_tiny_inclusive_0_0, shape_tiny_inclusive_0)

    Pochoir_Kernel_1D_Begin(k_tiny_inclusive_0_1, t, i)
#if APP_DEBUG
        printf("<k_tiny_inclusive_0_1> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.511 * a(t-1, i-1) + 0.551 * a(t-1, i) - 0.585 * a(t-1, i+1) + 0.81;
    Pochoir_Kernel_1D_End(k_tiny_inclusive_0_1, shape_tiny_inclusive_0)

    Pochoir_Kernel_1D_Begin(k_tiny_inclusive_0_2, t, i)
#if APP_DEBUG
        printf("<k_tiny_inclusive_0_2> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.512 * a(t-1, i-1) + 0.552 * a(t-1, i) - 0.585 * a(t-1, i+1) + 0.82;
    Pochoir_Kernel_1D_End(k_tiny_inclusive_0_2, shape_tiny_inclusive_0)

    Pochoir_Kernel_1D_Begin(k_tiny_inclusive_0_3, t, i)
#if APP_DEBUG
        printf("<k_tiny_inclusive_0_3> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.513 * a(t-1, i-1) + 0.553 * a(t-1, i) - 0.585 * a(t-1, i+1) + 0.83;
    Pochoir_Kernel_1D_End(k_tiny_inclusive_0_3, shape_tiny_inclusive_0)

    Pochoir_Kernel<1> tile_tiny_inclusive_0[2][2] = {{k_tiny_inclusive_0_0, k_tiny_inclusive_0_1}, {k_tiny_inclusive_0_2, k_tiny_inclusive_0_3}};

    Pochoir_Kernel_1D_Begin(k_tiny_inclusive_1_0, t, i)
#if APP_DEBUG
        printf("<k_tiny_inclusive_1_0> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.61 * a(t-1, i-1) + 0.65 * a(t-1, i) - 0.685 * a(t-1, i+1) + 0.8;
    Pochoir_Kernel_1D_End(k_tiny_inclusive_1_0, shape_tiny_inclusive_1)

    Pochoir_Kernel_1D_Begin(k_tiny_inclusive_1_1, t, i)
#if APP_DEBUG
        printf("<k_tiny_inclusive_1_1> : a(%d, %d)\n", t, i);
#endif
        a(t, i) = 0.616 * a(t-1, i-1) + 0.656 * a(t-1, i) - 0.685 * a(t-1, i+1) + 0.86;
    Pochoir_Kernel_1D_End(k_tiny_inclusive_1_1, shape_tiny_inclusive_1)

    Pochoir_Kernel<1> tile_tiny_inclusive_1[2] = {k_tiny_inclusive_1_0, k_tiny_inclusive_1_1};

    /* this is a 2D checkerboard style tiling of the entire rectangular region/domain */
    overlap.Register_Exclusive_Tile_Kernels(g_exclusive_0, tile_exclusive_0);
    overlap.Register_Exclusive_Tile_Kernels(g_exclusive_1, tile_exclusive_1);
    overlap.Register_Inclusive_Tile_Kernels(g_inclusive_0, tile_inclusive_0);
    overlap.Register_Inclusive_Tile_Kernels(g_inclusive_1, tile_inclusive_1);
    overlap.Register_Tiny_Inclusive_Tile_Kernels(g_tiny_inclusive_0, tile_tiny_inclusive_0);
    overlap.Register_Tiny_Inclusive_Tile_Kernels(g_tiny_inclusive_1, tile_tiny_inclusive_1);
    overlap.Register_Array(a);

    /* initialization */
    for (int i = 0; i < N; ++i) {
        a(0, i) = 1.0 * (rand() % BASE);
        b(0, i) = a(0, i);
    }

    Pochoir_Plan<1> & l_plan = overlap.Gen_Plan(T);
    sprintf(pochoir_plan_file_name, "pochoir_%d_%d.dat", N, T);
    overlap.Store_Plan(pochoir_plan_file_name, l_plan);
    Pochoir_Plan<1> & ll_plan = overlap.Load_Plan(pochoir_plan_file_name);
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        overlap.Run(ll_plan);
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    printf("Pochoir time = %.6f ms\n", min_tdiff);

    printf("\n--------------------------------------------------------------------------\n");
    min_tdiff = INF;
    /* cilk_for */
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        for (int t = 1; t < T + 1; ++t) {
            for (int i = 0; i < N; ++i) {
                if (g_exclusive_0(t-1, i)) {
                    /* k_exclusive_0 */
                    if ((t-1) % 2 == 0 && i % 2 == 0) {
#if APP_DEBUG
                        printf("<k_exclusive_0_0> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.11 * b(t-1, i-1) + 0.15 * b(t-1, i) + 0.185 * b(t-1, i+1) + 0.8;
                    } else if ((t-1) % 2 == 0 && i % 2 == 1) {
#if APP_DEBUG
                        printf("<k_exclusive_0_1> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.111 * b(t-1, i-1) + 0.151 * b(t-1, i) + 0.185 * b(t-1, i+1) + 0.81;
                    } else if ((t-1) % 2 == 1 && i % 2 == 0) {
#if APP_DEBUG
                        printf("<k_exclusive_0_2> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.112 * b(t-1, i-1) + 0.152 * b(t-1, i) + 0.185 * b(t-1, i+1) + 0.82;
                    } else if ((t-1) % 2 == 1 && i % 2 == 1) {
 #if APP_DEBUG
                        printf("<k_exclusive_0_3> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.113 * b(t-1, i-1) + 0.153 * b(t-1, i) + 0.185 * b(t-1, i+1) + 0.83;
                    }
                } else if (g_exclusive_1(t-1, i)) {
                    /* k_exclusive_1 */
                    if ((t-1) % 2 == 0) {
#if APP_DEBUG
                        printf("<k_exclusive_1_0> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.21 * b(t-1, i-1) + 0.25 * b(t-1, i) + 0.285 * b(t-1, i+1) + 0.8;
                    } else if ((t-1) % 2 == 1) {
#if APP_DEBUG
                        printf("<k_exclusive_1_1> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.212 * b(t-1, i-1) + 0.252 * b(t-1, i) + 0.285 * b(t-1, i+1) + 0.82;
                    }
                }
                if (g_inclusive_0(t-1, i)) {
                    /* k_inclusive_0 */
                    if ((t-1) % 2 == 0 && i % 2 == 0) {
#if APP_DEBUG
                        printf("<k_inclusive_0_0> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.31 * b(t-1, i-1) - 0.35 * b(t-1, i) + 0.385 * b(t-1, i+1) - 0.8;
                    } else if ((t-1) % 2 == 0 && i % 2 == 1) {
#if APP_DEBUG
                        printf("<k_inclusive_0_1> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.311 * b(t-1, i-1) - 0.351 * b(t-1, i) + 0.385 * b(t-1, i+1) - 0.81;
                    } else if ((t-1) % 2 == 1 && i % 2 == 0) {
#if APP_DEBUG
                        printf("<k_inclusive_0_2> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.312 * b(t-1, i-1) - 0.352 * b(t-1, i) + 0.385 * b(t-1, i+1) - 0.82;

                    } else if ((t-1) % 2 == 1 && i % 2 == 1) {
#if APP_DEBUG
                        printf("<k_inclusive_0_3> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.313 * b(t-1, i-1) - 0.353 * b(t-1, i) + 0.385 * b(t-1, i+1) - 0.83;
                    }
                }
                if (g_inclusive_1(t-1, i)) {
                    /* k_inclusive_1 */
                    if ((t-1) % 2 == 0) {
#if APP_DEBUG
                        printf("<k_inclusive_1_0> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.41 * b(t-1, i-1) - 0.45 * b(t-1, i) + 0.485 * b(t-1, i+1) - 0.8;
                    } else if ((t-1) % 2 == 1) {
#if APP_DEBUG
                        printf("<k_inclusive_1_1> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.414 * b(t-1, i-1) - 0.454 * b(t-1, i) + 0.485 * b(t-1, i+1) - 0.84;
                    }
                }
                if (g_tiny_inclusive_0(t-1, i)) {
                    /* k_tiny_inclusive_0 */
                    if ((t-1) % 2 == 0 && i % 2 == 0) {
#if APP_DEBUG
                        printf("<k_tiny_inclusive_0_0> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.51 * b(t-1, i-1) + 0.55 * b(t-1, i) - 0.585 * b(t-1, i+1) + 0.8;
                    } else if ((t-1) % 2 == 0 && i % 2 == 1) {
#if APP_DEBUG
                        printf("<k_tiny_inclusive_0_1> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.511 * b(t-1, i-1) + 0.551 * b(t-1, i) - 0.585 * b(t-1, i+1) + 0.81;
                    } else if ((t-1) % 2 == 1 && i % 2 == 0) {
#if APP_DEBUG
                        printf("<k_tiny_inclusive_0_2> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.512 * b(t-1, i-1) + 0.552 * b(t-1, i) - 0.585 * b(t-1, i+1) + 0.82;
                    } else if ((t-1) % 2 == 1 && i % 2 == 1) {
#if APP_DEBUG
                        printf("<k_tiny_inclusive_0_3> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.513 * b(t-1, i-1) + 0.553 * b(t-1, i) - 0.585 * b(t-1, i+1) + 0.83;
                    }
                }
                if (g_tiny_inclusive_1(t-1, i)) {
                    /* k_tiny_inclusive_1 */
                    if ((t-1) % 2 == 0) {
#if APP_DEBUG
                        printf("<k_tiny_inclusive_1_0> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.61 * b(t-1, i-1) + 0.65 * b(t-1, i) - 0.685 * b(t-1, i+1) + 0.8; 
                    } else if ((t-1) % 2 == 1) {
#if APP_DEBUG
                        printf("<k_tiny_inclusive_1_1> : b(%d, %d)\n", t, i);
#endif
                        b(t, i) = 0.616 * b(t-1, i-1) + 0.656 * b(t-1, i) - 0.685 * b(t-1, i+1) + 0.86; 
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

    return 0;
}

