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

/* test bench for 3D checkerboard style stencil in Pochoir
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

void check_result(int t, int i, int j, double a, double b)
{
    if (abs(a - b) < TOLERANCE) {
    //    printf("a(%d, %d, %d) == b(%d, %d, %d) == %f : passed!\n", t, i, j, t, i, j, a);
    } else {
        printf("a(%d, %d, %d) = %f, b(%d, %d, %d) = %f : FAILED!\n", t, i, j, a, t, i, j, b);
    }
}

Pochoir_Boundary_2D(periodic_2D, arr, t, i, j)
    const int arr_size_0 = arr.size(0);
    const int arr_size_1 = arr.size(1);

    int new_i = (i >= arr_size_1) ? (i - arr_size_1) : (i < 0 ? i + arr_size_1 : i);
    int new_j = (j >= arr_size_0) ? (j - arr_size_0) : (j < 0 ? j + arr_size_0 : j);

    return arr.get(t, new_i, new_j);
Pochoir_Boundary_End

Pochoir_Boundary_2D(aperiodic_2D, arr, t, i, j)
    return 0;
Pochoir_Boundary_End

int main(int argc, char * argv[])
{
    const int BASE = 1024798;
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
    Pochoir_Shape_2D twod_5pt[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_k0[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_k1[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_k2[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_k3[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_k4[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_k5[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_k6[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_k7[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Array_2D(double) a(N, N);
    Pochoir_Array_2D(double) b(N, N);
    Pochoir_2D leap_frog;
    a.Register_Boundary(aperiodic_2D);
    b.Register_Shape(twod_5pt);
    b.Register_Boundary(aperiodic_2D);

    Pochoir_Kernel_2D_Begin(k0, t, i, j)
        a(t, i, j) = 
            0.1 * a(t-1, i-1, j) + 0.15 * a(t-1, i, j) + 0.15 * a(t-1, i+1, j) 
          - 0.1 * a(t-1, i, j-1) - 0.15 * a(t-1, i, j) - 0.15 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k0, shape_k0)

    Pochoir_Kernel_2D_Begin(k1, t, i, j)
        a(t, i, j) = 
            0.2 * a(t-1, i-1, j) + 0.25 * a(t-1, i, j) + 0.25 * a(t-1, i+1, j) 
          - 0.2 * a(t-1, i, j-1) - 0.25 * a(t-1, i, j) - 0.25 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k1, shape_k1)

    Pochoir_Kernel_2D_Begin(k2, t, i, j)
        a(t, i, j) = 
            0.3 * a(t-1, i-1, j) + 0.35 * a(t-1, i, j) + 0.35 * a(t-1, i+1, j) 
          - 0.3 * a(t-1, i, j-1) - 0.35 * a(t-1, i, j) - 0.35 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k2, shape_k2)

    Pochoir_Kernel_2D_Begin(k3, t, i, j)
        a(t, i, j) = 
            0.4 * a(t-1, i-1, j) + 0.45 * a(t-1, i, j) + 0.45 * a(t-1, i+1, j) 
          - 0.4 * a(t-1, i, j-1) - 0.45 * a(t-1, i, j) - 0.45 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k3, shape_k3)

    Pochoir_Kernel_2D_Begin(k4, t, i, j)
        a(t, i, j) = 
            0.5 * a(t-1, i-1, j) + 0.55 * a(t-1, i, j) + 0.55 * a(t-1, i+1, j) 
          - 0.5 * a(t-1, i, j-1) - 0.55 * a(t-1, i, j) - 0.55 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k4, shape_k4)

    Pochoir_Kernel_2D_Begin(k5, t, i, j)
        a(t, i, j) = 
            0.6 * a(t-1, i-1, j) + 0.65 * a(t-1, i, j) + 0.65 * a(t-1, i+1, j) 
          - 0.6 * a(t-1, i, j-1) - 0.65 * a(t-1, i, j) - 0.65 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k5, shape_k5)

    Pochoir_Kernel_2D_Begin(k6, t, i, j)
        a(t, i, j) = 
            0.7 * a(t-1, i-1, j) + 0.75 * a(t-1, i, j) + 0.75 * a(t-1, i+1, j) 
          - 0.7 * a(t-1, i, j-1) - 0.75 * a(t-1, i, j) - 0.75 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k6, shape_k6)

    Pochoir_Kernel_2D_Begin(k7, t, i, j)
        a(t, i, j) = 
            0.8 * a(t-1, i-1, j) + 0.85 * a(t-1, i, j) + 0.85 * a(t-1, i+1, j) 
          - 0.8 * a(t-1, i, j-1) - 0.85 * a(t-1, i, j) - 0.85 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k7, shape_k7)

    printf("Pochoir_Kernel size = %d, sizeof(int) = %d\n", sizeof(k3), sizeof(int));
    /* this is a 3D checkerboard style tiling of the entire rectangular region/domain */
    Pochoir_Kernel<2> tile_3D_checkerboard[2][2][2] = {{{k0, k1}, {k2, k3}}, {{k4, k5}, {k6, k7}}};
    leap_frog.Register_Tile_Kernels(Default_Guard_2D, tile_3D_checkerboard);
    leap_frog.Register_Array(a);

    /* initialization */
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // a(0, i, j) = 1.0 * (rand() % BASE);
            a(0, i, j) = 1.0 * (rand() / BASE);
            b(0, i, j) = a(0, i, j);
            a(1, i, j) = 0;
            b(1, i, j) = 0;
        }
    }

    Pochoir_Plan<2> & l_plan = leap_frog.Gen_Plan(T);
    sprintf(pochoir_plan_file_name, "pochoir_%d_%d.dat\0", N, T);
    leap_frog.Store_Plan(pochoir_plan_file_name, l_plan);
    Pochoir_Plan<2> & ll_plan = leap_frog.Load_Plan(pochoir_plan_file_name);
    min_tdiff = INF;
//    leap_frog.Load_Plan();
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        leap_frog.Run(ll_plan);
//        leap_frog.Run(T);
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    printf("Pochoir time = %.6f ms\n", min_tdiff);
//    std::cout << "Pochoir time : " << min_tdiff << " ms" << std::endl;

    min_tdiff = INF;
    /* cilk_for */
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        for (int t = 1; t < T + 1; ++t) {
            for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if ((t-1) % 2 == 0 && i % 2 == 0 && j % 2 == 0) {
                /* k0 */
                b(t, i, j) = 
                    0.1 * b(t-1, i-1, j) + 0.15 * b(t-1, i, j) + 0.15 * b(t-1, i+1, j) 
                  - 0.1 * b(t-1, i, j-1) - 0.15 * b(t-1, i, j) - 0.15 * b(t-1, i, j+1);
            } else
            if ((t-1) % 2 == 0 && i % 2 == 0 && j % 2 == 1) {
                /* k1 */
                b(t, i, j) = 
                    0.2 * b(t-1, i-1, j) + 0.25 * b(t-1, i, j) + 0.25 * b(t-1, i+1, j) 
                  - 0.2 * b(t-1, i, j-1) - 0.25 * b(t-1, i, j) - 0.25 * b(t-1, i, j+1);
            } else
            if ((t-1) % 2 == 0 && i % 2 == 1 && j % 2 == 0) {
                /* k2 */
                b(t, i, j) = 
                    0.3 * b(t-1, i-1, j) + 0.35 * b(t-1, i, j) + 0.35 * b(t-1, i+1, j) 
                  - 0.3 * b(t-1, i, j-1) - 0.35 * b(t-1, i, j) - 0.35 * b(t-1, i, j+1);
            } else
            if ((t-1) % 2 == 0 && i % 2 == 1 && j % 2 == 1) {
                /* k3 */
                b(t, i, j) = 
                    0.4 * b(t-1, i-1, j) + 0.45 * b(t-1, i, j) + 0.45 * b(t-1, i+1, j) 
                  - 0.4 * b(t-1, i, j-1) - 0.45 * b(t-1, i, j) - 0.45 * b(t-1, i, j+1);
            } else
            if ((t-1) % 2 == 1 && i % 2 == 0 && j % 2 == 0) {
                /* k4 */
                b(t, i, j) = 
                    0.5 * b(t-1, i-1, j) + 0.55 * b(t-1, i, j) + 0.55 * b(t-1, i+1, j) 
                  - 0.5 * b(t-1, i, j-1) - 0.55 * b(t-1, i, j) - 0.55 * b(t-1, i, j+1);
            } else
            if ((t-1) % 2 == 1 && i % 2 == 0 && j % 2 == 1) {
                /* k5 */
                b(t, i, j) = 
                    0.6 * b(t-1, i-1, j) + 0.65 * b(t-1, i, j) + 0.65 * b(t-1, i+1, j) 
                  - 0.6 * b(t-1, i, j-1) - 0.65 * b(t-1, i, j) - 0.65 * b(t-1, i, j+1);
            } else
            if ((t-1) % 2 == 1 && i % 2 == 1 && j % 2 == 0) {
                /* k6 */
                b(t, i, j) = 
                    0.7 * b(t-1, i-1, j) + 0.75 * b(t-1, i, j) + 0.75 * b(t-1, i+1, j) 
                  - 0.7 * b(t-1, i, j-1) - 0.75 * b(t-1, i, j) - 0.75 * b(t-1, i, j+1);
            } else
            if ((t-1) % 2 == 1 && i % 2 == 1 && j % 2 == 1) {
                /* k7 */
                b(t, i, j) = 
                    0.8 * b(t-1, i-1, j) + 0.85 * b(t-1, i, j) + 0.85 * b(t-1, i+1, j) 
                  - 0.8 * b(t-1, i, j-1) - 0.85 * b(t-1, i, j) - 0.85 * b(t-1, i, j+1);
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
        for (int j = 0; j < N; ++j) {
            check_result(t, i, j, a(t, i, j), b(t, i, j));
        }
    } 

    return 0;
}

