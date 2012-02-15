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
#define APP_DEBUG 0
#define TIMES 1
#define TOLERANCE (1e-6)

double max_diff = 0;
double max_a = 0, max_b = 0;

void check_result(int t, int i, int j, double a, double b)
{
    double l_diff = pabs(a, b);
    if (l_diff < TOLERANCE) {
    //    printf("a(%d, %d, %d) == b(%d, %d, %d) == %f : passed!\n", t, i, j, t, i, j, a);
    } else {
        // printf("a(%d, %d, %d) = %f, b(%d, %d, %d) = %f : FAILED!\n", t, i, j, a, t, i, j, b);
        if (l_diff > max_diff) {
            max_diff = l_diff;
            max_a = a; max_b = b;
        }
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

#define N 65
#define T 65

int main(int argc, char * argv[])
{
    const int BASE = 1024798;
    int t;
    struct timeval start, end;
    double min_tdiff = INF;
    /* the 1D spatial dimension has 'N' points */
    // int N = 0, T = 0;
    char pochoir_plan_file_name[100];

    if (argc < 3) {
        printf("argc < 3, quit! \n");
        exit(1);
    }
    // N = StrToInt(argv[1]);
    // T = StrToInt(argv[2]);
    printf("N = %d, T = %d\n", N, T);
    Pochoir_Shape_2D twod_5pt[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_0_0[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_0_1[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_0_2[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_0_3[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_0_4[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_0_5[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_0_6[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_0_7[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_1_0[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_1_1[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_1_2[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_1_3[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_2_0[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_2_1[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_3_0[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_3_1[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_3_2[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_3_3[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_3_4[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_3_5[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_3_6[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_3_7[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_4_0[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_4_1[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_4_2[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_4_3[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_5_0[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Shape_2D shape_5_1[] = {{0, 0, 0}, {-1, 0, -1}, {-1, 0, 1}, {-1, 1, 0}, {-1, -1, 0}, {-1, 0, 0}};
    Pochoir_Array_2D(double) a(N, N);
    Pochoir_Array_2D(double) b(N, N);
    Pochoir_2D leap_frog;
    a.Register_Boundary(aperiodic_2D);
    b.Register_Shape(twod_5pt);
    b.Register_Boundary(aperiodic_2D);

    /* begin Pochoir_Guard functions */
    Pochoir_Guard_2D_Begin(g_exclusive_0, t, i, j)
        if (i < N/2 && j < N/2)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_exclusive_0)

    Pochoir_Guard_2D_Begin(g_exclusive_1, t, i, j)
        if (i >= N/2 && j >= N/2)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_exclusive_1)

    Pochoir_Guard_2D_Begin(g_inclusive_0, t, i, j)
        if (i < N/2 && j < N/2 && t < T/2)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_inclusive_0)

    Pochoir_Guard_2D_Begin(g_inclusive_1, t, i, j)
        if (i >= N/4 && i < N/3 && j >= N/4 && j < N/3 && t < T*2/3)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_inclusive_1);

    Pochoir_Guard_2D_Begin(g_tiny_inclusive_0, t, i, j)
        if (i >= 1 && i < 5 && j >= 1 && j < 5 && t > 1 && t <= 2)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_tiny_inclusive_0)

    Pochoir_Guard_2D_Begin(g_tiny_inclusive_1, t, i, j)
        if (i > 5 && i <= 10 && j > 5 && j <= 10 && t > 2 && t < 4)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_tiny_inclusive_1)
    /* end Pochoir_Guard functions */

    /* begin Pochoir_Kernel functions */
    Pochoir_Kernel_2D_Begin(k_0_0, t, i, j)
#if APP_DEBUG
        printf("<k_0_0> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.1 * a(t-1, i-1, j) + 0.15 * a(t-1, i, j) + 0.15 * a(t-1, i+1, j) 
          - 0.1 * a(t-1, i, j-1) - 0.15 * a(t-1, i, j) - 0.15 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_0_0, shape_0_0)

    Pochoir_Kernel_2D_Begin(k_0_1, t, i, j)
#if APP_DEBUG
        printf("<k_0_1> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.2 * a(t-1, i-1, j) + 0.25 * a(t-1, i, j) + 0.25 * a(t-1, i+1, j) 
          - 0.2 * a(t-1, i, j-1) - 0.25 * a(t-1, i, j) - 0.25 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_0_1, shape_0_1)

    Pochoir_Kernel_2D_Begin(k_0_2, t, i, j)
#if APP_DEBUG
        printf("<k_0_2> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.3 * a(t-1, i-1, j) + 0.35 * a(t-1, i, j) + 0.35 * a(t-1, i+1, j) 
          - 0.3 * a(t-1, i, j-1) - 0.35 * a(t-1, i, j) - 0.35 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_0_2, shape_0_2)

    Pochoir_Kernel_2D_Begin(k_0_3, t, i, j)
#if APP_DEBUG
        printf("<k_0_3> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.4 * a(t-1, i-1, j) + 0.45 * a(t-1, i, j) + 0.45 * a(t-1, i+1, j) 
          - 0.4 * a(t-1, i, j-1) - 0.45 * a(t-1, i, j) - 0.45 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_0_3, shape_0_3)

    Pochoir_Kernel_2D_Begin(k_0_4, t, i, j)
#if APP_DEBUG
        printf("<k_0_4> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.5 * a(t-1, i-1, j) + 0.55 * a(t-1, i, j) + 0.55 * a(t-1, i+1, j) 
          - 0.5 * a(t-1, i, j-1) - 0.55 * a(t-1, i, j) - 0.55 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_0_4, shape_0_4)

    Pochoir_Kernel_2D_Begin(k_0_5, t, i, j)
#if APP_DEBUG
        printf("<k_0_5> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.6 * a(t-1, i-1, j) + 0.65 * a(t-1, i, j) + 0.65 * a(t-1, i+1, j) 
          - 0.6 * a(t-1, i, j-1) - 0.65 * a(t-1, i, j) - 0.65 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_0_5, shape_0_5)

    Pochoir_Kernel_2D_Begin(k_0_6, t, i, j)
#if APP_DEBUG
        printf("<k_0_6> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.7 * a(t-1, i-1, j) + 0.75 * a(t-1, i, j) + 0.75 * a(t-1, i+1, j) 
          - 0.7 * a(t-1, i, j-1) - 0.75 * a(t-1, i, j) - 0.75 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_0_6, shape_0_6)

    Pochoir_Kernel_2D_Begin(k_0_7, t, i, j)
#if APP_DEBUG
        printf("<k_0_7> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.8 * a(t-1, i-1, j) + 0.85 * a(t-1, i, j) + 0.85 * a(t-1, i+1, j) 
          - 0.8 * a(t-1, i, j-1) - 0.85 * a(t-1, i, j) - 0.85 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_0_7, shape_0_7)

    /* this is a 3D checkerboard style tiling of the entire rectangular region/domain */
    Pochoir_Kernel<2> tile_3D_checkerboard_0[2][2][2] = {{{k_0_0, k_0_1}, {k_0_2, k_0_3}}, {{k_0_4, k_0_5}, {k_0_6, k_0_7}}};

    /**************************************************************************/
    Pochoir_Kernel_2D_Begin(k_1_0, t, i, j)
#if APP_DEBUG
        printf("<k_1_0> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.11 * a(t-1, i-1, j) + 0.151 * a(t-1, i, j) + 0.151 * a(t-1, i+1, j) 
          - 0.11 * a(t-1, i, j-1) - 0.151 * a(t-1, i, j) - 0.151 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_1_0, shape_1_0)

    Pochoir_Kernel_2D_Begin(k_1_1, t, i, j)
#if APP_DEBUG
        printf("<k_1_1> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.21 * a(t-1, i-1, j) + 0.251 * a(t-1, i, j) + 0.251 * a(t-1, i+1, j) 
          - 0.21 * a(t-1, i, j-1) - 0.251 * a(t-1, i, j) - 0.251 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_1_1, shape_1_1)

    Pochoir_Kernel_2D_Begin(k_1_2, t, i, j)
#if APP_DEBUG
        printf("<k_1_2> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.31 * a(t-1, i-1, j) + 0.351 * a(t-1, i, j) + 0.351 * a(t-1, i+1, j) 
          - 0.31 * a(t-1, i, j-1) - 0.351 * a(t-1, i, j) - 0.351 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_1_2, shape_1_2)

    Pochoir_Kernel_2D_Begin(k_1_3, t, i, j)
#if APP_DEBUG
        printf("<k_1_3> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.41 * a(t-1, i-1, j) + 0.451 * a(t-1, i, j) + 0.451 * a(t-1, i+1, j) 
          - 0.41 * a(t-1, i, j-1) - 0.451 * a(t-1, i, j) - 0.451 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_1_3, shape_1_3)

    /* this is a 2D checkerboard style tiling of the entire rectangular region/domain */
    Pochoir_Kernel<2> tile_2D_checkerboard_1[2][2] = {{k_1_0, k_1_1}, {k_1_2, k_1_3}};

    /**************************************************************************/
    Pochoir_Kernel_2D_Begin(k_2_0, t, i, j)
#if APP_DEBUG
        printf("<k_2_0> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.12 * a(t-1, i-1, j) + 0.152 * a(t-1, i, j) + 0.152 * a(t-1, i+1, j) 
          - 0.12 * a(t-1, i, j-1) - 0.152 * a(t-1, i, j) - 0.152 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_2_0, shape_2_0)

    Pochoir_Kernel_2D_Begin(k_2_1, t, i, j)
#if APP_DEBUG
        printf("<k_2_1> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.22 * a(t-1, i-1, j) + 0.252 * a(t-1, i, j) + 0.252 * a(t-1, i+1, j) 
          - 0.22 * a(t-1, i, j-1) - 0.252 * a(t-1, i, j) - 0.252 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_2_1, shape_2_1)

    /* this is a 2D checkerboard style tiling of the entire rectangular region/domain */
    Pochoir_Kernel<2> tile_1D_checkerboard_2[2] = {k_2_0, k_2_1};

    /**********************************************************************************/
    Pochoir_Kernel_2D_Begin(k_3_0, t, i, j)
#if APP_DEBUG
        printf("<k_3_0> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.13 * a(t-1, i-1, j) + 0.153 * a(t-1, i, j) + 0.153 * a(t-1, i+1, j) 
          - 0.13 * a(t-1, i, j-1) - 0.153 * a(t-1, i, j) - 0.153 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_3_0, shape_3_0)

    Pochoir_Kernel_2D_Begin(k_3_1, t, i, j)
#if APP_DEBUG
        printf("<k_3_1> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.23 * a(t-1, i-1, j) + 0.253 * a(t-1, i, j) + 0.253 * a(t-1, i+1, j) 
          - 0.23 * a(t-1, i, j-1) - 0.253 * a(t-1, i, j) - 0.253 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_3_1, shape_3_1)

    Pochoir_Kernel_2D_Begin(k_3_2, t, i, j)
#if APP_DEBUG
        printf("<k_3_2> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.33 * a(t-1, i-1, j) + 0.353 * a(t-1, i, j) + 0.353 * a(t-1, i+1, j) 
          - 0.33 * a(t-1, i, j-1) - 0.353 * a(t-1, i, j) - 0.353 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_3_2, shape_3_2)

    Pochoir_Kernel_2D_Begin(k_3_3, t, i, j)
#if APP_DEBUG
        printf("<k_3_3> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.43 * a(t-1, i-1, j) + 0.453 * a(t-1, i, j) + 0.453 * a(t-1, i+1, j) 
          - 0.43 * a(t-1, i, j-1) - 0.453 * a(t-1, i, j) - 0.453 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_3_3, shape_3_3)

    Pochoir_Kernel_2D_Begin(k_3_4, t, i, j)
#if APP_DEBUG
        printf("<k_3_4> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.53 * a(t-1, i-1, j) + 0.553 * a(t-1, i, j) + 0.553 * a(t-1, i+1, j) 
          - 0.53 * a(t-1, i, j-1) - 0.553 * a(t-1, i, j) - 0.553 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_3_4, shape_3_4)

    Pochoir_Kernel_2D_Begin(k_3_5, t, i, j)
#if APP_DEBUG
        printf("<k_3_5> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.63 * a(t-1, i-1, j) + 0.653 * a(t-1, i, j) + 0.653 * a(t-1, i+1, j) 
          - 0.63 * a(t-1, i, j-1) - 0.653 * a(t-1, i, j) - 0.653 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_3_5, shape_3_5)

    Pochoir_Kernel_2D_Begin(k_3_6, t, i, j)
#if APP_DEBUG
        printf("<k_3_6> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.73 * a(t-1, i-1, j) + 0.753 * a(t-1, i, j) + 0.753 * a(t-1, i+1, j) 
          - 0.73 * a(t-1, i, j-1) - 0.753 * a(t-1, i, j) - 0.753 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_3_6, shape_3_6)

    Pochoir_Kernel_2D_Begin(k_3_7, t, i, j)
#if APP_DEBUG
        printf("<k_3_7> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.83 * a(t-1, i-1, j) + 0.853 * a(t-1, i, j) + 0.853 * a(t-1, i+1, j) 
          - 0.83 * a(t-1, i, j-1) - 0.853 * a(t-1, i, j) - 0.853 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_3_7, shape_3_7)

    /* this is a 3D checkerboard style tiling of the entire rectangular region/domain */
    Pochoir_Kernel<2> tile_3D_checkerboard_3[2][2][2] = {{{k_3_0, k_3_1}, {k_3_2, k_3_3}}, {{k_3_4, k_3_5}, {k_3_6, k_3_7}}};

    /**************************************************************************/
    Pochoir_Kernel_2D_Begin(k_4_0, t, i, j)
#if APP_DEBUG
        printf("<k_4_0> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.14 * a(t-1, i-1, j) + 0.154 * a(t-1, i, j) + 0.154 * a(t-1, i+1, j) 
          - 0.14 * a(t-1, i, j-1) - 0.154 * a(t-1, i, j) - 0.154 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_4_0, shape_4_0)

    Pochoir_Kernel_2D_Begin(k_4_1, t, i, j)
#if APP_DEBUG
        printf("<k_4_1> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.24 * a(t-1, i-1, j) + 0.254 * a(t-1, i, j) + 0.254 * a(t-1, i+1, j) 
          - 0.24 * a(t-1, i, j-1) - 0.254 * a(t-1, i, j) - 0.254 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_4_1, shape_4_1)

    Pochoir_Kernel_2D_Begin(k_4_2, t, i, j)
#if APP_DEBUG
        printf("<k_4_2> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.34 * a(t-1, i-1, j) + 0.354 * a(t-1, i, j) + 0.354 * a(t-1, i+1, j) 
          - 0.34 * a(t-1, i, j-1) - 0.354 * a(t-1, i, j) - 0.354 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_4_2, shape_4_2)

    Pochoir_Kernel_2D_Begin(k_4_3, t, i, j)
#if APP_DEBUG
        printf("<k_4_3> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.44 * a(t-1, i-1, j) + 0.454 * a(t-1, i, j) + 0.454 * a(t-1, i+1, j) 
          - 0.44 * a(t-1, i, j-1) - 0.454 * a(t-1, i, j) - 0.454 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_4_3, shape_4_3)

    /* this is a 2D checkerboard style tiling of the entire rectangular region/domain */
    Pochoir_Kernel<2> tile_2D_checkerboard_4[2][2] = {{k_4_0, k_4_1}, {k_4_2, k_4_3}};

    /**************************************************************************/
    Pochoir_Kernel_2D_Begin(k_5_0, t, i, j)
#if APP_DEBUG
        printf("<k_5_0> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.15 * a(t-1, i-1, j) + 0.155 * a(t-1, i, j) + 0.155 * a(t-1, i+1, j) 
          - 0.15 * a(t-1, i, j-1) - 0.155 * a(t-1, i, j) - 0.155 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_5_0, shape_5_0)

    Pochoir_Kernel_2D_Begin(k_5_1, t, i, j)
#if APP_DEBUG
        printf("<k_5_1> : rec_a(%d, %d, %d)\n", t, i, j);
#endif
        a(t, i, j) = 
            0.25 * a(t-1, i-1, j) + 0.255 * a(t-1, i, j) + 0.255 * a(t-1, i+1, j) 
          - 0.25 * a(t-1, i, j-1) - 0.255 * a(t-1, i, j) - 0.255 * a(t-1, i, j+1);
    Pochoir_Kernel_2D_End(k_5_1, shape_5_1)

    /* this is a 2D checkerboard style tiling of the entire rectangular region/domain */
    Pochoir_Kernel<2> tile_1D_checkerboard_5[2] = {k_5_0, k_5_1};

    /* end Pochoir_Kernel functions */

    leap_frog.Register_Tile_Kernels(g_exclusive_0, tile_3D_checkerboard_0);
    leap_frog.Register_Tile_Kernels(g_exclusive_1, tile_2D_checkerboard_1);
    leap_frog.Register_Tile_Kernels(g_inclusive_0, tile_1D_checkerboard_2);
    leap_frog.Register_Tile_Kernels(g_inclusive_1, tile_3D_checkerboard_3);
    leap_frog.Register_Tile_Kernels(g_tiny_inclusive_0, tile_2D_checkerboard_4);
    leap_frog.Register_Tile_Kernels(g_tiny_inclusive_1, tile_1D_checkerboard_5);
    leap_frog.Register_Array(a);

    /* initialization */
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            auto tmp = 1.0 * (rand() % BASE);
            a(0, i, j) = tmp;
            b(0, i, j) = tmp;
            a(1, i, j) = 0;
            b(1, i, j) = 0;
        }
    }

    Pochoir_Plan<2> & l_plan = leap_frog.Gen_Plan(T);
    sprintf(pochoir_plan_file_name, "pochoir_%d_%d.dat", N, T);
    leap_frog.Store_Plan(pochoir_plan_file_name, l_plan);
    // Pochoir_Plan<2> & ll_plan = leap_frog.Load_Plan(pochoir_plan_file_name);
    min_tdiff = INF;
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        leap_frog.Run(l_plan);
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    fflush(stdout);
    printf("Pochoir time = %.6f ms\n", min_tdiff);
    leap_frog.Destroy_Plan(l_plan);
    // leap_frog.Destroy_Plan(ll_plan);

    min_tdiff = INF;
    /* cilk_for */
// #define b(t, i, j) b.interior(t, i, j)
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        for (int t = 1; t < T + 1; ++t) {
            cilk_for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (g_exclusive_0(t, i, j)) {
                if ((t) % 2 == 0 && i % 2 == 0 && j % 2 == 0) {
#if APP_DEBUG
                    printf("<k_0_0> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.1 * b(t-1, i-1, j) + 0.15 * b(t-1, i, j) + 0.15 * b(t-1, i+1, j)
                      - 0.1 * b(t-1, i, j-1) - 0.15 * b(t-1, i, j) - 0.15 * b(t-1, i, j+1);
                } else if ((t) % 2 == 0 && i % 2 == 0 && j % 2 == 1) {
#if APP_DEBUG
                    printf("<k_0_1> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.2 * b(t-1, i-1, j) + 0.25 * b(t-1, i, j) + 0.25 * b(t-1, i+1, j)
                      - 0.2 * b(t-1, i, j-1) - 0.25 * b(t-1, i, j) - 0.25 * b(t-1, i, j+1);
                } else if ((t) % 2 == 0 && i % 2 == 1 && j % 2 == 0) {
#if APP_DEBUG
                    printf("<k_0_2> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) = 
                        0.3 * b(t-1, i-1, j) + 
                        0.35 * b(t-1, i, j) + 
                        0.35 * b(t-1, i+1, j) - 
                        0.3 * b(t-1, i, j-1) - 
                        0.35 * b(t-1, i, j) - 
                        0.35 * b(t-1, i, j+1);
                } else if ((t) % 2 == 0 && i % 2 == 1 && j % 2 == 1) {
#if APP_DEBUG
                    printf("<k_0_3> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.4 * b(t-1, i-1, j) + 0.45 * b(t-1, i, j) + 0.45 * b(t-1, i+1, j)
                      - 0.4 * b(t-1, i, j-1) - 0.45 * b(t-1, i, j) - 0.45 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1 && i % 2 == 0 && j % 2 == 0) {
#if APP_DEBUG
                    printf("<k_0_4> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.5 * b(t-1, i-1, j) + 0.55 * b(t-1, i, j) + 0.55 * b(t-1, i+1, j)
                      - 0.5 * b(t-1, i, j-1) - 0.55 * b(t-1, i, j) - 0.55 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1 && i % 2 == 0 && j % 2 == 1) {
#if APP_DEBUG
                    printf("<k_0_5> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.6 * b(t-1, i-1, j) + 0.65 * b(t-1, i, j) + 0.65 * b(t-1, i+1, j)
                      - 0.6 * b(t-1, i, j-1) - 0.65 * b(t-1, i, j) - 0.65 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1 && i % 2 == 1 && j % 2 == 0) {
#if APP_DEBUG
                    printf("<k_0_6> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.7 * b(t-1, i-1, j) + 0.75 * b(t-1, i, j) + 0.75 * b(t-1, i+1, j)
                      - 0.7 * b(t-1, i, j-1) - 0.75 * b(t-1, i, j) - 0.75 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1 && i % 2 == 1 && j % 2 == 1) {
#if APP_DEBUG
                    printf("<k_0_7> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.8 * b(t-1, i-1, j) + 0.85 * b(t-1, i, j) + 0.85 * b(t-1, i+1, j)
                      - 0.8 * b(t-1, i, j-1) - 0.85 * b(t-1, i, j) - 0.85 * b(t-1, i, j+1);
                }
            } else if (g_exclusive_1(t, i, j)) {
                if ((t) % 2 == 0 && i % 2 == 0) {
#if APP_DEBUG
                    printf("<k_1_0> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.11 * b(t-1, i-1, j) + 0.151 * b(t-1, i, j) + 0.151 * b(t-1, i+1, j)
                      - 0.11 * b(t-1, i, j-1) - 0.151 * b(t-1, i, j) - 0.151 * b(t-1, i, j+1);
                } else if ((t) % 2 == 0 && i % 2 == 1) {
#if APP_DEBUG
                    printf("<k_1_1> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.21 * b(t-1, i-1, j) + 0.251 * b(t-1, i, j) + 0.251 * b(t-1, i+1, j)
                      - 0.21 * b(t-1, i, j-1) - 0.251 * b(t-1, i, j) - 0.251 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1 && i % 2 == 0) {
#if APP_DEBUG
                    printf("<k_1_2> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.31 * b(t-1, i-1, j) + 0.351 * b(t-1, i, j) + 0.351 * b(t-1, i+1, j)
                      - 0.31 * b(t-1, i, j-1) - 0.351 * b(t-1, i, j) - 0.351 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1 && i % 2 == 1) {
#if APP_DEBUG
                    printf("<k_1_3> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.41 * b(t-1, i-1, j) + 0.451 * b(t-1, i, j) + 0.451 * b(t-1, i+1, j)
                      - 0.41 * b(t-1, i, j-1) - 0.451 * b(t-1, i, j) - 0.451 * b(t-1, i, j+1);
                }
            }
            if (g_inclusive_0(t, i, j)) {
                if ((t) % 2 == 0) {
#if APP_DEBUG
                    printf("<k_2_0> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.12 * b(t-1, i-1, j) + 0.152 * b(t-1, i, j) + 0.152 * b(t-1, i+1, j)
                      - 0.12 * b(t-1, i, j-1) - 0.152 * b(t-1, i, j) - 0.152 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1) {
#if APP_DEBUG
                    printf("<k_2_1> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.22 * b(t-1, i-1, j) + 0.252 * b(t-1, i, j) + 0.252 * b(t-1, i+1, j)
                      - 0.22 * b(t-1, i, j-1) - 0.252 * b(t-1, i, j) - 0.252 * b(t-1, i, j+1);
                }
            }
            if (g_inclusive_1(t, i, j)) {
                if ((t) % 2 == 0 && i % 2 == 0 && j % 2 == 0) {
#if APP_DEBUG
                    printf("<k_3_0> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.13 * b(t-1, i-1, j) + 0.153 * b(t-1, i, j) + 0.153 * b(t-1, i+1, j)
                      - 0.13 * b(t-1, i, j-1) - 0.153 * b(t-1, i, j) - 0.153 * b(t-1, i, j+1);
                } else if ((t) % 2 == 0 && i % 2 == 0 && j % 2 == 1) {
#if APP_DEBUG
                    printf("<k_3_1> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.23 * b(t-1, i-1, j) + 0.253 * b(t-1, i, j) + 0.253 * b(t-1, i+1, j)
                      - 0.23 * b(t-1, i, j-1) - 0.253 * b(t-1, i, j) - 0.253 * b(t-1, i, j+1);
                } else if ((t) % 2 == 0 && i % 2 == 1 && j % 2 == 0) {
#if APP_DEBUG
                    printf("<k_3_2> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.33 * b(t-1, i-1, j) + 0.353 * b(t-1, i, j) + 0.353 * b(t-1, i+1, j)
                      - 0.33 * b(t-1, i, j-1) - 0.353 * b(t-1, i, j) - 0.353 * b(t-1, i, j+1);
                } else if ((t) % 2 == 0 && i % 2 == 1 && j % 2 == 1) {
#if APP_DEBUG
                    printf("<k_3_3> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.43 * b(t-1, i-1, j) + 0.453 * b(t-1, i, j) + 0.453 * b(t-1, i+1, j)
                      - 0.43 * b(t-1, i, j-1) - 0.453 * b(t-1, i, j) - 0.453 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1 && i % 2 == 0 && j % 2 == 0) {
#if APP_DEBUG
                    printf("<k_3_4> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.53 * b(t-1, i-1, j) + 0.553 * b(t-1, i, j) + 0.553 * b(t-1, i+1, j)
                      - 0.53 * b(t-1, i, j-1) - 0.553 * b(t-1, i, j) - 0.553 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1 && i % 2 == 0 && j % 2 == 1) {
#if APP_DEBUG
                    printf("<k_3_5> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.63 * b(t-1, i-1, j) + 0.653 * b(t-1, i, j) + 0.653 * b(t-1, i+1, j)
                      - 0.63 * b(t-1, i, j-1) - 0.653 * b(t-1, i, j) - 0.653 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1 && i % 2 == 1 && j % 2 == 0) {
#if APP_DEBUG
                    printf("<k_3_6> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.73 * b(t-1, i-1, j) + 0.753 * b(t-1, i, j) + 0.753 * b(t-1, i+1, j)
                      - 0.73 * b(t-1, i, j-1) - 0.753 * b(t-1, i, j) - 0.753 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1 && i % 2 == 1 && j % 2 == 1) {
#if APP_DEBUG
                    printf("<k_3_7> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.83 * b(t-1, i-1, j) + 0.853 * b(t-1, i, j) + 0.853 * b(t-1, i+1, j)
                      - 0.83 * b(t-1, i, j-1) - 0.853 * b(t-1, i, j) - 0.853 * b(t-1, i, j+1);
                }
            }
            if (g_tiny_inclusive_0(t, i, j)) {
                if ((t) % 2 == 0 && i % 2 == 0) {
#if APP_DEBUG
                    printf("<k_4_0> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.14 * b(t-1, i-1, j) + 0.154 * b(t-1, i, j) + 0.154 * b(t-1, i+1, j)
                      - 0.14 * b(t-1, i, j-1) - 0.154 * b(t-1, i, j) - 0.154 * b(t-1, i, j+1);
                } else if ((t) % 2 == 0 && i % 2 == 1) {
#if APP_DEBUG
                    printf("<k_4_1> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.24 * b(t-1, i-1, j) + 0.254 * b(t-1, i, j) + 0.254 * b(t-1, i+1, j)
                      - 0.24 * b(t-1, i, j-1) - 0.254 * b(t-1, i, j) - 0.254 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1 && i % 2 == 0) {
#if APP_DEBUG
                    printf("<k_4_2> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.34 * b(t-1, i-1, j) + 0.354 * b(t-1, i, j) + 0.354 * b(t-1, i+1, j)
                      - 0.34 * b(t-1, i, j-1) - 0.354 * b(t-1, i, j) - 0.354 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1 && i % 2 == 1) {
#if APP_DEBUG
                    printf("<k_4_3> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.44 * b(t-1, i-1, j) + 0.454 * b(t-1, i, j) + 0.454 * b(t-1, i+1, j)
                      - 0.44 * b(t-1, i, j-1) - 0.454 * b(t-1, i, j) - 0.454 * b(t-1, i, j+1);
                }
            }
            if (g_tiny_inclusive_1(t, i, j)) {
                if ((t) % 2 == 0) {
#if APP_DEBUG
                    printf("<k_5_0> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.15 * b(t-1, i-1, j) + 0.155 * b(t-1, i, j) + 0.155 * b(t-1, i+1, j)
                      - 0.15 * b(t-1, i, j-1) - 0.155 * b(t-1, i, j) - 0.155 * b(t-1, i, j+1);
                } else if ((t) % 2 == 1) {
#if APP_DEBUG
                    printf("<k_5_1> : loop_b(%d, %d, %d)\n", t, i, j);
#endif
                    b(t, i, j) =
                        0.25 * b(t-1, i-1, j) + 0.255 * b(t-1, i, j) + 0.255 * b(t-1, i+1, j)
                      - 0.25 * b(t-1, i, j-1) - 0.255 * b(t-1, i, j) - 0.255 * b(t-1, i, j+1);
                }
            }
        } 
            }
        }
// #undef b(t, i, j)
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
    printf("max_diff = %f, when a = %f, b = %f\n", max_diff, max_a, max_b);

    return 0;
}

