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

/*
 *    Description:  Animate a solution to the 1d wave equation with speed c,
 *                  discretized using a staggered-grid/leap-frog scheme
 *                  on M points from [0, L) with periodic boundary conditions,
 *                  running for a total time c*t = ct.
 *
 *                  Uses timestep c*dt = cdtdx * dx. For stability, cdtdx 
 *                  should be < 1.
 *
 *                  Initial matlab program got from Steven G. Johnson, 
 *                  ported to Pochoir by Yuan Tang.
 */

#include <cstdio>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>

#include <pochoir.hpp>

using namespace std;
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

typedef struct {
    double v;
    double u;
} stagger;

Pochoir_Boundary_1D(periodic_1D, arr, t, i)
    const int arr_size_0 = arr.size(0);

    int new_i = (i >= arr_size_0) ? (i - arr_size_0) : (i < 0 ? i + arr_size_0 : i);

    return arr.get(t, new_i);
Pochoir_Boundary_End

Pochoir_Boundary_1D(a1d, arr, t, i)
    return 0;
Pochoir_Boundary_End

Pochoir_Boundary_1D(b1d, arr, t, i)
    stagger tmp;
    tmp.v = 0; tmp.u = 0;
    return tmp;
Pochoir_Boundary_End

static inline double array_min_u(Pochoir_Array_1D(stagger) & arr, int t, int M) {
    double l_min = INF;
    for (int i = 0; i < M; ++i) {
        l_min = min(l_min, arr(t, i).u);
    }
    return l_min;
}

static inline double array_max_u(Pochoir_Array_1D(stagger) & arr, int t, int M) {
    double l_max = -INF;
    for (int i = 0; i < M; ++i) {
        l_max = max(l_max, arr(t, i).u);
    }
    return l_max;
}

int main(int argc, char * argv[])
{
    const int BASE = 1024;
    int t;
    struct timeval start, end;
    double min_tdiff = INF;
    /* c * dt = cdtdx * dx && cdtdx < 1 holds */
    const double dt = 0.2, dx = 0.5, cdtdx = 0.6, c = 1.5;
    /* ct = c * t = 1500 */
    const double ct = 1500;
    /* the 1D spatial dimension has 'M' points */
    int M = 0, nt = 0;
    double umin, umax;

    if (argc < 3) {
        printf("argc < 3, quit! \n");
        exit(1);
    }
    M = StrToInt(argv[1]);
    // nt = round(ct / (cdtdx * dx));
    nt = StrToInt(argv[2]);
    printf("M = %d, nt = %d\n", M, nt);
    Pochoir_Shape_1D oned_b[] = {{0, 0}, {-1, 0}, {-1, 1}, {-1, -1}};
    Pochoir_Shape_1D oned_a[] = {{0, 0}, {-1, -1}, {-1, 0}, {-1, 1}, {-2, 0}};
    Pochoir_Array_1D(double) a(M);
    Pochoir_Array_1D(stagger) b(M);
    Pochoir_1D leap_frog(oned_a);
    leap_frog.Register_Array(a);
    a.Register_Boundary(a1d);
    b.Register_Shape(oned_b);
    b.Register_Boundary(b1d);
    double * x = new double[M];

    /* we can define an arbitrary function here, and transfer it to the animwave
     * function as a templated input
     */
    auto f = [](double x) { return sin(x); };

    /* initialization */
    for (int i = 0; i < M; ++i) {
        x[i] = i * dx;
        /* a(0, i) corresponds to v */
        a(0, i) = -f(x[i] + 0.5 * dx * (1 + cdtdx));
        /* a(1, i) corresponds to u */
        a(1, i) = f(x[i]);
        b(0, i).v = a(0, i);
        b(0, i).u = a(1, i);
    }

#if 0
    /* it seems that umin and umax don't engage in the stencil computation */
    umin = array_min_u(a, 0, M) - (array_max_u(a, 0, M) - array_min_u(a, 0, M)) * 0.2;
    umax = array_max_u(a, 0, M) + (array_max_u(a, 0, M) - array_min_u(a, 0, M)) * 0.2;
#endif

    Pochoir_Kernel_1D(animwave_fn, t, i)
        if (t & 0x1) {
            /* update u */
            a(t, i) = a(t-2, i) + cdtdx * (a(t-1, i) - a(t-1, i-1));
        } else {
            /* update v */
            a(t, i) = a(t-2, i) + cdtdx * (a(t-1, i+1) - a(t-1, i));
        }
    Pochoir_Kernel_End

    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        leap_frog.Run(2*nt, animwave_fn);
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    std::cout << "Pochoir time : " << min_tdiff << " ms" << std::endl;

    min_tdiff = INF;

    /* cilk_for */
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        for (int t = 1; t <= nt; t += 1) {
            for (int i = 0; i < M; ++i) {
                b(t, i).v = b(t-1, i).v + cdtdx * (b(t-1, i+1).u - b(t-1, i).u);
                b(t, i).u = b(t-1, i).u + cdtdx * (b(t, i).v - b(t, i-1).v);
            }
        }
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    std::cout << "Parallel Loop time : " << min_tdiff << " ms" << std::endl;

    /* check results! */
    t = nt;
    for (int i = 0; i < M; ++i) {
        check_result(t, i, a(2*t, i), b(t, i).v);
        check_result(t, i, a(2*t+1, i), b(t, i).u);
    } 

    return 0;
}

