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


/*

Makefile:

fib : fib.cpp
#   Phase-II compilation
    ${CC} -o $@ ${OPT_FLAGS} -unroll-pointer $<
    ${ICC} -o fib_gdb ${ICC_DEBUG_FLAGS} fib_pochoir.cpp
#	Phase-I compilation with debugging aid
#	${CC} -o $@ ${POCHOIR_DEBUG_FLAGS} $<

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
//#define TOLERANCE (1e-6)

static double max_diff = 0;
static double diff_a = 0, diff_b = 0;
static double max_a = 0, max_b = 0;

//
// Correction comparation
//
void check_result(int t, int i, long a, long b)
{
    if (a == b) {
        printf("(a(%d, %d) = b(%d, %d)) = %ld : passed!\n", t, i, t, i, a);
    } else {
        printf("(a(%d, %d) =%ld) != (b(%d, %d) = %ld) : failed!\n", t, i, a, t, i, b);
    }
}

Pochoir_Boundary_1D(Pochoir_Aperiodic_1D, arr, t, i)
    return 0;
Pochoir_Boundary_End

//#define N 8000
//#define T 1000

int main(int argc, char * argv[])
{
    //const int BASE = 1024;
    int t;
    struct timeval start, end;
    double min_tdiff = INF;
    /* the 1D spatial dimension has 'N' points */
    /* We can NOT capture normal variables for now! 
     * - So to define the N, T as macros!
     */
    // int N = 0, T = 0;
    //double umin, umax;

    if (argc < 3) {
        printf("argc < 3, quit! \n");
        exit(1);
    }
    int N = StrToInt(argv[1]);
    int T = StrToInt(argv[2]);
    printf("N = %d, T = %d\n", N, T);

    Pochoir_Shape_1D oned_3pt[] = {{0, 0}, {-1, 0}, {-1, -1}, {-1, 1}, {-2, 0}}; // TODO: ÇÐ¸îÓÐÎÊÌâ¡£¡£
    Pochoir_Array_1D(long) a(1);
    Pochoir_Array_1D(long) b(1);
    Pochoir_1D fib;

    a.Register_Boundary(Pochoir_Aperiodic_1D);
    b.Register_Shape(oned_3pt);
    b.Register_Boundary(Pochoir_Aperiodic_1D);

    double coeff0 = 1.0, coeff1 = 2.0;
    Pochoir_Kernel_1D_Begin(fib_kernel, t, i)
#if APP_DEBUG
        printf("<k_exclusive_0_0> : a(%d, %d)\n", t, i);
#endif
        // a(t, i) = coeff0 * a(t-1, i) + coeff1 * a(t-2, i);
        a(t, i) = 1.0 * a(t-1, i) + 2.0 * a(t-2, i);
    Pochoir_Kernel_1D_End(fib_kernel, oned_3pt)


    /* this is a 2D checkerboard style tiling of the entire rectangular 
     * region/domain 
     */
    fib.Register_Tile_Kernels(Default_Guard_1D, fib_kernel);
    fib.Register_Array(a);

    /* initialization */
    a(0, 0) = 0;		// Use Pochoir
    a(1, 0) = 1;
    b(0, 0) = 0;		// Use normal loop
    b(1, 0) = 1;

    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        fib.Run(T);
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    printf("Pochoir time = %.6f ms\n", min_tdiff);
    printf("\n--------------------------------------------------------------------------\n");

    min_tdiff = INF;
    /* cilk_for */
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        for (int t = 2; t < T + 2; ++t) {
            b(t, 0) = coeff0 * b(t - 1, 0) + coeff1 * b(t - 2, 0);
        }
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    printf("Parallel Loop time = %.6f ms\n", min_tdiff);
//    std::cout << "Parallel Loop time : " << min_tdiff << " ms" << std::endl;

    /* check results! */
    t = T;
    check_result(t, 0, a(t, 0), b(t, 0));

    printf("max_diff = %f, when a = %f, b = %f\n", max_diff, diff_a, diff_b);
    printf("max_a = %f, max_b = %f\n", max_a, max_b);

    return 0;
}

