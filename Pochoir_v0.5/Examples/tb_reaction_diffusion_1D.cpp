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

/* Test bench - 1D reaction diffusion equation */
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
#define TOLERANCE (1e-12)

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

bool check_result(int t, int start, int end, Pochoir_Array_1D(double)* a, Pochoir_Array_1D(double)* b) 
{
    double cumulative_error = 0;
    double a_sum = 0;
    for (int i = start; i < end; ++i) {
        cumulative_error += abs(a->interior(t, i) - b->interior(t, i));
        a_sum += abs(a->interior(t, i));
    }

    #ifdef DEBUG
    printf("cumulative_error: %f\n", cumulative_error);
    printf("a_sum: %f\n", a_sum);
    #endif

    double relative_error = cumulative_error / a_sum;
    if (a_sum < TOLERANCE || relative_error < TOLERANCE) {
        return true;
    } else {
        for (int i = start; i < end; ++i) {
            check_result(t, i, a->interior(t, i), b->interior(t, i));
        }
        return false;
    }
}

bool check_result(int t, int start, int end, Pochoir_Array_1D(double)* au, Pochoir_Array_1D(double)* av, Pochoir_Array_2D(double)* b) 
{
    double cumulative_error = 0;
    double a_sum = 0;
    for (int i = start; i < end; ++i) {
        cumulative_error += abs(au->interior(t, i) - b->interior(t, 0, i));
        cumulative_error += abs(av->interior(t, i) - b->interior(t, 1, i));
        a_sum += abs(au->interior(t, i));
        a_sum += abs(av->interior(t, i));
    }

    #ifdef DEBUG
    printf("cumulative_error: %f\n", cumulative_error);
    printf("a_sum: %f\n", a_sum);
    #endif

    double relative_error = cumulative_error / a_sum;
    if (a_sum < TOLERANCE || relative_error < TOLERANCE) {
        return true;
    } else {
        for (int i = start; i < end; ++i) {
            check_result(t, i, au->interior(t, i), b->interior(t, 0, i));
            check_result(t, i, av->interior(t, i), b->interior(t, 1, i));
        }
        return false;
    }
}

Pochoir_Boundary_1D(reaction_diffusion_bv_1D, arr, t, i)
    return 0;
Pochoir_Boundary_End

Pochoir_Boundary_2D(reaction_diffusion_bv_2D, arr, t, s, i)
    return 0;
Pochoir_Boundary_End

int main(int argc, char * argv[])
{
    const int BASE = 1;
    int t;
    struct timeval start, end;
    double min_tdiff = INF;
    int N_SIZE = 0, T_SIZE = 0;

    // reaction-diffusion constants
    const double alpha = 0.1;
    const double beta = 0.9;
    const double gamma = 0.2;
    const double diffusion_rate = 10;

    if (argc < 3) {
        printf("Usage: %s N_SIZE T_SIZE\n", argv[0]);
        exit(1);
    }

    N_SIZE = StrToInt(argv[1]);
    T_SIZE = StrToInt(argv[2]);
    printf("N_SIZE = %d, T_SIZE = %d\n", N_SIZE, T_SIZE);
    /* data structure of Pochoir - row major */


    Pochoir_Shape_1D diffusion_shape_1D[] = {{1, 0}, {0, 1}, {0, -1}, {0, 0}}; // it's not really possible to express 
                                                                               // the reaction terms (since it involves separate arrays)
    Pochoir_Shape_2D diffusion_shape_2D[] = {{1, 0, 0}, // time, species, i
                                             {0, 0, 1}, {0, 0, -1}, {0, 0, 0}, // diffusion terms
                                             {0, -1, 0}, {0, 1, 0}}; // reaction terms
    Pochoir_1D reaction_diffusion_1D(diffusion_shape_1D);
    Pochoir_2D reaction_diffusion_2D(diffusion_shape_2D);

    Pochoir_Array_1D(double) au(N_SIZE), av(N_SIZE);
    Pochoir_Array_1D(double) bu(N_SIZE), bv(N_SIZE);
    Pochoir_Array_2D(double) c(2, N_SIZE); // species, i
    // Pochoir_Array_2D(double) c(N_SIZE, 2); // species, i

    cout << "a(T+1, J, I) = 0.125 * (a(T, J+1, I) - 2.0 * a(T, J, I) + a(T, J-1, I)) + 0.125 * (a(T, J, I+1) - 2.0 * a(T, J, I) + a(T, J, I-1)) + a(T, J, I)" << endl;
    Pochoir_Kernel_1D(reaction_diffusion_1D_fn, t, i)
       au(t+1, i) = 0.0125 * (au(t, i+1) - 2.0 * au(t, i) + au(t, i-1)) + gamma * (alpha - au(t, i) + (au(t, i) * au(t, i) * av(t, i)));
       av(t+1, i) = diffusion_rate * 0.0125 * (av(t, i+1) - 2.0 * av(t, i) + av(t, i-1)) + gamma * (beta - (au(t, i) * au(t, i) * av(t, i)));
    Pochoir_Kernel_End

    Pochoir_Kernel_2D(reaction_diffusion_2D_fn, t, s, i)
        if (s == 0) {
           c(t+1, s, i) = 0.0125 * (c(t, 0, i+1) - 2.0 * c(t, 0, i) + c(t, 0, i-1)) + gamma * (alpha - c(t, 0, i) + (c(t, 0, i) * c(t, 0, i) * c(t, 1, i)));
        } else if (s == 1) {
           c(t+1, s, i) = diffusion_rate * 0.0125 * (c(t, 1, i+1) - 2.0 * c(t, 1, i) + c(t, 1, i-1)) + gamma * (beta - (c(t, 0, i) * c(t, 0, i) * c(t, 1, i)));
        }
    Pochoir_Kernel_End

    au.Register_Boundary(reaction_diffusion_bv_1D);
    av.Register_Boundary(reaction_diffusion_bv_1D);
    reaction_diffusion_1D.Register_Array(au);
    reaction_diffusion_1D.Register_Array(av);

    bu.Register_Boundary(reaction_diffusion_bv_1D);
    bv.Register_Boundary(reaction_diffusion_bv_1D);
    bu.Register_Shape(diffusion_shape_1D);
    bv.Register_Shape(diffusion_shape_1D);

    c.Register_Boundary(reaction_diffusion_bv_2D);
    reaction_diffusion_2D.Register_Array(c);

    for (int i = 0; i < N_SIZE; ++i) {
        au(0, i) = 1.0 * (rand() % BASE); 
        au(1, i) = 0; 
        av(0, i) = 1.0 * (rand() % BASE);
        av(1, i) = 0;

        bu(0, i) = au(0, i);
        bu(1, i) = 0; 
        bv(0, i) = av(0, i);
        bv(1, i) = 0;

        c(0, 0, i) = au(0, i);
        c(1, 0, i) = 0;
        c(0, 1, i) = av(0, i);
        c(1, 1, i) = 0; 
    } 


#if 1
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        reaction_diffusion_1D.Run(T_SIZE, reaction_diffusion_1D_fn);
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    std::cout << "Pochoir Parallel-Arrays ET: consumed time :" << min_tdiff << "ms" << std::endl;
#endif

#if 1
    min_tdiff = INF;
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        reaction_diffusion_2D.Run(T_SIZE, reaction_diffusion_2D_fn);
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    std::cout << "Pochoir 2D-Array ET: consumed time :" << min_tdiff << "ms" << std::endl;
#endif

#if 1
    min_tdiff = INF;
    /* cilk_for + zero-padding */
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        for (int t = 0; t < T_SIZE; ++t) {
            cilk_for (int i = 0; i < N_SIZE; ++i) {
               bu(t+1, i) = 0.0125 * (bu(t, i+1) - 2.0 * bu(t, i) + bu(t, i-1)) + gamma * (alpha - bu(t, i) + (bu(t, i) * bu(t, i) * bv(t, i)));
               bv(t+1, i) = diffusion_rate * 0.0125 * (bv(t, i+1) - 2.0 * bv(t, i) + bv(t, i-1)) + gamma * (beta - (bu(t, i) * bu(t, i) * bv(t, i)));
            } 
        }
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    std::cout << "Naive Loop: consumed time :" << min_tdiff << "ms" << std::endl;
#endif

#if 1
    bool correct = true;
    correct &= check_result(T_SIZE, 0, N_SIZE, &au, &bu);
    correct &= check_result(T_SIZE, 0, N_SIZE, &av, &bv);
    correct &= check_result(T_SIZE, 0, N_SIZE, &au, &av, &c);
    printf("Correctness check: %s!\n", correct ? "PASSED": "FAILED");
#endif

    return 0;
}

