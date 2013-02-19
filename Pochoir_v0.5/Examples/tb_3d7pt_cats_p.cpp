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
 *   This program is mimicing the Berkeley Kaushik Datta's 3D 27 point stencil
 *
 *   Suggestsions:                  yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 *********************************************************************************
 */

/* It's an order 1, 3D 7 point stencil to match up with the 3d7pt stencil in paper
 * "NUMA aware iterative stencil computations on many-core systems"
 * -- periodic, call boundary function to handle boundary values
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <pochoir.hpp>

#define TIMES 1
#define TOLERANCE (1e-6)

using namespace std;

void check_result(int t, int i, int j, int k, double a, double b)
{
    if (abs(a - b) < TOLERANCE) {
//      printf("a(%d, %d, %d, %d) == b(%d, %d, %d, %d) == %f : passed!\n", t, i, j, k, t, i, j, k, a);
    } else {
        printf("a(%d, %d, %d, %d) = %f, b(%d, %d, %d, %d) = %f : FAILED!\n", t, i, j, k, a, t, i, j, k, b);
    }

}

Pochoir_Boundary_3D(aperiodic_3D, arr, t, i, j, k)
    return 0;
Pochoir_Boundary_End

Pochoir_Boundary_3D(periodic_3D, arr, t, i, j, k)
    const int arr_size_2 = arr.size(2);
    const int arr_size_1 = arr.size(1);
    const int arr_size_0 = arr.size(0);

    int new_i = (i >= arr_size_2) ? (i - arr_size_2) : (i < 0 ? i + arr_size_2 : i);
    int new_j = (j >= arr_size_1) ? (j - arr_size_1) : (j < 0 ? j + arr_size_1 : j);
    int new_k = (k >= arr_size_0) ? (k - arr_size_0) : (k < 0 ? k + arr_size_0 : k);

    return arr.get(t, new_i, new_j, new_k);
Pochoir_Boundary_End

int main(int argc, char *argv[])
{
    struct timeval start, end;
    const int BASE = 1024;
    const int ds = 1; /* this is the thickness for ghost cells */
    const double c0 = 0.0876, c1 = 0.0765, c2 = 0.0654, c3 = 0.0543;
    const double c4 = 0.0453, c5 = 0.0123, c6 = 0.0238;
    double min_tdiff = INF;
    int Nx, Ny, Nz, T;
    int t;

    if (argc > 3) {
      Nx = atoi(argv[1]);
      Ny = atoi(argv[2]);
      Nz = atoi(argv[3]);
    }
    /* T is time steps */
    if (argc > 4)
      T = atoi(argv[4]);

    printf("Order-%d 3D-Stencil (%d points) with space %dx%dx%d and time %d\n", 
       ds, 7, Nx, Ny, Nz, T);

    Pochoir_Shape_3D fd_shape_3D[] = {
        {0,0,0,0},
        {-1,0,0,0}, {-1,-1,0,0}, {-1,0,-1,0}, {-1,0,0,-1}, {-1,1,0,0}, {-1,0,1,0}, {-1,0,0,1}};

    Pochoir_Array_3D(double) pa(Nz, Ny, Nx), pb(Nz, Ny, Nx);

    Pochoir_3D fd_3D(fd_shape_3D);

    fd_3D.Register_Array(pa);
    pa.Register_Boundary(periodic_3D);
    pb.Register_Shape(fd_shape_3D);
    pb.Register_Boundary(periodic_3D);

    /* initialization! */
    for (int i = 0; i < Nz; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nx; ++ k) {
                pa(0, i, j, k) = 1.0 * (rand() % BASE);
                pa(1, i, j, k) = 0;
                pb(0, i, j, k) = pa(0, i, j, k);
                pb(1, i, j, k) = 0;
            }
        }
    }

    Pochoir_Kernel_3D(fd_3D_fn, t, i, j, k)
        pa(t, i, j, k) = c1 * pa(t-1, i-1, j, k) + c2 * pa(t-1, i, j-1, k) + 
                         c3 * pa(t-1, i, j, k-1) + c4 * pa(t-1, i+1, j, k) +
                         c5 * pa(t-1, i, j+1, k) + c6 * pa(t-1, i, j, k+1) +
                         c0 * pa(t-1, i, j, k);
    Pochoir_Kernel_End

    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        fd_3D.Run(T, fd_3D_fn);
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    std::cout << "Pochoir ET consumed time : " << min_tdiff << " ms " << std::endl;
    double gstencil_r = ((Nx)*(Ny)*(Nz)*(1e-6)*T)/(min_tdiff);
    std::cout << "GStencil/s : " << gstencil_r  << 
                 "\nGflop/s : " << gstencil_r * 13 << std::endl;

#if 0
    min_tdiff = INF;
    /* cilk_for + zero-padding */
    for (int times = 0; times < TIMES; ++times) {
    gettimeofday(&start, 0);
    for (int t = 1; t < T+1; ++t) {
    cilk_for (int i = 0; i < Nz; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nx; ++k) {
        pb(t, i, j, k) = c1 * pb(t-1, i-1, j, k) + c2 * pb(t-1, i, j-1, k) + 
                         c3 * pb(t-1, i, j, k-1) + c4 * pb(t-1, i+1, j, k) +
                         c5 * pb(t-1, i, j+1, k) + c6 * pb(t-1, i, j, k+1) +
                         c0 * pb(t-1, i, j, k);
    } } } }
    gettimeofday(&end, 0);
    min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    std::cout << "Naive Loop consumed time :" << min_tdiff << "ms" << std::endl;
    double gstencils_l = ((Nx)*(Ny)*(Nz)*(1e-6)*T)/(min_tdiff);
    std::cout << "GStencil/s : " << gstencils_l << 
                 "\nGflop/s : " << gstencils_l * 13 << std::endl;

    t = T+1;
    for (int i = ds; i < Nz-ds; ++i) {
    for (int j = ds; j < Ny-ds; ++j) {
    for (int k = ds; k < Nx-ds; ++k) {
        check_result(t, i, j, k, pa.interior(t, i, j, k), pb.interior(t, i, j, k));
    } } }
#endif

    return 0;
}
