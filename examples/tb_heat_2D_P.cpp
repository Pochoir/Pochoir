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
 *   Suggestions:                   yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 *********************************************************************************
 */

/* Test bench - 2D heat equation, Periodic version */
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

static double max_diff = 0;
static double diff_a = 0, diff_b = 0;
static double max_a = 0, max_b = 0;

void check_result(int t, int j, int i, double a, double b)
{
    double l_diff = abs(a - b);
	if (l_diff < TOLERANCE) {
//		printf("a(%d, %d, %d) == b(%d, %d, %d) == %f : passed!\n", t, j, i, t, j, i, a);
	} else {
        if (l_diff > max_diff) {
            max_diff = l_diff;
            diff_a = a; diff_b = b;
        }
	}
    if (a > max_a) max_a = a;
    if (b > max_b) max_b = b;
}

Pochoir_Boundary_2D(periodic_2D, arr, t, i, j)
    const int arr_size_1 = arr.size(1);
    const int arr_size_0 = arr.size(0);

    int new_i = (i >= arr_size_1) ? (i - arr_size_1) : (i < 0 ? i + arr_size_1 : i);
    int new_j = (j >= arr_size_0) ? (j - arr_size_0) : (j < 0 ? j + arr_size_0 : j);

    /* we use arr.get(...) instead of arr(...) to implement different boundary
     * checking strategy: In arr(...), if off-boundary access occurs, we call
     * boundary function to supply a value; in arr.get(...), if off-boundary 
     * access occurs, we will print the off-boundary access and quit!
     */
    return arr.get(t, new_i, new_j);
    // return arr.get(t, -1, -1);
Pochoir_Boundary_End

#define N 2000
#define T 200 

int main(int argc, char * argv[])
{
	const int BASE = 1024;
	int t;
	struct timeval start, end;
    // int N = 0, T = 0;
    char pochoir_plan_file_name[100];

    if (argc < 3) {
        printf("argc < 3, quit! \n");
        exit(1);
    }
    // N = StrToInt(argv[1]);
    // T = StrToInt(argv[2]);
    printf("N = %d, T = %d\n", N, T);
    Pochoir_Shape_2D heat_shape_2D[] = {{0, 0, 0}, {-1, 1, 0}, {-1, 0, 0}, {-1, -1, 0}, {-1, 0, -1}, {-1, 0, 1}};
    Pochoir_2D heat_2D;
	Pochoir_Array_2D(double) a(N, N), b(N, N);
    a.Register_Boundary(periodic_2D);
    b.Register_Shape(heat_shape_2D);
    b.Register_Boundary(periodic_2D);

	cout << "a(T+1, J, I) = 0.125 * (a(T, J+1, I) - 2.0 * a(T, J, I) + a(T, J-1, I)) + 0.125 * (a(T, J, I+1) - 2.0 * a(T, J, I) + a(T, J, I-1)) + a(T, J, I)" << endl;
    Pochoir_Kernel_2D_Begin(heat_2D_fn, t, i, j)
	    a(t, i, j) = 0.125 * (a(t-1, i+1, j) - 2.0 * a(t-1, i, j) + a(t-1, i-1, j)) + 0.125 * (a(t-1, i, j+1) - 2.0 * a(t-1, i, j) + a(t-1, i, j-1)) + a(t-1, i, j);
    Pochoir_Kernel_2D_End(heat_2D_fn, heat_shape_2D)

    heat_2D.Register_Tile_Kernels(Default_Guard_2D, heat_2D_fn);
    heat_2D.Register_Array(a);
    /* Now we can only access the Pochoir_Array after Register_Array,
     * or Register_Shape with the array, because we rely on the shape
     * to provide the depth of toggle array!!! 
     */
	for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j) {
        double tmp = 1.0 * (rand() / BASE); 
        a(0, i, j) = tmp; 
        b(0, i, j) = tmp;
        a(1, i, j) = tmp; 
        b(1, i, j) = tmp;
	} }

	gettimeofday(&start, 0);
    for (int times = 0; times < TIMES; ++times) {
        heat_2D.Run(T);
    }
	gettimeofday(&end, 0);
	std::cout << "Pochoir consumed time :" << 1.0e3 * tdiff(&end, &start)/TIMES << " ms" << std::endl;

	gettimeofday(&start, 0);
    for (int times = 0; times < TIMES; ++times) {
	for (int t = 0; t < T; ++t) {
    cilk_for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j) {
        b(t+1, i, j) = 0.125 * (b(t, i+1, j) - 2.0 * b(t, i, j) + b(t, i-1, j)) + 0.125 * (b(t, i, j+1) - 2.0 * b(t, i, j) + b(t, i, j-1)) + b(t, i, j); } } }
    }
	gettimeofday(&end, 0);
	std::cout << "Parallel Loop: consumed time :" << 1.0e3 * tdiff(&end, &start)/TIMES << " ms" << std::endl;

	t = T;
	for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j) {
		check_result(t, i, j, a.interior(t, i, j), b.interior(t, i, j));
	} } 
    printf("max_diff = %f, when a = %f, b = %f\n", max_diff, diff_a, diff_b);
    printf("max_a = %f, max_b = %f\n", max_a, max_b);

	return 0;
}
