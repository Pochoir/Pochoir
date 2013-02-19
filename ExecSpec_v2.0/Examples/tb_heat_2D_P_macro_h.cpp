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

#define N1 309
#define N2 307
#define T 2129


int main(int argc, char * argv[])
{
	const int BASE = 1024;
	int t;
	struct timeval start, end;
    char pochoir_plan_file_name[100];

	{
	Pochoir_Shape_2D heat_shape_2D[] = {{0, 0, 0}, {-1, 1, 0}, {-1, 0, 0}, {-1, -1, 0}, {-1, 0, -1}, {-1, 0, 1}};
	Pochoir_Array_2D(double) a(N1, N2), b(N1, N2);
	a.Register_Boundary(periodic_2D);
	b.Register_Shape(heat_shape_2D);
	b.Register_Boundary(periodic_2D);

	Pochoir_Kernel_2D_Begin(k1, t, i, j)
			a(t, i, j) = 0.125 * (a(t-1, i+1, j) - 2.0 * a(t-1, i, j) 
						+ a(t-1, i-1, j)) + 0.125 * (a(t-1, i, j+1) 
						- 2.0 * a(t-1, i, j) + a(t-1, i, j-1)) + 
						a(t-1, i, j);
	Pochoir_Kernel_2D_End(k1, heat_shape_2D)

	Pochoir_Kernel_2D_Begin(k2, t, i, j)
			a(t, i, j) = 0.25 * (a(t-1, i+1, j) - 2.0 * a(t-1, i, j) 
						+ a(t-1, i-1, j)) + 0.25 * (a(t-1, i, j+1) 
						- 2.0 * a(t-1, i, j) + a(t-1, i, j-1)) + 
						a(t-1, i, j);
	Pochoir_Kernel_2D_End(k2, heat_shape_2D)

	Pochoir_Kernel_2D_Begin(k3, t, i, j)
			a(t, i, j) = 0.375 * (a(t-1, i+1, j) - 2.0 * a(t-1, i, j) 
						+ a(t-1, i-1, j)) + 0.375 * (a(t-1, i, j+1) 
						- 2.0 * a(t-1, i, j) + a(t-1, i, j-1)) + 
						a(t-1, i, j);
	Pochoir_Kernel_2D_End(k3, heat_shape_2D)

	Pochoir_Kernel_2D_Begin(k4, t, i, j)
			a(t, i, j) = 0.5 * (a(t-1, i+1, j) - 2.0 * a(t-1, i, j) 
						+ a(t-1, i-1, j)) + 0.5 * (a(t-1, i, j+1) 
						- 2.0 * a(t-1, i, j) + a(t-1, i, j-1)) + 
						a(t-1, i, j);
	Pochoir_Kernel_2D_End(k4, heat_shape_2D)

	/* begin Pochoir_Guard functions */
    Pochoir_Guard_2D_Begin(g1, t, i, j)
        if (i < N1/2 && j < N2/2)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g1)

    Pochoir_Guard_2D_Begin(g2, t, i, j)
        if (i >= N1/2 && j >= N2/2)
            return true;
        else
            return false;
	Pochoir_Guard_2D_End(g2)

    Pochoir_Guard_2D_Begin(g3, t, i, j)
        if (i < N1/2 && j >= N2/2)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g3)

    Pochoir_Guard_2D_Begin(g4, t, i, j)
        if (i >= N1/2 && j < N2/2)
            return true;
        else
            return false;
	Pochoir_Guard_2D_End(g4)

	printf("N1 = %d, N2 = %d, T = %d\n", N1, N2, T);
	Pochoir_2D heat_2D;

	heat_2D.Register_Tile_Kernels(g1, k1);
	heat_2D.Register_Tile_Kernels(g2, k2);
	heat_2D.Register_Tile_Kernels(g3, k3);
	heat_2D.Register_Tile_Kernels(g4, k4);
	heat_2D.Register_Array(a);
	/* Now we can only access the Pochoir_Array after Register_Array,
	 * or Register_Shape with the array, because we rely on the shape
	 * to provide the depth of toggle array!!! 
	 */
	for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
		double tmp = 1.0 * (rand() / BASE); 
		a(0, i, j) = tmp; 
		b(0, i, j) = tmp;
		a(1, i, j) = tmp; 
		b(1, i, j) = tmp;
	} }

	Pochoir_Plan<2> & l_plan = heat_2D.Gen_Plan(T);
	//sprintf(pochoir_plan_file_name, "pochoir_%d_%d.dat\0", N, T);
	//heat_2D.Store_Plan(pochoir_plan_file_name, l_plan);
	// Pochoir_Plan<2> & ll_plan = heat_2D.Load_Plan(pochoir_plan_file_name);
	gettimeofday(&start, 0);
	for (int times = 0; times < TIMES; ++times) {
		heat_2D.Run(l_plan);
	}
	gettimeofday(&end, 0);
	std::cout << "Pochoir consumed time :" << 1.0e3 * tdiff(&end, &start)/TIMES << " ms" << std::endl;

	/* Do not check the correctness of the result */
	/*
	gettimeofday(&start, 0);
	for (int times = 0; times < TIMES; ++times) {
	for (int t = 0; t < T; ++t) {
	cilk_for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
		b(t+1, i, j) = 0.125 * (b(t, i+1, j) - 2.0 * b(t, i, j) + b(t, i-1, j)) + 0.125 * (b(t, i, j+1) - 2.0 * b(t, i, j) + b(t, i, j-1)) + b(t, i, j); } } }
	}
	gettimeofday(&end, 0);
	std::cout << "Parallel Loop: consumed time :" << 1.0e3 * tdiff(&end, &start)/TIMES << " ms" << std::endl;

	t = T;
	for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
		check_result(t, i, j, a.interior(t, i, j), b.interior(t, i, j));
	} } 
	printf("max_diff = %f, when a = %f, b = %f\n", max_diff, diff_a, diff_b);
	printf("max_a = %f, max_b = %f\n", max_a, max_b);
	*/
	}	


	//The equivalent of the above functionality using a single kernel
	{
	printf("N1 = %d, N2 = %d, T = %d\n", N1, N2, T);
	Pochoir_Shape_2D heat_shape_2D[] = {{0, 0, 0}, {-1, 1, 0}, {-1, 0, 0}, {-1, -1, 0}, {-1, 0, -1}, {-1, 0, 1}};
	Pochoir_2D heat_2D_single_kernel;
	Pochoir_Array_2D(double) a(N1, N2), b(N1, N2);
	a.Register_Boundary(periodic_2D);
	b.Register_Shape(heat_shape_2D);
	b.Register_Boundary(periodic_2D);
	Pochoir_Kernel_2D_Begin(heat_2D_fn, t, i, j)
		if (i < N1/2 && j < N2/2)
		{
			a(t, i, j) = 0.125 * (a(t-1, i+1, j) - 2.0 * a(t-1, i, j) 
						+ a(t-1, i-1, j)) + 0.125 * (a(t-1, i, j+1) 
						- 2.0 * a(t-1, i, j) + a(t-1, i, j-1)) + 
						a(t-1, i, j);
		}
		if (i >= N1/2 && j >= N2/2)
		{
			a(t, i, j) = 0.25 * (a(t-1, i+1, j) - 2.0 * a(t-1, i, j) 
						+ a(t-1, i-1, j)) + 0.25 * (a(t-1, i, j+1) 
						- 2.0 * a(t-1, i, j) + a(t-1, i, j-1)) + 
						a(t-1, i, j);
		}
		if (i < N1/2 && j >= N2/2)
		{
			a(t, i, j) = 0.375 * (a(t-1, i+1, j) - 2.0 * a(t-1, i, j) 
						+ a(t-1, i-1, j)) + 0.375 * (a(t-1, i, j+1) 
						- 2.0 * a(t-1, i, j) + a(t-1, i, j-1)) + 
						a(t-1, i, j);
		}
		if (i >= N1/2 && j < N2/2)
		{
			a(t, i, j) = 0.5 * (a(t-1, i+1, j) - 2.0 * a(t-1, i, j) 
						+ a(t-1, i-1, j)) + 0.5 * (a(t-1, i, j+1) 
						- 2.0 * a(t-1, i, j) + a(t-1, i, j-1)) + 
						a(t-1, i, j);
		}
	Pochoir_Kernel_2D_End(heat_2D_fn, heat_shape_2D)
	heat_2D_single_kernel.Register_Tile_Kernels(Default_Guard_2D, heat_2D_fn);
	heat_2D_single_kernel.Register_Array(a);
	/* Now we can only access the Pochoir_Array after Register_Array,
	 * or Register_Shape with the array, because we rely on the shape
	 * to provide the depth of toggle array!!! 
	 */
	for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
		double tmp = 1.0 * (rand() / BASE); 
		a(0, i, j) = tmp; 
		b(0, i, j) = tmp;
		a(1, i, j) = tmp; 
		b(1, i, j) = tmp;
	} }

	Pochoir_Plan<2> & l_plan = heat_2D_single_kernel.Gen_Plan(T);
	//sprintf(pochoir_plan_file_name, "pochoir_%d_%d.dat\0", N, T);
	//heat_2D_single_kernel.Store_Plan(pochoir_plan_file_name, l_plan);
	// Pochoir_Plan<2> & ll_plan = heat_2D_single_kernel.Load_Plan(pochoir_plan_file_name);
	gettimeofday(&start, 0);
	for (int times = 0; times < TIMES; ++times) {
		heat_2D_single_kernel.Run(l_plan);
	}
	gettimeofday(&end, 0);
	std::cout << "Pochoir consumed time :" << 1.0e3 * tdiff(&end, &start)/TIMES << " ms" << std::endl;

	/* Do not check the correctness of the result */
	/*
	gettimeofday(&start, 0);
	for (int times = 0; times < TIMES; ++times) {
	for (int t = 0; t < T; ++t) {
	cilk_for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
		b(t+1, i, j) = 0.125 * (b(t, i+1, j) - 2.0 * b(t, i, j) + b(t, i-1, j)) + 0.125 * (b(t, i, j+1) - 2.0 * b(t, i, j) + b(t, i, j-1)) + b(t, i, j); } } }
	}
	gettimeofday(&end, 0);
	std::cout << "Parallel Loop: consumed time :" << 1.0e3 * tdiff(&end, &start)/TIMES << " ms" << std::endl;

	t = T;
	for (int i = 0; i < N1; ++i) {
	for (int j = 0; j < N2; ++j) {
		check_result(t, i, j, a.interior(t, i, j), b.interior(t, i, j));
	} } 
	printf("max_diff = %f, when a = %f, b = %f\n", max_diff, diff_a, diff_b);
	printf("max_a = %f, max_b = %f\n", max_a, max_b);
	*/
	}

	return 0;
}
