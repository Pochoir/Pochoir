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

/* Test bench - employ Pochoir for dynamic programming - fibonacci number */
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

void check_result(int n, long a, long b)
{
	if (abs(a - b) < TOLERANCE) {
		printf("fib(%d) = %ld == %ld: Succeed!\n", n, a, b);
	} else {
		printf("fib(%d) = %ld != %ld: FAILED!\n", n, a, b);
	}

}

long lfib(long n) {
    if (n <= 1)
        return n;
    long a = cilk_spawn lfib(n-1);
    long b = lfib(n-2);
    cilk_sync;
    return (a+b);
}

int main(int argc, char * argv[])
{
	const int BASE = 1024;
	int t;
	struct timeval start, end;
    double min_tdiff = INF;
    int N = 0;

    if (argc < 2) {
        printf("argc < 3, quit! \n");
        exit(1);
    }
    N = StrToInt(argv[1]);
    printf("N = %d\n", N);
	/* data structure of Pochoir - row major */
    Pochoir_Shape_1D fib_shape[] = {{1, 0}, {0, 0}, {-1, 0}};
	Pochoir_Array_1D(long) a(1);
    Pochoir_1D fib(fib_shape);

	cout << "fib(" << N << ")" << endl;
    Pochoir_Kernel_1D(fib_fn, t, i)
	   a(t+1, 0) = a(t, 0) + a(t-1, 0);
    Pochoir_Kernel_End

    // a.Register_Boundary(aperiodic_1D);
    fib.Register_Array(a);

    a(0, 0) = 0; 
    a(1, 0) = 1; 

	gettimeofday(&start, 0);
    fib.Run(N, fib_fn);
	gettimeofday(&end, 0);
    min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
	std::cout << "Pochoir ET: consumed time :" << min_tdiff << "ms" << std::endl;

    min_tdiff = INF;
	gettimeofday(&start, 0);
    long b = lfib(N);
	gettimeofday(&end, 0);
    min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));

	std::cout << "Naive Loop: consumed time :" << min_tdiff << "ms" << std::endl;

	check_result(N, a.interior(N, 0), b);

	return 0;
}
