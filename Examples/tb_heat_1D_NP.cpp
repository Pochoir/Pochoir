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

/* Test bench - 1D heat equation, Non-periodic version */
#include <cstdio>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>
#include <vector>

#include <pochoir.hpp>

using namespace std;
#define TIMES 1
#define TOLERANCE (1e-6)
#define NAIVE

void pochoir_heat_1D_NP(double* seed, int N_SIZE, int T_SIZE, double* res);
void naive_pochoir_heat_1D_NP(double* seed, int N_SIZE, int T_SIZE, double* res);
void naive_heat_1D_NP(double* seed, int N_SIZE, int T_SIZE, double* res);

Pochoir_Boundary_1D(heat_bv_1D, arr, t, i)
  return 0;
Pochoir_Boundary_End

void check_result(int t, int i, double a, double b)
{
	if (a == b) {
//		printf("a(%d, %d) == b(%d, %d) == %f : passed!\n", t, i, t, i, a);
	} else {
		printf("a(%d, %d) = %f, b(%d, %d) = %f : FAILED!\n", t, i, a, t, i, b);
	}
}

int main(int argc, char * argv[])
{
	const int BASE = 1024;
	int t;
	struct timeval start, end;
  double min_tdiff = INF;
  int N_SIZE = 0, T_SIZE = 0;

  if (argc < 3) {
      printf("Usage: heat_1D_NP N_SIZE T_SIZE \n");
      exit(1);
  }
  
  N_SIZE = StrToInt(argv[1]);
  T_SIZE = StrToInt(argv[2]);
  printf("N_SIZE = %d, T_SIZE = %d\n", N_SIZE, T_SIZE);

	cout << 
  "a(T+1, I) = 0.125 * (a(T, I+1) - 2.0 * a(T, I) + a(T, I-1))" 
  << endl;
  
  // generate initial seed values
  double* seed = (double*) malloc(N_SIZE * sizeof(double));

	for (int i = 0; i < N_SIZE; ++i) {
    seed[i] = 1.0 * (rand() % BASE);
	}

  // compute various implementations
  double* naive = (double*) malloc(N_SIZE * sizeof(double));
  naive_heat_1D_NP(seed, N_SIZE, T_SIZE, naive);

  double* naive_pochoir = (double*) malloc(N_SIZE * sizeof(double));
  naive_pochoir_heat_1D_NP(seed, N_SIZE, T_SIZE, naive_pochoir);

  double* pochoir = (double*) malloc(N_SIZE * sizeof(double));
  pochoir_heat_1D_NP(seed, N_SIZE, T_SIZE, pochoir);
  
	for (int i = 0; i < N_SIZE; ++i) {
		check_result(T_SIZE, i, naive[i], naive_pochoir[i]);
	}

  for (int i = 0; i < N_SIZE; ++i) {
    check_result(T_SIZE, i, naive[i], pochoir[i]);
  }

	return 0;
}

/*
 * Pochoir implementation
 */
void pochoir_heat_1D_NP(double* seed, int N_SIZE, int T_SIZE, double* res) {
  struct timeval start, end;
  const int DIM = N_SIZE + 2;
  double min_tdiff = INF;

  Pochoir_Shape_1D heat_shape_1D[] = {{1, 0}, {0, 1}, {0, -1}, {0, 0}};
  Pochoir_1D heat_1D(heat_shape_1D);

  Pochoir_Array_1D(double) temp(N_SIZE);

  Pochoir_Kernel_1D(heat_1D_fn, t, i)
    temp(t+1, i) = 0.125 * (temp(t, i+1) - 2.0 * temp(t, i) + temp(t, i-1));
  Pochoir_Kernel_End

  temp.Register_Boundary(heat_bv_1D);
  heat_1D.Register_Array(temp);

  for (int i = 0; i < N_SIZE; ++i) {
    temp(0, i) = seed[i];
    temp(1, i) = 0;
  }

  for (int times = 0; times < TIMES; ++times) {
    gettimeofday(&start, 0);
    heat_1D.Run(T_SIZE, heat_1D_fn);
    gettimeofday(&end, 0);
    min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
  }
  std::cout << "Pochoir ET: consumed time :" << min_tdiff << "ms" << std::endl;

  for (int i = 0; i < N_SIZE; ++i) {
    res[i] = temp(T_SIZE, i);
  }
}

/*
 * Naive implementation using pochoir structures
 */
void naive_pochoir_heat_1D_NP(double* seed, int N_SIZE, int T_SIZE, double* res) {
  struct timeval start, end;
  const int DIM = N_SIZE + 2;
  double min_tdiff = INF;

  Pochoir_Shape_1D heat_shape_1D[] = {{1, 0}, {0, 1}, {0, -1}, {0, 0}};
  Pochoir_1D heat_1D(heat_shape_1D);

  Pochoir_Array_1D(double) temp(N_SIZE+2);
  temp.Register_Shape(heat_shape_1D);

  for (int i = 0; i < N_SIZE; ++i) {
    temp(0, i+1) = seed[i];
    temp(1, i+1) = 0;
  }

  /* cilk_for + zero-padding */
  for (int times = 0; times < TIMES; ++times) {
    gettimeofday(&start, 0);
    for (int t = 0; t < T_SIZE; ++t) {
      cilk_for (int i = 1; i < N_SIZE+1; ++i) {
        temp.interior(t+1, i) = 
          0.125 * (
            temp.interior(t, i+1) 
            - 2.0 * temp.interior(t, i) 
            + temp.interior(t, i-1)
          ); 
      } 
    }
    gettimeofday(&end, 0);
    min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
  }
  std::cout << "Naive Pochoir Loop: consumed time :" << min_tdiff << "ms" << std::endl;

  // copy over result
  for (int i = 0; i < N_SIZE; ++i) {
    res[i] = temp(T_SIZE, i+1);
  }
}

/*
 * Naive implementation using purely arrays (no pochoir references)
*/
void naive_heat_1D_NP(double* seed, int N_SIZE, int T_SIZE, double* res) {
  struct timeval start, end;
  const int DIM = N_SIZE + 2;
  double min_tdiff = INF;

  // initialize empty array
  double *temp = (double*) malloc(2 * DIM * sizeof(double));
  for (int i = 0; i < 2 * DIM; ++i) {
    temp[i] = 0;
  }

  // copy in seed elements
  for (int i = 0; i < N_SIZE; ++i) {
    temp[i+1] = seed[i];
  }

  // main loop, naive non-Pochoir implementation
  for (int times = 0; times < TIMES; ++times) {
    gettimeofday(&start, 0);
    for (int t = 0; t < T_SIZE; ++t) {
      cilk_for (int i = 1; i < N_SIZE+1; ++i) {
        temp[((t+1)%2) * DIM + i] = 
          0.125 * (
            temp[(t%2)*DIM + i + 1] 
            - 2.0 * temp[(t%2)*DIM + i]
            + temp[(t%2)*DIM + i - 1]
          );
      } 
    }
    gettimeofday(&end, 0);
    min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
  }
  std::cout << "Naive Loop: consumed time :" << min_tdiff << "ms" << std::endl;

  // copy over result
  for (int i = 0; i < N_SIZE; ++i) {
    res[i] = temp[(T_SIZE%2)*DIM + i + 1];
  }

  free(temp);
}