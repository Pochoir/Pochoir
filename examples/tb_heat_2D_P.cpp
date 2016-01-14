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

/* Test bench - 2D heat equation, Periodic version */
#include <cstdio>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>

//#define OPENTUNER
#ifdef OPENTUNER
  int opentuner_p[6];
#endif
#include <pochoir.hpp>

using namespace std;
#define TIMES 1
#define N_RANK 2
#define TOLERANCE (1e-6)
//#define CHECK_RESULT


void check_result(int t, int j, int i, double a, double b)
{
  if (abs(a - b) < TOLERANCE) {
    //printf("a(%d, %d, %d) == b(%d, %d, %d) == %f : passed!\n", t, j, i, t, j, i, a);
  } else {
    printf("a(%d, %d, %d) = %f, b(%d, %d, %d) = %f : FAILED!\n", t, j, i, a, t, j, i, b);
  }
}

Pochoir_Boundary_2D(aperiodic_2D, arr, t, i, j)
    return 0;
Pochoir_Boundary_End

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

int main(int argc, char * argv[])
{
  const int BASE = 1024;
  int t;
  //struct timeval start, end;
  struct timespec start, end;
  int N1 = 500, N2 = 100, T_SIZE = 731;

  if (argc < 4) {
    printf("argc < 4, quit! \n");
    exit(1);
  }
  N1 = StrToInt(argv[1]);
  N2 = StrToInt(argv[2]);
  T_SIZE = StrToInt(argv[3]);

#ifdef OPENTUNER
  if (argc < 10) {
    printf("heat_2d_p argc < 10, quit! \n");
    exit(1);
  }
  for (int i = 0 ; i < 6 ; i++) {
    opentuner_p [i] = StrToInt(argv[i + 4]);
  }
  printf("N1 = %d,N2 = %d, T_SIZE = %d, dt = %d, dx0 = %d, dx1 = %d, dt_b = %d, dx0_b = %d, dx1_b = %d \n", N1, N2, T_SIZE, opentuner_p [0], opentuner_p [1], opentuner_p [2], opentuner_p [3], opentuner_p [4], opentuner_p [5]);
#endif

  Pochoir_Shape_2D heat_shape_2D[] = {{0, 0, 0}, {-1, 1, 0}, {-1, 0, 0}, {-1, -1, 0}, {-1, 0, -1}, {-1, 0, 1}};
  //Pochoir_Shape_2D heat_shape_2D[] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 0}, {0, -1, 0}, {0, 0, -1}, {0, 0, 1}};
  Pochoir<N_RANK> heat_2D(heat_shape_2D);
  Pochoir_Array<double, N_RANK> a(N1, N2) ;
  a.Register_Boundary(periodic_2D);
  heat_2D.Register_Array(a);

#ifdef CHECK_RESULT
  Pochoir_Array<double, N_RANK> b(N1, N2);
  b.Register_Shape(heat_shape_2D);
  b.Register_Boundary(periodic_2D);
  double * c = new double [N1 * N2 * 2] ;
  long const area = N1 * N2 ;
#define c(t,i,j) c[t * area + i * N2 + j]
#endif

  /* Now we can only access the Pochoir_Array after Register_Array,
  * or Register_Shape with the array, because we rely on the shape
  * to provide the depth of toggle array!!! 
  */
  for (int i = 0; i < N1; ++i) {
    for (int j = 0; j < N2; ++j) {
      a(0, i, j) = 1.0 * (rand() % 2); 
      a(1, i, j) = 0; 
#ifdef CHECK_RESULT
      b(0, i, j) = a(0, i, j);
      b(1, i, j) = 0;
      c(0, i, j) = a(0, i, j);
      c(1, i, j) = 0;
#endif
  } }

  Pochoir_Kernel_2D(heat_2D_fn, t, i, j)
        a(t, i, j) = 0.125 * (a(t-1, i+1, j) - 2.0 * a(t-1, i, j) + a(t-1, i-1, j)) + 0.125 * (a(t-1, i, j+1) - 2.0 * a(t-1, i, j) + a(t-1, i, j-1)) + a(t-1, i, j);
        //a(t + 1, i, j) = 0.125 * (a(t, i+1, j) - 2.0 * a(t, i, j) + a(t, i-1, j)) + 0.125 * (a(t, i, j+1) - 2.0 * a(t, i, j) + a(t, i, j-1)) + a(t, i, j);
  Pochoir_Kernel_End

  char name [500] ;
  sprintf(name, "heat_2D_P") ;
  heat_2D.set_problem_name(name) ;
  clock_gettime(CLOCK_MONOTONIC, &start) ;
  for (int times = 0; times < TIMES; ++times) {
    heat_2D.Run(T_SIZE, heat_2D_fn);
  }
  clock_gettime(CLOCK_MONOTONIC, &end) ;
  std::cout << "Pochoir ET: consumed time :" << tdiff2(&end, &start)/TIMES * 1.0e-6 << "ms" << std::endl;

  printf("N1 = %d, N2 = %d, T_SIZE = %d\n", N1, N2, T_SIZE);
  cout << "a(T+1, J, I) = 0.125 * (a(T, J+1, I) - 2.0 * a(T, J, I) + a(T, J-1, I)) + 0.125 * (a(T, J, I+1) - 2.0 * a(T, J, I) + a(T, J, I-1)) + a(T, J, I)" << endl;
#ifdef CHECK_RESULT
  /* 
  clock_gettime(CLOCK_MONOTONIC, &start) ;
  for (int times = 0; times < TIMES; ++times) {
    for (int t = 0; t < T_SIZE; ++t) {
      for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
          b(t+1, i, j) = 0.125 * (b(t, i+1, j) - 2.0 * b(t, i, j) + b(t, i-1, j)) + 0.125 * (b(t, i, j+1) - 2.0 * b(t, i, j) + b(t, i, j-1)) + b(t, i, j); 
        } 
      } 
    }
  }
  clock_gettime(CLOCK_MONOTONIC, &end) ;
  //std::cout << "Naive Loop: consumed time :" << 1.0e3 * tdiff(&end, &start)/TIMES << "ms" << std::endl;
  std::cout << "Naive Loop: consumed time :" << tdiff2(&end, &start)/TIMES * 1.0e-6 << "ms" << std::endl;

  t = T_SIZE;
  for (int i = 0; i < N1; ++i) {
    for (int j = 0; j < N2; ++j) {
      check_result(t, i, j, a.interior(t, i, j), b.interior(t, i, j));
    } 
  }
  */
  
  int _M=N1-1 ;
  int _N=N2-1 ;
  int _T=T_SIZE ;

  clock_gettime(CLOCK_MONOTONIC, &start) ;
  for (int t = 0; t < _T; t++) {
    for (int i = 0; i <= _M; i++) {
#pragma ivdep 
      for (int j = 0; j <= _N; j++) {
        c(((t+1)%2),i,j) = 0.125 * (c((t%2),(i==_M?0:i+1),j) - 2.0 * c((t%2),i,j) + c((t%2),(i==0?_M:(i-1)),j)) + 0.125 * (c((t%2),i,(j==_N?0:j+1)) - 2.0 * c((t%2),i,j) + c((t%2),i,(j==0?_N:j-1))) + c((t%2),i,j) ;
      }
    }
  }
  clock_gettime(CLOCK_MONOTONIC, &end) ;
  std::cout << "non pochoir loop: consumed time :" << tdiff2(&end, &start)/TIMES * 1.0e-6 << "ms" << std::endl;
  delete [] c ;
#endif
  return 0;
}
