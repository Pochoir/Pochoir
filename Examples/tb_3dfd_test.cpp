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

/* It's order-4, 3D 15 point stencil, to match up with Matteo Frigo's
 * hand-optimized wave equation 
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

using namespace std;

int nthreads = 1;
const int ds = 4;
int Nx = 100;
int Ny = 100;
int Nz = 100;
int T = 40;
float ** A;
float coef[ds+1];
float *vsq;
int sx2, sx3, sx4;
int sxy2, sxy3, sxy4;

#define TOLERANCE (1e-6)
static inline double tdiff (struct timeval *a, struct timeval *b)
{
        return a->tv_sec - b->tv_sec + 1e-6 * (a->tv_usec - b->tv_usec);
}   

void print_summary(char *header, double interval) {
  /* print timing information */
  long total = (long)Nx * Ny * Nz;
//  int n_worker = cilk::current_worker_count();
  int n_worker = nthreads;
  printf("++++++++++++++++++++ %s ++++++++++++++++++++++\n", header);
  printf("first non-zero numbers\n");
  for(int i = 0; i < total; i++) {
    if(A[T%2][i] != 0) {
      printf("%d: %f\n", i, A[T%2][i]);
      break;
    }
  }
	
  long mul = (long)(Nx - 8) * (Ny  - 8) * (Nz - 8) * T;
  double perf = mul / (interval * 1e6);
  printf("time: %f\n", interval);
  printf("Perf: %f Mcells/sec (%f M-FAdd/s, %f M-FMul/s)\n", 
	 perf, 
	 perf * 26, 
	 perf * 7);
  printf("Perf per worker: %f Mcells/sec (%f M-FAdd/s, %f M-FMul/s)\n\n", 
	 perf / n_worker, 
	 perf * 26 / n_worker, 
	 perf * 7 / n_worker);
  //printf("count = %d\n\n", count);		
}

void check_result(int t, int i, int j, int k, double a, double b)
{
    if (abs(a - b) < TOLERANCE) {
//      printf("a(%d, %d, %d, %d) == b(%d, %d, %d, %d) == %f : passed!\n", t, i, j, k, t, i, j, k, a);
    } else {
        printf("a(%d, %d, %d, %d) = %f, b(%d, %d, %d, %d) = %f : FAILED!\n", t, i, j, k, a, t, i, j, k, b);
    }

}

int main(int argc, char *argv[])
{
  struct timeval start, end;
  if (argc > 3) {
    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    Nz = atoi(argv[3]);
  }
  /* T is time steps */
  if (argc > 4)
    T = atoi(argv[4]);

  printf("Order-%d 3D-Stencil (%d points) with space %dx%dx%d and time %d\n", 
	 ds, ds*2*3+1, Nx, Ny, Nz, T);
  const int Nxyz = Nx * Ny * Nz, Nxy = Nx * Ny;
  const int half_Nx = (Nx >> 1) << 3, half_Ny = Ny >> 1, half_Nz = Nz >> 1;
  const int half_Nxy = half_Nx * half_Ny;
  //initialization
  A = new float*[2];
  A[0] = new float[Nxyz];
  A[1] = new float[Nxyz];
  float * vsq = new float [Nxyz];
  float * pb = new float[2 * Nxyz];
  float * pc = new float[2 * Nxyz];
  // float vsq[Nz][Ny][Nx];
  // float pb[2][Nz/2][Ny/2][Nx/2][2][2][2];
  // float pc[2][Nz][Ny][Nx];
#define acc_pb(t, i, j, k) pb[((t)&0x1)*Nxyz + ((i)>>1) * half_Nxy+ ((j)>>1) * half_Nx + (((k)>>1) << 3) + (((i) & 0x1) << 2) + (((j) & 0x1) << 1) + (k) & 0x1]
#define acc_pbl(t, i, j, k, ii, jj, kk) pb[((t)&0x1)*Nxyz + (i) * half_Nxy + (j) * half_Nx + ((k) << 3) + ((ii) << 2) + ((jj) << 1) + (kk)]
#define acc_pc(t, i, j, k) pc[((t)&0x1)*Nxyz + (i) * Nxy + (j) * Nx + (k)]
#define acc_vsq(i, j, k) vsq[(i) * Nxy + (j) * Nx + (k)]
#define pabs(a) (a) < 0 ? -(a) : (a)
#define max(a, b) ((a) > (b) ? (a) : (b))
#define PB 0
#define PBL 1
#define PC 0

  sx2 = Nx * 2;
  sx3 = Nx * 3;
  sx4 = Nx * 4;
  sxy2 = Nxy * 2;
  sxy3 = Nxy * 3;
  sxy4 = Nxy * 4;

  coef[4] = -1.0f / 560.0f;
  coef[3] = 8.0f/315;
  coef[2] = -0.2f;
  coef[1] = 1.6f;
  coef[0] = -1435.0f/504 * 3;

  cilk_for (int z = 0; z < Nz; ++z)
    for (int y = 0; y < Ny; ++y) 
      for(int x = 0; x < Nx; ++x) {
	/* set initial values */
	float r = pabs((x - Nx/2 + y - Ny/2 + z - Nz/2) / 30);
	r = max(1 - r, 0.0f) + 1;
	
#if PB
	acc_pb(0, z, y, x) = r;
	acc_pb(1, z, y, x) = r;
#endif
#if PBL
	acc_pbl(0, z>>1, y>>1, x>>1, 0, 0, 0) = r;
	acc_pbl(0, z>>1, y>>1, x>>1, 0, 0, 1) = r;
	acc_pbl(0, z>>1, y>>1, x>>1, 0, 1, 0) = r;
	acc_pbl(0, z>>1, y>>1, x>>1, 0, 1, 1) = r;
	acc_pbl(0, z>>1, y>>1, x>>1, 1, 0, 0) = r;
	acc_pbl(0, z>>1, y>>1, x>>1, 1, 0, 1) = r;
	acc_pbl(0, z>>1, y>>1, x>>1, 1, 1, 0) = r;
	acc_pbl(0, z>>1, y>>1, x>>1, 1, 1, 1) = r;
	acc_pbl(1, z>>1, y>>1, x>>1, 0, 0, 0) = r;
	acc_pbl(1, z>>1, y>>1, x>>1, 0, 0, 1) = r;
	acc_pbl(1, z>>1, y>>1, x>>1, 0, 1, 0) = r;
	acc_pbl(1, z>>1, y>>1, x>>1, 0, 1, 1) = r;
	acc_pbl(1, z>>1, y>>1, x>>1, 1, 0, 0) = r;
	acc_pbl(1, z>>1, y>>1, x>>1, 1, 0, 1) = r;
	acc_pbl(1, z>>1, y>>1, x>>1, 1, 1, 0) = r;
	acc_pbl(1, z>>1, y>>1, x>>1, 1, 1, 1) = r;
#endif
#if PC
	acc_pc(0, z, y, x) = r;
	acc_pc(1, z, y, x) = r;
#endif
	acc_vsq(z, y, x) = 0.001f;
  }
  const int start_x = ds, end_x = (Nx - ds);
  const int start_y = ds, end_y = (Ny - ds);
  const int start_z = ds, end_z = (Nz - ds);
  const int half_start_x = ds/2, half_end_x = (Nx - ds)/2;
  const int half_start_y = ds/2, half_end_y = (Ny - ds)/2;
  const int half_start_z = ds/2, half_end_z = (Nz - ds)/2;
#if PB
  gettimeofday(&start, 0);
  for (int t = 0; t < T; ++t) {
     cilk_for (int i = start_z; i < end_z; i+=2) {
        for (int j = start_y; j < end_y; j+=2) {
#pragma ivdep
            for (int k = start_x; k < end_x; k+=2) {
          const float c0 = coef[0], c1 = coef[1], c2 = coef[2], c3 = coef[3], c4 = coef[4];
          float div = c0 * acc_pb(t, i, j, k) + 
                c1 * ((acc_pb(t, i, j, k+1) + acc_pb(t, i, j, k-1)) 
                    + (acc_pb(t, i, j+1, k) + acc_pb(t, i, j-1, k)) 
                    + (acc_pb(t, i+1, j, k) + acc_pb(t, i-1, j, k))) 
              + c2 * ((acc_pb(t, i, j, k+2) + acc_pb(t, i, j, k-2)) 
                    + (acc_pb(t, i, j+2, k) + acc_pb(t, i, j-2, k)) 
                    + (acc_pb(t, i+2, j, k) + acc_pb(t, i-2, j, k))) 
              + c3 * ((acc_pb(t, i, j, k+3) + acc_pb(t, i, j, k-3)) 
                    + (acc_pb(t, i, j+3, k) + acc_pb(t, i, j-3, k)) 
                    + (acc_pb(t, i+3, j, k) + acc_pb(t, i-3, j, k))) 
              + c4 * ((acc_pb(t, i, j, k+4) + acc_pb(t, i, j, k-4)) 
                    + (acc_pb(t, i, j+4, k) + acc_pb(t, i, j-4, k)) 
                    + (acc_pb(t, i+4, j, k) + acc_pb(t, i-4, j, k)));
     acc_pb(t+1, i, j, k) = 2 * acc_pb(t, i, j, k) - acc_pb(t+1, i, j, k) + acc_vsq(i, j, k) * div;
          div = c0 * acc_pb(t, i, j, k+1) + 
                c1 * ((acc_pb(t, i, j, k+1+1) + acc_pb(t, i, j, k-1+1)) 
                    + (acc_pb(t, i, j+1, k+1) + acc_pb(t, i, j-1, k+1)) 
                    + (acc_pb(t, i+1, j, k+1) + acc_pb(t, i-1, j, k+1))) 
              + c2 * ((acc_pb(t, i, j, k+2+1) + acc_pb(t, i, j, k-2+1)) 
                    + (acc_pb(t, i, j+2, k+1) + acc_pb(t, i, j-2, k+1)) 
                    + (acc_pb(t, i+2, j, k+1) + acc_pb(t, i-2, j, k+1))) 
              + c3 * ((acc_pb(t, i, j, k+3+1) + acc_pb(t, i, j, k-3+1)) 
                    + (acc_pb(t, i, j+3, k+1) + acc_pb(t, i, j-3, k+1)) 
                    + (acc_pb(t, i+3, j, k+1) + acc_pb(t, i-3, j, k+1))) 
              + c4 * ((acc_pb(t, i, j, k+4+1) + acc_pb(t, i, j, k-4+1)) 
                    + (acc_pb(t, i, j+4, k+1) + acc_pb(t, i, j-4, k+1)) 
                    + (acc_pb(t, i+4, j, k+1) + acc_pb(t, i-4, j, k+1)));
     acc_pb(t+1, i, j, k+1) = 2 * acc_pb(t, i, j, k+1) - acc_pb(t+1, i, j, k+1) + acc_vsq(i, j, k+1) * div;
          div = c0 * acc_pb(t, i, j+1, k) + 
                c1 * ((acc_pb(t, i, j+1, k+1) + acc_pb(t, i, j+1, k-1)) 
                    + (acc_pb(t, i, j+1+1, k) + acc_pb(t, i, j-1+1, k)) 
                    + (acc_pb(t, i+1, j+1, k) + acc_pb(t, i-1, j+1, k))) 
              + c2 * ((acc_pb(t, i, j+1, k+2) + acc_pb(t, i, j+1, k-2)) 
                    + (acc_pb(t, i, j+2+1, k) + acc_pb(t, i, j-2+1, k)) 
                    + (acc_pb(t, i+2, j+1, k) + acc_pb(t, i-2, j+1, k))) 
              + c3 * ((acc_pb(t, i, j+1, k+3) + acc_pb(t, i, j+1, k-3)) 
                    + (acc_pb(t, i, j+3+1, k) + acc_pb(t, i, j-3+1, k)) 
                    + (acc_pb(t, i+3, j+1, k) + acc_pb(t, i-3, j+1, k))) 
              + c4 * ((acc_pb(t, i, j+1, k+4) + acc_pb(t, i, j+1, k-4)) 
                    + (acc_pb(t, i, j+4+1, k) + acc_pb(t, i, j-4+1, k)) 
                    + (acc_pb(t, i+4, j+1, k) + acc_pb(t, i-4, j+1, k)));
     acc_pb(t+1, i, j+1, k) = 2 * acc_pb(t, i, j+1, k) - acc_pb(t+1, i, j+1, k) + acc_vsq(i, j+1, k) * div;
          div = c0 * acc_pb(t, i, j+1, k+1) + 
                c1 * ((acc_pb(t, i, j+1, k+1+1) + acc_pb(t, i, j+1, k-1+1)) 
                    + (acc_pb(t, i, j+1+1, k+1) + acc_pb(t, i, j-1+1, k+1)) 
                    + (acc_pb(t, i+1, j+1, k+1) + acc_pb(t, i-1, j+1, k+1))) 
              + c2 * ((acc_pb(t, i, j+1, k+2+1) + acc_pb(t, i, j+1, k-2+1)) 
                    + (acc_pb(t, i, j+2+1, k+1) + acc_pb(t, i, j-2+1, k+1)) 
                    + (acc_pb(t, i+2, j+1, k+1) + acc_pb(t, i-2, j+1, k+1))) 
              + c3 * ((acc_pb(t, i, j+1, k+3+1) + acc_pb(t, i, j+1, k-3+1)) 
                    + (acc_pb(t, i, j+3+1, k+1) + acc_pb(t, i, j-3+1, k+1)) 
                    + (acc_pb(t, i+3, j+1, k+1) + acc_pb(t, i-3, j+1, k+1))) 
              + c4 * ((acc_pb(t, i, j+1, k+4+1) + acc_pb(t, i, j+1, k-4+1)) 
                    + (acc_pb(t, i, j+4+1, k+1) + acc_pb(t, i, j-4+1, k+1)) 
                    + (acc_pb(t, i+4, j+1, k+1) + acc_pb(t, i-4, j+1, k+1)));
     acc_pb(t+1, i, j+1, k+1) = 2 * acc_pb(t, i, j+1, k+1) - acc_pb(t+1, i, j+1, k+1) + acc_vsq(i, j+1, k+1) * div;

          div = c0 * acc_pb(t, i+1, j, k) + 
                c1 * ((acc_pb(t, i+1, j, k+1) + acc_pb(t, i+1, j, k-1)) 
                    + (acc_pb(t, i+1, j+1, k) + acc_pb(t, i+1, j-1, k)) 
                    + (acc_pb(t, i+1+1, j, k) + acc_pb(t, i-1+1, j, k))) 
              + c2 * ((acc_pb(t, i+1, j, k+2) + acc_pb(t, i+1, j, k-2)) 
                    + (acc_pb(t, i+1, j+2, k) + acc_pb(t, i+1, j-2, k)) 
                    + (acc_pb(t, i+2+1, j, k) + acc_pb(t, i-2+1, j, k))) 
              + c3 * ((acc_pb(t, i+1, j, k+3) + acc_pb(t, i+1, j, k-3)) 
                    + (acc_pb(t, i+1, j+3, k) + acc_pb(t, i+1, j-3, k)) 
                    + (acc_pb(t, i+3+1, j, k) + acc_pb(t, i-3+1, j, k))) 
              + c4 * ((acc_pb(t, i+1, j, k+4) + acc_pb(t, i+1, j, k-4)) 
                    + (acc_pb(t, i+1, j+4, k) + acc_pb(t, i+1, j-4, k)) 
                    + (acc_pb(t, i+4+1, j, k) + acc_pb(t, i-4+1, j, k)));
     acc_pb(t+1, i+1, j, k) = 2 * acc_pb(t, i+1, j, k) - acc_pb(t+1, i+1, j, k) + acc_vsq(i+1, j, k) * div;
          div = c0 * acc_pb(t, i+1, j, k+1) + 
                c1 * ((acc_pb(t, i+1, j, k+1+1) + acc_pb(t, i+1, j, k-1+1)) 
                    + (acc_pb(t, i+1, j+1, k+1) + acc_pb(t, i+1, j-1, k+1)) 
                    + (acc_pb(t, i+1+1, j, k+1) + acc_pb(t, i-1+1, j, k+1))) 
              + c2 * ((acc_pb(t, i+1, j, k+2+1) + acc_pb(t, i+1, j, k-2+1)) 
                    + (acc_pb(t, i+1, j+2, k+1) + acc_pb(t, i+1, j-2, k+1)) 
                    + (acc_pb(t, i+2+1, j, k+1) + acc_pb(t, i-2+1, j, k+1))) 
              + c3 * ((acc_pb(t, i+1, j, k+3+1) + acc_pb(t, i+1, j, k-3+1)) 
                    + (acc_pb(t, i+1, j+3, k+1) + acc_pb(t, i+1, j-3, k+1)) 
                    + (acc_pb(t, i+3+1, j, k+1) + acc_pb(t, i-3+1, j, k+1))) 
              + c4 * ((acc_pb(t, i+1, j, k+4+1) + acc_pb(t, i+1, j, k-4+1)) 
                    + (acc_pb(t, i+1, j+4, k+1) + acc_pb(t, i+1, j-4, k+1)) 
                    + (acc_pb(t, i+4+1, j, k+1) + acc_pb(t, i-4+1, j, k+1)));
     acc_pb(t+1, i+1, j, k+1) = 2 * acc_pb(t, i+1, j, k+1) - acc_pb(t+1, i+1, j, k+1) + acc_vsq(i+1, j, k+1) * div;
          div = c0 * acc_pb(t, i+1, j+1, k) + 
                c1 * ((acc_pb(t, i+1, j+1, k+1) + acc_pb(t, i+1, j+1, k-1)) 
                    + (acc_pb(t, i+1, j+1+1, k) + acc_pb(t, i+1, j-1+1, k)) 
                    + (acc_pb(t, i+1+1, j+1, k) + acc_pb(t, i-1+1, j+1, k))) 
              + c2 * ((acc_pb(t, i+1, j+1, k+2) + acc_pb(t, i+1, j+1, k-2)) 
                    + (acc_pb(t, i+1, j+2+1, k) + acc_pb(t, i+1, j-2+1, k)) 
                    + (acc_pb(t, i+2+1, j+1, k) + acc_pb(t, i-2+1, j+1, k))) 
              + c3 * ((acc_pb(t, i+1, j+1, k+3) + acc_pb(t, i+1, j+1, k-3)) 
                    + (acc_pb(t, i+1, j+3+1, k) + acc_pb(t, i+1, j-3+1, k)) 
                    + (acc_pb(t, i+3+1, j+1, k) + acc_pb(t, i-3+1, j+1, k))) 
              + c4 * ((acc_pb(t, i+1, j+1, k+4) + acc_pb(t, i+1, j+1, k-4)) 
                    + (acc_pb(t, i+1, j+4+1, k) + acc_pb(t, i+1, j-4+1, k)) 
                    + (acc_pb(t, i+4+1, j+1, k) + acc_pb(t, i-4+1, j+1, k)));
     acc_pb(t+1, i+1, j+1, k) = 2 * acc_pb(t, i+1, j+1, k) - acc_pb(t+1, i+1, j+1, k) + acc_vsq(i+1, j+1, k) * div;
          div = c0 * acc_pb(t, i+1, j+1, k+1) + 
                c1 * ((acc_pb(t, i+1, j+1, k+1+1) + acc_pb(t, i+1, j+1, k-1+1)) 
                    + (acc_pb(t, i+1, j+1+1, k+1) + acc_pb(t, i+1, j-1+1, k+1)) 
                    + (acc_pb(t, i+1+1, j+1, k+1) + acc_pb(t, i-1+1, j+1, k+1))) 
              + c2 * ((acc_pb(t, i+1, j+1, k+2+1) + acc_pb(t, i+1, j+1, k-2+1)) 
                    + (acc_pb(t, i+1, j+2+1, k+1) + acc_pb(t, i+1, j-2+1, k+1)) 
                    + (acc_pb(t, i+2+1, j+1, k+1) + acc_pb(t, i-2+1, j+1, k+1))) 
              + c3 * ((acc_pb(t, i+1, j+1, k+3+1) + acc_pb(t, i+1, j+1, k-3+1)) 
                    + (acc_pb(t, i+1, j+3+1, k+1) + acc_pb(t, i+1, j-3+1, k+1)) 
                    + (acc_pb(t, i+3+1, j+1, k+1) + acc_pb(t, i-3+1, j+1, k+1))) 
              + c4 * ((acc_pb(t, i+1, j+1, k+4+1) + acc_pb(t, i+1, j+1, k-4+1)) 
                    + (acc_pb(t, i+1, j+4+1, k+1) + acc_pb(t, i+1, j-4+1, k+1)) 
                    + (acc_pb(t, i+4+1, j+1, k+1) + acc_pb(t, i-4+1, j+1, k+1)));
     acc_pb(t+1, i+1, j+1, k+1) = 2 * acc_pb(t, i+1, j+1, k+1) - acc_pb(t+1, i+1, j+1, k+1) + acc_vsq(i+1, j+1, k+1) * div;
            }}}}
  gettimeofday(&end, 0);
  print_summary("Data block + unrolling", tdiff(&end, &start));

#endif /* end PB */

#if PBL
  gettimeofday(&start, 0);
  for (int t = 0; t < T; ++t) {
     cilk_for (int i = half_start_z; i < half_end_z; ++i) {
        for (int j = half_start_y; j < half_end_y; ++j) {
#pragma ivdep
            for (int k = half_start_x; k < half_end_x; ++k) {
          const float c0 = coef[0], c1 = coef[1], c2 = coef[2], c3 = coef[3], c4 = coef[4];
          float div0 = c0 * acc_pbl(t, i, j, k, 0, 0, 0) + 
                c1 * ((acc_pbl(t, i, j, k+1, 0, 0, 0) + acc_pbl(t, i, j, k-1, 0, 0, 0)) 
                    + (acc_pbl(t, i, j+1, k, 0, 0, 0) + acc_pbl(t, i, j-1, k, 0, 0, 0)) 
                    + (acc_pbl(t, i+1, j, k, 0, 0, 0) + acc_pbl(t, i-1, j, k, 0, 0, 0))) 
              + c2 * ((acc_pbl(t, i, j, k+2, 0, 0, 0) + acc_pbl(t, i, j, k-2, 0, 0, 0)) 
                    + (acc_pbl(t, i, j+2, k, 0, 0, 0) + acc_pbl(t, i, j-2, k, 0, 0, 0)) 
                    + (acc_pbl(t, i+2, j, k, 0, 0, 0) + acc_pbl(t, i-2, j, k, 0, 0, 0))) 
              + c3 * ((acc_pbl(t, i, j, k+3, 0, 0, 0) + acc_pbl(t, i, j, k-3, 0, 0, 0)) 
                    + (acc_pbl(t, i, j+3, k, 0, 0, 0) + acc_pbl(t, i, j-3, k, 0, 0, 0)) 
                    + (acc_pbl(t, i+3, j, k, 0, 0, 0) + acc_pbl(t, i-3, j, k, 0, 0, 0))) 
              + c4 * ((acc_pbl(t, i, j, k+4, 0, 0, 0) + acc_pbl(t, i, j, k-4, 0, 0, 0)) 
                    + (acc_pbl(t, i, j+4, k, 0, 0, 0) + acc_pbl(t, i, j-4, k, 0, 0, 0)) 
                    + (acc_pbl(t, i+4, j, k, 0, 0, 0) + acc_pbl(t, i-4, j, k, 0, 0, 0)));
          float div1 = c0 * acc_pbl(t, i, j, k, 0, 0, 1) + 
                c1 * ((acc_pbl(t, i, j, k+1, 0, 0, 1) + acc_pbl(t, i, j, k-1, 0, 0, 1)) 
                    + (acc_pbl(t, i, j+1, k, 0, 0, 1) + acc_pbl(t, i, j-1, k, 0, 0, 1)) 
                    + (acc_pbl(t, i+1, j, k, 0, 0, 1) + acc_pbl(t, i-1, j, k, 0, 0, 1))) 
              + c2 * ((acc_pbl(t, i, j, k+2, 0, 0, 1) + acc_pbl(t, i, j, k-2, 0, 0, 1)) 
                    + (acc_pbl(t, i, j+2, k, 0, 0, 1) + acc_pbl(t, i, j-2, k, 0, 0, 1)) 
                    + (acc_pbl(t, i+2, j, k, 0, 0, 1) + acc_pbl(t, i-2, j, k, 0, 0, 1))) 
              + c3 * ((acc_pbl(t, i, j, k+3, 0, 0, 1) + acc_pbl(t, i, j, k-3, 0, 0, 1)) 
                    + (acc_pbl(t, i, j+3, k, 0, 0, 1) + acc_pbl(t, i, j-3, k, 0, 0, 1)) 
                    + (acc_pbl(t, i+3, j, k, 0, 0, 1) + acc_pbl(t, i-3, j, k, 0, 0, 1))) 
              + c4 * ((acc_pbl(t, i, j, k+4, 0, 0, 1) + acc_pbl(t, i, j, k-4, 0, 0, 1)) 
                    + (acc_pbl(t, i, j+4, k, 0, 0, 1) + acc_pbl(t, i, j-4, k, 0, 0, 1)) 
                    + (acc_pbl(t, i+4, j, k, 0, 0, 1) + acc_pbl(t, i-4, j, k, 0, 0, 1)));
          float div2 = c0 * acc_pbl(t, i, j, k, 0, 1, 0) + 
                c1 * ((acc_pbl(t, i, j, k+1, 0, 1, 0) + acc_pbl(t, i, j, k-1, 0, 1, 0)) 
                    + (acc_pbl(t, i, j+1, k, 0, 1, 0) + acc_pbl(t, i, j-1, k, 0, 1, 0)) 
                    + (acc_pbl(t, i+1, j, k, 0, 1, 0) + acc_pbl(t, i-1, j, k, 0, 1, 0))) 
              + c2 * ((acc_pbl(t, i, j, k+2, 0, 1, 0) + acc_pbl(t, i, j, k-2, 0, 1, 0)) 
                    + (acc_pbl(t, i, j+2, k, 0, 1, 0) + acc_pbl(t, i, j-2, k, 0, 1, 0)) 
                    + (acc_pbl(t, i+2, j, k, 0, 1, 0) + acc_pbl(t, i-2, j, k, 0, 1, 0))) 
              + c3 * ((acc_pbl(t, i, j, k+3, 0, 1, 0) + acc_pbl(t, i, j, k-3, 0, 1, 0)) 
                    + (acc_pbl(t, i, j+3, k, 0, 1, 0) + acc_pbl(t, i, j-3, k, 0, 1, 0)) 
                    + (acc_pbl(t, i+3, j, k, 0, 1, 0) + acc_pbl(t, i-3, j, k, 0, 1, 0))) 
              + c4 * ((acc_pbl(t, i, j, k+4, 0, 1, 0) + acc_pbl(t, i, j, k-4, 0, 1, 0)) 
                    + (acc_pbl(t, i, j+4, k, 0, 1, 0) + acc_pbl(t, i, j-4, k, 0, 1, 0)) 
                    + (acc_pbl(t, i+4, j, k, 0, 1, 0) + acc_pbl(t, i-4, j, k, 0, 1, 0)));
          float div3 = c0 * acc_pbl(t, i, j, k, 0, 1, 1) + 
                c1 * ((acc_pbl(t, i, j, k+1, 0, 1, 1) + acc_pbl(t, i, j, k-1, 0, 1, 1)) 
                    + (acc_pbl(t, i, j+1, k, 0, 1, 1) + acc_pbl(t, i, j-1, k, 0, 1, 1)) 
                    + (acc_pbl(t, i+1, j, k, 0, 1, 1) + acc_pbl(t, i-1, j, k, 0, 1, 1))) 
              + c2 * ((acc_pbl(t, i, j, k+2, 0, 1, 1) + acc_pbl(t, i, j, k-2, 0, 1, 1)) 
                    + (acc_pbl(t, i, j+2, k, 0, 1, 1) + acc_pbl(t, i, j-2, k, 0, 1, 1)) 
                    + (acc_pbl(t, i+2, j, k, 0, 1, 1) + acc_pbl(t, i-2, j, k, 0, 1, 1))) 
              + c3 * ((acc_pbl(t, i, j, k+3, 0, 1, 1) + acc_pbl(t, i, j, k-3, 0, 1, 1)) 
                    + (acc_pbl(t, i, j+3, k, 0, 1, 1) + acc_pbl(t, i, j-3, k, 0, 1, 1)) 
                    + (acc_pbl(t, i+3, j, k, 0, 1, 1) + acc_pbl(t, i-3, j, k, 0, 1, 1))) 
              + c4 * ((acc_pbl(t, i, j, k+4, 0, 1, 1) + acc_pbl(t, i, j, k-4, 0, 1, 1)) 
                    + (acc_pbl(t, i, j+4, k, 0, 1, 1) + acc_pbl(t, i, j-4, k, 0, 1, 1)) 
                    + (acc_pbl(t, i+4, j, k, 0, 1, 1) + acc_pbl(t, i-4, j, k, 0, 1, 1)));
          float div4 = c0 * acc_pbl(t, i, j, k, 1, 0, 0) + 
                c1 * ((acc_pbl(t, i, j, k+1, 1, 0, 0) + acc_pbl(t, i, j, k-1, 1, 0, 0)) 
                    + (acc_pbl(t, i, j+1, k, 1, 0, 0) + acc_pbl(t, i, j-1, k, 1, 0, 0)) 
                    + (acc_pbl(t, i+1, j, k, 1, 0, 0) + acc_pbl(t, i-1, j, k, 1, 0, 0))) 
              + c2 * ((acc_pbl(t, i, j, k+2, 1, 0, 0) + acc_pbl(t, i, j, k-2, 1, 0, 0)) 
                    + (acc_pbl(t, i, j+2, k, 1, 0, 0) + acc_pbl(t, i, j-2, k, 1, 0, 0)) 
                    + (acc_pbl(t, i+2, j, k, 1, 0, 0) + acc_pbl(t, i-2, j, k, 1, 0, 0))) 
              + c3 * ((acc_pbl(t, i, j, k+3, 1, 0, 0) + acc_pbl(t, i, j, k-3, 1, 0, 0)) 
                    + (acc_pbl(t, i, j+3, k, 1, 0, 0) + acc_pbl(t, i, j-3, k, 1, 0, 0)) 
                    + (acc_pbl(t, i+3, j, k, 1, 0, 0) + acc_pbl(t, i-3, j, k, 1, 0, 0))) 
              + c4 * ((acc_pbl(t, i, j, k+4, 1, 0, 0) + acc_pbl(t, i, j, k-4, 1, 0, 0)) 
                    + (acc_pbl(t, i, j+4, k, 1, 0, 0) + acc_pbl(t, i, j-4, k, 1, 0, 0)) 
                    + (acc_pbl(t, i+4, j, k, 1, 0, 0) + acc_pbl(t, i-4, j, k, 1, 0, 0)));
          float div5 = c0 * acc_pbl(t, i, j, k, 1, 0, 1) + 
                c1 * ((acc_pbl(t, i, j, k+1, 1, 0, 1) + acc_pbl(t, i, j, k-1, 1, 0, 1)) 
                    + (acc_pbl(t, i, j+1, k, 1, 0, 1) + acc_pbl(t, i, j-1, k, 1, 0, 1)) 
                    + (acc_pbl(t, i+1, j, k, 1, 0, 1) + acc_pbl(t, i-1, j, k, 1, 0, 1))) 
              + c2 * ((acc_pbl(t, i, j, k+2, 1, 0, 1) + acc_pbl(t, i, j, k-2, 1, 0, 1)) 
                    + (acc_pbl(t, i, j+2, k, 1, 0, 1) + acc_pbl(t, i, j-2, k, 1, 0, 1)) 
                    + (acc_pbl(t, i+2, j, k, 1, 0, 1) + acc_pbl(t, i-2, j, k, 1, 0, 1))) 
              + c3 * ((acc_pbl(t, i, j, k+3, 1, 0, 1) + acc_pbl(t, i, j, k-3, 1, 0, 1)) 
                    + (acc_pbl(t, i, j+3, k, 1, 0, 1) + acc_pbl(t, i, j-3, k, 1, 0, 1)) 
                    + (acc_pbl(t, i+3, j, k, 1, 0, 1) + acc_pbl(t, i-3, j, k, 1, 0, 1))) 
              + c4 * ((acc_pbl(t, i, j, k+4, 1, 0, 1) + acc_pbl(t, i, j, k-4, 1, 0, 1)) 
                    + (acc_pbl(t, i, j+4, k, 1, 0, 1) + acc_pbl(t, i, j-4, k, 1, 0, 1)) 
                    + (acc_pbl(t, i+4, j, k, 1, 0, 1) + acc_pbl(t, i-4, j, k, 1, 0, 1)));
          float div6 = c0 * acc_pbl(t, i, j, k, 1, 1, 0) + 
                c1 * ((acc_pbl(t, i, j, k+1, 1, 1, 0) + acc_pbl(t, i, j, k-1, 1, 1, 0)) 
                    + (acc_pbl(t, i, j+1, k, 1, 1, 0) + acc_pbl(t, i, j-1, k, 1, 1, 0)) 
                    + (acc_pbl(t, i+1, j, k, 1, 1, 0) + acc_pbl(t, i-1, j, k, 1, 1, 0))) 
              + c2 * ((acc_pbl(t, i, j, k+2, 1, 1, 0) + acc_pbl(t, i, j, k-2, 1, 1, 0)) 
                    + (acc_pbl(t, i, j+2, k, 1, 1, 0) + acc_pbl(t, i, j-2, k, 1, 1, 0)) 
                    + (acc_pbl(t, i+2, j, k, 1, 1, 0) + acc_pbl(t, i-2, j, k, 1, 1, 0))) 
              + c3 * ((acc_pbl(t, i, j, k+3, 1, 1, 0) + acc_pbl(t, i, j, k-3, 1, 1, 0)) 
                    + (acc_pbl(t, i, j+3, k, 1, 1, 0) + acc_pbl(t, i, j-3, k, 1, 1, 0)) 
                    + (acc_pbl(t, i+3, j, k, 1, 1, 0) + acc_pbl(t, i-3, j, k, 1, 1, 0))) 
              + c4 * ((acc_pbl(t, i, j, k+4, 1, 1, 0) + acc_pbl(t, i, j, k-4, 1, 1, 0)) 
                    + (acc_pbl(t, i, j+4, k, 1, 1, 0) + acc_pbl(t, i, j-4, k, 1, 1, 0)) 
                    + (acc_pbl(t, i+4, j, k, 1, 1, 0) + acc_pbl(t, i-4, j, k, 1, 1, 0)));
          float div7 = c0 * acc_pbl(t, i, j, k, 1, 1, 1) + 
                c1 * ((acc_pbl(t, i, j, k+1, 1, 1, 1) + acc_pbl(t, i, j, k-1, 1, 1, 1)) 
                    + (acc_pbl(t, i, j+1, k, 1, 1, 1) + acc_pbl(t, i, j-1, k, 1, 1, 1)) 
                    + (acc_pbl(t, i+1, j, k, 1, 1, 1) + acc_pbl(t, i-1, j, k, 1, 1, 1))) 
              + c2 * ((acc_pbl(t, i, j, k+2, 1, 1, 1) + acc_pbl(t, i, j, k-2, 1, 1, 1)) 
                    + (acc_pbl(t, i, j+2, k, 1, 1, 1) + acc_pbl(t, i, j-2, k, 1, 1, 1)) 
                    + (acc_pbl(t, i+2, j, k, 1, 1, 1) + acc_pbl(t, i-2, j, k, 1, 1, 1))) 
              + c3 * ((acc_pbl(t, i, j, k+3, 1, 1, 1) + acc_pbl(t, i, j, k-3, 1, 1, 1)) 
                    + (acc_pbl(t, i, j+3, k, 1, 1, 1) + acc_pbl(t, i, j-3, k, 1, 1, 1)) 
                    + (acc_pbl(t, i+3, j, k, 1, 1, 1) + acc_pbl(t, i-3, j, k, 1, 1, 1))) 
              + c4 * ((acc_pbl(t, i, j, k+4, 1, 1, 1) + acc_pbl(t, i, j, k-4, 1, 1, 1)) 
                    + (acc_pbl(t, i, j+4, k, 1, 1, 1) + acc_pbl(t, i, j-4, k, 1, 1, 1)) 
                    + (acc_pbl(t, i+4, j, k, 1, 1, 1) + acc_pbl(t, i-4, j, k, 1, 1, 1)));
     acc_pbl(t+1, i, j, k, 0, 0, 0) = 2 * acc_pbl(t, i, j, k, 0, 0, 0) - acc_pbl(t+1, i, j, k, 0, 0, 0) + acc_vsq(i<<1, j<<1, k<<1) * div0;
     acc_pbl(t+1, i, j, k, 0, 0, 1) = 2 * acc_pbl(t, i, j, k, 0, 0, 1) - acc_pbl(t+1, i, j, k, 0, 0, 1) + acc_vsq(i<<1, j<<1, (k<<1) + 1) * div1;
     acc_pbl(t+1, i, j, k, 0, 1, 0) = 2 * acc_pbl(t, i, j, k, 0, 1, 0) - acc_pbl(t+1, i, j, k, 0, 1, 0) + acc_vsq(i<<1, (j<<1)+1, k<<1) * div2;
     acc_pbl(t+1, i, j, k, 0, 1, 1) = 2 * acc_pbl(t, i, j, k, 0, 1, 1) - acc_pbl(t+1, i, j, k, 0, 1, 1) + acc_vsq(i<<1, (j<<1)+1, (k<<1)+1) * div3;
     acc_pbl(t+1, i, j, k, 1, 0, 0) = 2 * acc_pbl(t, i, j, k, 1, 0, 0) - acc_pbl(t+1, i, j, k, 1, 0, 0) + acc_vsq((i<<1)+1, j<<1, k<<1) * div4;
     acc_pbl(t+1, i, j, k, 1, 0, 1) = 2 * acc_pbl(t, i, j, k, 1, 0, 1) - acc_pbl(t+1, i, j, k, 1, 0, 1) + acc_vsq((i<<1)+1, j<<1, (k<<1)+1) * div5;
     acc_pbl(t+1, i, j, k, 1, 1, 0) = 2 * acc_pbl(t, i, j, k, 1, 1, 0) - acc_pbl(t+1, i, j, k, 1, 1, 0) + acc_vsq((i<<1)+1, (j<<1)+1, k) * div6;
     acc_pbl(t+1, i, j, k, 1, 1, 1) = 2 * acc_pbl(t, i, j, k, 1, 1, 1) - acc_pbl(t+1, i, j, k, 1, 1, 1) + acc_vsq((i<<1)+1, (j<<1)+1, (k<<1)+1) * div7;
            }}}}
  gettimeofday(&end, 0);
  print_summary("Data block + unrolling", tdiff(&end, &start));

#endif /* end PBL */

#if PC

  gettimeofday(&start, 0);
  for (int t = 0; t < T; ++t) {
     cilk_for (int i = start_z; i < end_z; ++i) {
        for (int j = start_y; j < end_y; ++j) {
#pragma ivdep
            for (int k = start_x; k < end_x; ++k) {
          const float c0 = coef[0], c1 = coef[1], c2 = coef[2], c3 = coef[3], c4 = coef[4];
          float div = c0 * acc_pc(t, i, j, k) + 
                c1 * ((acc_pc(t, i, j, k+1) + acc_pc(t, i, j, k-1)) 
                    + (acc_pc(t, i, j+1, k) + acc_pc(t, i, j-1, k)) 
                    + (acc_pc(t, i+1, j, k) + acc_pc(t, i-1, j, k))) 
              + c2 * ((acc_pc(t, i, j, k+2) + acc_pc(t, i, j, k-2)) 
                    + (acc_pc(t, i, j+2, k) + acc_pc(t, i, j-2, k)) 
                    + (acc_pc(t, i+2, j, k) + acc_pc(t, i-2, j, k))) 
              + c3 * ((acc_pc(t, i, j, k+3) + acc_pc(t, i, j, k-3)) 
                    + (acc_pc(t, i, j+3, k) + acc_pc(t, i, j-3, k)) 
                    + (acc_pc(t, i+3, j, k) + acc_pc(t, i-3, j, k))) 
              + c4 * ((acc_pc(t, i, j, k+4) + acc_pc(t, i, j, k-4)) 
                    + (acc_pc(t, i, j+4, k) + acc_pc(t, i, j-4, k)) 
                    + (acc_pc(t, i+4, j, k) + acc_pc(t, i-4, j, k)));
     acc_pc(t+1, i, j, k) = 2 * acc_pc(t, i, j, k) - acc_pc(t+1, i, j, k) + acc_vsq(i, j, k) * div;
            }}}}
  gettimeofday(&end, 0);
  print_summary("Naive", tdiff(&end, &start));
#endif /* end PC */

#if (PB && PC)
    int t = T;
    for (int i = start_z; i < end_z; ++i) {
    for (int j = start_y; j < end_y; ++j) {
    for (int k = start_x; k < end_x; ++k) {
        check_result(t, i, j, k, acc_pb(t, i, j, k), acc_pc(t, i, j, k));
    } } }
#endif


#if (PBL && PC)
    int t = T;
    for (int i = start_z; i < end_z; ++i) {
    for (int j = start_y; j < end_y; ++j) {
    for (int k = start_x; k < end_x; ++k) {
        check_result(t, i, j, k, acc_pbl(t, i>>1, j>>1, k>>1, i&0x1, j&0x1, k&0x1), acc_pc(t, i, j, k));
    } } }
    printf("check passed!\n");
#endif

  delete [] A[0];
  delete [] A[1];
  delete [] A;
  delete [] vsq;
  delete [] pc;
  delete [] pb;
  return 0;
}
