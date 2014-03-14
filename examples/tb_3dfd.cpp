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
#include <sys/time.h>

#include <pochoir.hpp>

using namespace std;

int nthreads = 1;
const int ds = 4;
int Nx = 100;
int Ny = 100;
int Nz = 100;
int T = 40;
static const int NPIECES = 2;
int N_CORES=1;
#if 1
static const int dt_threshold = 3;
static const int dx_threshold = 1000;
static const int dyz_threshold = 3;
#else
static const int dt_threshold = 5;
static const int dx_threshold = 150;
static const int dyz_threshold = 150;
#endif
float **A;

float coef[ds + 1];
float *vsq;

int N = 997;
int Nxy;
int sx2, sx3, sx4;
int sxy2, sxy3, sxy4;

void basecase(int t0, int t1, 
	      int x0, int dx0, int x1, int dx1,
	      int y0, int dy0, int y1, int dy1, 
	      int z0, int dz0, int z1, int dz1 )
{
  int _Nx = Nx;
  int Nxy = _Nx * Ny;
  int sx2 = _Nx * 2;
  int sx3 = _Nx * 3;
  int sx4 = _Nx * 4;
  int sxy2 = Nxy * 2;
  int sxy3 = Nxy * 3;
  int sxy4 = Nxy * 4;
  float c0 = coef[0], c1 = coef[1], c2 = coef[2], c3 = coef[3], c4 = coef[4];

  for(int t = t0; t < t1; ++t) {
    for(int z = z0; z < z1; ++z) {
      for(int y = y0; y < y1; ++y) {
	  float *A_cur = &A[t & 1][z * Nxy + y * _Nx];
	  float *A_next = &A[(t + 1) & 1][z * Nxy + y * _Nx];
	  float *vvv = &vsq[z * Nxy + y * _Nx];
#pragma ivdep
	for(int x = x0; x < x1; ++x) {
	  float div = c0 * A_cur[x] 
	    + c1 * ((A_cur[x + 1] + A_cur[x - 1])
		    + (A_cur[x + _Nx] + A_cur[x - _Nx])
		    + (A_cur[x + Nxy] + A_cur[x - Nxy]))
	    + c2 * ((A_cur[x + 2] + A_cur[x - 2])
		    + (A_cur[x + sx2] + A_cur[x - sx2])
		    + (A_cur[x + sxy2] + A_cur[x - sxy2]))
	    + c3 * ((A_cur[x + 3] + A_cur[x - 3])
		    + (A_cur[x + sx3] + A_cur[x - sx3])
		    + (A_cur[x + sxy3] + A_cur[x - sxy3]))
	    + c4 * ((A_cur[x + 4] + A_cur[x - 4])
		    + (A_cur[x + sx4] + A_cur[x - sx4])
		    + (A_cur[x + sxy4] + A_cur[x - sxy4]));
	  A_next[x] = 2 * A_cur[x] - A_next[x] + vvv[x] * div;
	}
      }
    }
    x0 += dx0; x1 += dx1;
    y0 += dy0; y1 += dy1;
    z0 += dz0; z1 += dz1;
  }
}

void basecase_loop(int t, 
                   int x0, int x1,
                   int y0, int y1, 
                   int z )
{
  int _Nx = Nx;
  int Nxy = _Nx * Ny;
  int sx2 = _Nx * 2;
  int sx3 = _Nx * 3;
  int sx4 = _Nx * 4;
  int sxy2 = Nxy * 2;
  int sxy3 = Nxy * 3;
  int sxy4 = Nxy * 4;
  float c0 = coef[0], c1 = coef[1], c2 = coef[2], c3 = coef[3], c4 = coef[4];

      for(int y = y0; y < y1; ++y) {
	  float *A_cur = &A[t & 1][z * Nxy + y * _Nx];
	  float *A_next = &A[(t + 1) & 1][z * Nxy + y * _Nx];
	  float *vvv = &vsq[z * Nxy + y * _Nx];
#pragma ivdep
	for(int x = x0; x < x1; ++x) {
	  float div = c0 * A_cur[x] 
	    + c1 * ((A_cur[x + 1] + A_cur[x - 1])
		    + (A_cur[x + _Nx] + A_cur[x - _Nx])
		    + (A_cur[x + Nxy] + A_cur[x - Nxy]))
	    + c2 * ((A_cur[x + 2] + A_cur[x - 2])
		    + (A_cur[x + sx2] + A_cur[x - sx2])
		    + (A_cur[x + sxy2] + A_cur[x - sxy2]))
	    + c3 * ((A_cur[x + 3] + A_cur[x - 3])
		    + (A_cur[x + sx3] + A_cur[x - sx3])
		    + (A_cur[x + sxy3] + A_cur[x - sxy3]))
	    + c4 * ((A_cur[x + 4] + A_cur[x - 4])
		    + (A_cur[x + sx4] + A_cur[x - sx4])
		    + (A_cur[x + sxy4] + A_cur[x - sxy4]));
	  A_next[x] = 2 * A_cur[x] - A_next[x] + vvv[x] * div;
	}
      }
}

/* map the triple (t, x, y) into a unique long long */
static inline long long encode(int t, int x, int y, int z)
{
  return N * (N * (N * (long long)t + x) + y) + z;
}

static inline float &aref(int t, int x, int y, int z)
{
  return A[t & 1][Nxy * z + Nx * y + x];
}

static inline float &aref(int t, int s) {
  return A[t & 1][s];
}

static inline float &vsqref(int x, int y, int z)
{
  return vsq[Nxy * z + Nx * y + x];
}

static inline float &vsqref(int s)
{
  return vsq[s];
}

//Kernel:
//	Addition: 26
//  Multiplication: 7
void loop_opt3(int t0, int t1, 
	       int x0, int x1,
	       int y0, int y1,
	       int z0, int z1)
{
  for(int t = t0; t < t1; ++t) {
      cilk_for (int z = z0; z < z1; ++ z) {
          basecase_loop(t, x0, x1, y0, y1, z);
      }
  }
}

void walk3(int t0, int t1, 
	   int x0, int dx0, int x1, int dx1,
	   int y0, int dy0, int y1, int dy1, 
	   int z0, int dz0, int z1, int dz1 )
{
  int dt = t1 - t0, dx = x1 - x0, dy = y1 - y0, dz = z1 - z0;
  int i;

  if (dx >= dx_threshold && dx >= dy && dx >= dz &&
      dt >= 1 && dx >= 2 * ds * dt * NPIECES) {
    int chunk = dx / NPIECES;

    for (i = 0; i < NPIECES - 1; ++i)
      cilk_spawn walk3(t0, t1,
		       x0 + i * chunk, ds, x0 + (i+1) * chunk, -ds,
		       y0, dy0, y1, dy1,
		       z0, dz0, z1, dz1);
    cilk_spawn walk3(t0, t1,
		     x0 + i * chunk, ds, x1, -ds,
		     y0, dy0, y1, dy1, 
		     z0, dz0, z1, dz1);
    cilk_sync;
    cilk_spawn walk3(t0, t1, 
		     x0, dx0, x0, ds,
		     y0, dy0, y1, dy1, 
		     z0, dz0, z1, dz1);
    for (i = 1; i < NPIECES; ++i)
      cilk_spawn walk3(t0, t1,
		       x0 + i * chunk, -ds, x0 + i * chunk, ds,
		       y0, dy0, y1, dy1, 
		       z0, dz0, z1, dz1);
    cilk_spawn walk3(t0, t1, 
		     x1, -ds, x1, dx1,
		     y0, dy0, y1, dy1, 
		     z0, dz0, z1, dz1);
  } else if (dy >= dyz_threshold && dy >= dz && dt >= 1 && dy >= 2 * ds * dt * NPIECES) {
    int chunk = dy / NPIECES;

    for (i = 0; i < NPIECES - 1; ++i)
      cilk_spawn walk3(t0, t1,
		       x0, dx0, x1, dx1,
		       y0 + i * chunk, ds, y0 + (i+1) * chunk, -ds, 
		       z0, dz0, z1, dz1);
    cilk_spawn walk3(t0, t1,
		     x0, dx0, x1, dx1,
		     y0 + i * chunk, ds, y1, -ds, 
		     z0, dz0, z1, dz1);
    cilk_sync;
    cilk_spawn walk3(t0, t1, 
		     x0, dx0, x1, dx1,
		     y0, dy0, y0, ds, 
		     z0, dz0, z1, dz1);
    for (i = 1; i < NPIECES; ++i)
      cilk_spawn walk3(t0, t1,
		       x0, dx0, x1, dx1,
		       y0 + i * chunk, -ds, y0 + i * chunk, ds, 
		       z0, dz0, z1, dz1);
    cilk_spawn walk3(t0, t1, 
		     x0, dx0, x1, dx1,
		     y1, -ds, y1, dy1, 
		     z0, dz0, z1, dz1);
  } else if (dz >= dyz_threshold && dt >= 1 && dz >= 2 * ds * dt * NPIECES) {
    int chunk = dz / NPIECES;

    for (i = 0; i < NPIECES - 1; ++i)
      cilk_spawn walk3(t0, t1,
		       x0, dx0, x1, dx1,
		       y0, dy0, y1, dy1,
		       z0 + i * chunk, ds, z0 + (i+1) * chunk, -ds);
    cilk_spawn walk3(t0, t1,
		     x0, dx0, x1, dx1,
		     y0, dy0, y1, dy1, 
		     z0 + i * chunk, ds, z1, -ds);
    cilk_sync;
    cilk_spawn walk3(t0, t1, 
		     x0, dx0, x1, dx1,
		     y0, dy0, y1, dy1,
		     z0, dz0, z0, ds);
    for (i = 1; i < NPIECES; ++i)
      cilk_spawn walk3(t0, t1,
		       x0, dx0, x1, dx1,
		       y0, dy0, y1, dy1,
		       z0 + i * chunk, -ds, z0 + i * chunk, ds);
    cilk_spawn walk3(t0, t1, 
		     x0, dx0, x1, dx1,
		     y0, dy0, y1, dy1,
		     z1, -ds, z1, dz1);
  }  else if (dt > dt_threshold) {
    int halfdt = dt / 2;
    walk3(t0, t0 + halfdt,
	  x0, dx0, x1, dx1,
	  y0, dy0, y1, dy1, 
	  z0, dz0, z1, dz1);
    walk3(t0 + halfdt, t1, 
	  x0 + dx0 * halfdt, dx0, x1 + dx1 * halfdt, dx1,
	  y0 + dy0 * halfdt, dy0, y1 + dy1 * halfdt, dy1, 
	  z0 + dz0 * halfdt, dz0, z1 + dz1 * halfdt, dz1);
  } else {
    basecase(t0, t1, 
	     x0, dx0, x1, dx1,
	     y0, dy0, y1, dy1,
	     z0, dz0, z1, dz1);
  } 
}

void init_variables() 
{
  int count = 0;
  Nxy = Nx * Ny;
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

  count = 0;

  for (int z = 0; z < Nz; ++z)
    for (int y = 0; y < Ny; ++y) 
      for(int x = 0; x < Nx; ++x) {
	/* set initial values */
	/*
	  aref(0, x, y, z) = encode(0, x, y, z);
	  aref(1, x, y, z) = encode(-1, x, y, z); // set to invalid
	*/
	float r = abs((float)(x - Nx/2 + y - Ny/2 + z - Nz/2) / 30);
	r = max(1 - r, 0.0f) + 1;
	
	aref(0, x, y, z) = r;
	aref(1, x, y, z) = r;
	vsqref(x, y, z) = 0.001f;
      }
    //N_CORES = max(2, __cilkrts_get_nworkers());
    N_CORES = max(1, __cilkrts_get_nworkers());
    printf("N_CORES = %d\n", N_CORES);
}

template <typename T_Array>
void init_pochoir_array(T_Array & arr) 
{
  int count = 0;
  Nxy = Nx * Ny;
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

  count = 0;

  for (int z = 0; z < Nz; ++z)
    for (int y = 0; y < Ny; ++y) 
      for(int x = 0; x < Nx; ++x) {
	/* set initial values */
	/*
	  aref(0, x, y, z) = encode(0, x, y, z);
	  aref(1, x, y, z) = encode(-1, x, y, z); // set to invalid
	*/
	float r = abs((float)(x - Nx/2 + y - Ny/2 + z - Nz/2) / 30);
	r = max(1 - r, 0.0f) + 1;
	
	arr(0, z, y, x) = r;
	arr(1, z, y, x) = r;
	vsqref(x, y, z) = 0.001f;
  }
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
  //printf("time: %f\n", interval);
  printf("time: %f ms \n", 1.0e3 * interval);
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

void print_y() {
  FILE *fout = fopen("y_points.txt", "w");
  int z = Nz/2;
  int x = Nx/2;
  for(int y = 0; y < Ny; y++) {
    fprintf(fout, "%f\n", aref(T, x, y, z));
  }
  fclose(fout);
  printf("Done writing output\n");
}

void dotest()
{
  //initialization
  A = new float*[2];
  A[0] = new float[Nx * Ny * Nz];
  A[1] = new float[Nx * Ny * Nz];
  vsq = new float[Nx * Ny * Nz];

  struct timeval start, end;
	
  ///////////////////////////////////////////////                                                                      
#if 0
  
  init_variables();
  gettimeofday(&start, 0);
  /* this is loop based version */
  loop_opt3(0, T,
            ds, Nx - ds, 
            ds, Ny - ds,
            ds, Nz - ds);
  gettimeofday(&end, 0);
  //basecase(0, T,
  //	    ds, 0, Nx - ds, 0, 
  //	    ds, 0, Ny - ds, 0, 
  //	    ds, 0, Nz - ds, 0);
  //copy_A_to_B();
  print_summary("base", tdiff(&end, &start));
  ///////////////////////////////////////////////
#endif
  
  init_variables();
  // verify_A_and_B();
  gettimeofday(&start, 0);
  /* this is the divide-and-conquer version in cilk++ */
  walk3(0, T,
	    ds, 0, Nx - ds, 0, 
		ds, 0, Ny - ds, 0, 
		ds, 0, Nz - ds, 0);
  gettimeofday(&end, 0);
  print_summary("COStencilTask", tdiff(&end, &start));

  // verify_A_and_B();
  //print_y();

}

Pochoir_Boundary_3D(fd_bv_3D, arr, t, i, j, k)
    return 0;
Pochoir_Boundary_End

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

  Pochoir_Shape_3D fd_shape_3D[26] = {{1, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 0, -1}, {0, 0, 1, 0}, {0, 0, -1, 0}, {0, 1, 0, 0}, {0, -1, 0, 0}, {0, 0, 0, 2}, {0, 0, 0, -2}, {0, 0, 2, 0}, {0, 0, -2, 0}, {0, 2, 0, 0}, {0, -2, 0, 0}, {0, 0, 0, 3}, {0, 0, 0, -3}, {0, 0, 3, 0}, {0, 0, -3, 0}, {0, 3, 0, 0}, {0, -3, 0, 0}, {0, 0, 0, 4}, {0, 0, 0, -4}, {0, 0, 4, 0}, {0, 0, -4, 0}, {0, 4, 0, 0}, {0, -4, 0, 0}};
  Pochoir_Array_3D(float) pa(Nz, Ny, Nx);
  Pochoir_3D fd_3D(fd_shape_3D);
  Pochoir_Domain I(0+ds, Nx-ds), J(0+ds, Ny-ds), K(0+ds, Nz-ds);

  pa.Register_Boundary(fd_bv_3D) ;
  fd_3D.Register_Array(pa);
  //fd_3D.Register_Domain(K, J, I);

  Pochoir_Kernel_3D(fd_3D_fn, t, i, j, k)
    float c0 = coef[0], c1 = coef[1], c2 = coef[2], c3 = coef[3], c4 = coef[4];
    float div = c0 * pa(t, i, j, k) + 
                c1 * ((pa(t, i, j, k+1) + pa(t, i, j, k-1)) 
                    + (pa(t, i, j+1, k) + pa(t, i, j-1, k)) 
                    + (pa(t, i+1, j, k) + pa(t, i-1, j, k))) 
              + c2 * ((pa(t, i, j, k+2) + pa(t, i, j, k-2)) 
                    + (pa(t, i, j+2, k) + pa(t, i, j-2, k)) 
                    + (pa(t, i+2, j, k) + pa(t, i-2, j, k))) 
              + c3 * ((pa(t, i, j, k+3) + pa(t, i, j, k-3)) 
                    + (pa(t, i, j+3, k) + pa(t, i, j-3, k)) 
                    + (pa(t, i+3, j, k) + pa(t, i-3, j, k))) 
              + c4 * ((pa(t, i, j, k+4) + pa(t, i, j, k-4)) 
                    + (pa(t, i, j+4, k) + pa(t, i, j-4, k)) 
                    + (pa(t, i+4, j, k) + pa(t, i-4, j, k)));
     pa(t+1, i, j, k) = 2 * pa(t, i, j, k) - pa(t+1, i, j, k) + vsq[i * Nxy + j * Nx + k] * div;
  Pochoir_Kernel_End

  dotest();

  char name [500] ;
  sprintf(name, "3dfd") ;
  fd_3D.set_problem_name(name) ;

  init_pochoir_array(pa);
  gettimeofday(&start, 0);
  fd_3D.Run(T, fd_3D_fn);
  gettimeofday(&end, 0);
  print_summary("Pochoir", tdiff(&end, &start));

  delete[] A;
  delete[] vsq;
  return 0;
}
