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
/* Simple implementation of a 2d scalar wave equation:

        \partial^2 u / \partial t^2 = b \nabla \cdot (a \nabla u)

   where a and b are material properties (usually varying in space,
   sometimes in time), with the phase velocity c given by c^2 = ab at each
   point in space.

   We will discretize this in an sx-by-sy domain, surrounded by
   perfectly matched layers (PML) with thickness dpml, and some given
   resolution.  For illustration purposes, we will choose b = 1 and
   a = quarter-circle radius-R "bend" of a waveguide (width w) formed
   by a = 0.1 surrounded by a = 1.0. We use a line-source term at one side 
   of the bend creating a gaussian pulse, periodic output of the non-PML
   region, and in the future a Fourier-transformed flux calculation at the
   other end of the bend (typical of a transmission calculation).

   To discretize the wave equation, we first write it in a first-order
   form:
          \partial u / \partial t = b \nabla \cdot \vec{v}
	  \partial \vec{v} / \partial t = a \nabla u
   where \vec{v} is a 2d vector field.  This is then discretized
   in a leap-frog fashion in time with a staggered (Yee) grid in space.

   That is, let u, vx, and vy by Nx by Ny arrays.  Then the mapping
   between the [ix,iy] entries of the arrays and the (x,y) points in
   space are:
                u[ix,iy] = u(ix*dx, iy*dy)
                vx[ix,iy] = vx((ix+0.5)*dx, iy*dy)
                vy[ix,iy] = vy(ix*dx, (iy+0.5)*dy)

   Our units of distance are chosen such that the "vacuum" (a = b = 1, c = 1)
   central wavelength of our Gaussian pulse is 1, with a frequency = 1.

   We will use zero initial conditions, but include a source term f(x,y,t)
   in the u equation (only nonzero in some small region), which becomes:
          \partial u / \partial t = b \nabla \cdot \vec{v} + f(x,y,t)

   We will use PML absorbing boundary layers as mentioned above, scaled
   so that the reflectivity R0 at infinite resolution is 10^{-10}, and
   a quadratic PML profile.  Given such absorbing layers, the actual boundary
   conditions are irrelevant (the fields are nearly zero at the boundaries),
   but for illustration purposes we will use periodic boundary conditions.
   In the PML regions adjacent to each boundary, auxiliary \psi arrays
   are needed, as described in section 4.2 of my notes:
               http://math.mit.edu/~stevenj/18.369/pml.pdf
   We only allocate these arrays for the boundary regions, so we need
   four of them (upper and lower boundaries in the x and y directions).

   All arrays are stored as X-by-Y in row-major order (Y contiguous).

   Initial implementation got from Steven G. Johnson
   ported to Pochoir by Yuan Tang
*/

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <cstddef>
#include <iostream>
#include <pochoir.hpp>

using namespace std;
#define TIMES 1
#define TOLERANCE (1e-6)

typedef struct {
    double u;
    double vx, vy;
} wave_2D_unit;

void check_result(int t, int i, int j, double a, double b)
{
    if (abs(a - b) < TOLERANCE) {
//        printf("a(%d, %d, %d) == b(%d, %d, %d) == %f : passed!\n", t, i, j, t, i, j, a);
    } else {
        printf("a(%d, %d, %d) = %f, b(%d, %d, %d) = %f : FAILED!\n", t, i, j, a, t, i, j, b);
    }
}

Pochoir_Boundary_2D(periodic_2D, arr, t, i, j)
    const int arr_size_0 = arr.size(0);
    const int arr_size_1 = arr.size(1);

    int new_i = (i >= arr_size_1) ? (i - arr_size_1) : (i < 0 ? i + arr_size_1 : i);
    int new_j = (j >= arr_size_0) ? (j - arr_size_0) : (j < 0 ? j + arr_size_0 : j);

    return arr.get(t, new_i, new_j);
Pochoir_Boundary_End

Pochoir_Boundary_2D(aperiodic_2D, arr, t, i, j)
    return 0;
Pochoir_Boundary_End

void zero(double *A, int N) // set array to zero
{
  for (int i = 0; i < N; ++i) A[i] = 0.0;
}

/* Output the Nx-by-Ny array A to a filename "fname-ttttt.csv" where
   ttttt is the timestep t, as comma-separated values.  (In a more
   realistic application we would probably use a standard binary
   format like HDF5.)  A(i,j) is stored at A[i*Sy + j], i.e. Sx is the
   stride in x.  Ny must be > 0. */
void output(const char *fname, int t, double *A, int Nx, int Ny, int Sx) {
    char *fn = new char[strlen(fname) + 1 + 5 + 4 + 1];
    sprintf(fn, "%s-%05d.csv", fname, t);
    FILE *f = fopen(fn, "w");
    if (!f) abort();
    for (int ix = 0; ix < Nx; ++ix) {
        for (int iy = 0; iy < Ny-1; ++iy)
            fprintf(f, "%g,", A[ix*Sx + iy]);
        fprintf(f, "%g\n", A[ix*Sx + Ny-1]);
    }
    fclose(f);
}

void output_p_u(const char * fname, int t, int P, Pochoir_Array_2D(wave_2D_unit) & p_uv, int Nx, int Ny) {
    char *fn = new char[strlen(fname) + 1 + 5 + 4 + 1];
    sprintf(fn, "%s-%05d.csv", fname, t);
    FILE *f = fopen(fn, "w");
    if (!f) abort();
    for (int i = P; i < Nx-1-P; ++i) {
        int j;
        for (j = P; j < Ny-1-P-1; ++j)
            fprintf(f, "%g,", p_uv(t, i, j).u);
        fprintf(f, "%g\n", p_uv(t, i, j).u);
    }
    fclose(f);
}

// some geometric parameters:
// We now get the domain size 'sx', 'sy' directly from user's command line argument,
// also it looks buggy to declare 'sx', 'sy' to be of type double
// const double sx = 5, sy = 5;
const double dpml = 0.5; // domain size
const double R = 3.0; // inner bend radius
const double w = 0.2; // waveguide width
const double aw = 0.1; // value of "a" inside waveguide
const double X0 = dpml + 1, Y0 = dpml + 1; // bend center

// Wave-equation coefficient function a(x,y)
double a(double x, double y) {
    x -= X0; y -= Y0;
    if (x > 0 && y > 0) {
        double r = sqrt(x*x + y*y);
        return r > R && r < R + w ? aw : 1.0;
    }
    else if (x <= 0)
        return y > R && y < R + w ? aw : 1.0;
    // else (y <= 0)
    return x > R && x < R + w ? aw : 1.0;
}

/* It's actual a stagger-stencil,
 * if (t-1 % 2 == 0) compute 'u'
 * else if (t-1 % 2 == 1) compute 'vx, vy'
 * But for each time step, it comprises of multiple different kernels in 
 * different sub-regions.
 */
int main(int argc, char * argv[]) {

    struct timeval start, end;
    double min_tdiff = INF;
    if (argc < 3) {
        printf("argc < 3, quit!\n");
        exit(1);
    }
    int sxy = StrToInt(argv[1]);
    int sx = sxy, sy = sxy;
    int T = StrToInt(argv[2]);
    printf("sxy = %d, T = %d\n", sxy, T);

    double resolution = 40; // resolution (pixels/distance)

    double dx = 1.0/resolution; // spatial grid spacing
    double dt = 0.7 * dx; // time-step size, from CFL condition for c = 1
    double dtdx = dt/dx; // ratio used for timestepping

    /* Number of pixels in x and y directions.  The +1 is
       for a halo pixel to implement the periodic boundaries. */
    int Nx = int((sx + 2 * dpml) * resolution) + 1;
    int Ny = int((sy + 2 * dpml) * resolution) + 1;

    int P = int(dpml * resolution); // PML thickness, in pixels

    // Allocate u and v = (vx,vy) arrays, initialized to 0:
    double *u = new double[Nx * Ny]; zero(u, Nx*Ny);
    double *vx = new double[Nx * Ny]; zero(u, Nx*Ny);
    double *vy = new double[Nx * Ny]; zero(u, Nx*Ny);
    /* this is a forged shape for temporary use only */
    Pochoir_Shape_2D shape_wave_2D_u[] = {{0, 0, 0}, {-1, 0, 0}, {-1, -1, 0}, {-1, 0, -1}, {-1, 1, 0}, {-1, 0, 1}};
    Pochoir_Shape_2D shape_wave_2D_v[] = {{1, 0, 0}, {0, 0, 0}, {0, -1, 0}, {0, 0, -1}, {0, 1, 0}, {0, 0, 1}};
    /* In Pochoir, we employ the boundary function to implement
     * the periodic boundary condition. So we should exclude that one point
     */
    Pochoir_Array_2D(wave_2D_unit) p_uv(Nx-1, Ny-1); 
    Pochoir_2D wave_2D;
    /* p_u is of periodic boundary condition */
    p_uv.Register_Boundary(periodic_2D);

    // Allocate auxiliary PML arrays for upper (U) and lower (L) boundary regions
    double *psiUx = new double[P * Ny]; zero(psiUx, P * Ny);
    double *psiLx = new double[P * Ny]; zero(psiLx, P * Ny);
    double *psiUy = new double[Nx * P]; zero(psiUy, Nx * P);
    double *psiLy = new double[Nx * P]; zero(psiLy, Nx * P);

    double *p_psiUx = new double[P * (Ny-1)]; zero(psiUx, P * (Ny-1));
    double *p_psiLx = new double[P * (Ny-1)]; zero(psiLx, P * (Ny-1));
    double *p_psiUy = new double[(Nx-1) * P]; zero(psiUy, (Nx-1) * P);
    double *p_psiLy = new double[(Nx-1) * P]; zero(psiLy, (Nx-1) * P);

    /* Set up the PML conductivity arrays sigma.  These only vary in the
       direction perpendicular to the PML so they can be 1d arrays, and
       we use the same conductivity profiles in the x and y direction and 
       for the upper and lower boundaries for simplicity.  However,
       because u and v are discretized on different (staggered) grids,
       we need to store sigma in an array of length 2P corresponding
       to twice the spatial resolution of the individual grids in order
       for it to have the correct values on both the u and v grids. */
    double R0 = 1e-10; // normal incidence reflectivity at infinite resolution
    // conductivity in PML is sigma(x) ~ sigma_coef * x^2
    double sigma_coef = -0.25*dx*dx * log(R0) / (4.0 * dpml*dpml*dpml / 3.0);
    double *sigma = new double[2*P];
    for (int i = 0; i < 2*P; ++i) sigma[i] = sigma_coef * (i+0.5) * (i+0.5);
    // It is also convenient to have 1/(1+sigma*dt/2) and (1-sigma*dt/2)
    // arrays for the timestepping equations below:
    double *sigmadt = new double[2*P], *sigmadtinv = new double[2*P];
    for (int i = 0; i < 2*P; ++i) {
        sigmadt[i] = 1 - sigma[i]*dt/2;
        sigmadtinv[i] = 1 / (1 + sigma[i]*dt/2);
    }
    // 2P x 2P arrays of 1/(1+sigma*dt/2) and (1-sigma*dt/2), for the corners
    double *sigmadt2 = new double[4*P*P], *sigmadtinv2 = new double[4*P*P];
    for (int ix = 0; ix < 2*P; ++ix)
        for (int iy = 0; iy < 2*P; ++iy) {
            int i = ix * (2*P) + iy;
            double s = sigma[ix] + sigma[iy];
            sigmadt2[i] = 1 - s*dt/2;
            sigmadtinv2[i] = 1 / (1 + s*dt/2);
    }

    /* Set up the coefficient array a.  (We won't allocate b since it = 1.)
       Since a is used in the vx and vy update equations, we need two arrays
       ax and ay, stored on the vx and vy grids. */
    double *ax = new double[Nx * Ny], *ay = new double[Nx * Ny];
    // First, set ax = a on u grid, for output:
    for (int ix = 0; ix < Nx; ++ix)
        for (int iy = 0; iy < Ny; ++iy) {
            int i = ix * Ny + iy;
            double x = ix * dx, y = iy * dx;
            ax[i] = a(x, y);
        }
    output("a", 0, ax + P*Ny + P, Nx - 2*P, Ny - 2*P, Ny); // output interior
    // Now, compute actual ax and ay arrays:
    for (int ix = 0; ix < Nx; ++ix)
        for (int iy = 0; iy < Ny; ++iy) {
            int i = ix * Ny + iy;
            double x = ix * dx, y = iy * dx;
            ax[i] = a(x + 0.5*dx, y);
            ay[i] = a(x, y + 0.5*dx);
        }
    output("ax", 0, ax, Nx, Ny, Ny);
    output("ay", 0, ay, Nx, Ny, Ny);

    /* Initialization of p_ax, p_ay */
    double *p_ax = new double[(Nx-1) * (Ny - 1)], *p_ay = new double[(Nx-1) * (Ny-1)];
    // First, set ax = a on u grid, for output:
    for (int ix = 0; ix < Nx-1; ++ix)
        for (int iy = 0; iy < Ny-1; ++iy) {
            int i = ix * (Ny-1) + iy;
            double x = ix * dx, y = iy * dx;
            /* a(x, y) is a function call to fill in the array of 'p_ax' */
            p_ax[i] = a(x, y);
        }
    output("p_a", 0, p_ax + P*(Ny-1) + P, (Nx-1) - 2*P, (Ny-1) - 2*P, (Ny-1)); // output interior
    // Now, compute actual p_ax and p_ay arrays:
    for (int ix = 0; ix < Nx-1; ++ix)
        for (int iy = 0; iy < Ny-1; ++iy) {
            int i = ix * (Ny-1) + iy;
            double x = ix * dx, y = iy * dx;
            p_ax[i] = a(x + 0.5*dx, y);
            p_ay[i] = a(x, y + 0.5*dx);
        }
    output("p_ax", 0, p_ax, Nx-1, Ny-1, Ny-1);
    output("p_ay", 0, p_ay, Nx-1, Ny-1, Ny-1);

    // Gaussian-pulse characteristics of source
    double fcen = 1.0, fwidth = 0.1; // center frequency and width
    const double twopi = atan(1.0)*8.0;
    // source(t) = cos(omega*t) * exp[- decay*(t-t0)^2]:
    double omega = twopi * fcen, decay = (twopi*twopi*fwidth*fwidth) * 0.5;
    double t0 = 10/fwidth; // start time of gaussian

    /* guard and kernel function for Pochoir starts */

    /* guard and kernel function for interior (non-PML) region */
    Pochoir_Guard_2D_Begin(g_interior, t, i, j)
        /* -1 because we don't have halo point */
        if (i >= P && i < Nx-1-P &&
                j >= P && j < Ny-1-P)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_interior)

    Pochoir_Kernel_2D_Begin(k_interior, t, i, j)
        p_uv(t, i, j).u = p_uv(t-1, i, j).u + dtdx * (p_uv(t, i, j).vx - p_uv(t, i-1, j).vx + p_uv(t, i, j).vy - p_uv(t, i, j-1).vy);
    Pochoir_Kernel_2D_End(k_interior, shape_wave_2D_u)

    /* guard and kernel functions for 8 PML regions */
    Pochoir_Guard_2D_Begin(g_lower_x, t, i, j)
        if (i >= 0 && i <= P-1 &&
                j >= P && j < Ny-1-P)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_lower_x)

    Pochoir_Kernel_2D_Begin(k_lower_x, t, i, j)
        int iSx = 2 * (P-i), iPx = (i-1) * (Ny-1) + j;
        double psi_old = p_psiLx[iPx];
        double dvy = p_uv(t, i, j).vy - p_uv(t, i, j-1).vy;
        p_psiLx[iPx] += dtdx * sigma[iSx] * dvy;
        p_uv(t, i, j).u = sigmadtinv[iSx] * (p_uv(t-1, i, j).u * sigmadt[iSx]
                            + dtdx * (p_uv(t, i, j).vx - p_uv(t, i-1, j).vx + dvy)
                            + (0.5 * dt) * (psi_old + p_psiLx[iPx]));
    Pochoir_Kernel_2D_End(k_lower_x, shape_wave_2D_u)

    Pochoir_Guard_2D_Begin(g_upper_x, t, i, j)
        if (i >= Nx-1-P && i < Nx-1 &&
                j >= P && j < Ny-1-P)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_upper_x)

    Pochoir_Kernel_2D_Begin(k_upper_x, t, i, j)
        int iSx = 2*(i-(Nx-1-P)), iPx = (i-(Nx-1-P))*(Ny-1) + j;
        double psi_old = p_psiUx[iPx];
        double dvy = p_uv(t, i, j).vy - p_uv(t, i, j-1).vy;
        p_psiUx[iPx] += dtdx * sigma[iSx] * dvy;
        p_uv(t, i, j).u = sigmadtinv[iSx] * (p_uv(t-1, i, j).u * sigmadt[iSx]
                      + dtdx * (p_uv(t, i, j).vx - p_uv(t, i-1, j).vx + dvy)
                      + (0.5 * dt) * (psi_old + p_psiUx[iPx]));
    Pochoir_Kernel_2D_End(k_upper_x, shape_wave_2D_u)

    Pochoir_Guard_2D_Begin(g_lower_y, t, i, j)
        if (i >= P && i < Nx-1-P &&
                j >= 0 && j <= P-1)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_lower_y)

    Pochoir_Kernel_2D_Begin(k_lower_y, t, i, j)
        int iSy = 2*(P-j), iPy = i*P + (j-1);
        double psi_old = p_psiLy[iPy];
        double dvx = p_uv(t, i, j).vx - p_uv(t, i-1, j).vx;
        p_psiLy[iPy] += dtdx * sigma[iSy] * dvx;
        p_uv(t, i, j).u = sigmadtinv[iSy] * (p_uv(t-1, i, j).u * sigmadt[iSy]
                      + dtdx * (dvx + p_uv(t, i, j).vy - p_uv(t, i, j-1).vy)
                      + (0.5 * dt) * (psi_old + p_psiLy[iPy]));
    Pochoir_Kernel_2D_End(k_lower_y, shape_wave_2D_u)

    Pochoir_Guard_2D_Begin(g_upper_y, t, i, j)
        if (i >= P && i < Nx-1-P &&
                j >= Ny-1-P && j < Ny-1)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_upper_y)

    Pochoir_Kernel_2D_Begin(k_upper_y, t, i, j)
        int iSy = 2*(j-(Ny-1-P)), iPy = i*P + (j-(Ny-1-P));
        double psi_old = p_psiUy[iPy];
        double dvx = p_uv(t, i, j).vx - p_uv(t, i-1, j).vx;
        p_psiUy[iPy] += dtdx * sigma[iSy] * dvx;
        p_uv(t, i, j).u = sigmadtinv[iSy] * (p_uv(t-1, i, j).u * sigmadt[iSy]
                      + dtdx * (dvx + p_uv(t, i, j).vy - p_uv(t, i, j-1).vy)
                      + (0.5 * dt) * (psi_old + p_psiUy[iPy]));
    Pochoir_Kernel_2D_End(k_upper_y, shape_wave_2D_u)

    Pochoir_Guard_2D_Begin(g_lower_x_lower_y, t, i, j)
        if (i >= 0 && i <= P-1 &&
                j >= 0 && j <= P-1)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_lower_x_lower_y)

    Pochoir_Kernel_2D_Begin(k_lower_x_lower_y, t, i, j)
        int iSx = 2*(P-i), iPx = (i-1)*(Ny-1) + j;
        int iSy = 2*(P-j), iPy = i*P + (j-1);
        int iS = iSx * (2*P) + iSy;
        double psi_old = p_psiLx[iPx] + p_psiLy[iPy];
        double dvx = p_uv(t, i, j).vx - p_uv(t, i-1, j).vx;
        double dvy = p_uv(t, i, j).vy - p_uv(t, i, j-1).vy;
        p_psiLx[iPx] += dtdx * sigma[iSx] * dvy;
        p_psiLy[iPy] += dtdx * sigma[iSy] * dvx;
        p_uv(t, i, j).u = sigmadtinv2[iS] * (p_uv(t-1, i, j).u * sigmadt2[iS]
                      + dtdx * (dvx + dvy)
                      + (0.5 * dt) * (psi_old +
                              p_psiLx[iPx] + p_psiLy[iPy]));
    Pochoir_Kernel_2D_End(k_lower_x_lower_y, shape_wave_2D_u)

    Pochoir_Guard_2D_Begin(g_lower_x_upper_y, t, i, j)
        if (i >= 0 && i <= P-1 &&
                j >= Ny-1-P && j < Ny-1)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_lower_x_upper_y)

    Pochoir_Kernel_2D_Begin(k_lower_x_upper_y, t, i, j)
        int iSx = 2*(P-i), iPx = (i-1)*(Ny-1) + j;
        int iSy = 2*(j-(Ny-1-P)), iPy = i*P + (j-(Ny-1-P));
        int iS = iSx * (2*P) + iSy;
        double psi_old = p_psiLx[iPx] + p_psiUy[iPy];
        double dvx = p_uv(t, i, j).vx - p_uv(t, i-1, j).vx;
        double dvy = p_uv(t, i, j).vy - p_uv(t, i, j-1).vy;
        p_psiLx[iPx] += dtdx * sigma[iSx] * dvy;
        p_psiUy[iPy] += dtdx * sigma[iSy] * dvx;
        p_uv(t, i, j).u = sigmadtinv2[iS] * (p_uv(t-1, i, j).u * sigmadt2[iS]
                      + dtdx * (dvx + dvy)
                      + (0.5 * dt) * (psi_old +
                              p_psiLx[iPx] + p_psiUy[iPy]));
    Pochoir_Kernel_2D_End(k_lower_x_upper_y, shape_wave_2D_u)

    Pochoir_Guard_2D_Begin(g_upper_x_lower_y, t, i, j)
        if (i >= Nx-1-P && i < Nx-1  &&
                j >= 0 && j <= P-1)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_upper_x_lower_y)

    Pochoir_Kernel_2D_Begin(k_upper_x_lower_y, t, i, j)
        int iSx = 2*(i-(Nx-1-P)), iPx = (i-(Nx-1-P))*(Ny-1) + j;
        int iSy = 2*(P-j), iPy = i*P + (j-1);
        int iS = iSx * (2*P) + iSy;
        double psi_old = p_psiUx[iPx] + p_psiLy[iPy];
        double dvx = p_uv(t, i, j).vx - p_uv(t, i-1, j).vx;
        double dvy = p_uv(t, i, j).vy - p_uv(t, i, j-1).vy;
        p_psiUx[iPx] += dtdx * sigma[iSx] * dvy;
        p_psiLy[iPy] += dtdx * sigma[iSy] * dvx;
        p_uv(t, i, j).u = sigmadtinv2[iS] * (p_uv(t-1, i, j).u * sigmadt2[iS]
                      + dtdx * (dvx + dvy)
                      + (0.5 * dt) * (psi_old +
                              p_psiUx[iPx] + p_psiLy[iPy]));
    Pochoir_Kernel_2D_End(k_upper_x_lower_y, shape_wave_2D_u)

    Pochoir_Guard_2D_Begin(g_upper_x_upper_y, t, i, j)
        if (i >= Nx-1-P && i < Nx-1 &&
                j >= Ny-1-P && j < Ny-1)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_upper_x_upper_y)

    Pochoir_Kernel_2D_Begin(k_upper_x_upper_y, t, i, j)
        int iSx = 2*(i-(Nx-1-P)), iPx = (i-(Nx-1-P))*(Ny-1) + j;
        int iSy = 2*(j-(Ny-1-P)), iPy = i*P + (j-(Ny-1-P));
        int iS = iSx * (2*P) + iSy;
        double psi_old = p_psiUx[iPx] + p_psiUy[iPy];
        double dvx = p_uv(t, i, j).vx - p_uv(t, i-1, j).vx;
        double dvy = p_uv(t, i, j).vy - p_uv(t, i, j-1).vy;
        p_psiUx[iPx] += dtdx * sigma[iSx] * dvy;
        p_psiUy[iPy] += dtdx * sigma[iSy] * dvx;
        p_uv(t, i, j).u = sigmadtinv2[iS] * (p_uv(t-1, i, j).u * sigmadt2[iS]
                      + dtdx * (dvx + dvy)
                      + (0.5 * dt) * (psi_old +
                              p_psiUx[iPx] + p_psiUy[iPy]));
    Pochoir_Kernel_2D_End(k_upper_x_upper_y, shape_wave_2D_u)

    Pochoir_Guard_2D_Begin(g_f_g, t, i, j)
        int ix0 = int((X0 + R) * resolution), ix1 = ix0 + int(w * resolution);
        if (t * dt < 2 * t0 && i >= ix0 && i <= ix1 &&
                j == P)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_f_g)

    Pochoir_Kernel_2D_Begin(k_f_g, t, i, j)
        double l_t = t * dt;
        double g = cos(omega * (l_t)) * exp(-(l_t-t0)*(l_t-t0)*decay);
        p_uv(t, i, j).u = p_uv(t-1, i, j).u + dt * g;
    Pochoir_Kernel_2D_End(k_f_g, shape_wave_2D_u)

    /* we have to stagger the computation between u and (vx, vy),
     * so it's if (t - 1 % 2 == 0) compute 'u' 
     *         else if (t - 1 % 2 == 1) compute "vx, vy" 
     */ 
    Pochoir_Guard_2D_Begin(g_vx_interior, t, i, j)
        if (i >= P-1 && i < Nx-1-P &&
                j >= 0 && j < Ny-1)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_vx_interior)

    Pochoir_Kernel_2D_Begin(k_vx_interior, t, i, j)
        int idx = i * (Ny - 1) + j;
        /* we employ the (t+1, ...) = f(t, ...) to avoid the dependency
         * on the same level -- a general graph dependency
         */
        p_uv(t+1, i, j).vx = p_uv(t, i, j).vx + dtdx * p_ax[idx] * (p_uv(t, i+1, j).u - p_uv(t, i, j).u);
    Pochoir_Kernel_2D_End(k_vx_interior, shape_wave_2D_v)

    Pochoir_Guard_2D_Begin(g_vx_lower_x, t, i, j)
        if (i >= 0 && i < P-1 &&
                j >= 0 && j <= Ny-1)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_vx_lower_x)

    Pochoir_Kernel_2D_Begin(k_vx_lower_x, t, i, j)
        int idx = i * (Ny-1) + j, iSx = 2 * (P-i) - 1;
        p_uv(t+1, i, j).vx = sigmadtinv[iSx] * (p_uv(t, i, j).vx * sigmadt[iSx]
                    + dtdx * p_ax[idx] * (p_uv(t, i+1, j).u - p_uv(t, i, j).u));
    Pochoir_Kernel_2D_End(k_vx_lower_x, shape_wave_2D_v)

    Pochoir_Guard_2D_Begin(g_vx_upper_x, t, i, j)
        if (i >= Nx-1-P && i < Nx-1 &&
                j >= 0 && j < Ny-1)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_vx_upper_x)

    Pochoir_Kernel_2D_Begin(k_vx_upper_x, t, i, j)
        int idx = i*Ny + j, iSx = 2*(i-(Nx-1-P)) + 1;
        p_uv(t+1, i, j).vx = sigmadtinv[iSx] * (p_uv(t, i, j).vx * sigmadt[iSx]
                   + dtdx * p_ax[idx] * (p_uv(t, i+1, j).u - p_uv(t, i, j).u));
    Pochoir_Kernel_2D_End(k_vx_upper_x, shape_wave_2D_v)

    Pochoir_Guard_2D_Begin(g_vy_interior, t, i, j)
        if (i >= 0 && i < Nx-1 &&
                j >= P-1 && j < Ny-1-P)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_vy_interior)

    Pochoir_Kernel_2D_Begin(k_vy_interior, t, i, j)
        int idx = i * (Ny-1) + j;
        p_uv(t+1, i, j).vy = p_uv(t, i, j).vy + dtdx * p_ay[idx] * (p_uv(t, i, j+1).u - p_uv(t, i, j).u); 
    Pochoir_Kernel_2D_End(k_vy_interior, shape_wave_2D_v)

    Pochoir_Guard_2D_Begin(g_vy_lower_y, t, i, j)
        if (i >= 0 && i < Nx-1 &&
                j >= 0 && j < P-1)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_vy_lower_y)

    Pochoir_Kernel_2D_Begin(k_vy_lower_y, t, i, j)
        int idx = i * (Ny-1) + j, iSy = 2*(P-j) - 1;
        p_uv(t+1, i, j).vy = sigmadtinv[iSy] * (p_uv(t, i, j).vy * sigmadt[iSy]
                   + dtdx * p_ay[idx] * (p_uv(t, i, j+1).u - p_uv(t, i, j).u));
    Pochoir_Kernel_2D_End(k_vy_lower_y, shape_wave_2D_v)

    Pochoir_Guard_2D_Begin(g_vy_upper_y, t, i, j)
        if (i >= 0 && i < Nx-1 &&
                j >= Ny-1-P && j < Ny-1)
            return true;
        else
            return false;
    Pochoir_Guard_2D_End(g_vy_upper_y)

    Pochoir_Kernel_2D_Begin(k_vy_upper_y, t, i, j)
        int idx = i*(Ny-1) + j, iSy = 2*(j-(Ny-1-P)) + 1;
        p_uv(t+1, i, j).vy = sigmadtinv[iSy] * (p_uv(t, i, j).vy * sigmadt[iSy]
                   + dtdx * p_ay[idx] * (p_uv(t, i, j+1).u - p_uv(t, i, j).u));
    Pochoir_Kernel_2D_End(k_vy_upper_y, shape_wave_2D_v)

    Pochoir_Guard_2D_Begin(g_output_u, t, i, j)
        if (t % 100 == 0)
            return true;
        else 
            return false;
    Pochoir_Guard_2D_End(g_output_u)

    Pochoir_Kernel_2D_Begin(k_output_u, t, i, j)
        printf("on time step %d / %d\n", t, int(T/dt)); fflush(stdout);
        if (t * dt > t0 - 1/fwidth)
            output_p_u("p_u", t, P, p_uv, Nx-1, Ny-1);
    Pochoir_Kernel_2D_End(k_output_u, shape_wave_2D_u)

#if 0
    wave_2D.Register_Exclusive_Tile_Kernels(g_interior, k_interior);
    wave_2D.Register_Exclusive_Tile_Kernels(g_lower_x, k_lower_x);
    wave_2D.Register_Exclusive_Tile_Kernels(g_upper_x, k_upper_x);
    wave_2D.Register_Exclusive_Tile_Kernels(g_lower_y, k_lower_y);
    wave_2D.Register_Exclusive_Tile_Kernels(g_upper_y, k_upper_y);
    wave_2D.Register_Exclusive_Tile_Kernels(g_lower_x_lower_y, k_lower_x_lower_y);
    wave_2D.Register_Exclusive_Tile_Kernels(g_lower_x_upper_y, k_lower_x_upper_y);
    wave_2D.Register_Exclusive_Tile_Kernels(g_upper_x_lower_y, k_upper_x_lower_y);
    wave_2D.Register_Exclusive_Tile_Kernels(g_upper_x_upper_y, k_upper_x_upper_y);
    // Original location
    // wave_2D.Register_Inclusive_Tile_Kernels(g_f_g, k_f_g);
    wave_2D.Register_Exclusive_Tile_Kernels(g_vx_interior, k_vx_interior);
    wave_2D.Register_Exclusive_Tile_Kernels(g_vx_lower_x, k_vx_lower_x);
    wave_2D.Register_Exclusive_Tile_Kernels(g_vx_upper_x, k_vx_upper_x);
    wave_2D.Register_Exclusive_Tile_Kernels(g_vy_interior, k_vy_interior);
    wave_2D.Register_Exclusive_Tile_Kernels(g_vy_lower_y, k_vy_lower_y);
    wave_2D.Register_Exclusive_Tile_Kernels(g_vy_upper_y, k_vy_upper_y);
    wave_2D.Register_Inclusive_Tile_Kernels(g_f_g, k_f_g);
    wave_2D.Register_Inclusive_Tile_Kernels(g_output_u, k_output_u);
#else
    wave_2D.Register_Exclusive_Tile_Kernels(g_interior, k_interior);
    wave_2D.Register_Exclusive_Tile_Kernels(g_lower_x, k_lower_x);
    wave_2D.Register_Exclusive_Tile_Kernels(g_upper_x, k_upper_x);
    wave_2D.Register_Exclusive_Tile_Kernels(g_lower_y, k_lower_y);
    wave_2D.Register_Exclusive_Tile_Kernels(g_upper_y, k_upper_y);
    wave_2D.Register_Exclusive_Tile_Kernels(g_lower_x_lower_y, k_lower_x_lower_y);
    wave_2D.Register_Exclusive_Tile_Kernels(g_lower_x_upper_y, k_lower_x_upper_y);
    wave_2D.Register_Exclusive_Tile_Kernels(g_upper_x_lower_y, k_upper_x_lower_y);
    wave_2D.Register_Exclusive_Tile_Kernels(g_upper_x_upper_y, k_upper_x_upper_y);
    // Original location
    wave_2D.Register_Inclusive_Tile_Kernels(g_f_g, k_f_g);
    wave_2D.Register_Inclusive_Tile_Kernels(g_vx_interior, k_vx_interior);
    wave_2D.Register_Inclusive_Tile_Kernels(g_vx_lower_x, k_vx_lower_x);
    wave_2D.Register_Inclusive_Tile_Kernels(g_vx_upper_x, k_vx_upper_x);
    wave_2D.Register_Inclusive_Tile_Kernels(g_vy_interior, k_vy_interior);
    wave_2D.Register_Inclusive_Tile_Kernels(g_vy_lower_y, k_vy_lower_y);
    wave_2D.Register_Inclusive_Tile_Kernels(g_vy_upper_y, k_vy_upper_y);
    // wave_2D.Register_Inclusive_Tile_Kernels(g_f_g, k_f_g);
    wave_2D.Register_Inclusive_Tile_Kernels(g_output_u, k_output_u);
#endif
    wave_2D.Register_Array(p_uv);
    for (int i = 0; i < Nx - 1; ++i)
        for (int j = 0; j < Ny - 1; ++j) {
            p_uv(1, i, j).u = p_uv(1, i, j).vx = p_uv(1, i, j).vy = 0;
            p_uv(0, i, j).u = p_uv(0, i, j).vx = p_uv(0, i, j).vy = 0;
        }

    int l_T = int(T/dt);
    Pochoir_Plan<2> & l_plan = wave_2D.Gen_Plan(l_T);
    min_tdiff = INF;
    for (int times = 0; times < TIMES; ++times) {
        gettimeofday(&start, 0);
        wave_2D.Run(l_plan);
        gettimeofday(&end, 0);
        min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    }
    cout << "Pochoir time : " << min_tdiff << " ms " << endl;

    /* guard and kernel function for Pochoir ends */

    /* we now get the total time step T directly from user's command line argument
     * Also note there looks a bug to define T to be of double
     * double T = 200; // total run time
     */
    /* start original version of stencil computation on 2D wave equation */
  
    gettimeofday(&start, 0);
    int it;
    for (int times = 0; times < TIMES; ++times) {
    for (it = 0; it * dt < T; ++it) {
        /////////////////////////////////////////////////////////////////////////
        // Update u in all 9 disjoint interior/boundary regions:
        
        // Update u in interior (non-PML) region: du = dt * div(v)
        for (int ix = 1 + P; ix < Nx - P; ++ix)
            for (int iy = 1 + P; iy < Ny - P; ++iy) {
                int i = ix*Ny + iy;
                u[i] += dtdx * (vx[i]-vx[i-Ny] + vy[i]-vy[i-1]);
            }

        // Update u in 8 PML regions (upper and lower x and y + 4 corners),
        // at the same time updating the (co-located) psi arrays.  (This
        // has to be done at the same time to preserve 2nd-order accuracy.)
        for (int ix = 1; ix <= P; ++ix) // lower x layer
            for (int iy = 1 + P; iy < Ny - P; ++iy) {
                int i = ix*Ny + iy, iSx = 2*(P-ix), iPx = (ix-1)*Ny + iy;
                double psi_old = psiLx[iPx];
                double dvy = vy[i]-vy[i-1];
                psiLx[iPx] += dtdx * sigma[iSx] * dvy;
                u[i] = sigmadtinv[iSx] * (u[i] * sigmadt[iSx]
                            + dtdx * (vx[i]-vx[i-Ny] + dvy)
                            + (0.5 * dt) * (psi_old + psiLx[iPx]));
            }
        for (int ix = Nx - P; ix < Nx; ++ix) // upper x layer
            for (int iy = 1 + P; iy < Ny - P; ++iy) {
                int i = ix*Ny + iy, iSx = 2*(ix-(Nx-P)), iPx = (ix-(Nx-P))*Ny + iy;
                double psi_old = psiUx[iPx];
                double dvy = vy[i]-vy[i-1];
                psiUx[iPx] += dtdx * sigma[iSx] * dvy;
                u[i] = sigmadtinv[iSx] * (u[i] * sigmadt[iSx]
                            + dtdx * (vx[i]-vx[i-Ny] + dvy)
                            + (0.5 * dt) * (psi_old + psiUx[iPx]));
            }
        for (int ix = 1 + P; ix < Nx - P; ++ix)
            for (int iy = 1; iy <= P; ++iy) { // lower y layer
                int i = ix*Ny + iy, iSy = 2*(P-iy), iPy = ix*P + (iy-1);
                double psi_old = psiLy[iPy];
                double dvx = vx[i]-vx[i-Ny];
                psiLy[iPy] += dtdx * sigma[iSy] * dvx;
                u[i] = sigmadtinv[iSy] * (u[i] * sigmadt[iSy]
                            + dtdx * (dvx + vy[i]-vy[i-1])
                            + (0.5 * dt) * (psi_old + psiLy[iPy]));
            }
        for (int ix = 1 + P; ix < Nx - P; ++ix)
            for (int iy = Ny - P; iy < Ny; ++iy) { // upper y layer
        int i = ix*Ny + iy, iSy = 2*(iy-(Ny-P)), iPy = ix*P + (iy-(Ny-P));
        double psi_old = psiUy[iPy];
        double dvx = vx[i]-vx[i-Ny];
        psiUy[iPy] += dtdx * sigma[iSy] * dvx;
        u[i] = sigmadtinv[iSy] * (u[i] * sigmadt[iSy]
                    + dtdx * (dvx + vy[i]-vy[i-1])
                    + (0.5 * dt) * (psi_old + psiUy[iPy]));
            }
        for (int ix = 1; ix <= P; ++ix) // lower x layer
            for (int iy = 1; iy <= P; ++iy) { // lower y layer
                int i = ix*Ny + iy;
                int iSx = 2*(P-ix), iPx = (ix-1)*Ny + iy;
                int iSy = 2*(P-iy), iPy = ix*P + (iy-1);
                int iS = iSx * (2*P) + iSy;
                double psi_old = psiLx[iPx] + psiLy[iPy];
                double dvx = vx[i]-vx[i-Ny], dvy = vy[i]-vy[i-1];
                psiLx[iPx] += dtdx * sigma[iSx] * dvy;
                psiLy[iPy] += dtdx * sigma[iSy] * dvx;
                u[i] = sigmadtinv2[iS] * (u[i] * sigmadt2[iS]
                            + dtdx * (dvx + dvy)
                            + (0.5 * dt) * (psi_old + 
                                    psiLx[iPx] + psiLy[iPy]));
            }
        for (int ix = 1; ix <= P; ++ix) // lower x layer
            for (int iy = Ny - P; iy < Ny; ++iy) { // upper y layer
                int i = ix*Ny + iy;
                int iSx = 2*(P-ix), iPx = (ix-1)*Ny + iy;
                int iSy = 2*(iy-(Ny-P)), iPy = ix*P + (iy-(Ny-P));
                int iS = iSx * (2*P) + iSy;
                double psi_old = psiLx[iPx] + psiUy[iPy];
                double dvx = vx[i]-vx[i-Ny], dvy = vy[i]-vy[i-1];
                psiLx[iPx] += dtdx * sigma[iSx] * dvy;
                psiUy[iPy] += dtdx * sigma[iSy] * dvx;
                u[i] = sigmadtinv2[iS] * (u[i] * sigmadt2[iS]
                            + dtdx * (dvx + dvy)
                            + (0.5 * dt) * (psi_old + 
                                    psiLx[iPx] + psiUy[iPy]));
            }

        for (int ix = Nx - P; ix < Nx; ++ix) // upper x layer
            for (int iy = 1; iy <= P; ++iy) { // lower y layer
                int i = ix*Ny + iy;
                int iSx = 2*(ix-(Nx-P)), iPx = (ix-(Nx-P))*Ny + iy;
                int iSy = 2*(P-iy), iPy = ix*P + (iy-1);
                int iS = iSx * (2*P) + iSy;
                double psi_old = psiUx[iPx] + psiLy[iPy];
                double dvx = vx[i]-vx[i-Ny], dvy = vy[i]-vy[i-1];
                psiUx[iPx] += dtdx * sigma[iSx] * dvy;
                psiLy[iPy] += dtdx * sigma[iSy] * dvx;
                u[i] = sigmadtinv2[iS] * (u[i] * sigmadt2[iS]
                            + dtdx * (dvx + dvy)
                            + (0.5 * dt) * (psi_old + 
                                    psiUx[iPx] + psiLy[iPy]));
            }
        for (int ix = Nx - P; ix < Nx; ++ix) // upper x layer
            for (int iy = Ny - P; iy < Ny; ++iy) { // upper y layer
                int i = ix*Ny + iy;
                int iSx = 2*(ix-(Nx-P)), iPx = (ix-(Nx-P))*Ny + iy;
                int iSy = 2*(iy-(Ny-P)), iPy = ix*P + (iy-(Ny-P));
                int iS = iSx * (2*P) + iSy;
                double psi_old = psiUx[iPx] + psiUy[iPy];
                double dvx = vx[i]-vx[i-Ny], dvy = vy[i]-vy[i-1];
                psiUx[iPx] += dtdx * sigma[iSx] * dvy;
                psiUy[iPy] += dtdx * sigma[iSy] * dvx;
                u[i] = sigmadtinv2[iS] * (u[i] * sigmadt2[iS]
                            + dtdx * (dvx + dvy)
                            + (0.5 * dt) * (psi_old + 
                                    psiUx[iPx] + psiUy[iPy]));
            }

        // Update halo of u from periodic boundary conditions.
        // (In u's case, update x=0 and y=0 edges from opposite edges.)
        for (int ix = 0; ix < Nx; ++ix) u[ix*Nx] = u[ix*Nx + Ny-1];
        for (int iy = 0; iy < Ny; ++iy) u[iy] = u[Ny*(Nx-1) + iy];

        /////////////////////////////////////////////////////////////////////////
        // Add source term f(x,y,t) in the u equation.  As is typical,
        // we will have a separatble f(x,y,t) = f(x,y) g(t), where in
        // this case g is a Gaussian pulse and f(x,y) is 1 in a small
        // region [a line segment crossing one end of the bend, designed
        // to excite waveguide mode(s)].

        double t = it * dt;
        if (t < 2*t0) {
            double g = cos(omega*t) * exp(-(t-t0)*(t-t0)*decay);
            int iy = P + 1; // right on the edge of the PML
            int ix0 = int((X0 + R) * resolution), ix1 = ix0 + int(w * resolution);
            for (int ix = ix0; ix <= ix1; ++ix)
                u[ix*Ny + iy] += dt * g;
        }
        else
            ; // nop: turn off source after gaussian has decayed to almost zero
        
        /////////////////////////////////////////////////////////////////////////
        // Update vx and vy.  (Note that these could be updated in parallel.)

        /* Update vx in 3 disjoint interior/boundary regions: 
             -- The vx equation is only modified by the PML in the x direction. 
           (Note that, because vx is discretized at (ix+0.5)*dx, the
            *upper* x edge is now the halo, so some loop bounds change.
            Also, for the same reason, the lower PML boundary shifts by 1.) */
        
        // Update vx in interior (non-PML) region: dvx = dt * ax * du/dx
        for (int ix = P; ix < Nx - P; ++ix)
            for (int iy = 1; iy < Ny; ++iy) {
                int i = ix*Ny + iy;
                vx[i] += dtdx * ax[i] * (u[i+Ny]-u[i]);
            }

        // Update vx in 2 PML regions (upper and lower x):
        for (int ix = 0; ix < P; ++ix) // lower x layer
            for (int iy = 1; iy < Ny; ++iy) {
                int i = ix*Ny + iy, iSx = 2*(P-ix) - 1;
                vx[i] = sigmadtinv[iSx] * (vx[i] * sigmadt[iSx]
                         + dtdx * ax[i] * (u[i+Ny]-u[i]));
            }
        for (int ix = Nx - P; ix < Nx - 1; ++ix) // upper x layer
            for (int iy = 1; iy < Ny; ++iy) {
                int i = ix*Ny + iy, iSx = 2*(ix-(Nx-P)) + 1;
                vx[i] = sigmadtinv[iSx] * (vx[i] * sigmadt[iSx]
                         + dtdx * ax[i] * (u[i+Ny]-u[i]));
            }

        // Update halo of vx from periodic boundary conditions.
        // (In vx's case, update *upper* x and *lower* y from opposite sides.)
        for (int ix = 0; ix < Nx; ++ix) vx[ix*Nx] = vx[ix*Nx + Ny-1];
        for (int iy = 0; iy < Ny; ++iy) vx[Ny*(Nx-1) + iy] = vx[iy];

        /* Update vy in 3 disjoint interior/boundary regions: 
             -- The vy equation is only modified by the PML in the y direction. 
           (Note that, because vy is discretized at (iy+0.5)*dy, the
            *upper* y edge is now the halo, so some loop bounds change.
            Also, for the same reason, the lower PML boundary shifts by 1.) */
        
        // Update vy in interior (non-PML) region: dvy = dt * ax * du/dx
        for (int ix = 1; ix < Nx; ++ix)
            for (int iy = P; iy < Ny - P; ++iy) {
                int i = ix*Ny + iy;
                vy[i] += dtdx * ay[i] * (u[i+1]-u[i]);
            }

        // Update vy in 2 PML regions (upper and lower y):
        for (int ix = 1; ix < Nx; ++ix)
            for (int iy = 0; iy < P; ++iy) { // lower y layer
                int i = ix*Ny + iy, iSy = 2*(P-iy) - 1;
                vy[i] = sigmadtinv[iSy] * (vy[i] * sigmadt[iSy]
                         + dtdx * ay[i] * (u[i+1]-u[i]));
            }
        for (int ix = 1; ix < Nx; ++ix)
            for (int iy = Ny - P; iy < Ny - 1; ++iy) { // upper y layer
                int i = ix*Ny + iy, iSy = 2*(iy-(Ny-P)) + 1;
                vy[i] = sigmadtinv[iSy] * (vy[i] * sigmadt[iSy]
                         + dtdx * ay[i] * (u[i+1]-u[i]));
            }

        // Update halo of vy from periodic boundary conditions.
        // (In vy's case, update *lower* x and *upper* y from opposite sides.)
        for (int ix = 0; ix < Nx; ++ix) vy[ix*Nx + Ny-1] = vy[ix*Nx];
        for (int iy = 0; iy < Ny; ++iy) vy[iy] = vy[Ny*(Nx-1) + iy];

        /////////////////////////////////////////////////////////////////////////
        // Output u in non-PML (interior) region every 100 time steps
        // after source has turned on to a reasonable magnitude:
        if (it % 100 == 0) {
            printf("on time step %d / %d\n", it, int(T/dt)); fflush(stdout);
            if (it * dt > t0 - 1/fwidth)
                output("u", it, u + P*Ny + P, Nx - 2*P, Ny - 2*P, Ny);
        }

        /////////////////////////////////////////////////////////////////////////

    } /* end iteration on it */
    } /* end iteration on times */
    gettimeofday(&end, 0);
    min_tdiff = min(min_tdiff, (1.0e3 * tdiff(&end, &start)));
    cout << "Serial Loop time : " << min_tdiff << " ms " << endl;

#if 1
    /* check results! */
    int t = it;
    for (int i = 0; i < Nx - 1; ++i)
        for (int j = 0; j < Ny - 1; ++j) {
            /* for array 'u', we shift '1' on both i and j dimension because
             * of the halo point
             */
            check_result(t, i, j, p_uv(t, i, j).u, u[(i+1) * Ny + (j+1)]);
        }
#endif
    return 0;
}
