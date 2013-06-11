/* Test case - 2D diffusion equation, Non-Periodic version */
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

int main(int argc, char * argv[])
{
	const int BASE = 1024;
	int t;
	struct timeval start, end;
    int Nx = 500, Ny = 100, T = 731;
	int resolution = 40 ;
	// some geometric parameters:
	const double sx = 5, sy = 5, dpml = 0.5; // domain size
	const double R = 3.0; // inner bend radius
	const double w = 0.2; // waveguide width
	const double aw = 0.1; // value of "a" inside waveguide
	const double X0 = dpml + 1, Y0 = dpml + 1; // bend center
	// Gaussian-pulse characteristics of source
	double fcen = 1.0, fwidth = 0.1; // center frequency and width
	const double twopi = atan(1.0)*8.0;
	// source(t) = cos(omega*t) * exp[- decay*(t-t0)^2]:
	//double omega = twopi * fcen, decay = (twopi*twopi*fwidth*fwidth) * 0.5;
	double omega = twopi * fcen, decay = (twopi*twopi*fwidth*fwidth) * 1e-15 ;
	double t0_ = 5/fwidth; // start time of gaussian
	if (argc < 2) {
        printf("argc < 2, quit! \n");
        exit(1);
    }
	resolution = StrToInt(argv[1]);

	int ix0 = int((X0 + R) * resolution), ix1 = ix0 + int(w * resolution); 
	double dx = 1.0/resolution; // spatial grid spacing
	double dt = 0.7 * dx; // time-step size, from CFL condition for c = 1
	//choose diffusion co-efficient Dm < 0.5 * dx * dx / dt
	double Dm = 0.5 * dx * dx / dt * 1e-2;
	cout << "Dm " << Dm << endl ;
	T = 200 / dt ;
	Nx = int((sx + 2 * dpml) * resolution) + 1;
	Ny = int((sy + 2 * dpml) * resolution) + 1;
	int P = int(dpml * resolution); // PML thickness, in piels
    printf("Resolution = %d, P = %d, Nx = %d, Ny = %d, T = %d\n", resolution, 
					P, Nx, Ny, T) ;
    Pochoir_Shape_2D heat_shape_2D[] = {{0, 0, 0}, {-1, 1, 0}, {-1, 0, 0}, {-1, -1, 0}, {-1, 0, -1}, {-1, 0, 1}};
    Pochoir<N_RANK> heat_2D(heat_shape_2D);
	Pochoir_Array<double, N_RANK> b(Nx, Ny);
	Pochoir_Array<double, N_RANK> a(Nx, Ny) ; 
    a.Register_Boundary(aperiodic_2D) ;
    heat_2D.Register_Array(a) ;
	heat_2D.set_resolution (resolution) ;
    b.Register_Shape(heat_shape_2D);
    b.Register_Boundary(aperiodic_2D);

	double * c = new double [2 * (P + 1)] ;

	for (int i = 0 ; i < 2 * P + 2 ; i++)
	{
		c [i] = (1 - i / (2 * P)) * (1 - i / (2 * P)) ;
	} 

    /* Now we can only access the Pochoir_Array after Register_Array,
     * or Register_Shape with the array, because we rely on the shape
     * to provide the depth of toggle array!!! 
     */
	for (int i = 0; i < Nx; ++i) {
	for (int j = 0; j < Ny; ++j) {
        a(0, i, j) = 1.0 * (rand() % BASE) * 1e-3; 
        a(1, i, j) = 0 ; 
        b(0, i, j) = a(0, i, j);
        b(1, i, j) = 0;
	} }

	cout << "a(T+1, J, I) = 0.125 * (a(T, J+1, I) - 2.0 * a(T, J, I) + a(T, J-1, I)) + 0.125 * (a(T, J, I+1) - 2.0 * a(T, J, I) + a(T, J, I-1)) + a(T, J, I)" << endl;
	cout << "t0 " << t0_ << " ix0 " << ix0 << " ix1 " << ix1 << endl ;
    Pochoir_Kernel_2D(heat_2D_fn, t, i, j)
	assert (i >= 0 && i < Nx) ;
	assert (j >= 0 && j < Ny) ;

	if (i >= P && i <= Nx - P && j >= P && j <= Ny - P)
	{
		//interior
	   	a(t, i, j) = Dm * (a(t-1, i+1, j) - 2.0 * a(t-1, i, j) + a(t-1, i-1, j)) + Dm * (a(t-1, i, j+1) - 2.0 * a(t-1, i, j) + a(t-1, i, j-1)) + a(t-1, i, j);

		if (t * dt < 2 * t0_ && j == P + 1 && i >= ix0 && i <= ix1)
		{
			//source
			double ts = t * dt ;
			double g = cos(omega * ts) * exp(- (ts - t0_) * (ts - t0_) * decay);
			a(t, i, j) += g * dt ;
		}
	}
	else if (j >= P && j <= Ny - P)
	{
		int k ;
		if (i < P)
		{
			assert (i >= 0) ;
			//west
			k = 2 * (P - i) ;
		}
		else
		{
			assert (i > Nx - P && i < Nx) ;
			//east
			k = 2 * (i - (Nx - P)) ;
		}
	   	a(t, i, j) = Dm * c [k] * 
			(c [k + 1] * (a(t-1, i+1, j) - a(t-1, i, j)) - 
			 c [k - 1] * (a(t-1, i, j) - a(t-1, i-1, j))) + 
			Dm * (a(t-1, i, j+1) - 2.0 * a(t-1, i, j) + a(t-1, i, j-1)) 
			+ a(t-1, i, j) ;
	}
	else if (i >= P && i <= Nx - P)
	{
		int k ;
		if (j < P)
		{
			assert (j >= 0) ;
			//south
			k = 2 * (P - j) ;
		}
		else
		{
			assert (j > Ny - P && j < Ny) ;
			//north
			k = 2 * (j - (Ny - P)) ;
		}
	   	a(t, i, j) = Dm * c [k] * 
			(c [k + 1] * (a(t-1, i, j + 1) - a(t-1, i, j)) - 
			 c [k - 1] * (a(t-1, i, j) - a(t-1, i, j-1))) + 
			Dm * (a(t-1, i+1, j) - 2.0 * a(t-1, i, j) + a(t-1, i-1, j)) 
			+ a(t-1, i, j) ; 
	}
	else
	{
		int kx, ky ;
		if (i < P)
		{
			assert (i >= 0) ; 
			//west
			kx = 2 * (P - i) ;
		}
		else
		{
			assert (i > Nx - P) ;
			//east
			kx = 2 * (i - (Nx - P)) ;
		}
		
		if (j < P)
		{
			assert (j >= 0) ;
			//south
			ky = 2 * (P - j) ;
		}
		else
		{
			assert(j > Ny - P) ;
			//north
			ky = 2 * (j - (Ny - P)) ;
		}
		a(t, i, j) = Dm * c [ky] * 
			(c [ky + 1] * (a(t-1, i, j + 1) - a(t-1, i, j)) - 
			 c [ky - 1] * (a(t-1, i, j) - a(t-1, i, j-1))) + 
			Dm * c [kx] *
			(c [kx + 1] * (a(t-1, i+1, j) - a(t-1, i, j)) -
			 c [kx - 1] * (a(t-1, i, j) - a(t-1, i-1, j))) +
			+ a(t-1, i, j) ; 
	}
	
    Pochoir_Kernel_End

	gettimeofday(&start, 0);
    for (int times = 0; times < TIMES; ++times) {
        heat_2D.Run(T, heat_2D_fn);
    }
	gettimeofday(&end, 0);
	std::cout << "Pochoir ET: consumed time :" << 1.0e3 * tdiff(&end, &start)/TIMES << "ms" << std::endl; 

#if 1	
	gettimeofday(&start, 0);
    for (int times = 0; times < TIMES; ++times) {
	for (int t = 0; t < T; ++t) {

	cilk_for (int i = P; i <= Nx - P; ++i) { //interior
	for (int j = P; j <= Ny - P; ++j) {
		 b(t+1, i, j) = Dm * (b(t, i+1, j) - 2.0 * b(t, i, j) + b(t, i-1, j)) + Dm * (b(t, i, j+1) - 2.0 * b(t, i, j) + b(t, i, j-1)) + b(t, i, j);
	}
	}

	cilk_for (int i = 0 ; i < P; ++i)  { // lower x layer
	for (int j = P; j <= Ny - P; ++j) {	
		int k = 2 * (P - i) ;
	   	b(t+1, i, j) = Dm * c [k] * 
			(c [k + 1] * (b(t, i+1, j) - b(t, i, j)) - 
			 c [k - 1] * (b(t, i, j) - b(t, i-1, j))) + 
			Dm * (b(t, i, j+1) - 2.0 * b(t, i, j) + b(t, i, j-1)) 
			+ b(t, i, j) ; 
	}
	}

    cilk_for (int i = Nx - P + 1 ; i < Nx; ++i) {// upper x layer
      for (int j = P; j <= Ny - P; ++j) {
		int k = 2 * (i - (Nx - P)) ;
	   	b(t+1, i, j) = Dm * c [k] * 
			(c [k + 1] * (b(t, i+1, j) - b(t, i, j)) - 
			 c [k - 1] * (b(t, i, j) - b(t, i-1, j))) + 
			Dm * (b(t, i, j+1) - 2.0 * b(t, i, j) + b(t, i, j-1)) 
			+ b(t, i, j) ;
	}
	}

    cilk_for (int i = P; i <= Nx - P; ++i) {
      for (int j = 0; j < P; ++j) { // lower y layer
		int k = 2 * (P - j) ;
	   	b(t+1, i, j) = Dm * c [k] * 
			(c [k + 1] * (b(t, i, j + 1) - b(t, i, j)) - 
			 c [k - 1] * (b(t, i, j) - b(t, i, j-1))) + 
			Dm * (b(t, i+1, j) - 2.0 * b(t, i, j) + b(t, i-1, j)) 
			+ b(t, i, j) ; 
	}
	}

	cilk_for (int i = P; i <= Nx - P; ++i) {
		for (int j = Ny - P + 1; j < Ny; ++j) { // upper y layer	
			int k = 2 * (j - (Ny - P)) ;
	   		b(t+1, i, j) = Dm * c [k] * 
			(c [k + 1] * (b(t, i, j + 1) - b(t, i, j)) - 
			 c [k - 1] * (b(t, i, j) - b(t, i, j-1))) + 
			Dm * (b(t, i+1, j) - 2.0 * b(t, i, j) + b(t, i-1, j)) 
			+ b(t, i, j) ; 
	}
	}

	cilk_for (int i = 0; i < P; ++i) { // lower x layer
		for (int j = 0; j < P; ++j) { // lower y layer
			int kx = 2 * (P - i) ;
			int ky = 2 * (P - j) ; 
			b(t+1, i, j) = Dm * c [ky] * 
				(c [ky + 1] * (b(t, i, j + 1) - b(t, i, j)) - 
				 c [ky - 1] * (b(t, i, j) - b(t, i, j-1))) + 
				Dm * c [kx] *
				(c [kx + 1] * (b(t, i+1, j) - b(t, i, j)) -
				 c [kx - 1] * (b(t, i, j) - b(t, i-1, j))) +
				+ b(t, i, j) ; 
	}
	}

	cilk_for (int i = 0; i < P; ++i) { // lower x layer
		for (int j = Ny - P + 1; j < Ny; ++j) { // upper y layer	
			int kx = 2 * (P - i) ;
			int ky = 2 * (j - (Ny - P)) ;
			b(t+1, i, j) = Dm * c [ky] * 
				(c [ky + 1] * (b(t, i, j + 1) - b(t, i, j)) - 
				 c [ky - 1] * (b(t, i, j) - b(t, i, j-1))) + 
				Dm * c [kx] *
				(c [kx + 1] * (b(t, i+1, j) - b(t, i, j)) -
				 c [kx - 1] * (b(t, i, j) - b(t, i-1, j))) +
				+ b(t, i, j) ; 
	}
	}

	cilk_for (int i = Nx - P + 1 ; i < Nx; ++i) { // upper x layer
		for (int j = 0; j < P; ++j) { // lower y layer
			int kx = 2 * (i - (Nx - P)) ;
			int ky = 2 * (P - j) ; 
			b(t+1, i, j) = Dm * c [ky] * 
				(c [ky + 1] * (b(t, i, j + 1) - b(t, i, j)) - 
				 c [ky - 1] * (b(t, i, j) - b(t, i, j-1))) + 
				Dm * c [kx] *
				(c [kx + 1] * (b(t, i+1, j) - b(t, i, j)) -
				 c [kx - 1] * (b(t, i, j) - b(t, i-1, j))) +
				+ b(t, i, j) ; 
	}
	}

	cilk_for (int i = Nx - P + 1 ; i < Nx; ++i) {// upper x layer
		for (int j = Ny - P + 1 ; j < Ny; ++j) { // upper y layer
			int kx = 2 * (i - (Nx - P)) ;
			int ky = 2 * (j - (Ny - P)) ;
			b(t+1, i, j) = Dm * c [ky] * 
				(c [ky + 1] * (b(t, i, j + 1) - b(t, i, j)) - 
				 c [ky - 1] * (b(t, i, j) - b(t, i, j-1))) + 
				Dm * c [kx] *
				(c [kx + 1] * (b(t, i+1, j) - b(t, i, j)) -
				 c [kx - 1] * (b(t, i, j) - b(t, i-1, j))) +
				+ b(t, i, j) ; 
	}
	}

	//source
	double ts = (t + 1) * dt ;
	if (ts < 2 * t0_) {
		int j = P + 1; // right on the edge of the PML
		double g = cos(omega * ts) * exp(- (ts - t0_) * (ts - t0_) * decay) ;
		cilk_for (int i = ix0; i <= ix1; ++i) {
			b(t+1, i, j) += g * dt ; 
		}
	}
	}
	}
	gettimeofday(&end, 0);
	std::cout << "Naive Loop: consumed time :" << 1.0e3 * tdiff(&end, &start)/TIMES << "ms" << std::endl;

	t = T;
	for (int i = 0; i < Nx; ++i) {
	for (int j = 0; j < Ny; ++j) {
		check_result(t, i, j, a.interior(t, i, j), b.interior(t, i, j));
	} }
#endif
	delete [] c ;
	return 0;
}
