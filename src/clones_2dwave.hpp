#ifndef CLONES_2D_WAVE_HPP
#define CLONES_2D_WAVE_HPP

template <int DIM>
class predicate
{
	int m_width [DIM] ;

	int resolution = 40 ;
	// some geometric parameters:
	const double sx = 5, sy = 5, dpml = 0.5; // domain size
	const double R = 3.0; // inner bend radius
	const double w = 0.2; // waveguide width
	const double aw = 0.1; // value of "a" inside waveguide
	const double X0 = dpml + 1, Y0 = dpml + 1; // bend center
	// Gaussian-pulse characteristics of source
	const double fcen = 1.0, fwidth = 0.1 ; // center frequency and width
	//const double t0_ = 5/fwidth; // start time of gaussian
	//const double twopi = atan(1.0)*8.0;
	// source(t) = cos(omega*t) * exp[- decay*(t-t0)^2]:
	//double omega = twopi * fcen, decay = (twopi*twopi*fwidth*fwidth) * 0.5;
	
	int ix0, ix1  ; 
	double dx ; // spatial grid spacing
	double dt ; // time-step size, from CFL condition for c = 1
	int Nx ;
    int Ny ;
    int P ; // PML thickness, in pixels

	void initialize_parameters()
	{
		// some geometric parameters:
		// Gaussian-pulse characteristics of source
		//fcen = 1.0, fwidth = 0.1; // center frequency and width
		//const double twopi = atan(1.0)*8.0;
		// source(t) = cos(omega*t) * exp[- decay*(t-t0)^2]:
		//double omega = twopi * fcen, decay = (twopi*twopi*fwidth*fwidth) * 0.5;
		
		ix0 = int((X0 + R) * resolution), ix1 = ix0 + int(w * resolution); 
		dx = 1.0/resolution; // spatial grid spacing
		dt = 0.7 * dx; // time-step size, from CFL condition for c = 1
		Nx = int((sx + 2 * dpml) * resolution) + 1;
		Ny = int((sy + 2 * dpml) * resolution) + 1;
		P = int(dpml * resolution); // PML thickness, in pixels
		cout << "resolution " << resolution << " Nx " << Nx << " Ny " << Ny <<
			" P " << P << " t0 " << t0_ << endl ;
	}
public :
	const double t0_ = 5/fwidth; // start time of gaussian
	predicate(const grid_info<DIM> & physical_grid)
	{
		for (int i = 0 ; i < DIM ; i++)
		{
			m_width [i] = physical_grid.x1 [i] - physical_grid.x0 [i] ;
		}
		initialize_parameters () ;
	}

	//return the index of the predicate corresponding to the spatial cell
	int operator() (int x) const ;
	int operator() (int x1, int x0) const ;

	void set_resolution(int r)
	{
		resolution = r ;
		initialize_parameters ();
	}
} ;

template<>
int predicate<1>::operator() (int x) const
{
	//divide the line into two segments
	int index = 0 ;
	if (x <=  m_width [0] / 2)
	{
		index = 0 ;
	}
	else
	{
		index = 1 ;
	}
	return index ;
}

template<>
int predicate<2>::operator() (int i, int j) const
{
	int index = 0 ;
	//cout << "( " << i << " , " << j << " ), " ;
	//if (t < 2 * t0_)
	//{
		//cout << t << " " << t0_ << endl ;
	if (j == P + 1 && i >= ix0 && i <= ix1)
	{
		//source
		index = 0 ;
	}
	//}
	else if (i > P && i < Nx - P && j > P && j < Ny - P)
	{
		//interior
		index = 1 ;
	}
	else if (i <= P && j > P && j < Ny - P)
	{
		assert (i >= 0) ;
		//lower x layer
		index = 2 ;
	}
	else if (i >= Nx - P && j > P && j < Ny - P)
	{
		assert (i < Nx) ;
		//upper x layer
		index = 3 ;
	}
	else if (i > P && i < Nx - P && j <= P)
	{
		assert (j >= 0) ;
		//lower y layer
		index = 4 ;
	}
	else if (i > P && i < Nx - P && j >= Ny - P)
	{
		assert (j < Ny) ;
		//upper y layer
		index = 5 ;
	}
	else if (i <= P && j <= P)
	{
		assert (i >= 0 && j >= 0) ; 
		//lower x layer, lower y layer
		index = 6 ;
	}
	else if (i <= P && j >= Ny - P)
	{
		assert (i >= 0 && j < Ny) ;
		//lower x layer, upper y layer
		index = 7 ;
	}
	else if (i >= Nx - P && j <= P)
	{
		assert (i < Nx && j >= 0) ;
		//upper x layer, lower y layer
		index = 8 ;
	}
	else if (i >= Nx - P && j >= Ny - P)
	{
		assert (i < Nx && j < Ny) ;
		//upper x layer, upper y layer
		index = 9 ;
	}
	return index ;
}

//The following are simple clones for 1D.
template <int DIM>
class pochoir_clone
{
public :
	pochoir_clone()
    {}

	virtual void operator() (int t, int i) const
	{
		cout << "base clone boundary" << endl ;
	}

	virtual void operator() (int t, int i, int j) const
	{
		cout << "base clone boundary 2" << endl ;
	}

	virtual void operator() (int t0, int t1, grid_info<DIM> const & grid) const
	{
		cout << "base clone interior" << endl ;
	}
};

template <int DIM>
class pochoir_clone_0 : public pochoir_clone<DIM>
{
	Pochoir_Array <double, DIM> & a ;
public:
	pochoir_clone_0(Pochoir_Array <double, DIM> & array):a(array)
	{}

	virtual void operator() (int t, int i)const
	{
		cout << "clone 0 boundary" << endl ;
	}

	virtual void operator() (int t, int i, int j)const
	{
		cout << "clone 0 boundary 2" << endl ;
	}

	virtual void operator() (int t0, int t1, grid_info<DIM> const & grid)const
	{
		cout << "clone 0 interior" << endl ;
	}
};

//template <int 1>
class pochoir_clone_1 : public pochoir_clone<1>
{
	Pochoir_Array <double, 1> & a ;
public:
	pochoir_clone_1(Pochoir_Array <double, 1> & array):a(array)
	{}

	virtual void operator() (int t, int i)const
	{
#define a(t, i) a.boundary(t, i)
		//cout << "clone 1  boundary" << endl ;
		a(t + 1, i) = 0.13 * (a(t, i + 1) - 2.0 * a(t, i) + 
							a(t, i - 1));
		//cout << "0.13 t " << t+1 << " i " << i << endl ;
#undef a(t, i)
	}

	virtual void operator() (int t0, int t1, grid_info<1> const & grid)const
	{
		//cout << "clone 1  interior" << endl ;
		grid_info<1> l_grid = grid;
		double * iter3;
		double * iter2;
		double * iter1;
		double * iter0;
		
		double * a_base = a.data();
		const int l_a_total_size = a.total_size();
		
		int gap_a_0;
		const int l_stride_a_0 = a.stride(0);

		for (int t = t0; t < t1; ++t) { 
		double * baseIter_1;
		double * baseIter_0;
		baseIter_0 = a_base + ((t + 1) & 0x1) * l_a_total_size + (l_grid.x0[0]) * l_stride_a_0;
		baseIter_1 = a_base + ((t) & 0x1) * l_a_total_size + (l_grid.x0[0]) * l_stride_a_0;
		iter0 = baseIter_0 + (0) * l_stride_a_0;
		iter1 = baseIter_1 + (1) * l_stride_a_0;
		iter2 = baseIter_1 + (0) * l_stride_a_0;
		iter3 = baseIter_1 + (-1) * l_stride_a_0;
		
		
		#pragma ivdep
		for (int i = l_grid.x0[0]; i < l_grid.x1[0]; ++i, 
		++iter0, 
		++iter1, 
		++iter2, 
		++iter3) {
			(*iter0) = 0.13 * ((*iter1) - 2.0 * (*iter2) + (*iter3));
		} /* end for (sub-trapezoid) */ 
		/* Adjust sub-trapezoid! */
		for (int i = 0; i < 1; ++i) {
			l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
		}
		} /* end for t */
	}
} ;

//template <int 1>
class pochoir_clone_2 : public pochoir_clone<1>
{
	Pochoir_Array <double, 1> & a ;
public:
	pochoir_clone_2(Pochoir_Array <double, 1> & array):a(array)
	{}

	virtual void operator() (int t, int i)const
	{
#define a(t, i) a.boundary(t, i)
		//cout << "clone 2  interior" << endl ;
		a(t + 1, i) = 0.125 * (a(t, i + 1) - 2.0 * a(t, i) + 
							a(t, i - 1));
		//cout << "0.125 t " << t+1 << " i " << i << endl ;
#undef a(t, i)
	}

	virtual void operator() (int t0, int t1, grid_info<1> const & grid)const
	{
		//cout << "clone 2 boundary" << endl ;
		grid_info<1> l_grid = grid;
		double * iter3;
		double * iter2;
		double * iter1;
		double * iter0;
		
		double * a_base = a.data();
		const int l_a_total_size = a.total_size();
		
		int gap_a_0;
		const int l_stride_a_0 = a.stride(0);

		for (int t = t0; t < t1; ++t) { 
		double * baseIter_1;
		double * baseIter_0;
		baseIter_0 = a_base + ((t + 1) & 0x1) * l_a_total_size + (l_grid.x0[0]) * l_stride_a_0;
		baseIter_1 = a_base + ((t) & 0x1) * l_a_total_size + (l_grid.x0[0]) * l_stride_a_0;
		iter0 = baseIter_0 + (0) * l_stride_a_0;
		iter1 = baseIter_1 + (1) * l_stride_a_0;
		iter2 = baseIter_1 + (0) * l_stride_a_0;
		iter3 = baseIter_1 + (-1) * l_stride_a_0;
		
		
		#pragma ivdep
		for (int i = l_grid.x0[0]; i < l_grid.x1[0]; ++i, 
		++iter0, 
		++iter1, 
		++iter2, 
		++iter3) {
			(*iter0) = 0.125 * ((*iter1) - 2.0 * (*iter2) + (*iter3));
		} /* end for (sub-trapezoid) */ 
		/* Adjust sub-trapezoid! */
		for (int i = 0; i < 1; ++i) {
			l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
		}
		} /* end for t */
	}
};

class pochoir_clone_2d_1 : public pochoir_clone<2>
{
	Pochoir_Array <double, 2> & a ;
public:
	pochoir_clone_2d_1(Pochoir_Array <double, 2> & array):a(array)
	{}

	#define a(t, i, j) a.boundary(t, i, j)
	virtual void operator() (int t, int i, int j) const
	{
		/*if (t < 2 * predicate<2>::t0_)
		{
			a(t, i, j) = 0.135 * (a(t - 1, i + 1, j) - 2.0 * a(t - 1, i, j) + a(t - 1, i - 1, j)) + 0.125 * (a(t - 1, i, j + 1) - 2.0 * a(t - 1, i, j) + a(t - 1, i, j - 1)) + a(t - 1, i, j);
		}
		else
		{*/
			a(t, i, j) = 0.125 * (a(t - 1, i + 1, j) - 2.0 * a(t - 1, i, j) + a(t - 1, i - 1, j)) + 0.125 * (a(t - 1, i, j + 1) - 2.0 * a(t - 1, i, j) + a(t - 1, i, j - 1)) + a(t - 1, i, j);
		//}
	}
	#undef a(t, i, j)

	virtual void operator() (int t0, int t1, grid_info<2> const & grid)const
	{
	grid_info<2> l_grid = grid;
	double * iter5;
	double * iter4;
	double * iter3;
	double * iter2;
	double * iter1;
	double * iter0;
	
	double * a_base = a.data();
	const int l_a_total_size = a.total_size();
	
	int gap_a_1, gap_a_0;
	const int l_stride_a_1 = a.stride(1), l_stride_a_0 = a.stride(0);

	for (int t = t0; t < t1; ++t) { 
	double * baseIter_1;
	double * baseIter_0;
	baseIter_0 = a_base + ((t) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
	baseIter_1 = a_base + ((t - 1) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
	iter0 = baseIter_0 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
	iter1 = baseIter_1 + (1) * l_stride_a_1 + (0) * l_stride_a_0;
	iter2 = baseIter_1 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
	iter3 = baseIter_1 + (-1) * l_stride_a_1 + (0) * l_stride_a_0;
	iter4 = baseIter_1 + (0) * l_stride_a_1 + (1) * l_stride_a_0;
	iter5 = baseIter_1 + (0) * l_stride_a_1 + (-1) * l_stride_a_0;
	
	gap_a_1 = l_stride_a_1 + (l_grid.x0[0] - l_grid.x1[0]) * l_stride_a_0;
	for (int i = l_grid.x0[1]; i < l_grid.x1[1]; ++i, 
	iter0 += gap_a_1, 
	iter1 += gap_a_1, 
	iter2 += gap_a_1, 
	iter3 += gap_a_1, 
	iter4 += gap_a_1, 
	iter5 += gap_a_1) {
	
	#pragma ivdep
	for (int j = l_grid.x0[0]; j < l_grid.x1[0]; ++j, 
	++iter0, 
	++iter1, 
	++iter2, 
	++iter3, 
	++iter4, 
	++iter5) {
		/*if (t0 < 2 * predicate<2>::t0_ ) {
			(*iter0) = 0.135 * ((*iter1) - 2.0 * (*iter2) + (*iter3)) + 0.125 * ((*iter4) - 2.0 * (*iter2) + (*iter5)) + (*iter2);
		}
		else 
		{*/
			(*iter0) = 0.125 * ((*iter1) - 2.0 * (*iter2) + (*iter3)) + 0.125 * ((*iter4) - 2.0 * (*iter2) + (*iter5)) + (*iter2);
		//}
	} } /* end for (sub-trapezoid) */ 
	/* Adjust sub-trapezoid! */
	for (int i = 0; i < 2; ++i) {
		l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
	}
	} /* end for t */
	}
} ;

class pochoir_clone_2d_2 : public pochoir_clone<2>
{
	Pochoir_Array <double, 2> & a ;
public:
	pochoir_clone_2d_2(Pochoir_Array <double, 2> & array):a(array)
	{}

	virtual void operator() (int t, int i, int j) const
	{
#define a(t, i, j) a.boundary(t, i, j)
		//cout << "2d - clone 2  boundary" << endl ;
		a(t, i, j) = 0.125 * (a(t - 1, i + 1, j) - 2.0 * a(t - 1, i, j) + a(t - 1, i - 1, j)) + 0.125 * (a(t - 1, i, j + 1) - 2.0 * a(t - 1, i, j) + a(t - 1, i, j - 1)) + a(t - 1, i, j);
		//cout << "0.126 t " << t << " i " << i << " j " << j << endl ;
#undef a(t, i, j)
	}

	virtual void operator() (int t0, int t1, grid_info<2> const & grid)const
	{
		//cout << "2d - clone 2  interior" << endl ;
		grid_info<2> l_grid = grid;
		double * iter5;
		double * iter4;
		double * iter3;
		double * iter2;
		double * iter1;
		double * iter0;
		
		double * a_base = a.data();
		const int l_a_total_size = a.total_size();
		
		int gap_a_1, gap_a_0;
		const int l_stride_a_1 = a.stride(1), l_stride_a_0 = a.stride(0);

		for (int t = t0; t < t1; ++t) { 
		double * baseIter_1;
		double * baseIter_0;
		baseIter_0 = a_base + ((t) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		baseIter_1 = a_base + ((t - 1) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		iter0 = baseIter_0 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter1 = baseIter_1 + (1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter2 = baseIter_1 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter3 = baseIter_1 + (-1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter4 = baseIter_1 + (0) * l_stride_a_1 + (1) * l_stride_a_0;
		iter5 = baseIter_1 + (0) * l_stride_a_1 + (-1) * l_stride_a_0;
		
		gap_a_1 = l_stride_a_1 + (l_grid.x0[0] - l_grid.x1[0]) * l_stride_a_0;
		for (int i = l_grid.x0[1]; i < l_grid.x1[1]; ++i, 
		iter0 += gap_a_1, 
		iter1 += gap_a_1, 
		iter2 += gap_a_1, 
		iter3 += gap_a_1, 
		iter4 += gap_a_1, 
		iter5 += gap_a_1) {
		
		#pragma ivdep
		for (int j = l_grid.x0[0]; j < l_grid.x1[0]; ++j, 
		++iter0, 
		++iter1, 
		++iter2, 
		++iter3, 
		++iter4, 
		++iter5) {
		
		(*iter0) = 0.125 * ((*iter1) - 2.0 * (*iter2) + (*iter3)) + 0.125 * ((*iter4) - 2.0 * (*iter2) + (*iter5)) + (*iter2);
		} } /* end for (sub-trapezoid) */ 
		/* Adjust sub-trapezoid! */
		for (int i = 0; i < 2; ++i) {
			l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
		}
		} /* end for t */
	//cout << "2d - clone 2 interior done" << endl ;
	}
} ;

class pochoir_clone_2d_3 : public pochoir_clone<2>
{
	Pochoir_Array <double, 2> & a ;
public:
	pochoir_clone_2d_3(Pochoir_Array <double, 2> & array):a(array)
	{}

	virtual void operator() (int t, int i, int j) const
	{
#define a(t, i, j) a.boundary(t, i, j)
		//cout << "clone 2  interior" << endl ;
		a(t, i, j) = 0.126 * (a(t - 1, i + 1, j) - 2.0 * a(t - 1, i, j) + a(t - 1, i - 1, j)) + 0.125 * (a(t - 1, i, j + 1) - 2.0 * a(t - 1, i, j) + a(t - 1, i, j - 1)) + a(t - 1, i, j);
		//cout << "0.126 t " << t+1 << " i " << i << endl ;
#undef a(t, i, j)
	}

	virtual void operator() (int t0, int t1, grid_info<2> const & grid)const
	{
		grid_info<2> l_grid = grid;
		double * iter5;
		double * iter4;
		double * iter3;
		double * iter2;
		double * iter1;
		double * iter0;
		
		double * a_base = a.data();
		const int l_a_total_size = a.total_size();
		
		int gap_a_1, gap_a_0;
		const int l_stride_a_1 = a.stride(1), l_stride_a_0 = a.stride(0);

		for (int t = t0; t < t1; ++t) { 
		double * baseIter_1;
		double * baseIter_0;
		baseIter_0 = a_base + ((t) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		baseIter_1 = a_base + ((t - 1) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		iter0 = baseIter_0 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter1 = baseIter_1 + (1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter2 = baseIter_1 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter3 = baseIter_1 + (-1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter4 = baseIter_1 + (0) * l_stride_a_1 + (1) * l_stride_a_0;
		iter5 = baseIter_1 + (0) * l_stride_a_1 + (-1) * l_stride_a_0;
		
		gap_a_1 = l_stride_a_1 + (l_grid.x0[0] - l_grid.x1[0]) * l_stride_a_0;
		for (int i = l_grid.x0[1]; i < l_grid.x1[1]; ++i, 
		iter0 += gap_a_1, 
		iter1 += gap_a_1, 
		iter2 += gap_a_1, 
		iter3 += gap_a_1, 
		iter4 += gap_a_1, 
		iter5 += gap_a_1) {
		
		#pragma ivdep
		for (int j = l_grid.x0[0]; j < l_grid.x1[0]; ++j, 
		++iter0, 
		++iter1, 
		++iter2, 
		++iter3, 
		++iter4, 
		++iter5) {
		
		(*iter0) = 0.126 * ((*iter1) - 2.0 * (*iter2) + (*iter3)) + 0.125 * ((*iter4) - 2.0 * (*iter2) + (*iter5)) + (*iter2);
		} } /* end for (sub-trapezoid) */ 
		/* Adjust sub-trapezoid! */
		for (int i = 0; i < 2; ++i) {
			l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
		}
		} /* end for t */
	}
} ;

class pochoir_clone_2d_4 : public pochoir_clone<2>
{
	Pochoir_Array <double, 2> & a ;
public:
	pochoir_clone_2d_4(Pochoir_Array <double, 2> & array):a(array)
	{}

	virtual void operator() (int t, int i, int j) const
	{
#define a(t, i, j) a.boundary(t, i, j)
		//cout << "2d clone 4  boundary" << endl ;
		a(t, i, j) = 0.127 * (a(t - 1, i + 1, j) - 2.0 * a(t - 1, i, j) + a(t - 1, i - 1, j)) + 0.125 * (a(t - 1, i, j + 1) - 2.0 * a(t - 1, i, j) + a(t - 1, i, j - 1)) + a(t - 1, i, j);
		//cout << "2d clone 4 boundary done t " << t << " i " << i << " j " << j << endl ;
#undef a(t, i, j)
	}

	virtual void operator() (int t0, int t1, grid_info<2> const & grid)const
	{
		grid_info<2> l_grid = grid;
		double * iter5;
		double * iter4;
		double * iter3;
		double * iter2;
		double * iter1;
		double * iter0;
		
		double * a_base = a.data();
		const int l_a_total_size = a.total_size();
		
		int gap_a_1, gap_a_0;
		const int l_stride_a_1 = a.stride(1), l_stride_a_0 = a.stride(0);

		for (int t = t0; t < t1; ++t) { 
		double * baseIter_1;
		double * baseIter_0;
		baseIter_0 = a_base + ((t) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		baseIter_1 = a_base + ((t - 1) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		iter0 = baseIter_0 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter1 = baseIter_1 + (1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter2 = baseIter_1 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter3 = baseIter_1 + (-1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter4 = baseIter_1 + (0) * l_stride_a_1 + (1) * l_stride_a_0;
		iter5 = baseIter_1 + (0) * l_stride_a_1 + (-1) * l_stride_a_0;
		
		gap_a_1 = l_stride_a_1 + (l_grid.x0[0] - l_grid.x1[0]) * l_stride_a_0;
		for (int i = l_grid.x0[1]; i < l_grid.x1[1]; ++i, 
		iter0 += gap_a_1, 
		iter1 += gap_a_1, 
		iter2 += gap_a_1, 
		iter3 += gap_a_1, 
		iter4 += gap_a_1, 
		iter5 += gap_a_1) {
		
		#pragma ivdep
		for (int j = l_grid.x0[0]; j < l_grid.x1[0]; ++j, 
		++iter0, 
		++iter1, 
		++iter2, 
		++iter3, 
		++iter4, 
		++iter5) {
		
		(*iter0) = 0.127 * ((*iter1) - 2.0 * (*iter2) + (*iter3)) + 0.125 * ((*iter4) - 2.0 * (*iter2) + (*iter5)) + (*iter2);
		} } /* end for (sub-trapezoid) */ 
		/* Adjust sub-trapezoid! */
		for (int i = 0; i < 2; ++i) {
			l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
		}
		} /* end for t */
	}
} ;

class pochoir_clone_2d_5 : public pochoir_clone<2>
{
	Pochoir_Array <double, 2> & a ;
public:
	pochoir_clone_2d_5(Pochoir_Array <double, 2> & array):a(array)
	{}

	virtual void operator() (int t, int i, int j) const
	{
#define a(t, i, j) a.boundary(t, i, j)
		//cout << "2d clone 4  boundary" << endl ;
		a(t, i, j) = 0.128 * (a(t - 1, i + 1, j) - 2.0 * a(t - 1, i, j) + a(t - 1, i - 1, j)) + 0.125 * (a(t - 1, i, j + 1) - 2.0 * a(t - 1, i, j) + a(t - 1, i, j - 1)) + a(t - 1, i, j);
		//cout << "2d clone 4 boundary done t " << t << " i " << i << " j " << j << endl ;
#undef a(t, i, j)
	}

	virtual void operator() (int t0, int t1, grid_info<2> const & grid)const
	{
		grid_info<2> l_grid = grid;
		double * iter5;
		double * iter4;
		double * iter3;
		double * iter2;
		double * iter1;
		double * iter0;
		
		double * a_base = a.data();
		const int l_a_total_size = a.total_size();
		
		int gap_a_1, gap_a_0;
		const int l_stride_a_1 = a.stride(1), l_stride_a_0 = a.stride(0);

		for (int t = t0; t < t1; ++t) { 
		double * baseIter_1;
		double * baseIter_0;
		baseIter_0 = a_base + ((t) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		baseIter_1 = a_base + ((t - 1) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		iter0 = baseIter_0 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter1 = baseIter_1 + (1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter2 = baseIter_1 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter3 = baseIter_1 + (-1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter4 = baseIter_1 + (0) * l_stride_a_1 + (1) * l_stride_a_0;
		iter5 = baseIter_1 + (0) * l_stride_a_1 + (-1) * l_stride_a_0;
		
		gap_a_1 = l_stride_a_1 + (l_grid.x0[0] - l_grid.x1[0]) * l_stride_a_0;
		for (int i = l_grid.x0[1]; i < l_grid.x1[1]; ++i, 
		iter0 += gap_a_1, 
		iter1 += gap_a_1, 
		iter2 += gap_a_1, 
		iter3 += gap_a_1, 
		iter4 += gap_a_1, 
		iter5 += gap_a_1) {
		
		#pragma ivdep
		for (int j = l_grid.x0[0]; j < l_grid.x1[0]; ++j, 
		++iter0, 
		++iter1, 
		++iter2, 
		++iter3, 
		++iter4, 
		++iter5) {
		
		(*iter0) = 0.128 * ((*iter1) - 2.0 * (*iter2) + (*iter3)) + 0.125 * ((*iter4) - 2.0 * (*iter2) + (*iter5)) + (*iter2);
		} } /* end for (sub-trapezoid) */ 
		/* Adjust sub-trapezoid! */
		for (int i = 0; i < 2; ++i) {
			l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
		}
		} /* end for t */
	}
} ;

class pochoir_clone_2d_6 : public pochoir_clone<2>
{
	Pochoir_Array <double, 2> & a ;
public:
	pochoir_clone_2d_6(Pochoir_Array <double, 2> & array):a(array)
	{}

	virtual void operator() (int t, int i, int j) const
	{
#define a(t, i, j) a.boundary(t, i, j)
		//cout << "2d clone 4  boundary" << endl ;
		a(t, i, j) = 0.129 * (a(t - 1, i + 1, j) - 2.0 * a(t - 1, i, j) + a(t - 1, i - 1, j)) + 0.125 * (a(t - 1, i, j + 1) - 2.0 * a(t - 1, i, j) + a(t - 1, i, j - 1)) + a(t - 1, i, j);
		//cout << "2d clone 4 boundary done t " << t << " i " << i << " j " << j << endl ;
#undef a(t, i, j)
	}

	virtual void operator() (int t0, int t1, grid_info<2> const & grid)const
	{
		grid_info<2> l_grid = grid;
		double * iter5;
		double * iter4;
		double * iter3;
		double * iter2;
		double * iter1;
		double * iter0;
		
		double * a_base = a.data();
		const int l_a_total_size = a.total_size();
		
		int gap_a_1, gap_a_0;
		const int l_stride_a_1 = a.stride(1), l_stride_a_0 = a.stride(0);

		for (int t = t0; t < t1; ++t) { 
		double * baseIter_1;
		double * baseIter_0;
		baseIter_0 = a_base + ((t) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		baseIter_1 = a_base + ((t - 1) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		iter0 = baseIter_0 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter1 = baseIter_1 + (1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter2 = baseIter_1 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter3 = baseIter_1 + (-1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter4 = baseIter_1 + (0) * l_stride_a_1 + (1) * l_stride_a_0;
		iter5 = baseIter_1 + (0) * l_stride_a_1 + (-1) * l_stride_a_0;
		
		gap_a_1 = l_stride_a_1 + (l_grid.x0[0] - l_grid.x1[0]) * l_stride_a_0;
		for (int i = l_grid.x0[1]; i < l_grid.x1[1]; ++i, 
		iter0 += gap_a_1, 
		iter1 += gap_a_1, 
		iter2 += gap_a_1, 
		iter3 += gap_a_1, 
		iter4 += gap_a_1, 
		iter5 += gap_a_1) {
		
		#pragma ivdep
		for (int j = l_grid.x0[0]; j < l_grid.x1[0]; ++j, 
		++iter0, 
		++iter1, 
		++iter2, 
		++iter3, 
		++iter4, 
		++iter5) {
		
		(*iter0) = 0.129 * ((*iter1) - 2.0 * (*iter2) + (*iter3)) + 0.125 * ((*iter4) - 2.0 * (*iter2) + (*iter5)) + (*iter2);
		} } /* end for (sub-trapezoid) */ 
		/* Adjust sub-trapezoid! */
		for (int i = 0; i < 2; ++i) {
			l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
		}
		} /* end for t */
	}
} ;

class pochoir_clone_2d_7 : public pochoir_clone<2>
{
	Pochoir_Array <double, 2> & a ;
public:
	pochoir_clone_2d_7(Pochoir_Array <double, 2> & array):a(array)
	{}

	virtual void operator() (int t, int i, int j) const
	{
#define a(t, i, j) a.boundary(t, i, j)
		//cout << "2d clone 4  boundary" << endl ;
		a(t, i, j) = 0.130 * (a(t - 1, i + 1, j) - 2.0 * a(t - 1, i, j) + a(t - 1, i - 1, j)) + 0.125 * (a(t - 1, i, j + 1) - 2.0 * a(t - 1, i, j) + a(t - 1, i, j - 1)) + a(t - 1, i, j);
		//cout << "2d clone 4 boundary done t " << t << " i " << i << " j " << j << endl ;
#undef a(t, i, j)
	}

	virtual void operator() (int t0, int t1, grid_info<2> const & grid)const
	{
		grid_info<2> l_grid = grid;
		double * iter5;
		double * iter4;
		double * iter3;
		double * iter2;
		double * iter1;
		double * iter0;
		
		double * a_base = a.data();
		const int l_a_total_size = a.total_size();
		
		int gap_a_1, gap_a_0;
		const int l_stride_a_1 = a.stride(1), l_stride_a_0 = a.stride(0);

		for (int t = t0; t < t1; ++t) { 
		double * baseIter_1;
		double * baseIter_0;
		baseIter_0 = a_base + ((t) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		baseIter_1 = a_base + ((t - 1) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		iter0 = baseIter_0 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter1 = baseIter_1 + (1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter2 = baseIter_1 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter3 = baseIter_1 + (-1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter4 = baseIter_1 + (0) * l_stride_a_1 + (1) * l_stride_a_0;
		iter5 = baseIter_1 + (0) * l_stride_a_1 + (-1) * l_stride_a_0;
		
		gap_a_1 = l_stride_a_1 + (l_grid.x0[0] - l_grid.x1[0]) * l_stride_a_0;
		for (int i = l_grid.x0[1]; i < l_grid.x1[1]; ++i, 
		iter0 += gap_a_1, 
		iter1 += gap_a_1, 
		iter2 += gap_a_1, 
		iter3 += gap_a_1, 
		iter4 += gap_a_1, 
		iter5 += gap_a_1) {
		
		#pragma ivdep
		for (int j = l_grid.x0[0]; j < l_grid.x1[0]; ++j, 
		++iter0, 
		++iter1, 
		++iter2, 
		++iter3, 
		++iter4, 
		++iter5) {
		
		(*iter0) = 0.130 * ((*iter1) - 2.0 * (*iter2) + (*iter3)) + 0.125 * ((*iter4) - 2.0 * (*iter2) + (*iter5)) + (*iter2);
		} } /* end for (sub-trapezoid) */ 
		/* Adjust sub-trapezoid! */
		for (int i = 0; i < 2; ++i) {
			l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
		}
		} /* end for t */
	}
} ;

class pochoir_clone_2d_8 : public pochoir_clone<2>
{
	Pochoir_Array <double, 2> & a ;
public:
	pochoir_clone_2d_8(Pochoir_Array <double, 2> & array):a(array)
	{}

	virtual void operator() (int t, int i, int j) const
	{
#define a(t, i, j) a.boundary(t, i, j)
		//cout << "2d clone 4  boundary" << endl ;
		a(t, i, j) = 0.131 * (a(t - 1, i + 1, j) - 2.0 * a(t - 1, i, j) + a(t - 1, i - 1, j)) + 0.125 * (a(t - 1, i, j + 1) - 2.0 * a(t - 1, i, j) + a(t - 1, i, j - 1)) + a(t - 1, i, j);
		//cout << "2d clone 4 boundary done t " << t << " i " << i << " j " << j << endl ;
#undef a(t, i, j)
	}

	virtual void operator() (int t0, int t1, grid_info<2> const & grid)const
	{
		grid_info<2> l_grid = grid;
		double * iter5;
		double * iter4;
		double * iter3;
		double * iter2;
		double * iter1;
		double * iter0;
		
		double * a_base = a.data();
		const int l_a_total_size = a.total_size();
		
		int gap_a_1, gap_a_0;
		const int l_stride_a_1 = a.stride(1), l_stride_a_0 = a.stride(0);

		for (int t = t0; t < t1; ++t) { 
		double * baseIter_1;
		double * baseIter_0;
		baseIter_0 = a_base + ((t) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		baseIter_1 = a_base + ((t - 1) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		iter0 = baseIter_0 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter1 = baseIter_1 + (1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter2 = baseIter_1 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter3 = baseIter_1 + (-1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter4 = baseIter_1 + (0) * l_stride_a_1 + (1) * l_stride_a_0;
		iter5 = baseIter_1 + (0) * l_stride_a_1 + (-1) * l_stride_a_0;
		
		gap_a_1 = l_stride_a_1 + (l_grid.x0[0] - l_grid.x1[0]) * l_stride_a_0;
		for (int i = l_grid.x0[1]; i < l_grid.x1[1]; ++i, 
		iter0 += gap_a_1, 
		iter1 += gap_a_1, 
		iter2 += gap_a_1, 
		iter3 += gap_a_1, 
		iter4 += gap_a_1, 
		iter5 += gap_a_1) {
		
		#pragma ivdep
		for (int j = l_grid.x0[0]; j < l_grid.x1[0]; ++j, 
		++iter0, 
		++iter1, 
		++iter2, 
		++iter3, 
		++iter4, 
		++iter5) {
		
		(*iter0) = 0.131 * ((*iter1) - 2.0 * (*iter2) + (*iter3)) + 0.125 * ((*iter4) - 2.0 * (*iter2) + (*iter5)) + (*iter2);
		} } /* end for (sub-trapezoid) */ 
		/* Adjust sub-trapezoid! */
		for (int i = 0; i < 2; ++i) {
			l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
		}
		} /* end for t */
	}
} ;

class pochoir_clone_2d_9 : public pochoir_clone<2>
{
	Pochoir_Array <double, 2> & a ;
public:
	pochoir_clone_2d_9(Pochoir_Array <double, 2> & array):a(array)
	{}

	virtual void operator() (int t, int i, int j) const
	{
#define a(t, i, j) a.boundary(t, i, j)
		//cout << "2d clone 4  boundary" << endl ;
		a(t, i, j) = 0.132 * (a(t - 1, i + 1, j) - 2.0 * a(t - 1, i, j) + a(t - 1, i - 1, j)) + 0.125 * (a(t - 1, i, j + 1) - 2.0 * a(t - 1, i, j) + a(t - 1, i, j - 1)) + a(t - 1, i, j);
		//cout << "2d clone 4 boundary done t " << t << " i " << i << " j " << j << endl ;
#undef a(t, i, j)
	}

	virtual void operator() (int t0, int t1, grid_info<2> const & grid)const
	{
		grid_info<2> l_grid = grid;
		double * iter5;
		double * iter4;
		double * iter3;
		double * iter2;
		double * iter1;
		double * iter0;
		
		double * a_base = a.data();
		const int l_a_total_size = a.total_size();
		
		int gap_a_1, gap_a_0;
		const int l_stride_a_1 = a.stride(1), l_stride_a_0 = a.stride(0);

		for (int t = t0; t < t1; ++t) { 
		double * baseIter_1;
		double * baseIter_0;
		baseIter_0 = a_base + ((t) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		baseIter_1 = a_base + ((t - 1) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		iter0 = baseIter_0 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter1 = baseIter_1 + (1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter2 = baseIter_1 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter3 = baseIter_1 + (-1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter4 = baseIter_1 + (0) * l_stride_a_1 + (1) * l_stride_a_0;
		iter5 = baseIter_1 + (0) * l_stride_a_1 + (-1) * l_stride_a_0;
		
		gap_a_1 = l_stride_a_1 + (l_grid.x0[0] - l_grid.x1[0]) * l_stride_a_0;
		for (int i = l_grid.x0[1]; i < l_grid.x1[1]; ++i, 
		iter0 += gap_a_1, 
		iter1 += gap_a_1, 
		iter2 += gap_a_1, 
		iter3 += gap_a_1, 
		iter4 += gap_a_1, 
		iter5 += gap_a_1) {
		
		#pragma ivdep
		for (int j = l_grid.x0[0]; j < l_grid.x1[0]; ++j, 
		++iter0, 
		++iter1, 
		++iter2, 
		++iter3, 
		++iter4, 
		++iter5) {
		
		(*iter0) = 0.132 * ((*iter1) - 2.0 * (*iter2) + (*iter3)) + 0.125 * ((*iter4) - 2.0 * (*iter2) + (*iter5)) + (*iter2);
		} } /* end for (sub-trapezoid) */ 
		/* Adjust sub-trapezoid! */
		for (int i = 0; i < 2; ++i) {
			l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
		}
		} /* end for t */
	}
} ;

class pochoir_clone_2d_10 : public pochoir_clone<2>
{
	Pochoir_Array <double, 2> & a ;
public:
	pochoir_clone_2d_10(Pochoir_Array <double, 2> & array):a(array)
	{}

	virtual void operator() (int t, int i, int j) const
	{
#define a(t, i, j) a.boundary(t, i, j)
		//cout << "2d clone 4  boundary" << endl ;
		a(t, i, j) = 0.133 * (a(t - 1, i + 1, j) - 2.0 * a(t - 1, i, j) + a(t - 1, i - 1, j)) + 0.125 * (a(t - 1, i, j + 1) - 2.0 * a(t - 1, i, j) + a(t - 1, i, j - 1)) + a(t - 1, i, j);
		//cout << "2d clone 4 boundary done t " << t << " i " << i << " j " << j << endl ;
#undef a(t, i, j)
	}

	virtual void operator() (int t0, int t1, grid_info<2> const & grid)const
	{
		grid_info<2> l_grid = grid;
		double * iter5;
		double * iter4;
		double * iter3;
		double * iter2;
		double * iter1;
		double * iter0;
		
		double * a_base = a.data();
		const int l_a_total_size = a.total_size();
		
		int gap_a_1, gap_a_0;
		const int l_stride_a_1 = a.stride(1), l_stride_a_0 = a.stride(0);

		for (int t = t0; t < t1; ++t) { 
		double * baseIter_1;
		double * baseIter_0;
		baseIter_0 = a_base + ((t) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		baseIter_1 = a_base + ((t - 1) & 0x1) * l_a_total_size + (l_grid.x0[1]) * l_stride_a_1 + (l_grid.x0[0]) * l_stride_a_0;
		iter0 = baseIter_0 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter1 = baseIter_1 + (1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter2 = baseIter_1 + (0) * l_stride_a_1 + (0) * l_stride_a_0;
		iter3 = baseIter_1 + (-1) * l_stride_a_1 + (0) * l_stride_a_0;
		iter4 = baseIter_1 + (0) * l_stride_a_1 + (1) * l_stride_a_0;
		iter5 = baseIter_1 + (0) * l_stride_a_1 + (-1) * l_stride_a_0;
		
		gap_a_1 = l_stride_a_1 + (l_grid.x0[0] - l_grid.x1[0]) * l_stride_a_0;
		for (int i = l_grid.x0[1]; i < l_grid.x1[1]; ++i, 
		iter0 += gap_a_1, 
		iter1 += gap_a_1, 
		iter2 += gap_a_1, 
		iter3 += gap_a_1, 
		iter4 += gap_a_1, 
		iter5 += gap_a_1) {
		
		#pragma ivdep
		for (int j = l_grid.x0[0]; j < l_grid.x1[0]; ++j, 
		++iter0, 
		++iter1, 
		++iter2, 
		++iter3, 
		++iter4, 
		++iter5) {
		
		(*iter0) = 0.133 * ((*iter1) - 2.0 * (*iter2) + (*iter3)) + 0.125 * ((*iter4) - 2.0 * (*iter2) + (*iter5)) + (*iter2);
		} } /* end for (sub-trapezoid) */ 
		/* Adjust sub-trapezoid! */
		for (int i = 0; i < 2; ++i) {
			l_grid.x0[i] += l_grid.dx0[i]; l_grid.x1[i] += l_grid.dx1[i];
		}
		} /* end for t */
	}
} ;

template <int DIM>
class pochoir_clone_array
{
	vector<pochoir_clone<DIM> *> clones ;
public:
	pochoir_clone_array(Pochoir_Array <double, DIM> & array) ;
	/*{
		clones.push_back(new pochoir_clone_0<DIM> (array)) ;
		clones.push_back(new pochoir_clone_1<DIM> (array)) ;
		clones.push_back(new pochoir_clone_2<DIM> (array)) ;
	}*/

	pochoir_clone <DIM> & operator [] (int i)
	{
		return *(clones [i]) ;
	}

	~pochoir_clone_array()
	{
		for (int i = 0 ; i < clones.size() ; i++)
		{
			delete clones [i] ;
		} 
	}
};
template<>
pochoir_clone_array<1>::pochoir_clone_array(Pochoir_Array <double, 1> & array)
{
		clones.push_back(new pochoir_clone_0<1> (array)) ;
		clones.push_back(new pochoir_clone_1 (array)) ;
		clones.push_back(new pochoir_clone_2 (array)) ;
}

template<>
pochoir_clone_array<2>::pochoir_clone_array(Pochoir_Array <double, 2> & array)
{
		clones.push_back(new pochoir_clone_0<2> (array)) ;
		clones.push_back(new pochoir_clone_2d_1 (array)) ;
		clones.push_back(new pochoir_clone_2d_2 (array)) ;
		clones.push_back(new pochoir_clone_2d_3 (array)) ;
		clones.push_back(new pochoir_clone_2d_4 (array)) ;
		clones.push_back(new pochoir_clone_2d_5 (array)) ;
		clones.push_back(new pochoir_clone_2d_6 (array)) ;
		clones.push_back(new pochoir_clone_2d_7 (array)) ;
		clones.push_back(new pochoir_clone_2d_8 (array)) ;
		clones.push_back(new pochoir_clone_2d_9 (array)) ;
		clones.push_back(new pochoir_clone_2d_10 (array)) ;
}

#endif
