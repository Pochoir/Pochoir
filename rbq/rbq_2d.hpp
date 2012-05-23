/*
 * ============================================================================
 *		 Filename:  rbq_2d.hpp
 *    Description:  Range bit query interface for a 2D grid
 *                  using inclusion/exclusion.
 *                  Preprocessing time: Theta(Nq)
 *                  where N is the number of points in the 2D grid
 *                  q is the number of bits per point
 *                  Storage: Theta(Nq) words
 *                  Query time : Theta(q). 
 *					To answer a 2D octagonal range, 4 parallelograms
 *					and 2 rectangles are queried, and results are combined.
 *		  Created:  4/26/2012
 *
 *         Author:  Eka Palamadai, epn@mit.edu
 * ============================================================================
 */
#ifndef RBQ_2D_HPP
#define RBQ_2D_HPP

#include <cassert>
#include <ctime>
#include <climits>
#include <cmath>
#include <new>
#include <iostream>
#include "ktiming.h"
#include "rbq_common.h"
using namespace std ;

template <class word_type>
class rbq_2d
{
public :
	//constructor
	rbq_2d(word_type * row_input, unsigned long n1, unsigned long n2,
		   unsigned int q, unsigned int dxdt, unsigned int dydt) ;

	//returns 1 if successfully preprocessed, 0 otherwise
	int preprocess() ;


	//query without preprocessing. Used to check the correctness of query 
	//results of a octagonal query.
	void query_without_preprocessing(unsigned long x1, 
			unsigned long x2, unsigned long x3, unsigned long x4,
			unsigned long y1, unsigned long y2, unsigned long y3, 
			unsigned long y4, int * result) ;


	/*A bit in result is the result of applying the 'and/or' 
	  operation on all bits in the query range */
	void query(unsigned long x1, unsigned long x2, 
				unsigned long x3, unsigned long x4, 
				unsigned long y1, unsigned long y2,
				unsigned long y3, unsigned long y4, int * result) ;

	//free internal data structures
	virtual ~rbq_2d() ;
	
	//operator to compute a + b
	virtual inline word_type op(const word_type & a, const word_type & b) = 0 ; 

	static const word_type WORD_MAX = ULONG_MAX ;
	static const int WORD_SIZE = sizeof(word_type) * 8 ;
	static const int W_SROOT ;
protected :
	void print_data() ;

	void compute_bit_count(word_type * C, word_type * input, unsigned long n1, 
						   unsigned long n2) ;
	void skew_input() ;

	virtual void initialize_skew_data() = 0 ;

	virtual void query(word_type * C, unsigned long n1, unsigned long n2, 
				unsigned long i, unsigned long j, unsigned long k, 
				unsigned long l, unsigned long N, int * result) = 0 ;

	//Used to check the correctness of query results of a rectangular query.
	virtual void query_without_preprocessing(unsigned long i, unsigned long j, 
							unsigned long k, unsigned long l, int * result) = 0;

	/* Used to check the correctness of query results */
	inline void query_without_preprocessing(unsigned long x1,
                unsigned long x2, unsigned long x3, unsigned long x4, 
				unsigned long y1, unsigned long y2, unsigned long y3, 
				unsigned long y4, word_type * data_input, 
				unsigned long num_words, int * result) ;

	/* Query without preprocessing. Used to check the correctness of query 
	   results */
	virtual int query_1d_without_preprocessing(unsigned long begin, 
					unsigned long end, word_type * p_input) = 0 ;

	//initialize a variable to a default value.
	virtual void init(word_type * a) = 0 ;


	//pointer to input array in row major order.
	word_type * input ; //size n1 * n2 * q bits
	unsigned long n1 ; //# of rows 
	unsigned long n2 ; //# of cols
	unsigned int q ; //number of bits per input point.
	unsigned long N ; //volume of the grid, N = n1 * n2
	unsigned int dxdt ; //slope w.r.t dimension x
	unsigned int dydt ; //slope w.r.t dimension y
	unsigned long ns ; //# of rows/columns in the skewed input
	unsigned long Ns ; //volume of skewed grid. Ns = ns * ns
	//pointer to the skewed input array.
	word_type * skewed_input ; //size ns * ns * q bits
	//pointer to the bit count array of input
	word_type * C ; //size n1 * n2 * q words
	//pointer to the bit count array of skewed input
	word_type * Cs ; //size ns * ns * q words
} ;

template <class word_type>
const int rbq_2d<word_type>::W_SROOT = sqrt(WORD_SIZE) ;

template <class word_type>
rbq_2d<word_type>::rbq_2d(word_type * p_input, unsigned long p_n1, 
						unsigned long p_n2, unsigned int p_q,
						unsigned int p_dxdt, unsigned int p_dydt)
{
    assert(p_input) ;
    assert(p_n1 > 0) ;
    assert(p_n2 > 0) ;
    assert(p_q > 0) ;
	assert(p_dxdt >= 1) ;
	assert(p_dydt >= 1) ;
	C = 0 ;
    unsigned long num_words = p_n1 * p_n2 ;
	ns = p_dydt * (p_n2 - 1) + p_dxdt * (p_n1 - 1) + 1 ; 
	unsigned long num_words_skewed_input = (ns + WORD_SIZE - 1) / WORD_SIZE * 
											ns ; 
	try
	{
    	C = new word_type [num_words * p_q] ;
		skewed_input = new word_type [num_words_skewed_input * p_q] ;
		Cs = new word_type [ns * ns * p_q] ;
	}
	catch (bad_alloc)
	{
		delete [] C ;
		delete [] Cs ;
		delete [] skewed_input ;
		C = 0 ;
		Cs = 0 ;
		skewed_input = 0 ;
		ns = 0 ;
		throw ;
	}
	input = p_input ;
	n1 = p_n1 ;
	n2 = p_n2 ;
	q = p_q ;
	N = n1 * n2 ;
	Ns = ns * ns ;
	dxdt = p_dxdt ;
	dydt = p_dydt ;
}

/* preprocess the input */
template <class word_type>
int rbq_2d<word_type>::preprocess()
{
    assert(input) ;
	assert(C) ;

    //compute the bit count array
	compute_bit_count(C, input, n1, n2) ;
	initialize_skew_data() ;
	skew_input() ;
	compute_bit_count(Cs, skewed_input, ns, ns) ;
#ifndef NDEBUG
	print_data() ;
#endif
    return 1 ;
}


template <class word_type>
void rbq_2d<word_type>::skew_input()
{
	word_type * s_input = skewed_input ;
	word_type * in = input ;
	unsigned long num_words_row_skewed = (ns + WORD_SIZE - 1) / WORD_SIZE ; 
	unsigned long num_words_per_bit = num_words_row_skewed * ns ;
	unsigned long num_words_row = (n2 + WORD_SIZE - 1) / WORD_SIZE ;

	//copy the input into skewed input
	for (unsigned int k = 0 ; k < q ; k++)
	{
		for (unsigned long i = 0 ; i < n1 ; i++)
		{
			unsigned long offset = i * dxdt ;
			for (unsigned long j = 0 ; j < n2 ; j++)
			{
				//Find the row index and column index in the skewed grid.
				unsigned long row = (n2 - 1 - j) * dydt + offset ;
				unsigned long col = j * dydt + offset ;
				set_bit(s_input + row * num_words_row_skewed + col / WORD_SIZE, 
						col % WORD_SIZE, 
						get_bit(in + j / WORD_SIZE, j % WORD_SIZE)) ;
			}
			in += num_words_row ;
		}
		s_input += num_words_per_bit ;
	}
}

template <class word_type>
void rbq_2d<word_type>::compute_bit_count(word_type * p_C, word_type * p_input,
										 unsigned long p_n1, unsigned long p_n2)
{
	assert (p_C) ;
	assert (p_input) ;
	word_type * c_curr = p_C, 
			  * in = p_input ;
	int num_words_row = (p_n2 + WORD_SIZE - 1) / WORD_SIZE ;
	for (unsigned int k = 0 ; k < q ; k++)
	{
		//initialize row 0
		unsigned long r = 0 ;
		for (unsigned long j = 0 ; j < p_n2 ; j++)
		{
			r += get_bit(in + j / WORD_SIZE, j % WORD_SIZE) ;
			c_curr [j] = r ;
		}
		in += num_words_row ;
		word_type * c_prev = c_curr ;
		c_curr += p_n2 ; 
		//count bits from rows 1:(p_n1 - 1)
		for (unsigned long i = 1 ; i < p_n1 ; i++)
		{
			unsigned long p = 0 ;
			for (unsigned long j = 0 ; j < p_n2 ; j++)
			{
				//get the bit count on current row
				p += get_bit(in + j / WORD_SIZE, j % WORD_SIZE) ;
				//C [i, j, k] = C [i - 1, j, k] + p
				c_curr [j] = c_prev [j] + p ;
			}
			in += num_words_row ;
			c_prev = c_curr ;
			c_curr += p_n2 ;
		}
	}
}


/*The query co-ordinates correspond to those of a octagon

y4       --------------
       /               \
      /     3           \
     /___________________\
y3  |                     |
    |       2             |
    |_____________________|
y2   \                   /
      \     1           /
       \               /
y1       -------------
    x1  x2           x3   x4
*/

//A bit in result is set to 1 if homogeneous and 0 otherwise
template <class word_type>
void rbq_2d<word_type>::query(unsigned long x1, unsigned long x2, 
							unsigned long x3, unsigned long x4, 
							unsigned long y1, unsigned long y2,
							unsigned long y3, unsigned long y4, int * result)
{
	assert (result) ;
    assert (x1 <= x2 && x2 <= x3 && x3 <=x4) ;
    assert (y1 <= y2 && y2 <= y3 && y3 <=y4) ;
    assert (x4 < n2) ;
    assert (y4 < n1) ;
	//initialize the result to one
	for (unsigned int i = 0 ; i < q ; i++)
	{
		init((word_type *) (result + i)); 
	}
	//query region ([y1, y4] * [x2, x3])
	query(C, n1, n2, y1, y4, x2, x3, N, result) ;
	
	//query region ([y2, y3] * [x1, x4])
	query(C, n1, n2, y2, y3, x1, x4, N, result) ;

	unsigned long a1, a2, b1, b2 ;
	//query the four parallelograms
	//parallelogram 1 at top left
	if (x1 != x2 && y3 != y4)
	{
		a2 = (n2 - 1 - x1) * dydt + y3 * dxdt ;
		b1 = x1 * dydt + y3 * dxdt ;
		a1 = (n2 - 1 - x2) * dydt + y3 * dxdt ;
		b2 = x2 * dydt + y4 * dxdt ;
		query(Cs, ns, ns, a1, a2, b1, b2, Ns, result) ;
	}
	//parallelogram 2 at top right
	if (y3 != y4 && x3 != x4)
	{
		a1 = (n2 - 1 - x4) * dydt + y3 * dxdt ;
		b2 = x4 * dydt + y3 * dxdt ;
		a2 = (n2 - 1 - x3) * dydt + y4 * dxdt ;
		b1 = x3 * dydt + y3 * dxdt ;
		query(Cs, ns, ns, a1, a2, b1, b2, Ns, result) ;
	}

	//parallelogram 3 at bottom left
	if (y1 != y2 && x1 != x2)
	{
		a2 = (n2 - 1 - x1) * dydt + y2 * dxdt ;
		b1 = x1 * dydt + y2 * dxdt ;
		a1 = (n2 - 1 - x2) * dydt + y1 * dxdt ;
		b2 = x2 * dydt + y2 * dxdt ;
		query(Cs, ns, ns, a1, a2, b1, b2, Ns, result);
	}

	//parallelogram 4 at bottom right
	if (y1 != y2 && x3 != x4)
	{
		a1 = (n2 - 1 - x4) * dydt + y2 * dxdt ;
		b2 = x4 * dydt + y2 * dxdt ;
		a2 = (n2 - 1 - x3) * dydt + y2 * dxdt ;
		b1 = x3 * dydt + y1 * dxdt ;
		query(Cs, ns, ns, a1, a2, b1, b2, Ns, result);
	}
}

/* Free the memory */
template <class word_type>
rbq_2d<word_type>::~rbq_2d()
{
	assert(input) ;
	assert(C) ;
	assert(Cs) ;
	delete [] C ;
	delete [] Cs ;
	delete [] skewed_input ;
	input = 0 ;
}


/* The interface to check the correctness of query results */
template <class word_type>
void rbq_2d<word_type>::query_without_preprocessing(unsigned long x1, 
			unsigned long x2, unsigned long x3, unsigned long x4,
			unsigned long y1, unsigned long y2, unsigned long y3, 
			unsigned long y4, int * result)
{
    assert (result) ;
    assert (x1 <= x2 && x2 <= x3 && x3 <=x4) ;
    assert (y1 <= y2 && y2 <= y3 && y3 <=y4) ;
    assert (x4 < n2) ;
    assert (y4 < n1) ;
	unsigned long num_words_row = (n2 + WORD_SIZE - 1) / WORD_SIZE ;
    query_without_preprocessing(x1, x2, x3, x4, y1, y2, y3, y4, input,
								num_words_row, result) ;
}


/* This routine is used to check the correctness of query results */
template <class word_type>
inline void rbq_2d<word_type>::query_without_preprocessing(unsigned long x1,
                unsigned long x2, unsigned long x3, unsigned long x4, 
				unsigned long y1, unsigned long y2, unsigned long y3, 
				unsigned long y4, word_type * data_input, 
				unsigned long num_words, int * result)
{
    for (unsigned int bit = 0 ; bit < q ; bit++)
    {
        word_type answer ;
		init(&answer) ;
		double m1 = 0., m2 = 0. ;
        if (y2 > y1)
        {
            m1 = -1. * (double) (x2 - x1) /  (y2 - y1) ;
            m2 = (double) (x4 - x3) /  (y2 - y1) ;
		}
        //region 1
        word_type * in = data_input + y1 * num_words ;
        for (unsigned long i = y1 ; i < y2 ; i++)
        {
            unsigned long begin = (long) ceil(m1 * (i - y1) + x2) ;
            unsigned long end = (long) floor(m2 * (i - y1) + x3) ;
            //make a 1d query
            answer = op(answer, 
				   (word_type) query_1d_without_preprocessing(begin, end, in)) ;
            in += num_words ;
        }
        //region 2
        in = data_input + y2 * num_words ;
        for (unsigned long i = y2 ; i < y3 ; i++)
        {
            answer = op(answer, 
					(word_type) query_1d_without_preprocessing(x1, x4, in)) ;
            in += num_words ;
        }
		m1 = 0. ;
		m2 = 0. ;
        if (y4 > y3)
        {
            m1 = (double) (x2 - x1) / (y4 - y3) ;
            m2 = -1. * (double) (x4 - x3) / (y4 - y3) ;
		}
        //region 3
        in = data_input + y3 * num_words ;
        for (unsigned long i = y3 ; i <= y4 ; i++)
        {
            unsigned long begin = (long) ceil(m1 * (i - y3) + x1) ;
            unsigned long end = (long) floor(m2 * (i - y3) + x4) ;
            //make a 1d query
            answer = op(answer, 
				(word_type) query_1d_without_preprocessing(begin, end, in)) ;
            in += num_words ;
        }
        *result++ = answer ;
		data_input += n1 * num_words ;
    }
}


template <class word_type>
void rbq_2d<word_type>::print_data()
{
	std::cout << "n1 : " << n1 << " n2 : " << n2 << " ns : " << 
				ns << std::endl ;
	word_type * in = input ;
	unsigned long num_words_row = (n2 + WORD_SIZE - 1) / WORD_SIZE ; 
	std::cout << "Input " << std::endl ;
	for (unsigned int k = 0 ; k < q ; k++)
	{
		for (unsigned long i = 0 ; i < n1 ; i++)
		{
			std::cout << "row : " << i <<std::endl ;
			print_bits(in, 0ul, n2, WORD_SIZE) ;
			in += num_words_row ;
		}
	}
	word_type * c = C ;
	std::cout << "Count " << std::endl ;
	for (unsigned int k = 0 ; k < q ; k++)
	{
		for (unsigned long i = 0 ; i < n1 ; i++)
		{
			for (unsigned long j = 0 ; j < n2 ; j++)
			{
				std::cout << c [j] << " " ;
			}
			std::cout << std::endl ;
			c += n2 ;
		}
	}
	word_type * s_in = skewed_input ;
	unsigned long num_words_row_skewed = (ns + WORD_SIZE - 1) / WORD_SIZE ; 
	std::cout << "Skewed Input " << std::endl ;
	for (unsigned int k = 0 ; k < q ; k++)
	{
		for (unsigned long i = 0 ; i < ns ; i++)
		{
			print_bits(s_in, 0ul, ns, WORD_SIZE) ;
			s_in += num_words_row_skewed ;
		}
	}
	word_type * cs = Cs ;
	std::cout << "Skewed Count " << std::endl ;
	for (unsigned int k = 0 ; k < q ; k++)
	{
		for (unsigned long i = 0 ; i < ns ; i++)
		{
			for (unsigned long j = 0 ; j < ns ; j++)
			{
				std::cout << cs [j] << " " ;
			}
			std::cout << std::endl ;
			cs += ns ;
		}
	}
}
#endif //RBQ_2D_HPP
