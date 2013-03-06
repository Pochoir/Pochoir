/*
 * ============================================================================
 *		 Filename:  rbq_2d_simple.hpp
 *    Description:  Range bit query implementation to find the "and"
 *                  of the points in a 2D grid using inclusion/exclusion.
 *                  Preprocessing time: Theta(Nq)
 *                  where N is the number of points in the 2D grid
 *                  q is the number of bits per point
 *                  Storage: Theta(Nq) words
 *                  Query time : Theta(1). 
 *					To answer a 2D octagonal range, 4 parallelograms
 *					and 2 rectangles are queried, and results are combined.
 *		  Created:  4/26/2012
 *
 *         Author:  Eka Palamadai, epn@mit.edu
 * ============================================================================
 */
#ifndef RBQ_2D_SIMPLE_HPP
#define RBQ_2D_SIMPLE_HPP

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
	//returns 1 if successfully initialized, 0 otherwise
	rbq_2d(word_type * row_input, unsigned long n1, unsigned long n2,
		   unsigned long q, unsigned int dxdt, unsigned int dydt) ;

	//returns 1 if successfully preprocessed, 0 otherwise
	int preprocess() ;


	//query without preprocessing. Used to check the correctness of query 
	//results of a octagonal query.
	void query_without_preprocessing(unsigned long x1, 
			unsigned long x2, unsigned long x3, unsigned long x4,
			unsigned long y1, unsigned long y2, unsigned long y3, 
			unsigned long y4, int * result) ;


	//A bit in result is the and of all bits in the query range
	void query(unsigned long x1, unsigned long x2, 
				unsigned long x3, unsigned long x4, 
				unsigned long y1, unsigned long y2,
				unsigned long y3, unsigned long y4, int * result) ;

	//free internal data structures
	~rbq_2d() ;
	
	static const word_type WORD_MAX = ULONG_MAX ;
	static const int WORD_SIZE = sizeof(word_type) * 8 ;
	static const int W_SROOT ;
private :
	void print_data() ;

	void compute_bit_count(word_type * C, word_type * input, unsigned long n1, 
						   unsigned long n2) ;
	void skew_input() ;

	void query(word_type * C, unsigned long n1, unsigned long n2, 
				unsigned long i, unsigned long j, unsigned long k, 
				unsigned long l, unsigned long N, int * result) ;

	//Used to check the correctness of query results of a rectangular query.
	void query_without_preprocessing(unsigned long i, unsigned long j, 
							unsigned long k, unsigned long l, int * result) ;

	/* Used to check the correctness of query results */
	inline void query_without_preprocessing(unsigned long x1,
                unsigned long x2, unsigned long x3, unsigned long x4, 
				unsigned long y1, unsigned long y2, unsigned long y3, 
				unsigned long y4, word_type * data_input, 
				unsigned long num_words, int * result) ;

	/* Query without preprocessing. Used to check the correctness of query 
	   results */
	inline int query_1d_without_preprocessing(unsigned long begin, 
					unsigned long end, word_type * p_input) ;

	//pointer to input array in row major order.
	word_type * input ; //size n1 * n2 * q bits
	unsigned long n1 ; //# of rows 
	unsigned long n2 ; //# of cols
	unsigned long q ; //number of bits per input point.
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
						unsigned long p_n2, unsigned long p_q,
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
	unsigned long num_words = num_words_per_bit * q ;
	unsigned long num_words_row = (n2 + WORD_SIZE - 1) / WORD_SIZE ;

	//initialize the skewed input structure to all ones.
	for (unsigned long k = 0 ; k < num_words ; k++)
	{
		s_input [k] = WORD_MAX ;
	}
	//copy the input into skewed input
	for (int k = 0 ; k < q ; k++)
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
	for (int k = 0 ; k < q ; k++)
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
	for (int i = 0 ; i < q ; i++)
	{
		result [i] = 1 ; 
	}
	//query region ([y1, y4] * [x2, x3])
	query(C, n1, n2, y1, y4, x2, x3, N, result) ;
	
	//query region ([y2, y3] * [x1, x4])
	query(C, n1, n2, y2, y3, x1, x4, N, result) ;

	//calculate the offsets for the parallelogram queries
	unsigned long yoffset = y3 * dxdt ;
	unsigned long yoffset2 = y2 * dxdt ;
	unsigned long xoffset = (n2 - 1 - x1) * dydt ;
	unsigned long xoffset2 = x1 * dydt ;
	unsigned long xm = 2 * x2 - x1 ;
	unsigned long xoffset3 = (n2 - 1 - xm) * dydt ;
	unsigned long xoffset4 = xm * dydt ;
	unsigned long xoffset5 = (n2 - 1 - x4) * dydt ;
	unsigned long xoffset6 = x4 * dydt ;
	unsigned long xn = 2 * x3 - x4 ;
	unsigned long xoffset7 = (n2 - 1 - xn) * dydt ;
	unsigned long xoffset8 = xn * dydt ;

	//query the four parallelograms
	//parallelogram 1
	unsigned long a2 = xoffset + yoffset ;
	unsigned long b1 = xoffset2 + yoffset ;
	unsigned long a1 = xoffset3 + yoffset ;
	unsigned long b2 = xoffset4 + yoffset ;
	query(Cs, ns, ns, a1, a2, b1, b2, Ns, result) ;
	
	//parallelogram 2
	a1 = xoffset5 + yoffset ;
	b2 = xoffset6 + yoffset ;
	a2 = xoffset7 + yoffset ;
	b1 = xoffset8 + yoffset ;
	query(Cs, ns, ns, a1, a2, b1, b2, Ns, result) ;
	
	//parallelogram 3
	a2 = xoffset + yoffset2 ;
	b1 = xoffset2 + yoffset2 ;
	a1 = xoffset3 + yoffset2 ;
	b2 = xoffset4 + yoffset2 ;
	query(Cs, ns, ns, a1, a2, b1, b2, Ns, result);
	
	//parallelogram 4
	a1 = xoffset5 + yoffset2 ;
	b2 = xoffset6 + yoffset2 ;
	a2 = xoffset7 + yoffset2 ;
	b1 = xoffset8 + yoffset2 ;
	query(Cs, ns, ns, a1, a2, b1, b2, Ns, result);
}





/* 
   --------------------------
  |						     |
j |		(j,k)	    (j,l)	 |
  |						     |
  |						     |
i |		(i,k)		(i,l)    |
  |						     |
   --------------------------
         k          l
*/
//A bit in result is set to 1 if homogeneous and 0 otherwise
template <class word_type>
void rbq_2d<word_type>::query(word_type * p_C, unsigned long p_n1,
							unsigned long p_n2, unsigned long i, 
							unsigned long j, unsigned long k, unsigned long l,
							unsigned long p_N, int * result)
{
	assert (i <= j && j <= p_n1) ;
	assert (k <= l && l <= p_n2) ;
	unsigned long size = (j - i + 1) * (l - k + 1) ;
	word_type * c2 = p_C + j * p_n2  ;
	//std::cout << "x1 " << k << " x2 " << l << " y1 " << i << " y2 " << j << 
	//			std::endl ;
	if (i > 0 && k > 0)
	{
		word_type * c1 = p_C + (i - 1) * p_n2 ;
		for (int p = 0 ; p < q ; p++)
		{
			result [p] = result [p] & 
					  (c2 [l] - c2 [k - 1] - c1 [l] + c1 [k - 1] == size) ;
			c1 += p_N ;
			c2 += p_N ;
		}
	}
	else if (i > 0)  //k = 0
	{
		word_type * c1 = p_C + (i - 1) * p_n2 ;
		for (int p = 0 ; p < q ; p++)
		{
			result [p] = result [p] & (c2 [l] - c1 [l] == size) ;
			c1 += p_N ;
			c2 += p_N ;
		}
	}
	else if (k >0) //i = 0
	{
		for (int p = 0 ; p < q ; p++)
		{
			result [p] = result [p] & (c2 [l] - c2 [k - 1] == size) ;
			c2 += p_N ;
		}
	}
	else //i = 0 and k = 0
	{
		for (int p = 0 ; p < q ; p++)
		{
			result [p] = result [p] & (c2 [l] == size) ;
			c2 += p_N ;
		}
	}
	//std::cout << "result [0] " << result [0] << std::endl ;
	//std::cout << "result [1] " << result [1] << std::endl ;
}

//query without preprocessing. 
//Used to check the correctness of query results of a rectangular query.
template <class word_type>
void rbq_2d<word_type>::query_without_preprocessing(unsigned long i, 
							unsigned long j, unsigned long k, unsigned long l,
							int * result)
{
	assert (i <= j && j <= n1) ;
	assert (k <= l && l <= n2) ;
	word_type * in = input, * ref = input ;
	int num_words_row = (n2 + WORD_SIZE - 1) / WORD_SIZE ;
	unsigned long num_words = num_words_row * n1 ;
	for (int p = 0 ; p < q ; p++)
	{
		in = ref + num_words_row * i ; //start of row i
		int answer = 1 ;
		for (unsigned long m = i ; answer && m <= j ; m++)
		{
			for (unsigned long n = k ; answer && n <= l ; n++)
			{
				answer = get_bit(in + n / WORD_SIZE, n % WORD_SIZE) ;
			}
			in += num_words_row ;
		}
		*result++ = answer ;
		ref += num_words ;
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
    unsigned long bit, i ;
    for (bit = 0 ; bit < q ; bit++)
    {
        int answer = 1 ;
        if (y2 > y1)
        {
            //region 1
            double m1 = -1. * (double) (x2 - x1) /  (y2 - y1) ;
            double m2 = (double) (x4 - x3) /  (y2 - y1) ;
            word_type * in = data_input + y1 * num_words ;
            for (i = y1 ; answer && i < y2 ; i++)
            {
                unsigned long begin = (long) ceil(m1 * (i - y1) + x2) ;
                unsigned long end = (long) floor(m2 * (i - y1) + x3) ;
                //make a 1d query
                answer = query_1d_without_preprocessing(begin, end, in) ;
                in += num_words ;
            }
        }
        //region 2
        word_type * in = data_input + y2 * num_words ;
        for (i = y2 ; answer && i < y3 ; i++)
        {
            answer = query_1d_without_preprocessing(x1, x4, in) ;
            in += num_words ;
        }

        if (y4 > y3)
        {
            //region 3
            double m3 = (double) (x2 - x1) / (y4 - y3) ;
            double m4 = -1. * (double) (x4 - x3) / (y4 - y3) ;
            word_type * in = data_input + y3 * num_words ;
            for (i = y3 ; answer && i <= y4 ; i++)
            {
                unsigned long begin = (long) ceil(m3 * (i - y3) + x1) ;
                unsigned long end = (long) floor(m4 * (i - y3) + x4) ;
                //make a 1d query
                answer = query_1d_without_preprocessing(begin, end, in) ;
                in += num_words ;
            }
        }
        *result++ = answer ;
		data_input += n1 * num_words ;
    }
}


/* Query without preprocessing. Used to check the correctness of query results.
 */
template <class word_type>
inline int rbq_2d<word_type>::query_1d_without_preprocessing(
                    unsigned long begin, unsigned long end,
                    word_type * p_input)
{
    long c = begin / WORD_SIZE ;
    long d = end / WORD_SIZE ;

    int i = begin % WORD_SIZE ;
    int j = end % WORD_SIZE ;
#ifndef NDEBUG
	std::cout << "begin " << begin << " end " << end << std::endl ;
	print_bits(p_input, begin, end, WORD_SIZE) ;
#endif
    if (c + 1 <= d - 1)
    {
        //query spans atleast three blocks
        int answer = 1 ;
        unsigned long k ;
        for (k = c + 1 ; answer && k <= d - 1 ; k++)
        {
            answer = (p_input [k] == WORD_MAX) ;
        }
        return answer &&
            ! geneity(i, WORD_SIZE - 1, p_input [c], WORD_SIZE) &&
            ! geneity(0, j, p_input [d], WORD_SIZE) &&
            get_bit(p_input + c, i) == 1 && get_bit(p_input + d, j) == 1 ;
    }
    else
    {
        if (c < d)
        {
            //query spans two adjacent blocks
            return ! geneity(i, WORD_SIZE - 1, p_input [c], WORD_SIZE) &&
                ! geneity(0, j, p_input [d], WORD_SIZE) &&
                get_bit(p_input + c, i) == 1 && get_bit(p_input + d, j) == 1;
        }
        else
        {
            //query lies within a block
            assert(c == d) ;
            return ! geneity(i, j, p_input [c], WORD_SIZE) && 
					get_bit(p_input + c, i) == 1;
        }
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
	for (int k = 0 ; k < q ; k++)
	{
		for (unsigned long i = 0 ; i < n1 ; i++)
		{
			print_bits(in, 0ul, n2, WORD_SIZE) ;
			in += num_words_row ;
		}
	}
	word_type * c = C ;
	std::cout << "Count " << std::endl ;
	for (int k = 0 ; k < q ; k++)
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
	for (int k = 0 ; k < q ; k++)
	{
		for (unsigned long i = 0 ; i < ns ; i++)
		{
			print_bits(s_in, 0ul, ns, WORD_SIZE) ;
			s_in += num_words_row_skewed ;
		}
	}
	word_type * cs = Cs ;
	std::cout << "Skewed Count " << std::endl ;
	for (int k = 0 ; k < q ; k++)
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
#endif //RBQ_2D_SIMPLE_H
