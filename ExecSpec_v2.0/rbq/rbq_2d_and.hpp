/*
 * ============================================================================
 *		 Filename:  rbq_2d_and.hpp
 *    Description:  Range bit query implementation to find the "and"
 *                  of the points in a 2D grid using inclusion/exclusion.
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
#ifndef RBQ_2D_AND_HPP
#define RBQ_2D_AND_HPP

#include "rbq_2d.hpp"

template <class word_type>
class rbq_2d_and : public rbq_2d <word_type>
{
public :
	//constructor
	rbq_2d_and(word_type * row_input, unsigned long n1, unsigned long n2,
		   unsigned int q, unsigned int dxdt, unsigned int dydt) :
			rbq_2d<word_type> (row_input, n1, n2, q, dxdt, dydt)
	{}

	//free internal data structures
	~rbq_2d_and() {}
	
	//operator to compute a + b
	inline virtual word_type op(const word_type & a, const word_type & b)
	{
		return a & b ;
	}
protected :
	virtual void initialize_skew_data() ;

	virtual void query(word_type * C, unsigned long n1, unsigned long n2, 
				unsigned long i, unsigned long j, unsigned long k, 
				unsigned long l, unsigned long N, int * result) ;

	//Used to check the correctness of query results of a rectangular query.
	virtual void query_without_preprocessing(unsigned long i, unsigned long j, 
							unsigned long k, unsigned long l, int * result) ;

	/* Query without preprocessing. Used to check the correctness of query 
	   results */
	inline virtual int query_1d_without_preprocessing(unsigned long begin, 
					unsigned long end, word_type * p_input) ;

	//initialize a variable to a default value.
	virtual void init(word_type * a) {*a = 1 ;}
} ;


template <class word_type>
void rbq_2d_and<word_type>::initialize_skew_data()
{
	word_type * s_input = this->skewed_input ;
	unsigned long num_words_row_skewed = 
				(this->ns + this->WORD_SIZE - 1) / this->WORD_SIZE ; 
	unsigned long num_words = num_words_row_skewed * this->ns * this->q ;
	//initialize the skewed input structure to all ones.
	for (unsigned long k = 0 ; k < num_words ; k++)
	{
		s_input [k] = this->WORD_MAX ;
	}
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

 A bit in result is the bitwise and of the corresponding bits in the query range
*/
template <class word_type>
void rbq_2d_and<word_type>::query(word_type * p_C, unsigned long p_n1,
							unsigned long p_n2, unsigned long i, 
							unsigned long j, unsigned long k, unsigned long l,
							unsigned long p_N, int * result)
{
	assert (i <= j && j <= p_n1) ;
	assert (k <= l && l <= p_n2) ;
	unsigned long size = (j - i + 1) * (l - k + 1) ;
	word_type * c2 = p_C + j * p_n2  ;
	if (i > 0 && k > 0)
	{
		word_type * c1 = p_C + (i - 1) * p_n2 ;
		for (unsigned int p = 0 ; p < this->q ; p++)
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
		for (unsigned int p = 0 ; p < this->q ; p++)
		{
			result [p] = result [p] & (c2 [l] - c1 [l] == size) ;
			c1 += p_N ;
			c2 += p_N ;
		}
	}
	else if (k >0) //i = 0
	{
		for (unsigned int p = 0 ; p < this->q ; p++)
		{
			result [p] = result [p] & (c2 [l] - c2 [k - 1] == size) ;
			c2 += p_N ;
		}
	}
	else //i = 0 and k = 0
	{
		for (unsigned int p = 0 ; p < this->q ; p++)
		{
			result [p] = result [p] & (c2 [l] == size) ;
			c2 += p_N ;
		}
	}
}

//query without preprocessing. 
//Used to check the correctness of query results of a rectangular query.
template <class word_type>
void rbq_2d_and<word_type>::query_without_preprocessing(unsigned long i, 
							unsigned long j, unsigned long k, unsigned long l,
							int * result)
{
	assert (i <= j && j <= this->n1) ;
	assert (k <= l && l <= this->n2) ;
	word_type * in = this->input, * ref = this->input ;
	int num_words_row = (this->n2 + this->WORD_SIZE - 1) / this->WORD_SIZE ;
	unsigned long num_words = num_words_row * this->n1 ;
	for (unsigned int p = 0 ; p < this->q ; p++)
	{
		in = ref + num_words_row * i ; //start of row i
		int answer = 1 ;
		for (unsigned long m = i ; answer && m <= j ; m++)
		{
			for (unsigned long n = k ; answer && n <= l ; n++)
			{
#if 0
				answer = get_bit(in + n / this->WORD_SIZE, n % this->WORD_SIZE);
#else
				answer &= get_bit(in + n / this->WORD_SIZE, n % this->WORD_SIZE);
#endif
			}
			in += num_words_row ;
		}
		*result++ = answer ;
		ref += num_words ;
	}
}

/* Query without preprocessing. Used to check the correctness of query results.
 */
template <class word_type>
inline int rbq_2d_and<word_type>::query_1d_without_preprocessing(
                    unsigned long begin, unsigned long end,
                    word_type * p_input)
{
    long c = begin / this->WORD_SIZE ;
    long d = end / this->WORD_SIZE ;

    int i = begin % this->WORD_SIZE ;
    int j = end % this->WORD_SIZE ;
#ifndef NDEBUG
	std::cout << "begin " << begin << " end " << end << std::endl ;
	print_bits(p_input, begin, end + 1, this->WORD_SIZE) ;
#endif
    if (c + 1 <= d - 1)
    {
        //query spans atleast three blocks
        int answer = 1 ;
        for (long k = c + 1 ; answer && k <= d - 1 ; k++)
        {
            answer = (p_input [k] == this->WORD_MAX) ;
        }
        return answer &&
            ! geneity(i, this->WORD_SIZE - 1, p_input [c], this->WORD_SIZE) &&
            ! geneity(0, j, p_input [d], this->WORD_SIZE) &&
            get_bit(p_input + c, i) == 1 && get_bit(p_input + d, j) == 1 ;
    }
    else
    {
        if (c < d)
        {
            //query spans two adjacent blocks
            return ! geneity(i, this->WORD_SIZE - 1, p_input [c], 
							 this->WORD_SIZE) &&
                ! geneity(0, j, p_input [d], this->WORD_SIZE) &&
                get_bit(p_input + c, i) == 1 && get_bit(p_input + d, j) == 1;
        }
        else
        {
            //query lies within a block
            assert(c == d) ;
            return ! geneity(i, j, p_input [c], this->WORD_SIZE) && 
					get_bit(p_input + c, i) == 1;
        }
    }
}

#endif //RBQ_2D_AND_HPP
