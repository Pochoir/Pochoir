/*
 * ============================================================================
 *       Filename: test_rbq_2d.cpp
 *    Description: Test driver for the 2D range bit query 
 *		  Created: 4/26/2012	
 *         Author: Eka Palamadai, epn@mit.edu
 * ============================================================================
 */

#include "rbq_2d_simple.hpp"
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>

typedef unsigned long word_type ;
#define NUM_QUERIES 100

/* This routine simulates a divide and conquer on a 2D grid
   invoking a rectangular query at each step. 
   This is only simulating a single time step.
   To do : Modify to add time cuts in addition to space cuts */
void simulate(unsigned long x1, unsigned long x2, unsigned long y1, 
			 unsigned long y2, rbq_2d<word_type> * rbq, unsigned long q,
			 double * result)
{
	unsigned long dx = x2 - x1 + 1 ;
	unsigned long dy = y2 - y1 + 1 ;
    clockmark_t time1, time2 ;
	int query_result [q] ;
	int check_result [q] ;
	double qtime = 0., q2time = 0. ;
	unsigned long j, k ;
	for (k = 0 ; k < NUM_QUERIES ; k++)
	{
		time1 = ktiming_getmark( );
		rbq->query(x1, x1, x2, x2, y1, y1, y2, y2, query_result) ;
		time2 = ktiming_getmark( );
		qtime += ktiming_diff_sec( &time1, &time2 );

		int error = 0 ;
		time1 = ktiming_getmark( );
		rbq->query_without_preprocessing(x1, x1, x2, x2, y1, y1, y2, y2, 
										check_result) ;
		time2 = ktiming_getmark( );
		q2time += ktiming_diff_sec( &time1, &time2 );

		//validate results for correctness
		for (j = 0 ; j < q ; j++)
		{
			if (query_result [j] != check_result [j])
			{
				error = 1 ;
				std::cout << "bit " << j << " result :" << query_result [j] << 
				 "  result no preprocess: " << check_result [j] << std::endl ;
			}
		}
		
		if (error)
		{
			printf("Error. Results are different!\n") ;
			break ;
		}
	}
	//printf("%lu\t %lu\t %lu\t %lu\t %lu\t, %g\t, %g\t\n", 
	//		x1, x2, y1, y2, dx * dy, qtime, q2time) ;
	printf("%lu\t\t %g\t\t\n", 
			dx * dy, qtime) ;
	
	if (dx * dy <= 128)
	{
		int j = dx * dy ;
		int i ;
		for(i = 0 ; i < j ; i++)
		{
			double a = sqrt(i * qtime);
			double b = sqrt(i * q2time) ;
			*result = a * b ;
		}
	}
	else if (dx > dy)
	{
		//cut along x
		unsigned long mid = x1 + (x2 - x1) / 2 ;
		simulate(x1, mid, y1, y2, rbq, q, result) ;
		simulate(mid + 1, x2, y1, y2, rbq, q, result) ;
	}
	else
	{
		//cut along y
		unsigned long mid = y1 + (y2 - y1) / 2 ;
		simulate(x1, x2, y1, mid, rbq, q, result) ;
		simulate(x1, x2, mid + 1, y2, rbq, q, result) ;
	}
}

int main(int argc, char *argv[])
{
    unsigned long n1 = 0 ;
    unsigned long n2 = 0 ;
    unsigned long N = 0 ;
    unsigned long q = 1 ;
	unsigned long i, j, k ;
	unsigned long max_q = 1 ;

	char * row_input_bits = 0 ; 
	//test input for testing a small case.
	/*char * test_input = "11111" \
						"11111" \
						"11111" \
						"11111" \
						"11111" \
						"11111" \
						"11111" \
						"11111" \
						"00000" \
						"00000" ; */
	char * test_input = "1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"1111111111111111111111" \
						"0000111100000000001100" \
						"0000111100000000001100" \
						"0000111100000000001100" \
						"0000111100000000001100" \
						"0000000000000000000000" \
						"0000111100000000000000" \
						"0000111100000000001100" \
						"0000111100000000001100" \
						"0000111100000000001100" \
						"0000111100000000001100" \
						"0000000000000000000000" \
						"0000111100000000000000" \
						"0000111100000000001100" \
						"0000111100000000001100" \
						"0000111100000000001100" \
						"0000111100000000001100" \
						"0000000000000000000000" \
						"0000111100000000000000" \
						"0000111100000000001100" \
						"0000111100000000001100" \
						"0000111100000000001100" \
						"0000111100000000001100" \
						"0000000000000000000000" \
						"0000111100000000000000" ; 
	char c ;
	while ((c = getopt(argc, argv, "q:r:m:n:u")) != EOF) 
	{
    	switch (c) 
		{
			case 'm': 
        		n1 = atol(optarg);
        		break;
			case 'n': 
        		n2 = atol(optarg);
        		break;
			case 'q': 
        		q = atol(optarg);
        		break;
			case 'r': 
        		max_q = atol(optarg);
        		break;
			case 'u':
				printf("Usage: ./rmq [-m # of rows] [-n # of cols] [-q number of bits at each point] [-r max number of bits at each point] \n") ;
				exit(0) ;
      		default:
				row_input_bits = test_input ;
    	}
	}
	//use test input only when the user has not specified
	//any valid command line input
	if (n1 == 0 || n2 == 0)
	{
		row_input_bits = test_input ;
		n1 = 24 ;
		n2 = 22 ;
		//n1 = 10 ;
		//n2 = 5 ;
		//q = 1 ;
	}
	if (max_q < q)
	{
		max_q = q ;
	}
	N = n1 * n2 ;
	//int num_q = __builtin_ctz (max_q) + 1 ;
    printf("N %lu \n", N) ;
	//W_SROOT = (int) sqrt((double)WORD_SIZE) ;
	unsigned long q1, Q = q, end = max_q - q + 1 ;
    double itime[end], ptime[end], qtime[end], q2time[end] ;
	int w_sroot = rbq_2d<word_type>::W_SROOT ;
	int word_size = rbq_2d<word_type>::WORD_SIZE ;
	for (q1 = 0 ; q1 < end ; q1++)
	{
    	printf("Number of bits %lu \n", q) ;

		unsigned long num_block_rows = (n1 + w_sroot - 1) / w_sroot ;
    	unsigned long num_block_cols = (n2 + w_sroot - 1) / w_sroot ;
    	unsigned long num_words = num_block_rows * num_block_cols ;
    	word_type * block_input = (word_type * )
				malloc(num_words * q * sizeof(word_type)) ;
    	assert(block_input) ;
		unsigned long num_words_row = (n2 + word_size - 1) / word_size ;
	    word_type * row_input = (word_type *)
				malloc(n1 * num_words_row * q * sizeof(word_type)) ;
		assert(row_input) ;
		unsigned long num_words_col = (n1 + word_size - 1) / word_size ;
    	word_type * col_input = (word_type *)
				malloc(num_words_col * n2 * q * sizeof(word_type)) ;
		assert(col_input) ;
		if (block_input == (word_type *) 0 || row_input == (word_type *) 0 ||
			col_input == (word_type *) 0)
		{
			printf("Error allocating memory\n") ;
			free(block_input) ;
			free(row_input) ;
			free(col_input) ;
			return 1 ;
		}
		//Generate Input either automatically or from the sample test input
		if (row_input_bits == (char *) 0)
		{
			printf("generating input..\n") ;
    		srand (time(NULL));
			word_type * in = row_input ;
			for (j = 0 ; j < q ; j++)
			{
				for (i = 0 ; i < n1 ; i++)
				{
					for (k = 0 ; k < num_words_row ; k++)
					{
						//in [k] = (rand() % 2 == 0 ? 0 : WORD_MAX) ;
						in [k] = rbq_2d<word_type>::WORD_MAX ;
					}
					in += num_words_row ;
				}
			}

			//generate column input
			word_type * col_in = col_input ;
			in = row_input ;
            for (j = 0 ; j < q ; j++)
			{
                for (i = 0 ; i < n1 ; i++)
                {
					int word_index = i / word_size ;
					int bit_index = i % word_size ;
					for (k = 0 ; k < n2 ; k++)
					{
						int bit = get_bit(in + k / word_size, k % word_size) ;
                    	set_bit(col_in + k * num_words_col + word_index, 
								bit_index, bit) ;
					}
					in += num_words_row ;
                }
				col_in += num_words_col * n2 ;
            }

			//generate block input
			word_type * block_in = block_input ;
			in = row_input ;
			for (j = 0 ; j < q ; j++)
			{
				for (i = 0 ; i < n1 ; i++)
				{
					int row = i / w_sroot ;
					int row_index = i % w_sroot ;
					block_in = block_input + row * num_block_cols ;
					for (k = 0 ; k < n2 ; k++)
					{
						int col = k / w_sroot ;
						int col_index = k % w_sroot ;
						int bit = get_bit(in + k / word_size, k % word_size) ;
						set_bit(block_in + col, 
								row_index * w_sroot + col_index, bit) ;
					}
					in += num_words_row ;
				}
				block_in += num_block_rows * num_block_cols ;
			}
		}
		else
		{
			//use the test case input
			//generate col input
			char * c = row_input_bits ;
			word_type * col_in = col_input ;
            for (j = 0 ; j < q ; j++)
			{
                for (i = 0 ; i < n1 ; i++)
                {
					int word_index = i / word_size ;
					int bit_index = i % word_size ;
					for (k = 0 ; k < n2 ; k++)
					{
                    	set_bit(col_in + k * num_words_col + word_index, 
								bit_index, (*c++ == '1')) ;
					}
                }
				col_in += num_words_col * n2 ;
            }

			//generate row input
			word_type * row_in = row_input ;
			c = row_input_bits ;
            for (j = 0 ; j < q ; j++)
            {
                for (i = 0 ; i < n1 ; i++)
                {
					for (k = 0 ; k < n2 ; k++)
					{
                    	set_bit(row_in + k / word_size, k % word_size,
                        	    (*c++ == '1')) ;
					}
                	row_in += num_words_row ;
                }
            }
			
			//generate block input
			c = row_input_bits ;
			word_type * in = block_input ;
			for (j = 0 ; j < q ; j++)
			{
				for (i = 0 ; i < n1 ; i++)
				{
					int row = i / w_sroot ;
					int row_index = i % w_sroot ;
					in = block_input + row * num_block_cols ;
					for (k = 0 ; k < n2 ; k++)
					{
						int col = k / w_sroot ;
						int col_index = k % w_sroot ;
						set_bit(in + col, row_index * w_sroot + col_index, 
								(*c++ == '1')) ;
					}
				}
				in += num_block_rows * num_block_cols ;
			}
		}
		//print input for small sizes
		/*
		if (N < 600)
		{
			printf("\nblock_input : \n") ;
			word_type * in = block_input ;
			for (i = 0 ; i < q ; i++)
			{
				for (j = 0 ; j < num_block_rows ; j++)
				{
					printf("row %lu\n", j) ;
					for (k = 0 ; k < num_block_cols ; k++)
					{
    					print_bits(in + k, 0, (word_type) word_size, word_size) ;
					}
					in += num_block_cols ;
				}
				printf("\n") ;
			}
			for (i = 0 ; i < num_words * q ; i++)
			{
				printf("Input [%lu] = %lu\n", i, block_input[i]) ;
			}
			
			printf("col input...\n") ;
			word_type * col_in = col_input ;
            for (j = 0 ; j < q ; j++)
			{
                for (i = 0 ; i < n2 ; i++)
                {
					printf("col %lu\n", i) ;
					for (k = 0 ; k < num_words_col ; k++)
					{
						print_bits(col_in + k, 0, (word_type) word_size, 
									word_size) ;
					}
                	col_in += num_words_col ;
                }
            }
			printf("row input...\n") ;
			word_type * row_in = row_input ;
            for (j = 0 ; j < q ; j++)
			{
                for (i = 0 ; i < n1 ; i++)
                {
					printf("row %lu\n", i) ;
					for (k = 0 ; k < num_words_row ; k++)
					{
						print_bits(row_in + k, 0, (word_type) word_size,
									word_size) ;
					}
                	row_in += num_words_row ;
					printf("\n") ;
                }
            } 
		}*/	
    	clockmark_t time1, time2 ;

		itime [q1] = 0. ;
    	time1 = ktiming_getmark( );
		rbq_2d<word_type> * rbq = new rbq_2d<word_type>(row_input, 
														n1, n2, q, 1, 1) ;	
    	time2 = ktiming_getmark( );
    	itime [q1] = ktiming_diff_sec( &time1, &time2 );
	
		ptime [q1] = 0. ;
    	time1 = ktiming_getmark( );
		int result = rbq->preprocess() ;
    	time2 = ktiming_getmark( );
    	ptime [q1] = ktiming_diff_sec( &time1, &time2 );
		if (! result)
		{
			printf("Error preprocessing\n") ;
			return 1 ;
		}

		double temp ;
		printf("Number of queries  %d\n\n", NUM_QUERIES) ;
		printf("Number of\t  Query \t \n" \ 
				"points\t\t  time, s \t \n") ;
		simulate(0, n2 - 1, 0, n1 - 1, rbq, q, &temp) ;
		free(block_input);
		free(col_input);
		free(row_input);
		delete rbq ;
		q++ ;
	}

	q = Q ;
	printf("\nq = [") ;
	for (q1 = 0 ; q1 < end ; q1++)
	{
		printf("%lu ", q++) ;
	}
	printf ("]\n") ;

	printf("initializion_time, s = [") ;
	for (q1 = 0 ; q1 < end ; q1++)
	{
		printf("%g, ", itime [q1]) ;
	}
	printf ("]\n") ;

	printf("preprocessing_time, s = [") ;
	for (q1 = 0 ; q1 < end ; q1++)
	{
		printf("%g, ", ptime [q1]) ;
	}
	printf ("]\n") ;
}
