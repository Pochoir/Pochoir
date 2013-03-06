/*
 * ============================================================================
 *	     Filename: rbq_common.h
 *    Description: General helper functions in rbq implementation
 *	      Created: 4/26/2012
 *         Author: Eka Palamadai, epn@mit.edu
 * ============================================================================
 */
#ifndef RBQ_COMMON_H
#define RBQ_COMMON_H

#include <stdio.h>

//helper functions
/* returns 0 if word[begin, end] is homogeneous
   and 1 otherwise */
template <class word_type>
inline unsigned int geneity(int begin, int end, word_type word, int word_size)
{
    assert (begin <= end && begin < word_size && end < word_size) ;
    unsigned int offset = word_size - (end - begin + 1) ;
    word = word >> begin ;
    return ((word_type) (word << offset) | 0 &&
            (word_type) (~word << offset) | 0) ;
}

template <class word_type>
inline void set_bit(word_type * word, int pos, int value)
{
    *word = (*word & ~(1ul << pos)) | ((word_type)(value != 0) << pos);
}

template <class word_type>
inline int get_bit(word_type * word, int pos)
{
    return 1 & (*word >> pos);
}

template <class word_type>
inline void print_bits(word_type * bits, word_type begin, word_type end, 
						int word_size)
{
    word_type i ;
    for (i = begin ; i < end ; i++)
    {
        printf(" %d", get_bit(bits + i / word_size, i % word_size)) ;
    }
    printf("\n") ;
}

template <class word_type>
inline void merge_subwords(word_type * word1, word_type * word2, int pos,
                        int word_size, word_type * result)
{
    // extract word1[pos, word_size - 1]  and word2[0, pos - 1]
    *result = *word1 >> pos | *word2 << (word_size - pos) ;
}

template <class word_type>
inline word_type get_subword(int begin, int end, word_type word, int word_size)
{
	//extract the subword word[begin, end]
    assert (begin <= end && begin < word_size && end < word_size) ;
    unsigned int offset = word_size - (end - begin + 1) ;
    word = word >> begin ;
    return (word_type) (word << offset) ;
}

template <class word_type>
inline void copy_subword(word_type * source, word_type * dest, int src_begin,
                    int src_end, int dest_begin, int dest_end,
                    int word_size)
{
    //copy source[src_begin...src_end] to dest[dest_begin...dest_end]
    assert(src_end - src_begin = dest_begin - dest_end) ;
    int offset = src_begin + word_size - src_end - 1 ;
    word_type t = *source >> src_begin ;
    t = t << offset ;
    t = t >> (offset - dest_begin) ;
    *dest = *dest | t ;
}

#endif //RBQ_COMMON_H
