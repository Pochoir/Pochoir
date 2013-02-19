/*
 * ============================================================================
 *
 *       Filename:  tb_tmp.cpp
 *
 *    Description:  misc test bench 
 *
 *        Version:  1.0
 *        Created:  06/15/2011 09:11:22 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Tang (), yuantang@csail.mit.edu
 *        Company:  CSAIL, M.I.T.
 *
 * ============================================================================
 */
#include <cstdio>
#include <cstdlib>

int main (void) {
    int a = 5, b = 2;
    printf("%d/%d = %d, %d/%d = %d\n", a, b, a/b, a-1, b, (a-1)/b);
    return 0;
}

