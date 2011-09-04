/*
 * ============================================================================
 *
 *       Filename:  tb_lcm.cpp
 *
 *    Description:  test bench for lowest common multiple
 *
 *        Version:  1.0
 *        Created:  07/10/2011 09:31:50 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Tang (), yuantang@csail.mit.edu
 *        Company:  CSAIL, M.I.T.
 *
 * ============================================================================
 */

#include <cstdio>
#include <cstddef>
#include <pochoir.hpp>

int main(int argc, char * argv[])
{
    int a = 0, b = 0;

    if (argc < 3) {
        printf("argc < 3, quit! \n");
        exit(1);
    }
    a = StrToInt(argv[1]);
    b = StrToInt(argv[2]);
    printf("a = %d, b = %d, gcd(%d, %d) = %d, lcm(%d, %d) = %d\n", 
            a, b, a, b, gcd(a, b), a, b, lcm(a, b));
    return 0;
}


