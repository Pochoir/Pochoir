/*
 * ============================================================================
 *
 *       Filename:  tb_cpp.cpp
 *
 *    Description:  test template lambda function in c++0x
 *
 *        Version:  1.0
 *        Created:  05/31/2011 09:37:07 AM
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
/* template function can only reside outside any specific function! */
    template <typename T>
    auto func(T x) -> T { return 2*x;}
int main(void) {
    auto func1 = [](int x) -> int { return 2+x;};
    printf("result = %d\n", func1(3));
    return 0;
}

