/*
 * ============================================================================
 *
 *       Filename:  tb_function_object.cpp
 *
 *    Description:  test bench for function object
 *
 *        Version:  1.0
 *        Created:  11/12/2010 11:05:27 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Tang (), yuantang@csail.mit.edu
 *        Company:  CSAIL, M.I.T.
 *
 * ============================================================================
 */

#include <functional>
#include <cstdio>
using namespace std;
int main(void) {
    function<int (int, int)> f2 = [] (int x, int y) { return x + y; };
    printf("result = %d\n", f2(1, 2));
    return 0;
}
