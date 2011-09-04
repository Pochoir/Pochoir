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
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;
template <int N_RANK> 
struct Pochoir_Kernel {
    typedef function<void (int, int)> T;
};
template <> 
struct Pochoir_Kernel<1> { 
    typedef function<void (int, int)> T;
};
template <> 
struct Pochoir_Kernel<2> {
    typedef function<void (int, int, int)> T;
};
template <>
struct Pochoir_Kernel<3> {
    typedef function<void (int, int, int, int)> T;
};

int main(void) {
//    typedef function<int (int, int)> T;
    typedef function<void (int, int)> TP;
//    Pochoir_Kernel<1>::T f2;
    TP f1;
    Pochoir_Kernel<1>::T f2;
    f2 = [] (int x, int y) { printf("In f2, x+y = %d\n",  x + y); };
    f1 = [] (int x, int y) { printf("In f1, x+y = %d\n", x+y); };
    // printf("result = %d\n", f2(1, 2));
    f1(3,4);
    f2(4, 5);
    return 0;
}
