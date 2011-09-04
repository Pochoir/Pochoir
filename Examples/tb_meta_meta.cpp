/*
 * ============================================================================
 *
 *       Filename:  tb_meta_meta.cpp
 *
 *    Description:  test bench for meta-meta templates
 *
 *        Version:  1.0
 *        Created:  06/28/2011 02:55:11 PM
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
#include <functional>

using namespace std;

template <bool Cond, class THEN, class ELSE>
struct IF {
    template <bool Condition> 
    struct selector {
        typedef THEN SELECT_CLASS;
    };
    typedef typename selector<Cond>::SELECT_CLASS RESULT;
};

template<bool Cond, class THEN, class ELSE> template <>
IF::selector <false>{
    typedef ELSE SELECT_CLASS;
};

struct THEN {
    static int func() {
        printf("inside THEN!\n");
        return 42;
    }
};

struct ELSE {
    static int func() {
        printf("inside ELSE!\n");
        return 3;
    }
};

int main(void) {
    int result = IF<4 == sizeof(int), THEN, ELSE>::RESULT::func();
    printf("returning : %d\n", result);
    return 0;
}
