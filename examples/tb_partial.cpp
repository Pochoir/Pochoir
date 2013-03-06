/*
 * ============================================================================
 *
 *       Filename:  tb_partial.cpp
 *
 *    Description:  test bench for C++ template partial specialization
 *
 *        Version:  1.0
 *        Created:  08/20/2011 10:56:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Tang (), yuantang@csail.mit.edu
 *        Company:  CSAIL, M.I.T.
 *
 * ============================================================================
 */

#include <iostream>
using namespace std;

template <class T, class U, int I> struct X {
    void f() {
        cout << "Primary template" << endl;
    }
};

template <class T, int I> struct X<T, T*, I> {
    void f() {
        cout << "Partial specialization 1" << endl;
    }
};

template <class T, class U, int I> struct X<T*, U, I> {
    void f() {
        cout << "Partial specialization 2" << endl;
    }
};

template <class T> struct X<int, T*, 10> {
    void f() {
        cout << "Partial specialization 3" << endl;
    }
};

template <class T, class U, int I> struct X<T, U*, I> {
    void f() {
        cout << "Partial specialization 4" << endl;
    }
};

int main(void) {
    X<int, int, 10> a;
    X<int, int*, 5> b;
    X<int*, float, 10> c;
    X<int, char*, 10> d;
    X<float, int*, 10> e;
    a.f(); b.f(); c.f(); d.f(); e.f();
}
