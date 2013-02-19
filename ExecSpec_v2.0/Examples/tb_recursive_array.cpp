/*
 * ============================================================================
 *
 *       Filename:  tb_recursive_array.cpp
 *
 *    Description:  test bench for recursive array in C++
 *
 *        Version:  1.0
 *        Created:  08/22/2011 10:45:41 PM
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
#include <iostream>

#if 1
template <int N_RANK>
struct Pochoir_Tile {
    __TILE__<N_RANK-1, N_RANK> tile_[];
    int size_;
    friend std::ostream & operator << (std::ostream & os, Pochoir_Tile<N_RANK> const & tile) {
        int i;
        os << "{";
        for (i = 0; i < tile.size_ - 1; ++i) {
            os << tile.tile_[i];
            os << ", ";
        }
        os << tile.tile_[i];
        os << "}";
        return os;
    }
};

template <int N1, int N2>
struct __TILE__ {
    __TILE__<N1-1, N2> tile_[];
    int size_;
    friend std::ostream & operator << (std::ostream & os, __TILE__<N1, N2> const & tile) {
        int i;
        os << "{";
        for (i = 0; i < tile.size_ - 1; ++i) {
            os << tile.tile_[i];
            os << ", ";
        }
        os << tile.tile_[i];
        os << "}";
        return os;
    }
};

template <int N2>
struct __TILE__<0, N2> {
    int tile_[];
    int size_;
    friend std::ostream & operator << (std::ostream & os, __TILE__<0, N2> const & tile) {
        int i;
        os << "{";
        for (int i = 0; i < tile.size_ - 1; ++i) {
            os << tile.tile_[i];
            os << ", ";
        }
        os << tile.tile_[i];
        os << "}";
        return os;
    }
};
#endif

int main(void) {
    int tile[2][2] = {{0, 1}, {2, 3}};
    int tile1[] = {1,2, 3,4};
    std::cout << tile[1][0] << std::endl;
    std::cout << tile1 << std::endl;
    return 0;
}
