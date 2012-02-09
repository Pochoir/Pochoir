/*
 **********************************************************************************
 *  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 *  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 *                           Charles E. Leiserson <cel@mit.edu>
 *   
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 *   This helper class 'proxy' was originally written by Dahua Lin@csail.mit.edu
 *   adapted to Pochoir by Yuan Tang
 * 
 *
 *   Suggestsions:                  yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 ********************************************************************************/

#ifndef POCHOIR_PROXY_H
#define POCHOIR_PROXY_H

/* if the type T is of a struct, the user has to employ the arrow operator '->'
 * to reference to the member, since the dot operator '.' is non-overload-able in C++
 */
template<typename T>
class Pochoir_Proxy
{
public:
    explicit Pochoir_Proxy(T * v) : val_(*v), ref_(v) { }
    explicit Pochoir_Proxy(T v) : val_(v), ref_(&val_) { }

    operator T() const { // the implicit conversion makes a proxy just like the value itself
	    return (*ref_);
    }
    operator T& () {
        return (*ref_);
    }
    T * operator->() {
        return ref_;
    }

    Pochoir_Proxy<T> & operator= (T const & rhs) {
        (*ref_) = rhs;
        return (*this);
    }

    Pochoir_Proxy<T> & operator= (Pochoir_Proxy<T> & rhs) {
        T const & l_rhs = T(rhs);
        (*ref_) = l_rhs;
        return (*this);
    }
private:
    T val_;
    T * ref_;
};

#endif /* POCHOIR_PROXY_H */
