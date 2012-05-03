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
    explicit Pochoir_Proxy(T * v) : ref_(v) { 
        val_ = *v;
    }
    explicit Pochoir_Proxy(T v) : val_(v) { 
        /* if this constructor is used later in assignment operator, it will cause a
         * segmentation fault!!!
         */
        ref_ = NULL;
    }

    Pochoir_Proxy(Pochoir_Proxy<T> const & rhs) { 
        ref_ = rhs.get_ref(); val_ = rhs.get_val();
    }
    /* the implicit conversion makes a proxy just like the value itself
     */
    operator T () const { return val_; }
    operator T () { return val_; }

    T * operator->() { return ref_; }
    T * get_ref() { return ref_; }
    T get_val() { return val_; }
    Pochoir_Proxy<T> & operator= (T const & rhs) {
        (*ref_) = rhs;
        // val_ = rhs;
        return (*this);
    }

    Pochoir_Proxy<T> & operator= (Pochoir_Proxy<T> & rhs) {
        T const & l_rhs = T(rhs);
        (*ref_) = l_rhs;
        // val_ = l_rhs;
        return (*this);
    }
private:
    T val_;
    T * ref_;
};

#endif /* POCHOIR_PROXY_H */
