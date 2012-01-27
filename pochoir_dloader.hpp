/*
 **********************************************************************************
 *  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 *  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 *                           Charles E. Leiserson <cel@mit.edu>
 *   
 *   This program is delete software: you can redistribute it and/or modify
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
 *   Suggestsions:                  yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 *********************************************************************************
 */

#ifndef POCHOIR_DLOADER_H
#define POCHOIR_DLOADER_H

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <functional>
#include <dlfcn.h>

class DynamicLoader
{
public:
    DynamicLoader(const char * filename)
    {
        char cmd[200];
        char cpp_filename[100], so_filename[100];

        sprintf(cpp_filename, "%s.cpp", filename);
        sprintf(so_filename, "%s.so", filename);
        sprintf(cmd, "icpc -o %s -shared -nostartfiles -fPIC -O0 -g -std=c++0x -I${POCHOIR_LIB_PATH} %s\0", so_filename, cpp_filename);

        printf("%s\n", cmd);
        int ret = system(cmd);
        if (ret == -1) {
            fprintf(stderr, "system() call failed!\n");
            exit(EXIT_FAILURE);
        }

        m_handle = dlopen(so_filename, RTLD_LAZY);
        if (!m_handle) {
            fprintf(stderr, "can't load library named %s\n", so_filename);
        }
        dlerror();
    }

    template<class T>
    std::function<T> load(const char * functionName)
    {
        void* result = dlsym(m_handle, functionName);
        if ((error = dlerror()) != NULL) {
            fprintf(stderr, "can't find symbol named %s\n", functionName);
            exit(EXIT_FAILURE);
        }
        return reinterpret_cast<T*>(result);
    }

    void close(void) { dlclose(m_handle); }
private:
    void * m_handle;
    char * error;
};

template <typename T>
std::function<T> dloader(void * handle, const char * func_name) {
    void * result = dlsym(handle, func_name);
    char * error;

    if ((error = dlerror()) != NULL) {
        fprintf(stderr, "can't find symbol named %s\n", func_name);
        exit(EXIT_FAILURE);
    }
    return reinterpret_cast<T*>(result);
};

#endif
