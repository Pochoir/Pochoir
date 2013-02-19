/*
 **********************************************************************************
 *  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 *  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 * 		                     Charles E. Leiserson <cel@mit.edu>
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
 *   Suggestsions:                  yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 *********************************************************************************
 */

#ifndef POCHOIR_ERRMSG_H
#define POCHOIR_ERRMSG_H

#include <cstdio>
#include <cstdlib>
#include <cassert>

#define INF 100000000

#define LOGGER

#ifdef LOGGER

#define LOG_LEVEL 1

#define LOG(log_level, fmt) \
if (log_level >= LOG_LEVEL) { \
    fprintf(stderr, "<%s:%d> " fmt "\n", __FUNCTION__, __LINE__); \
}

#define LOG_ARGS(log_level, fmt, ...) \
if (log_level >= LOG_LEVEL) { \
    fprintf(stderr, "<%s:%d> " fmt "\n", __FUNCTION__, __LINE__, __VA_ARGS__); \
}

#else 

#define LOG(log_level, fmt) ((void) 0) 

#define LOG_ARGS(log_level, fmt, ...) ((void) 0)

#endif /* defined LOGGER */

#define ERROR(fmt) \
    do { \
        fprintf(stderr, "<%s:%d> " fmt "\nExit!\n", __FUNCTION__, __LINE__); \
        exit(EXIT_FAILURE); \
    } while(0)

#define ERROR_ARGS(fmt, ...) \
    do { \
        fprintf(stderr, "<%s:%d> " fmt "\nExit!\n", __FUNCTION__, __LINE__, __VA_ARGS__); \
        exit(EXIT_FAILURE); \
    } while(0)

#define WARNING(fmt) \
    do { \
        fprintf(stderr, "<%s:%d> " fmt ".\n", __FUNCTION__, __LINE__); \
    } while(0)


#define WARNING_ARGS(fmt, ...) \
    do { \
        fprintf(stderr, "<%s:%d> " fmt ".\n", __FUNCTION__, __LINE__, __VA_ARGS__); \
    } while(0)

#endif /* POCHOIR_ERRMSG_H */
