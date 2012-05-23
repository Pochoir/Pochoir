/**
 Copyright (c) 2010 by 6.172 Staff <6.172-staff@mit.edu>

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
**/


/**
 * Linux kernel-assisted timing library -- provides high-precision time
 * measurements for the execution time of your algorithms.
 *
 * It also provides timing functionality on other platforms such as Cygwin and
 * Darwin for portability, but the meaning of the time reported may vary.  For
 * example, on Linux we try to measure CPU time used.  On Linux, timing
 * sleep(1) will report close to zero time elapsed, while on Darwin and Cygwin
 * it will report the wall time, which is about 1 second.
 *
 * You shouldn't need to modify this file. More importantly, you should not
 * depend on any modifications you make here, as we will replace it with a
 * fresh copy when we test your code.
 **/

// We need _POSIX_C_SOURCES to pick up 'struct timespec' and clock_gettime.
#define _POSIX_C_SOURCE 200112L

#include "ktiming.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include "CoreServices/CoreServices.h"
#include "mach/mach.h"
#include "mach/mach_time.h"
#else
#include <time.h>
#endif

#ifdef __CYGWIN__
#define KTIMING_CLOCK_ID CLOCK_REALTIME
#else
//#define KTIMING_CLOCK_ID CLOCK_MONOTONIC
#define KTIMING_CLOCK_ID CLOCK_PROCESS_CPUTIME_ID
//#define KTIMING_CLOCK_ID CLOCK_THREAD_CPUTIME_ID
//#define KTIMING_CLOCK_ID CLOCK_REALTIME
#endif

clockmark_t ktiming_getmark(void)
{
#ifdef __APPLE__
  uint64_t abs_time = mach_absolute_time();
  Nanoseconds nanos = AbsoluteToNanoseconds(*(AbsoluteTime*)&abs_time);
  return *(uint64_t*)&nanos;
#else
  struct timespec temp;
  uint64_t nanos;

  /* Try the default clock first, if it fails,
     try again with CLOCK_MONOTONIC */
  int stat = clock_gettime(KTIMING_CLOCK_ID, &temp);
  if (stat != 0) {
#ifndef __CYGWIN__
    stat = clock_gettime(CLOCK_MONOTONIC , &temp);
#endif
    if (stat != 0) {
      perror("ktiming_getmark()");
      exit(-1);
    }
  }
  nanos = temp.tv_nsec;
  nanos += ((uint64_t)temp.tv_sec) * 1000 * 1000 * 1000;
  return nanos;
#endif
}

uint64_t
ktiming_diff_nanosec(const clockmark_t* const start, const clockmark_t* const end)
{
  return *end - *start;
}

float
ktiming_diff_sec(const clockmark_t* const start, const clockmark_t* const end)
{
  return (float)ktiming_diff_nanosec(start, end) / 1000000000.0f;
}
