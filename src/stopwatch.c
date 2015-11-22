//#include "bug.h"
//#include "internal/stopwatch.h"
#include "stopwatch.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_TIMING
#   include <time.h>
#else

/* Indices into perf events. */
#define RAW_COUNT_INDEX 0
#define TIME_ENABLED_INDEX 1
#define TIME_RUNNING_INDEX 2

#endif

inline void stopwatches_setup() {
#ifndef USE_TIMING
    /* (kxing): Initialize system for profiling off hardware counters. */
    int return_value = pfm_initialize();
    if (return_value != PFM_SUCCESS) {
        exit(EXIT_FAILURE);
    }
#endif
}

inline void stopwatches_teardown() {
#ifndef USE_TIMING
    /* (kxing): Shut down the hardware counter measurements for profiling. */
    pfm_terminate();
#endif
}

inline void stopwatch_init(stopwatch *s) {
#ifdef KXING_DEBUG_STOPWATCH
    printf("starting: %p\n", s);
#endif
    s->initialized = true;

#ifndef USE_TIMING
    memset(&(s->perf_attributes), 0, sizeof(struct perf_event_attr));
    s->perf_attributes.type = PERF_TYPE_HARDWARE;
    s->perf_attributes.size = sizeof(struct perf_event_attr);
    s->perf_attributes.config = PERF_COUNT_HW_INSTRUCTIONS;
    s->perf_attributes.disabled = 1;
    s->perf_attributes.exclude_kernel = 1;
    s->perf_attributes.exclude_hv = 1;
    s->perf_attributes.read_format = PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;

    s->file_descriptor = perf_event_open(&(s->perf_attributes), 0, -1, -1, 0);
    if (s->file_descriptor == -1) {
        fprintf(stderr, "Could not initialize stopwatch.\n");
        exit(EXIT_FAILURE);
    }
#endif

    s->running = false;
    s->measurements_valid = false;
	s->num_calls = 0 ;
}

inline void stopwatch_destroy(stopwatch *s) {
#ifndef USE_TIMING
    close(s->file_descriptor);
#endif
}

inline void stopwatch_start(stopwatch *s) {
#ifdef KXING_DEBUG_STOPWATCH
    printf("stopping: %p\n", s);
#endif
    assert(s->initialized);
    assert(!s->running);

    //CILK_ASSERT(s->initialized);
    //CILK_ASSERT(!s->running);

#ifdef USE_TIMING
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &(s->start));
#else
    ioctl(s->file_descriptor, PERF_EVENT_IOC_RESET, 0);
    ioctl(s->file_descriptor, PERF_EVENT_IOC_ENABLE, 0);
#endif
    s->running = true;
    s->measurements_valid = false;
}

inline void stopwatch_stop(stopwatch *s) {
    assert(s->initialized);
    assert(s->running);

    //CILK_ASSERT(s->initialized);
    //CILK_ASSERT(s->running);

#ifdef USE_TIMING
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &(s->end));
    s->elapsed_time =
            (long long)(s->end.tv_sec - s->start.tv_sec) * 1000000000ll +
            s->end.tv_nsec - s->start.tv_nsec;
	assert (s->elapsed_time >= 0) ;
#else
    ioctl(s->file_descriptor, PERF_EVENT_IOC_DISABLE, 0);

    /*
      See:

      http://sourceforge.net/p/perfmon2/libpfm4/ci/master/tree/perf_examples/perf_util.h

      and

      http://web.eece.maine.edu/~vweaver/projects/perf_events/perf_event_open.html

      for details.
    */
    uint64_t values[3];

    read(s->file_descriptor, values, sizeof(values));

    uint64_t raw_count = values[RAW_COUNT_INDEX];
    uint64_t time_enabled = values[TIME_ENABLED_INDEX];
    uint64_t time_running = values[TIME_RUNNING_INDEX];

    if (time_enabled < time_running) {
        assert(false);
        //CILK_ASSERT(false);
    }

    if (time_running == 0) {
        assert(raw_count == 0);
        //CILK_ASSERT(raw_count == 0);
        s->elapsed_time = 0;
    } else {
        /* Scale the values. */
        s->elapsed_time = (long long)((double)(raw_count) * time_enabled / time_running);
    }

    assert(s->elapsed_time >= 0);
    //CILK_ASSERT(s->elapsed_time >= 0);
#endif

    s->running = false;
    s->measurements_valid = true;
}

inline bool stopwatch_is_running(stopwatch *s) {
    return s->running;
}

inline long long stopwatch_get_elapsed_time(stopwatch *s) {
    assert(s->initialized);
    //CILK_ASSERT(s->initialized);
    assert(s->measurements_valid);
    //CILK_ASSERT(s->measurements_valid);
    return s->elapsed_time;
}

inline long long stopwatch_stop_and_start(stopwatch *s) {
#ifdef KXING_DEBUG_STOPWATCH
    printf("stopping: %p\n", s);
    printf("starting: %p\n", s);
#endif

#if defined(USE_TIMING) && defined(KXING_OPTIMIZED_STOPWATCH)
    assert(s->running);
    //CILK_ASSERT(s->running);

    struct timespec current_time;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &current_time);
    long long elapsed_time =
        (long long)(current_time.tv_sec - s->start.tv_sec) * 1000000000ll +
                    current_time.tv_nsec - s->start.tv_nsec;
    s->start = current_time;
    return elapsed_time;
#else
    stopwatch_stop(s);
    long long elapsed_time = stopwatch_get_elapsed_time(s);
	s->num_calls++ ;
    stopwatch_start(s);
#ifdef USE_TIMING
    /* Avoid imprecision from doing two timing calls. */
    s->start = s->end;
#endif
    return elapsed_time;
#endif
}

inline void stopwatch_time_to_string(char *buffer, long long time) {
#ifdef USE_TIMING
    sprintf(buffer,
            "%lld.%09lld seconds",
            time / 1000000000ll,
            time % 1000000000ll);
#else
    sprintf(buffer, "%lld instructions", time);
#endif
}

inline void stopwatch_reset_num_calls(stopwatch *s)
{
	s->num_calls = 0 ;
}

