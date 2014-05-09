#ifndef _STOPWATCH_H_
#define _STOPWATCH_H_

#include <stdbool.h>
#include <assert.h>

#define USE_TIMING

#ifdef USE_TIMING
#   include <time.h>
#else
#   include <perfmon/pfmlib_perf_event.h>
#endif

typedef struct stopwatch stopwatch;

struct stopwatch {
    bool initialized;

    bool running;
    bool measurements_valid;

    long long elapsed_time;

#ifdef USE_TIMING
    struct timespec start;
    struct timespec end;
#else
    struct perf_event_attr perf_attributes;
    int file_descriptor;
#endif
};

/* Class functions. */
void stopwatches_setup();

void stopwatches_teardown();

/* Stopwatch functions. */
void stopwatch_init(stopwatch *s);

void stopwatch_destroy(stopwatch *s);

void stopwatch_start(stopwatch *s);

void stopwatch_stop(stopwatch *s);

bool stopwatch_is_running(stopwatch *s);

long long stopwatch_get_elapsed_time(stopwatch *s);

/* Returns the elapsed time of the stopwatch. */
long long stopwatch_stop_and_start(stopwatch *s);

/* Utility function. */
void stopwatch_time_to_string(char *buffer, long long time);

#endif  /* _STOPWATCH_H_ */
