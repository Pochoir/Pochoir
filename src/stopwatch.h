#ifndef _STOPWATCH_H_
#define _STOPWATCH_H_

#include <stdbool.h>
#include <assert.h>

#define USE_TIMING
//#define CLOCK_FLAG CLOCK_THREAD_CPUTIME_ID
#define CLOCK_FLAG CLOCK_MONOTONIC

#ifdef USE_TIMING
#   include <time.h>
#else
#   include <perfmon/pfmlib_perf_event.h>

/* Indices into perf events. */
#define RAW_COUNT_INDEX 0
#define TIME_ENABLED_INDEX 1
#define TIME_RUNNING_INDEX 2

#endif

#define NANO_SEC 1000000000ll

typedef struct stopwatch stopwatch;

struct stopwatch {
#ifndef NDEBUG
    bool initialized;
    bool running;
    bool measurements_valid;
#endif

    long long elapsed_time;
    long long measure_time;

	unsigned long long num_calls ;

#ifdef USE_TIMING
    struct timespec start;
    struct timespec end;
#else
    struct perf_event_attr perf_attributes;
    int file_descriptor;
#endif
};

/* Class functions. */
inline void stopwatches_setup();

inline void stopwatches_teardown();

/* Stopwatch functions. */
inline void stopwatch_init(stopwatch *s);

inline void stopwatch_destroy(stopwatch *s);

#ifndef NDEBUG
inline void stopwatch_start(stopwatch *s);
inline void stopwatch_stop(stopwatch *s);
inline bool stopwatch_is_running(stopwatch *s);
#endif

/* Utility function. */
inline void stopwatch_time_to_string(char *buffer, long long time);
inline double stopwatch_time_to_double(long long time) ;

inline void stopwatch_reset_num_calls(stopwatch *s) ;

inline long long stopwatch_compute_measurement_time() ;

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
#ifndef NDEBUG
    s->initialized = true;
#endif

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

#ifndef NDEBUG
    s->running = false;
    s->measurements_valid = false;
#endif
	s->num_calls = 0 ;
	s->measure_time = stopwatch_compute_measurement_time () ;
}

inline void stopwatch_destroy(stopwatch *s) {
#ifndef USE_TIMING
    close(s->file_descriptor);
#endif
}

#ifndef NDEBUG
inline void stopwatch_start(stopwatch *s) {
#ifdef KXING_DEBUG_STOPWATCH
    printf("stopping: %p\n", s);
#endif
#ifndef NDEBUG
    assert(s->initialized);
    assert(!s->running);
    //CILK_ASSERT(s->initialized);
    //CILK_ASSERT(!s->running);
#endif

#ifdef USE_TIMING
    clock_gettime(CLOCK_FLAG, &(s->start));
#else
    ioctl(s->file_descriptor, PERF_EVENT_IOC_RESET, 0);
    ioctl(s->file_descriptor, PERF_EVENT_IOC_ENABLE, 0);
#endif
#ifndef NDEBUG
    s->running = true;
    s->measurements_valid = false;
#endif
}

#else

#ifdef USE_TIMING
#define stopwatch_start(s) {\
    clock_gettime(CLOCK_FLAG, &(s->start)); \
}
#else 
#define stopwatch_start(s) { \
    ioctl(s->file_descriptor, PERF_EVENT_IOC_RESET, 0); \
    ioctl(s->file_descriptor, PERF_EVENT_IOC_ENABLE, 0); \
}
#endif
#endif

#ifndef NDEBUG
inline void stopwatch_stop(stopwatch *s) {
#ifndef NDEBUG
    assert(s->initialized);
    assert(s->running);
    //CILK_ASSERT(s->initialized);
    //CILK_ASSERT(s->running);
#endif

#ifdef USE_TIMING
    clock_gettime(CLOCK_FLAG, &(s->end));
    s->elapsed_time =
            (long long)(s->end.tv_sec - s->start.tv_sec) * NANO_SEC +
            s->end.tv_nsec - s->start.tv_nsec;
	assert (s->elapsed_time >= 0) ;
#else
    ioctl(s->file_descriptor, PERF_EVENT_IOC_DISABLE, 0);

    // See:
    //  http://sourceforge.net/p/perfmon2/libpfm4/ci/master/tree/perf_examples/perf_util.h
    // and
    //  http://web.eece.maine.edu/~vweaver/projects/perf_events/perf_event_open.html
    // for details.
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
        // Scale the values. 
        s->elapsed_time = (long long)((double)(raw_count) * time_enabled / time_running);
    }
    assert(s->elapsed_time >= 0);
    //CILK_ASSERT(s->elapsed_time >= 0);
#endif

#ifndef NDEBUG
    s->running = false;
    s->measurements_valid = true;
#endif
}
#else
#ifdef USE_TIMING
#define stopwatch_stop(s) { \
    clock_gettime(CLOCK_FLAG, &(s->end)); \
    s->elapsed_time = \
            (long long)(s->end.tv_sec - s->start.tv_sec) * NANO_SEC + \
            s->end.tv_nsec - s->start.tv_nsec; \
}
#else
#define stopwatch_stop(s) { \
    ioctl(s->file_descriptor, PERF_EVENT_IOC_DISABLE, 0); \
    /*
      See:
      http://sourceforge.net/p/perfmon2/libpfm4/ci/master/tree/perf_examples/perf_util.h
      and
      http://web.eece.maine.edu/~vweaver/projects/perf_events/perf_event_open.html
      for details.
    */\
    uint64_t values[3]; \
    read(s->file_descriptor, values, sizeof(values)); \
    uint64_t raw_count = values[RAW_COUNT_INDEX]; \
    uint64_t time_enabled = values[TIME_ENABLED_INDEX]; \
    uint64_t time_running = values[TIME_RUNNING_INDEX]; \
    if (time_enabled < time_running) { \
        assert(false); \
    } \
    if (time_running == 0) { \
        assert(raw_count == 0); \
        s->elapsed_time = 0; \
    } else { \
        /* Scale the values. */ \
        s->elapsed_time = (long long)((double)(raw_count) * time_enabled / time_running); \
    } \
}
#endif
#endif

inline bool stopwatch_is_running(stopwatch *s) {
#ifndef NDEBUG
    return s->running;
#else
    return false ; //can't say if the clock is running
#endif
}

#define stopwatch_get_elapsed_time(s, t) { \
    assert(s->initialized); \
    assert(s->measurements_valid); \
    t = s->elapsed_time; \
}

#define stopwatch_stop_and_start(s, t) { \
    stopwatch_stop(s); \
	stopwatch_get_elapsed_time(s, t); \
	s->num_calls++ ; \
    stopwatch_start(s); \
}

inline double stopwatch_time_to_double(long long time) {
	assert (time >= 0) ;
	double t ;
#ifdef USE_TIMING
	t = (double) time / (double) 1000000 ; //in milliseconds
#else
	t = time ;
#endif
	return t ;
}

inline void stopwatch_time_to_string(char *buffer, long long time) {
#ifdef USE_TIMING
    sprintf(buffer,
            "%lld.%09lld seconds",
            time / NANO_SEC,
            time % NANO_SEC);
#else
    sprintf(buffer, "%lld instructions", time);
#endif
}

inline void stopwatch_reset_num_calls(stopwatch *s)
{
	s->num_calls = 0 ;
}

#define stopwatch_diff(start, end, elapsed_time) \
{ \
	elapsed_time = (long long) (end.tv_sec - start.tv_sec) * NANO_SEC + end.tv_nsec - start.tv_nsec ; \
}

inline long long stopwatch_compute_measurement_time()
{
	struct timespec start1, start2 ;
	struct timespec end1, end2 ;
	clock_gettime(CLOCK_FLAG, &start1) ;
	int TIMES = 25000 ;
	for (int i = 0 ; i < TIMES ; i++)
	{
		//1
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//2
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//3
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//4
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//5
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//6
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//7
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//8
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//9
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//10
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	
	
		//11
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//12
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//13
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//14
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//15
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//16
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//17
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//18
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//19
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	

		//20
		clock_gettime(CLOCK_FLAG, &start2) ;
		clock_gettime(CLOCK_FLAG, &end2) ;	
	}
	clock_gettime(CLOCK_FLAG, &end1) ;
	long long elapsed_time = 0 ;
	stopwatch_diff (start1, end1, elapsed_time) ;
	elapsed_time /= (40 * TIMES) ;
	char c [100] ;
	stopwatch_time_to_string(c, elapsed_time) ;
	//cout << "measurement time " << c << endl ;
	
	return elapsed_time ;
}

#define stopwatch_stop_and_start2(s, t) { \
	clock_gettime(CLOCK_FLAG, &(s->end)); \
	stopwatch_diff (s->start, s->end, s->elapsed_time) ; \
	assert (s->elapsed_time >= 0) ; \
    t = s->elapsed_time ; \
	s->num_calls++ ; \
    clock_gettime(CLOCK_FLAG, &(s->start)); \
}


#define stopwatch_stop_and_start_sub(s, t) { \
	clock_gettime(CLOCK_FLAG, &(s->end)); \
	stopwatch_diff (s->start, s->end, s->elapsed_time) ; \
	s->elapsed_time -= s->measure_time ; \
	s->elapsed_time = max(s->elapsed_time, (long long) 0) ; \
	assert (s->elapsed_time >= 0) ; \
    t = s->elapsed_time ; \
	s->num_calls++ ; \
    clock_gettime(CLOCK_FLAG, &(s->start)); \
}
#endif  /* _STOPWATCH_H_ */
