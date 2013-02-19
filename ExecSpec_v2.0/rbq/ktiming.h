#ifndef _KTIMING_H_
#define _KTIMING_H_

#include <stdint.h>

typedef uint64_t clockmark_t;

uint64_t ktiming_diff_nanosec(const clockmark_t* const start, const clockmark_t* const end);
float ktiming_diff_sec(const clockmark_t* const start, const clockmark_t* const end);
clockmark_t ktiming_getmark(void);

#endif
