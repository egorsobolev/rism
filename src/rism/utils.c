#include <rism/utils.h>

#include <time.h>

walltime_t walltime()
{
	struct timespec tms;
	int64_t micros;
	if (clock_gettime(CLOCK_MONOTONIC, &tms)) {
		return -1;
	}
	micros = tms.tv_sec * 1000000L;
	micros += tms.tv_nsec / 1000 + (int64_t) (tms.tv_nsec % 1000 >= 500);
	return micros;
}
