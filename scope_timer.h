#ifndef SCOPE_TIMER_H
#define SCOPE_TIMER_H

#include <chrono>

class scope_timer
{
public:
	using clock      = std::chrono::high_resolution_clock;
	using time_point = clock::time_point;
	using duration   = clock::duration;
	
private:
	duration  &timer;
	time_point start;
	
public:
	inline scope_timer(duration &timer) : timer(timer)
	{
		 start = clock::now();
	}
	
	inline ~scope_timer()
	{
		timer += clock::now() - start;
	}
};

#endif
