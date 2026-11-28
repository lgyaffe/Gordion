#ifndef PARALLAL_H
#define PARALLAL_H

#ifdef PARALLEL
#include <tbb/task_arena.h>
#include <tbb/parallel_for_each.h>

#ifdef DELPHI
using intrange = tbb::blocked_range<int> ;
using sizrange = tbb::blocked_range<size_t> ;
#endif

#define FOR_EACH tbb::parallel_for_each
#define TASK_ARENA(n,cap,code) \
    if (n && n > 1) \
	{ tbb::task_arena arena(n); arena.execute([&cap]{code;}); } \
    else code ;
#else
#include <execution>
#define FOR_EACH std::for_each
#define TASK_ARENA(n,cap,code) code
#endif

#endif
