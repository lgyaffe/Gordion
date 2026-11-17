#ifndef PARALLAL_H
#define PARALLAL_H

#ifdef PARALLEL
#include <tbb/task_arena.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for_each.h>

using tbb::task_arena ;
using tbb::blocked_range ;
using tbb::parallel_for ;
using tbb::parallel_for_each ;
using intrange = blocked_range<int> ;
using sizrange = blocked_range<size_t> ;

#define UNSEQ std::execution::unseq,
#define FOR_EACH parallel_for_each
#define TASK_ARENA(n,cap,code) \
    if (n && n > 1) { task_arena arena(n); arena.execute([&cap]{code;}); } \
    else code ;
#else
#include <execution>

#define UNSEQ
#define FOR_EACH std::for_each
#define TASK_ARENA(n,cap,code) code
#endif

#endif
