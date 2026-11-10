#ifndef OS_H
#define OS_H

#if defined(__APPLE__) || defined(__MACH__)
#include <sys/sysctl.h>
#define PAGESIZE 16 * 1024
inline size_t memtot()
    {
    uint64_t mem { 0 } ;
    size_t   len { sizeof(mem) } ;
    return sysctlbyname ("hw.memsize", &mem, &len, NULL, 0) ? 0 : mem ;
    }
inline size_t memused()
    {
    struct rusage usage ;
    getrusage (RUSAGE_SELF, &usage) ;
    return usage.ru_maxrss ;
    }
inline double cpumsec()
    {
    struct rusage usage ;
    getrusage (RUSAGE_SELF, &usage) ;
    return usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1.e-3 ;
    }

#elif defined (__linux__)
#include <unistd.h>
#include <sys/resource.h>
#define PAGESIZE 4096
inline size_t memtot()
    {
    long pages	   { sysconf (_SC_PHYS_PAGES) } ;
    long page_size { sysconf (_SC_PAGE_SIZE) } ;
    return pages * page_size ;
    }
inline size_t memused()
    {
    struct rusage usage ;
    getrusage (RUSAGE_SELF, &usage) ;
    return usage.ru_maxrss * 1024 ;
    }
inline double cpumsec()
    {
    struct rusage usage ;
    getrusage (RUSAGE_SELF, &usage) ;
    return usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1.e-3 ;
    }

#elif defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#include <psapi.h>
#define PAGESIZE 4096
inline size_t memtot()
    {
    ULONGLONG mem_siz { 0 } ;
    if (GetPhysicallyInstalledSystemMemory (&mem_kb)) mem_siz *= 1024 ;
    return mem_siz ;
    }
inline size_t memused()
    {
    PROCESS_MEMORY_COUNTERS pmc ;
    GetProcessMemoryInfo (GetCurrentProcess(), &pmc, sizeof(pmc)) ;
    return pmc.PeakWorkingSetSize ;
    }
inline double cpumsec()
    {
    HANDLE	hproc { GetCurrentProcess() } ;
    FILETIME	creat, exit, kernel, user ;
    GetProcessTimes (hproc, &creat, &exit, &kern, &user) ;
    ULONGLONG uli ;
    uli.LowPart  = user.dwLowDateTime  ;
    uli.HighPart = user.dwHighDateTime ;
    return static_cast<double>(uli) * 1.e-4 ;
    }
#endif

#endif
