#ifndef GRIPE_H
#define GRIPE_H
#include <string>

class BadInput : public std::runtime_error
    {
    public:
    BadInput (const char* msg) : std::runtime_error (msg) {}
    } ;

class Abort : public std::runtime_error
    {
    public:
    Abort (const char* msg) : std::runtime_error (msg) {}
    } ;

class IOerror : public std::runtime_error
    {
    public:
    IOerror (const char* msg) : std::runtime_error (msg) {}
    } ;

class Fatal : public std::logic_error
    {
    public:
    Fatal (const char* msg) : std::logic_error (msg) {}
    } ;

[[noreturn]] inline void gripe (const std::string msg)
    {
    std::cout << std::flush ;
    throw BadInput (msg.data()) ;
    }

[[noreturn]] inline void abort (const std::string msg)
    {
    std::cout << std::flush ;
    throw Abort (msg.data()) ;
    }

[[noreturn]] inline void ioerror (const std::string msg)
    {
    std::cout << std::flush ;
    throw IOerror (msg.data()) ;
    }

[[noreturn]] inline void fatal (const std::string msg)
    {
    std::cout << std::flush ;
    throw Fatal (msg.data()) ;
    }

#endif
