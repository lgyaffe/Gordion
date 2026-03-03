#ifndef PARSE_H
#define PARSE_H
#include "Gordion.h"
#include <sstream>

namespace Parse
    {
    void read_input  (istream&&, bool) ;	// Read input line
    void read_input  (istream&,  bool) ;	// Read input line
    void parse_line  (const string&)  ;		// Parse input line
    void parse_cmd   (const string&)  ;		// Parse command word
    bool parse_set   (istringstream&) ;		// Parse "set" commands
    bool parse_reset (istringstream&) ;		// Parse "reset" commands
    bool parse_call  (istringstream&) ;		// Parse "call" commands
    bool parse_build (istringstream&) ;		// Parse "build" commands
    bool parse_eval  (istringstream&) ;		// Parse "evalaute" commands
    bool parse_do    (istringstream&) ;		// Parse "do" commands
    bool parse_add   (istringstream&) ;		// Parse "add" command
    bool parse_test  (istringstream&) ;		// Parse "test" commands
    bool parse_print (istringstream&) ;		// Parse "print" commands
    bool parse_read  (istringstream&) ;		// Parse "read" command
    bool parse_write (istringstream&) ;		// Parse "write" command
    bool parse_save  (istringstream&) ;		// Parse "save" command
    bool parse_load  (istringstream&) ;		// Parse "load" command
    bool isstar	     (istringstream&) ;		// Is next word == "*"?
    bool eos	     (istringstream&) ;		// End of string?

    void	 print_help  () ;		// Print command help
    void	 sig_catch (int) ;		// Catch interrupts
    [[noreturn]] void quit (int) ;		// Exit program

    template <typename... Args>
    bool parse_args (istringstream& line, Args&... args) // Parse command args
	{
	auto pos { line.tellg() } ;
	if (pos != -1)
	    {
	    if ((line >> ... >> args) && eos(line)) return true ;
	    line.clear() ;
	    line.seekg (pos) ;
	    }
	return false;
	}

    inline bool eos (istringstream& line)	// End of string?
	{
	return line.peek() == EOF ;
	}

    inline bool isword(const string& w, const string& cmd, int min=1)
	//
	// Abbreviating word comparison
	//
	{
	return (w.size() >= min) && 0 == cmd.compare (0, w.size(), w) ;
	}
    } ;

#endif
