#ifndef PARSE_H
#define PARSE_H
#include "Gordion.h"
#include <sstream>

namespace Parse
    {
    void parse_help  () ;			// Print command help
    void parse_line  (const string&)  ;		// Parse input line
    void parse_cmd   (const string&)  ;		// Parse command word
    bool parse_set   (istringstream&) ;		// Parse "set" commands
    bool parse_reset (istringstream&) ;		// Parse "reset" commands
    bool parse_build (istringstream&) ;		// Parse "build" commands
    bool parse_call  (istringstream&) ;		// Parse "call" commands
    bool parse_do    (istringstream&) ;		// Parse "do" commands
    bool parse_eval  (istringstream&) ;		// Parse "evalaute" commands
    bool parse_gen   (istringstream&) ;		// Parse "generator" commands
    bool parse_add   (istringstream&) ;		// Parse "add" command
    bool parse_print (istringstream&) ;		// Parse "print" commands
    bool parse_purge (istringstream&) ;		// Parse "purge" command
    bool parse_test  (istringstream&) ;		// Parse "test" commands
    bool parse_read  (istringstream&) ;		// Parse "read" command
    bool parse_write (istringstream&) ;		// Parse "write" command
    bool parse_save  (istringstream&) ;		// Parse "save" command
    bool parse_load  (istringstream&) ;		// Parse "load" command
    void read_input  (istream&,  bool) ;	// Read input line

    inline void read_input  (istream&& in, bool prompt)
	{
	read_input (in, prompt) ;
	}
    inline bool eos (istringstream& line)	// End of string?
	{
	return line.peek() == EOF ;
	}

    inline bool isword (const string& w, const string& cmd, int min=1)
	{					// Abbreviating comparison
	return (w.size() >= min) && 0 == cmd.compare (0, w.size(), w) ;
	}

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

    inline bool isstar (istringstream& line)	// Next word == "*"?
	{
	string	word ;
	auto	pos { line.tellg() } ;
	if (parse_args (line,word) && word == "*") return true ;
	line.clear() ;
	line.seekg (pos) ;
	return false ;
	}

    inline static bool	  echo	   { false } ;		// Echo commands?
    inline static bool	  timing   { true  } ;		// Report command times?
    inline static bool	  awaiting { false } ;		// Awaiting user input?
    inline static ostream myout	   { cout.rdbuf() } ;	// Initial output stream
    } ;

#endif
