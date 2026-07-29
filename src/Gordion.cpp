#include "Gordion.h"
#include "Init.h"
#include "Parse.h"
#include "Global.h"
#include "Save.h"
#include <fstream>
#include <filesystem>
#include <csignal>
#include <getopt.h>

static auto	starttime { high_resolution_clock::now() } ;	// Starting time
static string	dotfile  { "./.gordion" } ;		// Default start-up file
static string	infile ;				// Command input file

void process_flags (int argc, char** argv)	// Get command line flags
    {
    program = *argv ;
    int c ;
    while ((c = getopt (argc, argv, "f:?")) != -1)
	{
	switch (c)
	    {
	    case 'f': infile = optarg ; break ;
	    case '?': 
	    default : cout << program << cmdargs ; std::exit (1) ;
	    }
	}
    cout << std::boolalpha ;
    }

string nicetime (float secs)
    {
    constexpr sv fmt	{ "{:8.2f} {}" };
    float	 min	{ 60 } ;
    float	 hr	{ 60 * min } ;
    float	 day	{ 24 * hr  } ;
    if      (secs > day)	return format (fmt, secs / day, "day") ;
    else if (secs > hr)		return format (fmt, secs / hr,  "hr")  ;
    else if (secs > min)	return format (fmt, secs / min, "min") ;
    else			return format (fmt, secs,       "sec") ;
    }

string nicemem (float rss)
    {
    constexpr sv  fmt	{ "{:8.2f} {}" };
    constexpr int KB	{ 1024 } ;
    constexpr int MB	{ 1024 * KB } ;
    constexpr int GB	{ 1024 * MB } ;
    if      (rss > GB)	return format (fmt, rss / GB, "GB") ;
    else if (rss > MB)	return format (fmt, rss / MB, "MB")  ;
    else if (rss > KB)	return format (fmt, rss / KB, "KB") ;
    else		return format (fmt, rss,      "bytes") ;
    }

[[noreturn]] void quit (int code)		// Exit program
    {
    auto	end	{ high_resolution_clock::now() } ;	// Ending time
    float	delta	{ std::chrono::duration<float>(end - starttime).count() } ;

    if (Parse::timing && delta > 1)
	{
	Parse::myout << "Clock time: " << nicetime (delta)     << "\n" ;
	Parse::myout << "CPU time:   " << nicetime (cputime()) << "\n" ;
	Parse::myout << "Max memory: " << nicemem  (maxmem())  << "\n" ;
	}
    exit (code) ;
    }

void sig_catch (int signum)			// Catch interrupts
    {
    global.interrupt = true ;
    std::cerr << "\nCaught Interrupt!\n" ;
    if (Parse::awaiting) quit(1) ;
    }

int main (int argc, char** argv)		// Main program
    {
    signal (SIGINT, sig_catch) ;
    initialize() ;
    process_flags (argc, argv) ;
    for (int i (optind) ; i < argc ; i++)
	{
	try { Save::load_save (-1,argv[i]) ; }	// Load command line file
	catch (const BadInput& e) { cout << e.what() << "\n" ; }
	}
    if (infile == "" && std::filesystem::exists (dotfile))
	infile = dotfile ;
	
    if (infile != "" && infile != "-")
	{
	ifstream in { infile } ;
	if (in.good())
	    {
	    try { Parse::read_input (in, false) ; }
	    catch (const exception& e)
		{
		cout << "Aborting input from " << infile << ": " << e.what() << "\n" ;
		}
	    }
	else cout << "Cannot read input from file " << infile << "\n" ;
	}
    Parse::read_input (std::cin, true) ;	// Read commands from stdin
    quit (global.interrupt) ;
    }
