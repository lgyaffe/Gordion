#include "Main.h"
#include "Global.h"
#include "Init.h"
#include "Parse.h"
#include "Save.h"
#include "OS.h"
#include <fstream>
#include <filesystem>
#include <csignal>
#include <getopt.h>
#include <unistd.h>

static auto	starttime { high_resolution_clock::now() } ;	// Starting time
static string	dotfile  { "./." + progname } ;		// Default start-up file
static string	infile ;				// Command input file

void process_flags (int argc, char** argv)	// Get command line flags
    {
    cmdname = *argv ;
    int c ;
    while ((c = getopt (argc, argv, "f:?")) != -1)
	{
	switch (c)
	    {
	    case 'f': infile = optarg ; break ;
	    case '?': 
	    default : cout << progname << cmdargs ; std::exit (1) ;
	    }
	}
    cout << std::boolalpha ;
    }

[[noreturn]] void quit (int code)		// Exit program
    {
    auto end  { high_resolution_clock::now() } ;	// Ending time
    auto tics { duration_cast<milliseconds>(end - starttime).count() } ;

    if (Parse::timing && tics > 10)
	{
	Parse::myout << "Clock time:    " << nicetime (tics)      << "\n"
		     << "CPU time:      " << nicetime (cpumsec()) << "\n"
		     << "Memory used:   " << nicemem  (memused()) << "\n"
		     << "System memory: " << nicemem  (memtot())  << "\n" ;
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
	try { Save::load_save (-1,argv[i]) ; }	// Load file from command line
	catch (const exception& e) { cout << e.what() << "\n" ; infile = "" ; }
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
