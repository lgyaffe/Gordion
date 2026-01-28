#include "Init.h"
#include "Parse.h"
#include "Global.h"
#include "Save.h"
#include <fstream>
#include <csignal>
#include <getopt.h>
#include <unistd.h>

static	string	dfltstart  { "./.gordion" } ;	// Default start-up file
static	string	infile ;

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

int main (int argc, char** argv)		// Main program
    {
    signal (SIGINT, Parse::sig_catch) ;
    initialize() ;
    process_flags (argc, argv) ;
    for (int i (optind) ; i < argc ; i++)
	{
	try {
	    Save::load_save (-1,argv[i]) ;	// Load save file from command line
	    }
	catch (const BadInput& e)
	    {
	    cout << e.what() << "\n" ;
	    }
	}
    if (infile == "" && std::filesystem::exists (dfltstart))
	infile = dfltstart ;
	
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
    Parse::quit (global.interrupt) ;
    }
