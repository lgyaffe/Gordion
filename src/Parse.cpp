#include "Parse.h"
#include "Build.h"
#include "Commute.h"
#include "Rep.h"
#include "Numerics.h"
#include "Test.h"
#include "Print.h"
#include "Save.h"
#include "Blab.h"
#include "Gripe.h"
#include <fstream>
#include <chrono>
#include <regex>

using hires = std::chrono::high_resolution_clock ;
using msec  = std::chrono::milliseconds ;
using sec   = std::chrono::seconds ;
using std::chrono::duration_cast ;

static bool	echo	  { false } ;
static bool	timing	  { true  } ;
static bool	awaiting  { false } ;
static auto	starttime { hires::now() } ;
static ostream	myout	  { cout.rdbuf() } ;

void Parse::read_input (istream& in, bool prompt)	// Read input
    {
    string buf ;
    while (!global.interrupt)
	{
	myout << flush ;
	if (prompt) myout << ": " ;
	awaiting = true ;
	bool gotline { getline (in,buf) } ;
	awaiting = false ;
	if (!gotline) break ;
	if (echo) cout << ": " << buf << "\n" << flush ;
	try { parse_line (buf) ; }
	catch (const BadInput& e)
	    {
	    if (&in != &std::cin) throw e ;
	    else myout << e.what() << "\n" ;
	    }
	catch (const IOerror& e)
	    {
	    if (&in != &std::cin) throw e ;
	    else myout << e.what() << "\n" ;
	    }
	catch (const Fatal& e)
	    {
	    myout << "Fatal error: " << e.what() << "\n" ;
	    quit (1) ;
	    }
	}
    }

void Parse::read_input (istream&& in, bool prompt)	// Read input
    {
    read_input (in, prompt) ;
    }

void Parse::sig_catch (int signum)			// Catch interrupts
    {
    global.interrupt = true ;
    std::cerr << "\nCaught Interrupt!\n" ;
    if (awaiting) quit(1) ;
    }

[[noreturn]] void Parse::quit (int code)		// Exit program
    {
    if (global.interrupt) FINALIZE ;
    auto end  { hires::now() } ;
    auto secs { duration_cast<sec>(end - starttime).count() } ;
    if (timing && secs > 10) myout << "Total time: " << secs << " sec\n" ;
    exit (code) ;
    }

void Parse::parse_line (const string& buf)		// Parse input line
    {
    int	start(0) ;
    while (buf[start] == ':') ++start ;
    while (start < buf.size() && !global.interrupt)
	{
	while (isspace (buf[start])) ++start ;
	if (buf[start] == '#') break ;
	auto end { buf.find(';', start) } ;
	if (end == buf.npos) end = buf.size() ;
	auto next { end + 1 } ;
	while (end > start && isspace (buf[end-1])) --end ;
	parse_cmd (buf.substr (start, end - start)) ;
	start = next ;
	}
    }

void Parse::parse_cmd (const string& buf)		// Parse command
    {
    istringstream	line  { buf } ;
    bool		valid { true } ;
    string		cmd ;

    line >> cmd >> std::boolalpha ;			// get command

    auto start { hires::now() } ;

    if (cmd[0] == '#')			return ;	// skip cmoment
    else if (cmd == "?")		print_help() ;
    else if (isword(cmd,"help"))	print_help() ;
    else if (isword(cmd,"set"))		valid = parse_set   (line) ;
    else if (isword(cmd,"reset"))	valid = parse_reset (line) ;
    else if (isword(cmd,"call"))	valid = parse_call  (line) ;
    else if (isword(cmd,"add"))		valid = parse_add   (line) ;
    else if (isword(cmd,"build"))	valid = parse_build (line) ;
    else if (isword(cmd,"evaluate"))	valid = parse_eval  (line) ;
    else if (isword(cmd,"do"))		valid = parse_do    (line) ;
    else if (isword(cmd,"test"))	valid = parse_test  (line) ;
    else if (isword(cmd,"print"))	valid = parse_print (line) ;
    else if (isword(cmd,"read"))	valid = parse_read  (line) ;
    else if (isword(cmd,"write"))	valid = parse_write (line) ;
    else if (isword(cmd,"save"))	valid = parse_save  (line) ;
    else if (isword(cmd,"load"))	valid = parse_load  (line) ;
    else if (isword(cmd,"quit"))	quit (1) ;
    else valid = false ;
    if (!valid && cmd.size()) gripe ("Don't understand: " + buf) ;

    auto end  { hires::now() } ;
    auto tics { duration_cast<msec>(end - start).count() } ;
    if (timing && tics > 10) cout << tics << " msec\n" ;
    }

bool Parse::parse_read (istringstream& line)		// Parse "read" command
    {
    string word ;
    if (parse_args (line,word))
	{
	if (ifstream input {word})
	    {
	    read_input (std::move(input), false) ;
	    return true ;
	    }
	else gripe ("Cannot read input from file " + word) ;
	}
    return false ;
    }

bool Parse::parse_write (istringstream& line)		// Parse "write" command
    {
    static std::streambuf*	stdoutbuf ;
    static ofstream		fileout ;
    string			word ;
    if (parse_args (line,word))
	{
	if (word != "-")
	    {
	    cout.flush () ;
	    if (fileout.is_open()) fileout.close() ;
	    fileout.open (word) ;
	    if (fileout.is_open())
		{
		if (!stdoutbuf) stdoutbuf = std::cout.rdbuf() ;
		cout.rdbuf (fileout.rdbuf()) ;
		}
	    else gripe ("Cannot open " + word) ;
	    }
	else if (stdoutbuf)
	    {
	    if (fileout.is_open()) fileout.close() ;
	    cout.rdbuf (stdoutbuf) ;
	    }
	return true ;
	}
    return false ;
    }

bool Parse::parse_load (istringstream& line)		// Parse "load" command
    {
    int		set (-1) ;
    string	word ;

    if	    (parse_args (line,set))		Save::load_vev  (set) ;
    else if (parse_args (line,set,word))	Save::load_save (set,word) ;
    else if (parse_args (line,word))		Save::load_save (-1, word) ;
    else gripe ("Load what info?") ;
    return true ;
    }

bool Parse::parse_save (istringstream& line)		// Parse "save" command
    {
    bool	valid { true } ;
    string	word, word2 ;

    if (eos(line)) gripe ("Save which information?") ;
    else if (parse_args (line,word,word2) || parse_args(line,word))
	{
	if      (isword(word,"sys"))
	    Save::save_sys (word2.empty() ? global.sysfile : word2) ;
	else if (isword(word,"vev"))
	    Save::save_vev (word2.empty() ? global.vevfile : word2) ;
	else valid = false ;
	}
    else valid = false ;
    return valid ;
    }

bool Parse::parse_set (istringstream& line)		// Parse "set" commands
    {
    uint	i ;
    bool	flag ;
    bool	valid { true } ;
    doub	value ;
    string	word ;
    if (line >> word)
	{
	if (isword(word,"representation") && parse_args (line,word))
	    {
	    try { global.repnum = Rep::known (word) ; }
	    catch (const exception& e)
		{
		gripe ("Unknown representation " + word) ;
		}
	    }
	else if (isword(word,"stage") && parse_args (line,word))
	    {
	    if      (isword(word,"gauge"))	Global::stageinit (0) ;
	    else if (isword(word,"fermi"))	Global::stageinit (1) ;
	    else valid = false ;
	    }
	else if (isword(word,"approx") && parse_args (line,flag))
	    {
	    global.approx = ObsList::obs.approx = flag ;
	    }
	else if (isword(word,"approx") && parse_args (line,i))
	    {
	    global.approx = ObsList::obs.approx = i ;
	    }
	else if (isword(word,"autoEgens") && parse_args (line,flag))
	    {
	    Gen::autoEgens = flag ;
	    }
	else if (isword(word,"autosave") && parse_args (line,flag))
	    {
	    global.autosave = flag ;
	    }
	else if (isword(word,"blab") && parse_args (line,word,i))
	    {
	    Blab::setblab (word, i) ;
	    }
	else if (isword(word,"checkobs") && parse_args (line,flag))
	    {
	    Obs::check = flag ;
	    }
	else if (isword(word,"dots") && parse_args (line,flag))
	    {
	    SymbStr::dots = flag ;
	    }
	else if (isword(word,"echo") && parse_args (line,flag))
	    {
	    echo = flag ;
	    }
	else if (isword(word,"gennorm") && parse_args (line,flag))
	    {
	    if (flag != Gen::gennorm)
		{
		Gen::gennorm = flag ;
		Gen::normalize () ;
		Global::clearbuild (false) ;
		}
	    }
	else if (isword(word,"geoswap") && parse_args (line,flag))
	    {
	    global.geoswap = flag ;
	    if (global.geoswap) global.autosave = true ;
	    }
	else if (isword(word,"maxthread") && parse_args (line,i))
	    {
	    global.maxthread = i ;
	    }
	else if (isword(word,"minlim") && parse_args (line,i))
	    {
	    numerics.minlim = i ;
	    }
	else if (isword(word,"minmax") && parse_args (line,i))
	    {
	    numerics.minmax = i ;
	    }
	else if (isword(word,"mintol") && parse_args (line,value))
	    {
	    numerics.mintol = value ;
	    }
	else if (isword(word,"odemax") && parse_args (line,i))
	    {
	    numerics.odemax = i ;
	    }
	else if (isword(word,"odetol") && parse_args (line,value))
	    {
	    numerics.odetol = value ;
	    }
	else if (isword(word,"oknegeig") && parse_args (line,flag))
	    {
	    global.oknegeig = flag ;
	    }
	else if (isword(word,"symcurv") && parse_args (line,flag))
	    {
	    global.symcurv = flag ;
	    }
	else if (isword(word,"rkmethod") && parse_args (line,word))
	    {
	    auto beg { RKdef::list.cbegin() } ;
	    auto end { RKdef::list.cend() } ;
	    auto ptr { std::find_if (beg, end,
		    [word](const RKdef& rk) { return isword(word,rk.name,4) ;})} ;
	    if (ptr != end) numerics.rk = *ptr ;
	    else gripe ("Unknown RK method " + word) ;
	    }
	else if (isword(word,"savedir") && parse_args (line,word))
	    {
	    global.savedir = word ;
	    global.sysstream.close() ;
	    global.vevstream.close() ;
	    }
	else if (isword(word,"speclim") && parse_args (line,i))
	    {
	    numerics.speclim = i ;
	    }
	else if (isword(word,"svdlim") && parse_args (line,value))
	    {
	    numerics.svdlim = value ;
	    }
	else if (isword(word,"sysfile") && parse_args (line,word))
	    {
	    auto pos { word.find_last_of ('/') } ;
	    if (pos != word.npos)
		{
		global.savedir.assign (word, 0, pos) ;
		global.sysfile.assign (word, pos+1) ;
		}
	    else global.sysfile = word ;
	    global.sysstream.close() ;
	    }
	else if (isword(word,"timing") && parse_args (line,flag))
	    {
	    timing = flag ;
	    }
	else if (isword(word,"vevappend") && parse_args (line,flag))
	    {
	    global.vevappend = flag ;
	    }
	else if (isword(word,"vevfile") && parse_args (line,word))
	    {
	    auto pos { word.find_last_of ('/') } ;
	    if (pos != word.npos)
		{
		global.savedir.assign  (word, 0, pos) ;
		global.vevfile.assign (word, pos+1) ;
		}
	    else global.vevfile = word ;
	    global.vevstream.close() ;
	    }
	else if (isword(word,"MMAappend") && parse_args (line,flag))
	    {
	    global.MMAappend = flag ;
	    }
	else if (isword(word,"MMAdir") && parse_args (line,word))
	    {
	    global.MMAdir = word ;
	    global.MMAstream.close() ;
	    }
	else if (isword(word,"MMAfile") && parse_args (line,word))
	    {
	    auto pos { word.find_last_of ('/') } ;
	    if (pos != word.npos)
		{
		global.MMAdir.assign  (word, 0, pos) ;
		global.MMAfile.assign (word, pos+1) ;
		}
	    else global.MMAfile = word ;
	    global.MMAstream.close() ;
	    }
	else if (isword(word,"MMAlimit") && parse_args (line,i))
	    {
	    global.info().MMAlimit = i ;
	    }
	else if (isword(word,"MMAlist") && line >> word)
	    {
	    do  {
		Obs o { word } ; o.canon() ;
		uint indx { ObsList::obs.find (o) } ;
		if (indx != UINT_MAX)
		    {
		    global.info().MMAlist.insert (indx) ;
		    }
		else gripe ("Unknown Obs " + word) ;
		}
	    while (line >> word) ;
	    }
	else
	    {
	    int indx { Coupling::indx (word) } ;
	    if (indx >= 0 && parse_args (line,value))
		{
		doub oldv { Coupling::list[indx].value } ; 
		Coupling::list[indx].value = value ;
		if (oldv != value && global.stage && theory.euclid)
		    {
		    char8 mass	{ "mass" } ;
		    int   mindx	{ Coupling::indx (mass) } ;
		    if (indx == mindx && global.massreinit)
			cout << "Reinitializing fermion vev's\n" ;
			Numerics::numericsinit() ;
		    }
		}
	    else valid = false ;
	    }
	}
    else valid = false ;
    return valid ;
    }

bool Parse::parse_reset (istringstream& line)		// Parse "reset" commands
    {
    bool	valid { true } ;
    string	word ;
    if (parse_args (line,word))
	{
	if (isword(word,"representation"))	global.repnum = 0 ;
	else if (isword(word,"maxthread"))	global.maxthread = 0 ;
	else if (isword(word,"MMAlist"))	global.info().MMAlist.clear() ;
	else if (isword(word,"stage"))		Global::stageinit (0) ;
	else if (isword(word,"blab"))		Blab::resetblab () ;
	else if (isword(word,"mintol"))		numerics.mintol = numerics.dflttol ;
	else if (isword(word,"odetol"))		numerics.odetol = numerics.dflttol ;
	else if (isword(word,"svdlim"))		numerics.svdlim = numerics.dfltlim ;
	else if (isword(word,"rkmethod"))	numerics.rk = RKdef::list.back() ;
	else if (isword(word,"numerics"))	Numerics::numericsinit() ;
	else if (isword(word,"sysfile"))
	    {
	    global.sysfile.clear() ;
	    global.sysstream.close() ;
	    }
	else if (isword(word,"vevfile"))
	    {
	    global.vevfile.clear() ;
	    global.vevstream.close() ;
	    }
	else if (isword(word,"MMAfile"))
	    {
	    global.MMAfile.clear() ;
	    global.MMAstream.close() ;
	    }
	else valid = false ;
	}
    else valid = false ;
    return valid ;
    }

bool Parse::parse_build (istringstream& line)		// Parse "build" commands
    {
    int		i ;
    bool	valid { true } ;
    bool	isH   { !theory.euclid } ;
    auto&	rep   { global.repnum } ;
    string	word ;
    if (line >> word)
	{
	if (isword(word,"observable") && parse_args (line,i))
	    {
	    if (i < global.maxord()) Global::clearbuild (true) ;
	    Build::mk_obs (i) ;
	    }
	else if (isword(word,"all") && parse_args (line,i))
	    {
	    if (i < global.maxord()) Global::clearbuild (true) ;
	    Build::mk_obs (i) ;
	    Build::mk_grad () ;
	    for (int j(0) ; j < Rep::list.size() ; ++j)
		{
		Build::mk_curv (j) ;
		if (isH) Build::mk_lagr (j) ;
		}
	    Build::mk_geos () ;				// Must be last for autosave
	    }
	else if (isword(word,"geodesics") && eos(line))	Build::mk_geos () ;
	else if (isword(word,"gradient")  && eos(line))	Build::mk_grad () ;
	else if (isword(word,"curvature") && eos(line))	Build::mk_curv (rep) ;
	else if (isword(word,"curvature") && parse_args(line,word))
	    {
	    if (isword(word,"all"))
		{
		for (int j(0) ; j < Rep::list.size() ; ++j)	Build::mk_curv (j) ;
		}
	    else
		{
		try {
		    Build::mk_curv (rep = Rep::known (word)) ;
		    }
		catch (const exception& e)
		    {
		    gripe ("Unknown representation " + word) ;
		    }
		}
	    }
	else if (isword(word,"lagrange") && eos(line) && isH) Build::mk_lagr (rep) ;
	else if (isword(word,"lagrange") && parse_args(line,word) && isH)
	    {
	    if (isword(word,"all"))
		{
		for (int j(0) ; j < Rep::list.size() ; ++j)	Build::mk_lagr (j) ;
		}
	    else
		{
		try {
		    Build::mk_lagr (rep = Rep::known (word)) ;
		    }
		catch (const exception& e)
		    {
		    gripe ("Unknown representation " + word) ;
		    }
		}
	    }
	else valid = false ;
	}
    else valid = false ;
    return valid ;
    }

bool Parse::parse_eval (istringstream& line)		// Parse "evaluate" commands
    {
    int		i ;
    bool	valid { true } ;
    bool	isH   { !theory.euclid } ;
    auto&	rep   { global.repnum } ;
    string	word ;
    if (line >> word)
	{
	try {
	    if (isword(word,"freeenergy")	&& eos(line) && !isH)
		{
		numerics.eval_H (true) ;
		}
	    else if (isword(word,"hamiltonian")	&& eos(line) && isH)
		{
		numerics.eval_H (true) ;
		}
	    else if (isword(word,"gradient")    && eos(line))
		{
		numerics.eval_grad (true) ;
		}
	    else if (isword(word,"curvature"))
		{
		if (eos(line))		numerics.eval_curv (rep, true, 1) ;
		else if (isstar (line))	numerics.eval_curv (rep, true, 2) ;
		else valid = false ;
		}
	    else if (isword(word,"lagrange")    && eos(line) && isH)
		{
		numerics.eval_lagr (rep, true) ;
		}
	    else if (isword(word,"delta")       && eos(line))
		{
		numerics.eval_delta (true) ;
		}
	    else if (isword(word,"spectrum")    && eos(line) && isH)
		{
		numerics.eval_spectra (rep, true) ;
		}
	    else if (isword(word,"spectrum") && parse_args(line,word) && isH)
		{
		if (isword(word,"all"))	
		    {
		    for (int j(0) ; j < Rep::list.size() ; ++j)
			{
			numerics.eval_spectra (j, true) ;
			}
		    }
		else
		    {
		    try {
			numerics.eval_spectra (rep = Rep::known (word), true) ;
			}
		    catch (const exception& e)
			{
			gripe ("Unknown representation " + word) ;
			}
		    }
		}
	    else if (isword(word,"geodesics"))
		{
		if (eos(line)) gripe ("Print how many geodesic values?") ;
		else if (isstar (line))		numerics.eval_geos (0) ;
		else if (parse_args (line,i))	numerics.eval_geos (i) ;
		else valid = false ;
		}
	    else if (isword(word,"geostats") && eos(line))
		{
		Build::do_geostats() ;
		Print::print_geostats() ;
		}
	    else valid = false ;
	    }
	catch (const Abort& e) { cout << e.what() << "\n" ; }
	}
    else valid = false ;
    return valid ;
    }

bool Parse::parse_do (istringstream& line)		// Parse "do" commands
    {
    bool	valid { true } ;
    doub	v0, v1, inc ;
    string	word ;
    if (line >> word)
	{
	try {
	    if (isword(word,"step") && eos(line))
		{
		numerics.do_step () ;
		}
	    else if (isword(word,"minimize") && eos(line))
		{
		numerics.do_minimize () ;
		}
	    else 
		{
		int indx { Coupling::indx (word) } ;
		if (indx >= 0 && parse_args (line,v0,v1,inc))
		    {
		    numerics.do_flow (indx, v0, v1, inc) ;
		    }
		else valid = false ;
		}
	    }
	catch (const Abort& e) { cout << e.what() << "\n" ; }
	}
    else valid = false ;
    return valid ;
    }

bool Parse::parse_add (istringstream& line)		// Parse "add" command
    {
    bool	valid { true } ;
    string	word ;
    if (line >> word)
	{
	if (isword(word,"generator"))
	    {
	    doub	coef  (0) ;
	    short	order (0) ;
	    char	left, right ;
	    bool	gotcoef { line >> coef } ;
	    bool	gotword { false } ;

	    if (!gotcoef)
		{
		line.clear() ;
		if (line >> word)
		    {
		    std::regex	pattern { "\\((\\d+)\\)" } ;
		    std::smatch match ;
		    if (std::regex_match (word, match, pattern))
			{
			order = stoi (match[1].str()) ;
			gotword = false ;
			}
		    else
			{
			gotword = true ;
			coef = 1.0 ;
			}
		    }
		}

	    OpSum sum ;
	    do  {
		if (!gotcoef && !gotword && !(line >> coef))
		    {
		    line.clear() ;
		    coef = 1.0 ;
		    }
		gotcoef = false ;
		if (gotword || (line >> word))
		    {
		    try {
			Obs o { word } ; o.canon() ;
			uint indx { ObsList::obs.find (o) } ;
			if (indx != UINT_MAX)
			    {
			    short cord { ObsList::obs(indx).corder } ;
			    if (!order) order = cord ;
			    else if (order != cord)
				gripe ("Invalid generator: inconsistent orders") ;
			    }
			}
		    catch (const BadInput&) {}
		    if (order)
			{
			Op op { word, order } ;
			sum.emplace_back (Op::store(op), coef) ;
			}
		    else gripe ("unknown order for " + word) ;
		    gotword = false ;
		    }
		} while (!eos(line)) ;

	    if (sum.size())
		{
		int n { Gen::addgen (std::move(sum)) } ;
		if (n)	cout << "  " << n ;
		else	cout << "  No" ;
		cout << " generators added\n" ;
		}
	    else valid = false ;
	    }
	else valid = false ;
	}
    else valid = false ;
    return valid ;
    }

bool Parse::parse_print (istringstream& line)		// Parse "print" commands
    {
    using namespace Print ;

    bool	valid { true } ;
    bool	isH   { !theory.euclid } ;
    uint	i, j, k ;
    string	word ;
    if (line >> word)
	{
	if  (isword(word,"observables"))
	    {
	    if (eos(line)) gripe ("Print which observables?") ;
	    else if (isstar (line))		print_obs () ;
	    else if (parse_args (line,i))	print_obs (i) ;
	    else if (parse_args (line,i,j))	print_obs (i,j) ;
	    else if (parse_args (line,word))
		{
		if (word[0] == '(')		print_obs_select (word) ;
		else				print_obs (word) ;
		}
	    else valid = false ;
	    }
	else if (isword(word,"operator"))
	    {
	    if (eos(line)) gripe ("Print which operators?") ;
	    else if (isstar (line))		print_op () ;
	    else if (parse_args (line,i))	print_op (i) ;
	    else if (parse_args (line,i,j))	print_op (i,j) ;
	    else valid = false ;
	    }
	else if (isword(word,"symmetry"))
	    {
	    if (eos(line)) gripe ("Print which symmetries?") ;
	    else if (isstar (line))		print_symm () ;
	    else if (parse_args (line,word))	print_symm (word) ;
	    else valid = false ;
	    }
	else if (isword(word,"representation"))
	    {
	    if (eos(line))			print_rep () ;
	    else if (parse_args (line,word))	print_rep (word) ;
	    else valid = false ;
	    }
	else if (isword(word,"generator"))
	    {
	    if (eos(line))			print_gen (false) ;
	    else if (isstar (line))		print_gen (true) ;
	    else if (parse_args (line,i))	print_gen (i,true) ;
	    else valid = false ;
	    }
	else if (isword(word,"gradient"))
	    {
	    if (eos(line))  gripe ("Print which gradient elements?") ;
	    else if (isstar (line))		print_grad () ;
	    else if (parse_args (line,i))	print_grad (i) ;
	    else if (parse_args (line,i,j))	print_grad (i,j) ;
	    else valid = false ;
	    }
	else if (isword(word,"curvature"))
	    {
	    if (eos(line)) gripe ("Print which curvature elements?") ;
	    else if (isstar (line))		print_curv () ;
	    else if (parse_args (line,i))	print_curv (i) ;
	    else if (parse_args (line,i,j))	print_curv (i,j) ;
	    else if (parse_args (line,i,j,k))	print_curv (i,j,k) ;
	    else valid = false ;
	    }
	else if (isword(word,"geodesics"))
	    {
	    if (eos(line)) gripe ("Print which geodesic elements?") ;
	    else if (isstar (line))		print_geodesic () ;
	    else if (parse_args (line,i))	print_geodesic (i) ;
	    else if (parse_args (line,i,j))	print_geodesic (i,j) ;
	    else valid = false ;
	    }
	else if (isword(word,"lagrange") && isH)
	    {
	    if (eos(line)) gripe ("Print which lagrange elements?") ;
	    else if (isstar (line))		print_lagrange () ;
	    else if (parse_args (line,i))	print_lagrange (i) ;
	    else if (parse_args (line,i,j))	print_lagrange (i,j) ;
	    else valid = false ;
	    }
	else if (isword(word,"mode") && isH)
	    {
	    if (eos(line))  gripe ("Print which eigenmodes?") ;
	    else if (isstar (line))		 print_mode () ;
	    else if (parse_args (line,i))	 print_mode (i) ;
	    else valid = false ;
	    }
	else if (isword(word,"freeenergy")  && eos(line) && !isH) print_freeenergy  () ;
	else if (isword(word,"hamiltonian") && eos(line) &&  isH) print_hamiltonian () ;
	else if (isword(word,"spectrum")    && eos(line) &&  isH) print_spectrum    () ;
	else if (isword(word,"baseobs")     && eos(line))	print_base      () ;
	else if (isword(word,"blablevels")  && eos(line))	print_blab      () ;
	else if (isword(word,"bcktlist")    && eos(line))	print_bcktlist  () ;
	else if (isword(word,"cache")       && eos(line))	print_cache     () ;
	else if (isword(word,"couplings")   && eos(line))	print_couplings () ;
	else if (isword(word,"fermiinit")   && eos(line))	print_fermiinit () ;
	else if (isword(word,"geostats")    && eos(line))	print_geostats  () ;
	else if (isword(word,"obsstats")    && eos(line))	print_obsstats  () ;
	else if (isword(word,"primary")     && eos(line))	print_primary   () ;
	else if (isword(word,"rkmethods")   && eos(line))	print_rkmethods () ;
	else if (isword(word,"state")       && eos(line))	print_state     () ;
	else if (isword(word,"stats")       && eos(line))	print_stats     () ;
	else if (isword(word,"symmsets")    && eos(line))	print_symmsets  () ;
	else if (isword(word,"sysindex")    && eos(line))	print_sysindex  () ;
	else if (isword(word,"vevindex")    && eos(line))	print_vevindex  () ;
	else valid = false ;
	}
    else valid = false ;
    return valid ;
    }

bool Parse::parse_test (istringstream& line)		// Parse "test" commands
    {
    uint	i(0) ;
    uint	nobs  { ObsList::obs.nobs() } ;
    bool	valid { true } ;
    string	word, word2 ;
    if (line >> word)
	{
	if (isword(word,"jacobi"))
	    {
	    if (eos(line))
		{
		Test::jacobi () ;
		}
	    else if (parse_args (line,i))
		{
		if (i < nobs) Test::jacobi (i) ;
		else gripe (format("Invalid observable number {}", i)) ;
		}
	    else if (parse_args (line,word,word2,i))
		{
		if (i < nobs) Test::jacobi (word,word2,i) ;
		else gripe (format("Invalid observable number {}", i)) ;
		}
	    else valid = false ;
	    }
	else if (isword(word,"irreps") && eos(line))
	    {
	    Test::irreps () ;
	    }
	else valid = false ;
	}
    else valid = false ;
    return valid ;
    }

bool Parse::parse_call (istringstream& line)		// Parse "call" commands
    {
    string	word, word2 ;
    uint	i, j ;
    if (line >> word)
	{
	if (isword(word,"canon") && parse_args (line,word))
	    {
	    Obs		obs { word } ;
	    const char*	sgn { obs.canon() < 0 ? "-" : "" } ;
	    cout << word << " -> " << sgn << obs << "\n" ;
	    }
	else if (isword(word,"transform") && parse_args (line,word,word2))
	    {
	    auto	indx	{ Symm::known(std::move(word)).item } ;
	    const auto&	s	{ Symm::list[indx] } ;
	    Obs		oa	{ word2 } ;
	    Obs		ob	{ oa } ;
	    bool	Cflip	{ ob.is_fermi() || ob.is_Loop() } ;
	    int		start	( Cflip && s.isCodd() ? ob.size() - 1 : 0 ) ;
	    int		sgn	{ ob.trans(s,start) } ;
	    char	c	{ sgn == 1 ? '+' : '-' } ;
	    cout << s.name << "(" << oa << ") = " << c << " " << ob << "\n" ;
	    }
	else if (isword(word,"u1bound") && parse_args (line,word))
	    {
	    Obs		obs	{ word } ;
	    short	xord	{ obs.u1bound(0) } ;
	    cout << obs << ": xorder >= " << xord << "\n" ;
	    }
	else if (isword(word,"noEbound") && parse_args (line,word))
	    {
	    Obs		obs	{ word } ;
	    short	xord	{ obs.noEbound(ObsList::obs) } ;
	    cout << obs << ": xorder >= " << xord << "\n" ;
	    }
	else if (isword(word,"classify") && parse_args (line,word))
	    {
	    Obs		obs	{ word } ;
	    const char*	sgn	{ obs.canon() < 0 ? "-" : "" } ;
	    cout << word << " -> " << sgn << obs << "\n" ;
	    bool ok { obs.classify (ObsList::obs) } ;
	    if (ok) cout << "classified: " << obs
			 << " xorder " << obs.xorder << "\n" ;
	    else    cout << "classify failed\n" ;
	    }
	else if (isword(word,"commute") && parse_args (line,i,j))
	    {
	    auto&	gens  { global.info().gens[global.repnum] } ;
	    uint	ngens ( gens.size() ) ;
	    uint	nobs  ( ObsList::obs.size() ) ;

	    if (i < ngens && j < nobs)
		{
		Obs	obs	{ ObsList::obs(j) } ;
		ObsPoly	poly	{ obs, ObsList::obs } ;
		ObsList	tmplist	{ "ParseTemp" } ;
		PolyMap	ans	{ tmplist } ;

		Commute::commute_poly (gens[i], poly, ans) ;
		if (obs.imag() && gens[i].imag) ans.negate() ;
		cout << "[" << gens[i] << ", " << obs << "] = \n\t" ;
		cout << ans << "\n" ;
		}
	    }
	else if (isword(word,"commute") && parse_args (line,word,word2))
	    {
	    Op  	op	{ word, -1 } ;
	    Obs 	obs	{ word2 } ;
	    PolyTerm	factor	{ Polyindx(), 1 } ;
	    ObsList	tmplist	{ "ParseTemp" } ;
	    PolyMap	ans	{ tmplist } ;

	    Commute::do_commute (op, obs, factor, tmplist, ans) ;
	    cout << "[" << op << ", " << obs << "] = " ;
	    cout << (obs.imag() ? "i (\n\t" : "(\n\t") ;
	    cout << ans << " )\n" ;
	    }
	else if (isword(word,"do_inner") && parse_args (line,word))
	    {
	    Obs 	obs (word) ;
	    if (obs.is_Eloop())
		{
		ObsList	tmplist { "ParseTemp" } ;
		PolyTerm	factor { Polyindx(), 1 } ;
		PolyMap	ans { tmplist } ;

		Commute::do_inner (obs, factor, tmplist, ans) ;
		cout << "do_inner [" << obs << "] = \n\t" ;
		cout << ans << "\n" ;
		}
	    else cout << "Bad Obs type\n" ;
	    }
	else return false ;
	}
    return true ;
    }

void Parse::print_help ()					// Print command help
    {
    cout << "Usage: " << program << cmdargs << "\n" ;
    cout << R"(Commands:
add		generator	[(order)] [<coeff>] <Op> [...]

build		observables	<maxorder>
		geodesics
		gradient
		curvature	[all | <repname>]
		lagrange (+)	[all | <repname>]
		all		<maxorder>

call		canon		<obs>
		u1bound		<obs>
		noEbound	<obs>
		classify	<obs>
		commute		<op> <obs>
		do_inner	<Eloop_obs>
		transform	<symmetry> <obs>

do		step
		minimize
		<coupling>	<value0> <value1> <inc>

evaluate	hamiltonian (+)
		freeenergy (++)
		gradient
		curvature	[*]
		delta
		lagrange (+)
		spectrum (+)	[all | <repname>]
		geodesics	* | <print_limit>
		geostats

load		<sys_infofile>
		[<vev_set_#>] <vev_datafile>	
		<vev_set_#>

print		baseobs
		blablevels
		bucketlist
		cache
		couplings
		curvature	* | <term_#> [<gen_#> [<gen_#>]]
		fermiinit
		freeenergy (++)
		generator	[* | <gen_#>]
		geodesics	* | <obs_#> [<gen_#>]
		geostats
		gradient	* | <term_#> [<gen_#>]
		hamiltonian (+)
		lagrange (+)	* | <gen_#> [<gen_#>]
		mode (+)	* | <mode_#>
		observable	(<corder>,<xorder>)
		observable	* | <obs_#> [<obs_#>] | <obs>
		obsstats
		operator	* | <op_#>  [<op_#>]
		primary
		representation  [* | <repname>]
		rkmethods
		state
		stats
		spectrum (+)
		symmetry	* | <symmname>
		symmsets
		sysindex
		vevindex

reset		blab
		maxthread
		mintol
		numerics
		odetol
		representation
		rkmethod
		stage
		svdlim
		sysfile
		vevfile
		MMAfile
		MMAlist

save		sys		[<filename>]
		vev		[<filename>]

set		stage		gauge | fermi
		approx		false | true
		autoEgens	true | false
		autosave	false | true |
		blab		<source_file> <value>
		checkobs	false | true
		dots		false | true
		echo		false | true
		gennorm		false | true
		geoswap		false | true
		oknegeig	false | true
		symcurv		true | false
		massreinit	true | false
		maxthread	<integer>
		minlim		<integer>
		minmax		<integer>
		odemax		<integer>
		mintol		<value>
		odetol		<value>
		representation  <repname>
		rkmethod	<rkname>
		savedir		<directory>
		speclim		<integer>
		svdlim		<value>
		sysfile		<filename>
		timing		true | false
		vevfile		<filename>
		MMAappend	true | false
		MMAdir		<directory>
		MMAfile		<filename>
		MMAlimit	<integer>
		MMAlist		<obs> [<obs> [...]]
		<coupling>	<value>

test	 	irreps
		jacobi
		jacobi		<obs_#>
		jacobi		<op1> <op2> <obs_#>

read		<filename>
write		<filename> | -

#		<comment>
quit
help
?

N.B.:
<...> denotes specified user input
[...] denotes optional input
a | b denotes alternatives a or b
(+)  indicates command only valid for Hamiltonian theories
(++) indicates command only valid for Euclidean theories
Semicolons may spearate multiple commands on a single input line.
)" << "\n" ;
    }

bool Parse::isstar (istringstream& line)		// Next word == "*"?
    {
    string	word ;
    auto	pos { line.tellg() } ;
    if (parse_args (line,word) && word == "*") return true ;
    line.clear() ;
    line.seekg (pos) ;
    return false ;
    }

bool Parse::eos (istringstream& line)			// End of string?
    {
    return line.peek() == EOF ;
    }
