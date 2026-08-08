#include "Save.h"
#include "Canon.h"
#include "Gen.h"
#include "Numerics.h"
#include "Rep.h"
#include "Blab.h"
#include "Gripe.h"
#include <filesystem>

#ifndef PAGESIZE					// used below in alignas
#if __APPLE__ && __arm64__ 
#define PAGESIZE 16 * 1024
#else
#define PAGESIZE 4096
#endif
#endif

using std::ios_base ;

static constexpr int	indxsiz { sizeof (SysIndex) } ;
static constexpr int	elemsiz { sizeof (Element)  } ;

union Buffer						// Alignment buffer
    {
    Element		elem ;
    array<symb,elemsiz>	byte ;
    Buffer() {} ;
    } ;
static Buffer		buffer ;

void Save::save_sys ()					// Save sys info
    {
    const auto&	savedir		{ global.savedir } ;
    auto&	syspath		{ global.info().sysfile.path } ;
    auto&	sysstream	{ global.info().sysfile.stream } ;
    string	file		{ global.mk_filename("sys") } ;
    string	path		{ global.addsubdir (savedir) + file } ;
    auto	mode		{ ios::in | ios::out | ios::trunc | ios::binary } ;

    if (sysstream.is_open() && path == syspath)
	{
	const auto& geopolys	{ global.data().geos[0] } ;
	if (geopolys.entry().reclen && !geopolys.size())
	    gripe ("Rewriting sys-info file will lose swapped-out geodesic equations!") ;
	}
    sysstream.close() ;
    sysstream.open ( path, mode ) ;
    if (sysstream.is_open())
	{
	cout << "Writing sys-info file " << path << "\n" ;
	syspath = std::move (path) ;
	write_header (sysstream) ;
	}
    else gripe ("Cannot write sys-info file " + path) ;

    write_op	() ;
    write_obs	() ;
    write_gen 	() ;
    write_ham	() ;
    write_grad	() ;
    write_curv	() ;
    write_lagr	() ;
    write_stat	() ;
    write_geos	() ;
    write_sysindex () ;
    }

void Save::save_vev ()					// Save vev data
    {
    const auto&	blab		{ Blab::level(Blab::SAVE) } ;
    const auto&	nvev		{ numerics.nvev() } ;
    const auto& savedir		{ global.savedir } ;
    auto&	vevpath		{ global.info().vevfile.path } ;
    auto&	vevstream	{ global.info().vevfile.stream } ;
    auto	vevappend	{ global.info().vevfile.append } ;
    string	file		{ global.mk_filename("vev") } ;
    string	path		{ global.addsubdir (savedir) + file } ;
    uint	ncoup		{ Coupling::ncoup() } ;

    if (path != vevpath) vevstream.close() ;
    if (!vevstream.is_open())				// open vevfile
	{
	alignas (PAGESIZE) static char mybuf [512*1024] ;
	bool	exists	{ std::filesystem::exists (path) } ;
	auto	mode	{ ios::in | ios::out | ios::binary } ;
	fstream	stream ;
	stream.rdbuf()->pubsetbuf (mybuf, sizeof mybuf) ;

	if (!exists)					// touch into existence
	    {
	    stream.open (path, ios::out | ios::binary) ;
	    stream.close() ;
	    }
	else mode |= (vevappend ? ios::app : ios::trunc) ;

	stream.open (path, mode) ;
	if (stream.is_open())
	    {
	    if (exists && vevappend)
		{
		read_header (stream, path) ;
		if (!filehdr.is_vevfile())
		    abort ("Inconsistent vev data file " + path) ;
		cout << "Appending vev data to " << path << "\n" ;
		}
	    else
		{
		cout << "Writing vev data to " << path << "\n" ;
		write_header (stream, ncoup, nvev) ;
		}
	    vevpath	= std::move (path) ;
	    vevstream	= std::move (stream) ;
	    }
	else gripe ("Cannot write vev data file " + path) ;
	}
    if (blab) cout << "Writing data set to " << vevpath << "\n" ;
    write_coup () ;
    write_vev  () ;
    }

void Save::write_header (fstream& stream, uint ncoup, uint nvev) // Write file header
    {
    auto hashG { ncoup || nvev ? numerics.hashG : global.info(0).obshash } ;
    auto hashF { ncoup || nvev ? numerics.hashF : global.info(1).obshash } ;
    filehdr.version	= global.version ;
    filehdr.name	= global.stage ? theory.name : theory.parent() ;
    filehdr.hashF	= global.stage ? hashF : 0 ;
    filehdr.hashG	= hashG ;
    filehdr.ncoup	= ncoup ;
    filehdr.nvev	= nvev ;
    stream.write (cast_to<char*>(&filehdr), sizeof filehdr) ;
    if (stream.fail()) ioerror ("write_header: I/O error!") ;
    }

void Save::write_op ()						// Write Op record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	record	{ global.data().op  } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    const auto&	oplist	( global.info().ops ) ;
    uint	nop	( oplist.size() ) ;

    if (record.entry().id == RecordID::Op) record.clear() ;
    else fatal ("write_op: bad record ID!") ;
    for (int indx(0) ; indx < nop ; ++indx)
	{
	const Op& op	{ oplist[indx] } ;
	Element	  e	{ OpHdr ( op.order, (char) op.type, op.primary ) } ;

	record.push_back (e) ;
	auto	symbptr { op.data() } ;
	auto	elemptr { cast_to<const Element*> (symbptr) } ;
	auto	n       { op.size() / elemsiz } ;

	record.insert (record.end(), elemptr, elemptr + n) ;
	if (auto rem { op.size() % elemsiz })
	    {
	    buffer.byte.fill (X) ;
	    std::memcpy (buffer.byte.data(), symbptr + n * elemsiz, rem) ;
	    record.push_back (buffer.elem) ;
	    ++n ;
	    }
	(&record.back() - n)->hdr.len = n ;
	}
    record.entry().nelem = nop ;
    record.writerec (stream) ;
    record.free() ;
    if (stream.fail()) ioerror ("write_op: I/O error!") ;
    if (blab > 1) cout << "Saved Op\n" << flush ;
    }

void Save::write_obs ()					// Write Obs record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    const auto&	nobs	{ global.info().nobs } ;
    auto&	record	{ global.data().obs  } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    long	start	{ global.stage ? global.info(0).nobs : 0 } ;

    if (global.obs.swapped) reload_obs() ;
    if (record.entry().id == RecordID::Obs) record.clear() ;
    else fatal ("write_obs: bad record ID!") ;
    for (long indx(start) ; indx < start + nobs ; ++indx)
	{
	const Obs& obs	{ global.obs(indx) } ;
	Element	   e	{ ObsHdr ( obs.corder, obs.xorder, (char) obs.type ) } ;

	record.push_back (e) ;
	auto	symbptr	{ obs.data() } ;
	auto	elemptr	{ cast_to<const Element*> (symbptr) } ;
	auto	n       { obs.size() / elemsiz } ;

	record.insert (record.end(), elemptr, elemptr + n) ;
	if (auto rem { obs.size() % elemsiz })
	    {
	    buffer.byte.fill (X) ;
	    std::memcpy (buffer.byte.data(), symbptr + n * elemsiz, rem) ;
	    record.push_back (buffer.elem) ;
	    ++n ;
	    }
	(&record.back() - n)->hdr.len = n ;
	}
    record.entry().nelem = nobs ;
    record.writerec (stream) ;
    record.free() ;
    if (stream.fail()) ioerror ("write_obs: I/O error!") ;
    if (blab > 1) cout << "Saved Obs\n" << flush ;
    }

void Save::write_gen ()					// Write Gen record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	record	{ global.data().gen } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    int		optsiz	{ sizeof (OpTerm) / sizeof (Element) } ;
    uint	ngen	( 0 ) ;
    
    if (record.entry().id == RecordID::Gen) record.clear() ;
    else fatal ("write_gen: bad record ID!") ;
    for (short rep(0) ; rep < theory.nrep ; ++rep)
	{
	for (const auto& gen : global.info().gens[rep])
	    {
	    Element	e { GenHdr ( rep, gen.order, (char) gen.type, gen.T_odd ) } ;
	    auto	n { 1 + gen.size() * optsiz } ;
	    record.push_back (e) ;
	    e.coeff = gen.coeff ;
	    record.push_back (e) ;
	    auto elemptr { cast_to<const Element*> (gen.data()) } ;
	    record.insert (record.end(), elemptr, elemptr + n - 1) ;
	    (&record.back() - n)->hdr.len = n ;
	    ++ngen ;
	    }
	}
    record.entry().nelem = ngen ;
    record.writerec (stream) ;
    record.free() ;
    if (stream.fail()) ioerror ("write_gen: I/O error!") ;
    if (blab > 1) cout << "Saved Gen\n" << flush ;
    }

void Save::write_ham ()					// Write Ham record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	record	{ global.data().ham } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto&	Hterms	{ global.info().Hterms } ;
    ushort	nterms	( Hterms.size() ) ;

    if (record.entry().id == RecordID::Ham) record.clear() ;
    else fatal ("write_ham: bad record ID!") ;
    for (ushort i(0) ; i < nterms ; ++i)
	{
	record.add (Hterms[i].cpoly) ;
	}
    record.entry().nelem = nterms ;
    record.writerec (stream) ;
    record.free() ;
    if (stream.fail()) ioerror ("write_ham: I/O error!") ;
    if (blab > 1) cout << "Saved Ham\n" << flush ;
    }

void Save::write_grad ()				// Write Grad record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto&	record	{ global.data().grad } ;

    if (record.entry().id != RecordID::Grad)
	fatal ("write_grad: bad record ID!") ;

    record.writerec (stream) ;
    if (stream.fail()) ioerror ("write_grad: I/O error!") ;
    if (record.size() && blab > 1)
	cout << "Saved Grad\n" << flush ;
    }

void Save::write_curv ()				// Write Curv records
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;

    for (int rep(0) ; rep < global.data().curv.size() ; ++rep)
	{
	auto&	record	{ global.data().curv[rep] } ;
	auto&	repname { Rep::list[rep].name } ;

	if (record.entry().id != RecordID::Curv)
	    fatal (format("write_curv {}: bad record ID!",repname)) ;

	record.writerec (stream) ;
	if (stream.fail())
	    ioerror (format("write_curv {}: I/O error!", repname)) ;
	else if (record.size() && blab > 1)
	    cout << "Saved Curv " << repname << "\n" << flush ;
	}
    }

void Save::write_lagr ()				// Write Lagr records
    {
    if (theory.euclid) return ;
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;

    for (int rep(0) ; rep < global.data().lagr.size() ; ++rep)
	{
	auto&	record	{ global.data().lagr[rep] } ;
	auto&	repname { Rep::list[rep].name } ;

	if (record.entry().id != RecordID::Lagr)
	    fatal (format("write_lagr {}: bad record ID!",repname)) ;

	record.writerec (stream) ;
	if (stream.fail())
	    ioerror (format("write_lagr {}: I/O error!",repname)) ;
	else if (record.size() && blab > 1)
	    cout << "Saved Lagr " << Rep::list[rep].name << "\n" << flush ;
	}
    }

void Save::write_geos ()				// Write Geo records
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto	nbckts	{ global.info().bckt.size() } ;
    int		bcktnum (0) ;

    for (; bcktnum < nbckts ; ++bcktnum)
	{
	auto&	stream	{ global.info().sysfile.stream } ;
	auto&	record	{ global.data().geos[bcktnum] } ;

	if (record.entry().id != RecordID::Geos)
	    fatal ("write_geos: bad record ID!") ;

	record.writerec (stream) ;
	if (stream.fail())
	    ioerror ("write_geos: I/O error!") ;
	else if (record.size() && blab > 1)
	    cout << "Saved Geos [" << bcktnum << "/" << nbckts << "]\n" << flush ;
	}
    for (; bcktnum < global.data().geos.size() ; ++bcktnum)
	{
	auto& record { global.data().geos[bcktnum] } ;
	record.entry().reclen = 0 ;
	}
    }

void Save::write_geo_bckt (int bcktnum)			// Write single Geo bucket
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto&	record	{ global.data().geos[bcktnum] } ;
    auto	nbckts	{ global.info().bckt.size() } ;

    if (global.interrupt) return ;
    if (record.entry().id != RecordID::Geos)
	fatal (format("write_geos [{}]: bad record ID!",bcktnum)) ;

    std::lock_guard<std::mutex> lock (savemutex) ;
    stream.seekg (-indxsiz, ios_base::end) ;
    record.writerec (stream) ;
    if (stream.fail())
	ioerror (format("write_geo_bckt [{}]: I/O error!",bcktnum)) ;
    if (blab > 1)
	cout << "Saved Geos [" << bcktnum << "/" << nbckts << "]\n" << flush ;
    write_sysindex () ;
    rewrite_stat () ;
    if (global.geoswap) record.free() ;
    }

void Save::write_stat ()				// Write Stat record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto&	record	{ global.data().stat } ;
    auto&	count	{ global.info().count } ;
    uint	nelem	( sizeof (count) / sizeof (Counter) ) ;

    if (record.entry().id != RecordID::Stat)
	fatal ("write_stat: bad record ID!") ;

    record.entry().nelem   = nelem ;
    record.entry().reclen  = nelem * sizeof (Counter) / sizeof (Element) ;
    record.entry().filepos = static_cast<std::streamoff>(stream.tellp()) ;

    stream.write (cast_to<const char*>(&count), nelem * sizeof (Counter)) ;
    if (stream.fail()) ioerror ("write_stat: I/O error!") ;
    if (blab > 1) cout << "Saved Stat\n" << flush ;
    }

void Save::rewrite_stat ()				// Rewrite Stat record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto&	record	{ global.data().stat } ;
    auto&	count	{ global.info().count } ;
    uint	nelem	( sizeof (count) / sizeof (Counter) ) ;

    if (record.entry().id != RecordID::Stat)
	fatal ("rewrite_stat: bad record ID!") ;

    stream.seekg (record.entry().filepos, ios_base::beg) ;
    stream.write (cast_to<const char*>(&count), nelem * sizeof (Counter)) ;
    if (stream.fail()) ioerror ("rewrite_stat: I/O error!") ;
    if (blab > 1) cout << "Rewrote Stat\n" ;
    }

void Save::write_sysindex ()				// Write SysIndex
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto&	syspath	{ global.info().sysfile.path } ;
    auto&	index	{ global.info().sysindex } ;
    auto&	last	{ index.back() } ;

    last.id	 = RecordID::Indx ;
    last.nelem   = 1 ;
    last.reclen  = sizeof (index) / sizeof (Element) ;
    last.filepos = static_cast<std::streamoff>(stream.tellp()) ;
    
    stream.write (cast_to<const char*>(index.data()), sizeof index) ;
    stream.flush() ;
    if (stream.fail()) ioerror ("write_sysindex: I/O error!") ;
    std::filesystem::resize_file (syspath,stream.tellp()) ;
    if (blab > 1) cout << "Wrote SysIndex\n" ;
    }

void Save::write_coup ()				// Write Coupling's
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().vevfile.stream } ;
    auto	ptr	{ Coupling::list.data() } ;

    stream.seekp (0, ios::end) ;
    stream.write (cast_to<const char*>(ptr), coupsize()) ;
    if (stream.fail()) ioerror ("write_coup: I/O error!") ;
    if (blab > 1) cout << "Wrote Couplings\n" << flush ;
    }

void Save::write_vev ()					// Write Vev's
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    long	start	{ global.stage ? numerics.nvevG : 0 } ;
    auto&	stream	{ global.info().vevfile.stream } ;
    auto	ptr	{ numerics.vev.memptr() + start } ;

    stream.seekp (0, ios::end) ;
    stream.write (cast_to<const char*>(ptr), vevsize()) ;
    if (stream.fail()) ioerror ("write_vev: I/O error!") ;
    if (blab > 1) cout << "Wrote Vevs\n" << flush ;
    }

void Save::load_save (int set, string file)		// Load save file
    {
    string	path	{ file } ;
    const auto& savedir { global.savedir } ;

    if (!std::filesystem::exists (path))
	path = global.addsubdir (savedir) + file ;

    if (!std::filesystem::exists (path))
	path = global.addsubdir (savedir, !global.stage) + file ;

    if (!std::filesystem::exists (path))
	gripe ("Cannot find save file " + file) ;

    fstream stream { path, ios::in | ios::out | ios::binary } ;
    if (stream.is_open())
	{
	global.stageinit (read_header (stream, path)) ;
	if (filehdr.is_sysfile())			// sys-info file
	    {
	    if (set >= 0) gripe ("No set numbers in sys-info files") ;

	    cout << "Loading sys info from " << path << "\n" ;
	    global.info().sysfile.path	 = std::move (path) ;
	    global.info().sysfile.stream = std::move (stream) ;

	    read_sysindex () ;
	    read_op	() ;
	    read_obs	() ;
	    read_gen	() ;
	    read_ham	() ;
	    read_grad	() ;
	    read_curv	() ;
	    read_lagr	() ;
	    read_stat	() ;
	    read_geos	() ;
	    }
	else /* filehdr.is_vevfile() */		// vev data file
	    {
	    cout << "Loading vev data from " << path ;
	    if (set < 0)	cout << "\n" ;
	    else		cout << ", set " << set << "\n" ;

	    global.info().vevfile.path   = std::move (path) ;
	    global.info().vevfile.stream = std::move (stream) ;
	    load_vev (set) ;
	    }
	}
    else gripe ("Cannot open save file " + path) ;
    }

void Save::load_vev (int set)			// Load vev data set
    {
    if (global.info(global.stage).vevfile.stream.is_open())
	{
	if (read_coup (set,true)) read_vev (set) ;
	else gripe (format("Cannot load vev set {}", set)) ;
	Print::print_couplings() ;
	}
    else gripe ("No open vev data file") ;
    }

int Save::read_header (fstream& stream, const string& path) // Read save file header
    {
    const auto& ncoupG	 { Coupling::ncoup(0) } ;
    const auto& ncoupF	 { Coupling::ncoup(1) } ;
    const auto& nvevG	 { numerics.nvevG } ;
    const auto& nvevF	 { numerics.nvevF } ;
    const auto& vevhashG { numerics.hashG } ;
    const auto& vevhashF { numerics.hashF } ;
    const auto& obshashG { global.info(0).obshash } ;
    auto&	hdr	 { filehdr } ;

    stream.read (cast_to<char*>(&hdr), sizeof hdr) ;
    if (stream.fail()) ioerror ("read_header: I/O error!") ;

    string	hdrthy	{ hdr.name.data(), sizeof (hdr.name) } ;
    string	mythy	{ theory.name.data(), sizeof (theory.name) } ;
    string	mypar	{ theory.parent().data(), sizeof (theory.name) } ;

    static constexpr sv newer  { "Warning: file {} has newer program version {}" } ;
    static constexpr sv corrpt { "File {} has inconsistent Obs set" } ;
    static constexpr sv badnum { "File {} has wrong real num type" } ;
    static constexpr sv badset { "File {} has hash code mismatch, {} != {}" } ;
    static constexpr sv numcup { "File {} has wrong # couplings, {} != {}" } ;
    static constexpr sv numobs { "File {} has wrong # vev's, {} != {}" } ;
    static constexpr sv	badpar { "File {} theory {} != {} or {}" } ;
    static constexpr sv	badthy { "File {} theory {} != {}" } ;

    if (hdr.version.incompat()) gripe (format (badnum, path)) ;
    if (hdr.version.newer())	cout << format(newer, path, hdr.version.print()) ;

    if (theory.parent() == hdr.name)				// YM save file
	{
	if (hdr.hashF)
	    gripe (format (corrpt, path)) ;
	if (hdr.ncoup == 0 && hdr.nvev == 0)
	    return 0 ;						// YM sys file
	if (hdr.ncoup != ncoupG)
	    gripe (format (numcup, path, hdr.ncoup, ncoupG)) ;
	if (hdr.nvev  != nvevG)
	    gripe (format (numobs, path, hdr.nvev,  nvevG))  ;
	if (hdr.hashG != vevhashG)
	    gripe (format (badset, path, hdr.hashG, vevhashG)) ;
	return 0 ;						// YM vev file
	}
    else if (theory.name == hdr.name && theory.nf)		// QCD save file
	{
	if (hdr.hashG != obshashG)
	    gripe (format (badset, path, hdr.hashG, obshashG));
	if (hdr.ncoup == 0 && hdr.nvev == 0)
	    return 1 ;						// QCD sys file
	if (hdr.ncoup != ncoupF)
	    gripe (format (numcup, path, hdr.ncoup, ncoupF));
	if (hdr.nvev  != nvevF)
	    gripe (format (numobs, path, hdr.nvev,  nvevF)) ;
	if (hdr.hashF != vevhashF)
	    gripe (format (badset, path, hdr.hashF, vevhashF));
	return 1 ;						// QCD vev file
	}
    else if (theory.name != theory.parent())
	gripe (format(badpar, path, hdrthy, mypar, mythy)) ;
    else 
	gripe (format(badthy, path, hdrthy, mythy)) ;
    }

void Save::read_sysindex ()					// Load SysIndex
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto&	index	{ global.info().sysindex } ;
    auto	ptr	{ index.data() } ;

    stream.seekg (-sizeof index, ios_base::end) ;
    stream.read (cast_to<char*>(ptr), sizeof index) ;
    if (stream.fail()) ioerror ("read_sysindex: I/O error!") ;
    if (blab > 1) cout << "Read Index\n" ;
    }

void Save::read_op ()						// Read Op record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto&	record	{ global.data().op } ;

    if (record.entry().id != RecordID::Op)
	abort ("read_op: bad record ID!") ;
    record.readrec (stream) ;
    if (stream.fail()) ioerror ("read_op: I/O error!") ;

    uint	indx	{ 0 } ;
    auto	nop	{ record.entry().nelem } ;
    Element*	elemptr	{ record.data() } ;
    Element*	recend	{ elemptr + record.size() } ;
    OpList&	oplist	{ global.info().ops } ;

    oplist.clear () ;
    oplist.reserve (nop) ;
    while (elemptr < recend)
	{
	auto	len	{ elemptr->len() } ;
	OpHdr	hdr	{ elemptr++->hdr.op } ;
	symb*	ptr	{ cast_to<symb*> (elemptr) } ;
	symb*	end	{ ptr + len * elemsiz  } ; while (*--end == X) ;
	Str	s	{ ptr, end + 1 } ;
	Op	o	{ s, hdr } ;

	if (oplist.store(o) != indx++)
	    abort ("read_op: Inconsistent OpList") ;

	elemptr += len ;
	}
    record.free() ;
    if (indx != nop)
	abort ("read_op: Inconsistent save record!") ;

    if (blab > 1) cout << "Loaded Op\n" << flush ;
    }

void Save::read_obs ()					// Read Obs record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto&	record	{ global.data().obs } ;
    auto&	obslist	{ global.obs } ;
    auto	nobs	{ record.entry().nelem } ;

    if (record.entry().id != RecordID::Obs)
	abort ("read_obs: bad record ID!") ;

    if (global.obsswap) { global.obs.swapped = true ; return ; }

    record.readrec (stream) ;
    if (stream.fail()) ioerror ("read_obs: I/O error!") ;

    long	start	{ global.stage ? global.info(0).nobs : 0 } ;
    long	indx    { start } ;
    Element*	elemptr	{ record.data() } ;
    Element*	recend	{ elemptr + record.size() } ;

    obslist.purge   (start) ;
    obslist.reserve (start + nobs) ;
    while (elemptr < recend)
	{
	auto	len	{ elemptr->len() } ;
	ObsHdr	hdr	{ elemptr++->hdr.obs } ;
	symb*	ptr	{ cast_to<symb*> (elemptr) } ;
	symb*	end	{ ptr + len * elemsiz } ; while (*--end == X) ;
	Str	s	{ ptr, end + 1 } ;
	Obs	o	{ s, hdr } ;

	if (obslist.store (o) != indx++)
	    abort ("read_obs: Inconsistent global.obs") ;

	if (global.info().maxord < o.order() && o.corder == o.xorder)
	    global.info().maxord = o.order() ;

	elemptr += len ;
	}
    record.free() ;
    if (indx - start != nobs)
	abort ("read_obs: Inconsistent save record!") ;
    if (global.info(0).nobs + global.info(1).nobs != obslist.size())
	abort ("read_obs: Inconsistent ObsList size") ;

    if (blab > 1) cout << "Loaded Obs\n" << flush ;
    global.mk_bcktlist () ;
    if (global.stage) global.obs.do_fermiinit() ;
    Canon::cache.reload () ;
    numerics.hashG = global.info(0).obshash = filehdr.hashG ;
    numerics.hashF = global.info(1).obshash = filehdr.hashF ;
    numerics.nvevG = global.info(0).nobs ;
    numerics.nvevF = global.info(1).nobs ;
    numerics.initialize () ;
    }

void Save::reload_obs ()			// Reload Obs record
    {
    auto& obslist { global.obs } ;
    obslist.clear () ;
    for (int stage(0) ; stage < 2 - !theory.nf ; ++stage)
	{
	auto&	stream	{ global.info(stage).sysfile.stream } ;
	auto&	record	{ global.data(stage).obs } ;

	if (record.entry().id != RecordID::Obs)
	    abort ("reload_obs: bad record ID!") ;
	record.readrec (stream) ;
	if (stream.fail())
	    ioerror ("reload_obs: I/O error!") ;

	cout << "Reloading Obs: " << flush ;
	long		start	{ stage ? global.info(0).nobs : 0 } ;
	long		indx    { start } ;
	long		nobs	( record.entry().nelem ) ;
	Element*	elemptr	{ record.data() } ;
	Element*	recend	{ elemptr + record.size() } ;

	obslist.reserve (start + nobs) ;
	while (elemptr < recend)
	    {
	    auto	len	{ elemptr->len() } ;
	    ObsHdr	hdr	{ elemptr++->hdr.obs } ;
	    symb*	ptr	{ cast_to<symb*> (elemptr) } ;
	    symb*	end	{ ptr + len * elemsiz  } ; while (*--end == X) ;
	    Str	s	{ ptr, end + 1 } ;
	    Obs	o	{ s, hdr } ;

	    if (obslist.store (o) != indx++)
		abort ("reload_obs: Inconsistent global.obs") ;

	    elemptr += len ;
	    }
	record.free() ;
	if (indx - start != nobs)
	    abort ("reload_obs: Inconsistent save record!") ;
	}
    if (global.info(0).nobs + global.info(1).nobs != obslist.size())
	abort ("reload_obs: Inconsistent ObsList size") ;

    global.obs.swapped = false ;
    cout << "done\n" << flush ;
    }

void Save::read_gen ()					// Read Gen record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto&	record	{ global.data().gen } ;
    int		optsiz	{ sizeof (OpTerm) / sizeof (Element) } ;

    if (record.entry().id != RecordID::Gen)
	abort ("read_gen: bad record ID!") ;
    record.readrec (stream) ;
    if (stream.fail())
	ioerror ("read_gen: I/O error!") ;

    uint	nelem   ( 0 ) ;
    ulong	ngen	( record.entry().items() ) ;
    Element*	elemptr	{ record.data() } ;
    Element*	recend	{ elemptr + record.size() } ;

    for (auto& gens : global.info().gens)  gens.clear() ;
    for (auto& even : global.info().neven) even = 0 ;

    while (elemptr < recend)
	{
	Element	e	{ *elemptr++ } ;
	real	coeff	{ elemptr++->coeff } ;
	auto	len	{ e.len() } ;
	auto	rep	{ e.hdr.gen.rep } ;
	auto	nterm	{ (len - 1) / optsiz } ;
	OpTerm*	ptr	{ cast_to<OpTerm*> (elemptr) } ;
	OpSum	s	( ptr, ptr + nterm, global.info().ops ) ;
	Gen	gen	{ s, e, coeff } ;

	global.info().gens[rep].push_back (gen) ;

	if (!gen.T_odd) ++global.info().neven[rep] ;
	if (global.info().maxgen < gen.order)
	    global.info().maxgen = gen.order ;

	elemptr += len - 1 ;
	++nelem ;
	}
    record.free() ;

    if (nelem != record.entry().nelem)
	abort ("read_gen: Inconsistent save record!") ;

    if (blab > 1) cout << "Loaded Gen\n" << flush ;
    }

void Save::read_ham ()				// Read Ham record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto& 	record	{ global.data().ham } ;
    auto&	Hterms	{ global.info().Hterms } ;
    auto	nterms	{ Hterms.size() } ;
    int		i(0) ;

    if (record.entry().id != RecordID::Ham)
	abort ("read_ham: bad record ID!") ;
    record.readrec (stream) ;
    if (stream.fail())
	ioerror ("read_ham: I/O error!") ;

    if (Hterms.size() != record.entry().nelem)
	abort ("read_ham: Inconsistent save record!") ;

    for (const auto& poly : record)
	{
	if (i >= nterms) abort ("read_ham: Inconsistent save record!") ;
	auto& Hterm { Hterms[i++] } ;
	Hterm.cpoly.clear() ;
	Hterm.cpoly.add (poly) ;
	}
    if (i != nterms) abort ("read_ham: Inconsistent save record!") ;
    record.free() ;
    if (blab > 1) cout << "Loaded Ham\n" << flush ;
    }

void Save::read_grad ()				// Read Grad record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto& 	record	{ global.data().grad } ;

    if (record.entry().id != RecordID::Grad)
	abort ("read_grad: bad record ID!") ;
    record.readrec (stream) ;
    if (stream.fail())
	ioerror ("read_grad: I/O error!") ;
    else if (record.size() && blab > 1)
	cout << "Loaded Grad\n" << flush ;
    }

void Save::read_curv ()				// Read Curv records
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;

    for (int rep(0) ; rep < theory.nrep ; ++rep)
	{
	auto&	record	{ global.data().curv[rep] } ;
	auto&	repname	{ Rep::list[rep].name } ;

	if (record.entry().id != RecordID::Curv)
	    abort (format("read_curv {}: bad record ID!", repname)) ;
	record.readrec (stream) ;
	if (stream.fail())
	    ioerror (format("read_curv {}: I/O error!", repname)) ;
	else if (record.size() && blab > 1)
	    cout << "Loaded Curv " << repname << "\n" << flush ;
	}
    }

void Save::read_lagr ()					// Read Lagr records
    {
    if (theory.euclid) return ;

    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;

    for (int rep(0) ; rep < theory.nrep ; ++rep)
	{
	auto&	record	{ global.data().lagr[rep] } ;
	auto&	repname	{ Rep::list[rep].name } ;

	if (record.entry().id != RecordID::Lagr)
	    abort (format("read_lagr {}: bad record ID!", repname)) ;
	record.readrec (stream) ;
	if (stream.fail())
	    ioerror (format("read_lagr {}: I/O error!", repname)) ;
	else if (record.size() && blab > 1)
	    cout << "Loaded Lagr " << repname << "\n" << flush ;
	}
    }

void Save::read_geos ()					// Read Geo records
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto	nbckts	{ global.info().bckt.size() } ;

    if (global.geoswap) return ;
    for (int bcktnum(0) ; bcktnum < nbckts ; ++bcktnum)
	{
	auto&	record	{ global.data().geos[bcktnum] } ;

	if (record.entry().id != RecordID::Geos)
	    abort (format("read_geos [{}]: bad record ID!",bcktnum)) ;
	record.readrec (stream) ;
	if (stream.fail())
	    ioerror (format("read_geos [{}]: I/O error!",bcktnum)) ;
	else if (record.size() && blab > 1)
	    cout << "Loaded Geos [" << bcktnum << "/" << nbckts << "]\n" << flush ;
	}
    }

void Save::read_geo_bckt (int stage, int bcktnum)	// Read Geo bucket
    {
    const auto&			blab	 { Blab::level(Blab::SAVE) } ;
    const auto&			syspath  { global.info(stage).sysfile.path } ;
    thread_local string		mypath   { syspath } ;
    thread_local bool		checkout { false } ;
    thread_local fstream	mystream ;
    thread_local vector<Element> myvec ;

    if (bcktnum >= 0)
	{
	auto	 	nbckts	{ global.info(stage).bckt.size() } ;
	auto& 	 	record	{ global.data(stage).geos[bcktnum] } ;

	if (!global.info().sysfile.stream.is_open())
	    gripe ("Must write or load save file!") ;

	if (mypath != syspath) mystream.close() ;

	if (!mystream.is_open())
	    {
	    alignas (PAGESIZE) thread_local char mybuf [512*1024] ;
	    mystream.rdbuf()->pubsetbuf (mybuf, sizeof mybuf) ;
	    mystream.open (mypath = syspath, ios::in | ios::binary) ;
	    }

	if (record.entry().id != RecordID::Geos)
	    abort (format("read_geo_bckt [{}]: bad record ID!",bcktnum)) ;
	if (record.size())
	    abort (format("read_geo_bckt [{}]: non-empty record!",bcktnum)) ;
	if (!record.entry().reclen)
	    abort (format("Empty geodesic equation bucket {}!",bcktnum)) ;

	if (!checkout)
	    {
	    record.swap (myvec) ;
	    checkout = true ;
	    }
	else abort (format("read_geo_bckt [{}]: bad check out!",bcktnum)) ;

	record.readrec (mystream) ;
	if (mystream.fail())
	    ioerror (format("read_geo_bckt [{}]: I/O error!", bcktnum)) ;
	if (blab > 1)
	    cout << "Loaded Geos [" << bcktnum << "/" << nbckts << "]\n" << flush ;
	}
    else if (checkout)
	{
	myvec.swap (global.data(stage).geos[-bcktnum-1]) ;
	checkout = false ;
	}
    else abort (format("read_geo_bckt [{}]: bad check in!",bcktnum)) ;
    }

void Save::read_stat ()					// Read Stat record
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().sysfile.stream } ;
    auto&	record	{ global.data().stat } ;
    auto	nelem	{ record.entry().nelem } ;
    auto	ptr	{ &global.info().count } ;

    if (record.entry().id != RecordID::Stat)
	abort ("read_stat: bad record ID!") ;
    stream.seekg (record.entry().filepos, ios_base::beg) ;
    stream.read (cast_to<char*>(ptr), nelem * sizeof (Counter)) ;
    if (stream.fail()) ioerror ("read_stat: I/O error!") ;
    if (blab > 1) cout << "Loaded Stat\n" << flush ;
    }

ulong Save::vevsize ()				// Vev record bytes
    {
    return numerics.nvev() * sizeof (real) ;
    }

ulong Save::coupsize ()				// Coupling record bytes
    {
    return Coupling::ncoup() * sizeof (Coupling) ;
    }

ulong Save::datasetsize ()			// coup + vev record size
    {
    return coupsize() + vevsize() ;
    }

Couplings* Save::read_coup (int set, bool update)	// Read Coupling set
    {
    static Couplings	list ;
    const auto&		blab	{ Blab::level(Blab::SAVE) } ;
    auto&		stream	{ global.info().vevfile.stream } ;
    ulong		offset	{ set * datasetsize() } ;

    if (set < 0) stream.seekg (offset, ios_base::end) ;
    else	 stream.seekg (offset + sizeof filehdr, ios_base::beg) ;

    list.resize (Coupling::ncoup()) ;
    stream.read (cast_to<char*>(list.data()), coupsize()) ;
    if (stream.eof())
	{
	stream.clear() ;
	return nullptr ;
	}
    if (stream.fail()) ioerror ("read_coup: I/O error!") ;
    if (update)
	{
	if (Coupling::update (list))
	    {
	    if (blab > 1) cout << "Loaded Coup\n" << flush ;
	    return &Coupling::list ;
	    }
	else gripe ("Incompatible coupling sets in specified data set") ;
	}
    return &list ;
    }

void Save::read_vev (int set)				// Read Vev data
    {
    const auto&	blab	{ Blab::level(Blab::SAVE) } ;
    auto&	stream	{ global.info().vevfile.stream } ;
    long	start	{ global.stage ? numerics.nvevG : 0 } ;
    ulong	offset	{ coupsize() + set * datasetsize() } ;
    auto	ptr	{ numerics.vev.memptr() + start } ;

    if (set < 0) stream.seekg (offset, ios_base::end) ;
    else	 stream.seekg (offset + sizeof filehdr, ios_base::beg) ;

    stream.read (cast_to<char*>(ptr), filehdr.nvev * sizeof (real)) ;
    if (stream.fail()) ioerror ("read_vev: I/O error!") ;
    if (blab > 1) cout << "Loaded Vev\n" << flush ;
    }
