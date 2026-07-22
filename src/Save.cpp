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
static constexpr int	rhdrsiz { sizeof (RecHdr) } ;
alignas (rhdrsiz)
static symb		symbbuf [rhdrsiz] ;		// alignment buffer

void Save::save_sys ()					// Save sys info
    {
    const auto& savedir		{ global.savedir } ;
    auto&	syspath		{ global.info().syspath } ;
    auto&	sysstream	{ global.info().sysstream } ;
    string	file		{ global.mk_filename("sys") } ;
    string	path		{ global.addsubdir (savedir) + file } ;
    auto	mode		{ ios::out | ios::trunc | ios::binary } ;

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

    write_op	   () ;
    write_obs	   () ;
    write_gen 	   () ;
    write_ham	   () ;
    write_grad	   () ;
    write_curv	   () ;
    write_lagr	   () ;
    write_stat	   () ;
    write_geos	   () ;
    write_sysindex () ;
    }

void Save::save_vev ()					// Save vev data
    {
    const auto&	nvev		{ global.info().nobs } ;
    const auto& savedir		{ global.savedir } ;
    auto&	vevpath		{ global.info().vevpath } ;
    auto&	vevstream	{ global.info().vevstream } ;
    string	file		{ global.mk_filename("vev") } ;
    string	path		{ global.addsubdir (savedir) + file } ;
    uint	ncoup		{ Coupling::ncoup() } ;
    uint	blab		{ Blab::level(Blab::SAVE) } ;

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
	else mode |= (global.vevappend ? ios::app : ios::trunc) ;

	stream.open (path, mode) ;			// open for append
	if (stream.is_open())
	    {
	    if (exists && global.vevappend)
		{
		read_header (stream) ;
		if (!filehdr.is_vevfile())
		    gripe (format("Inconsistent vev-data file {}", path)) ;
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
	else gripe ("Cannot write vev-data file " + path) ;
	}
    write_coup () ;
    write_vev  () ;
    }

void Save::write_header (fstream& stream, uint ncoup, uint nvev) // Write file header
    {
    filehdr.version	= global.version ;
    filehdr.name	= global.stage ? theory.name : theory.parent() ;
    filehdr.hashF	= global.stage ? global.info(1).obshash : 0 ;
    filehdr.hashG	= global.info(0).obshash ;
    filehdr.ncoup	= ncoup ;
    filehdr.nvev	= nvev ;
    stream.write (cast_to<char*>(&filehdr), sizeof filehdr) ;
    if (stream.fail()) ioerror ("write_header: I/O error!") ;
    }

void Save::write_op ()						// Write Op record
    {
    auto&	record	{ global.data().op  } ;
    auto&	stream	{ global.info().sysstream } ;
    const auto&	oplist	( global.info().ops ) ;
    uint	nelem	( oplist.size() ) ;

    record.clear() ;
    for (int indx(0) ; indx < nelem ; ++indx)
	{
	const Op& op	{ oplist[indx] } ;
	OpHdr hdr {op.order, (char) op.type, op.primary, 0, (ushort) op.size()} ;

	record.push_back (hdr) ;
	if (hdr.len)
	    {
	    auto	symbptr { op.data() } ;
	    auto	recptr  { cast_to<const RecHdr*> (symbptr) } ;
	    int		n       { hdr.len / rhdrsiz } ;

	    record.insert (record.end(), recptr, recptr + n) ;
	    if (int rem { hdr.len % rhdrsiz })
		{
		std::memcpy (symbbuf, symbptr + n * rhdrsiz, rem) ;
		record.push_back (*cast_to<const RecHdr*>(symbbuf)) ;
		}
	    }
	}
    if (record.entry().id != RecordID::Op)
	fatal ("write_op: bad record ID!") ;

    record.entry().nelem = nelem ;
    record.writerec (stream) ;
    record.free() ;
    if (stream.fail()) ioerror ("write_op: I/O error!") ;
    if (Blab::level(Blab::SAVE)) cout << "Saved Op\n" << flush ;
    }

void Save::write_obs ()					// Write Obs record
    {
    const auto&	nelem	{ global.info().nobs } ;
    auto&	stream	{ global.info().sysstream } ;
    auto&	record	{ global.data().obs  } ;
    uint	start	{ global.stage ? global.info(0).nobs : 0 } ;

    record.clear() ;
    for (uint indx(start) ; indx < start + nelem ; ++indx)
	{
	const Obs& obs { ObsList::obs(indx) } ;
	ObsHdr hdr {obs.corder, obs.xorder, (char) obs.type, 0, (ushort) obs.size()} ;
	record.push_back (hdr) ;
	if (hdr.len)
	    {
	    auto	symbptr { obs.data() } ;
	    auto	recptr  { cast_to<const RecHdr*> (symbptr) } ;
	    int		n       { hdr.len / rhdrsiz } ;

	    record.insert (record.end(), recptr, recptr + n) ;
	    if (int rem { hdr.len % rhdrsiz })
		{
		std::memcpy (symbbuf, symbptr + n * rhdrsiz, rem) ;
		record.push_back (*cast_to<const RecHdr*>(symbbuf)) ;
		}
	    }
	}
    if (record.entry().id != RecordID::Obs)
	fatal ("write_obs: bad record ID!") ;

    record.entry().nelem  = nelem ;
    record.writerec (stream) ;
    record.free() ;
    if (stream.fail()) ioerror ("write_obs: I/O error!") ;
    if (Blab::level(Blab::SAVE)) cout << "Saved Obs\n" << flush ;
    }

void Save::write_gen ()					// Write Gen record
    {
    auto&	record	{ global.data().gen } ;
    auto&	stream	{ global.info().sysstream } ;
    int		optsiz	{ sizeof (OpTerm) / sizeof (RecHdr) } ;
    uint	nelem	( 0 ) ;
    
    record.clear() ;
    for (short rep(0) ; rep < theory.nrep ; ++rep)
	{
	for (const auto& gen : global.info().gens[rep])
	    {
	    GenHdr hdr {rep, gen.order, (char) gen.type, gen.T_odd, (ushort) gen.size()} ;
	    record.push_back (hdr) ;
	    record.push_back (*cast_to<const RecHdr*>(&gen.coeff)) ;
	    auto genptr { cast_to<const RecHdr*> (gen.data()) } ;
	    record.insert (record.end(), genptr, genptr + gen.size() * optsiz) ;
	    ++nelem ;
	    }
	}
    if (record.entry().id != RecordID::Gen)
	fatal ("write_gen: bad record ID!") ;

    record.entry().nelem = nelem ;
    record.writerec (stream) ;
    record.free() ;
    if (stream.fail()) ioerror ("write_gen: I/O error!") ;
    if (Blab::level(Blab::SAVE)) cout << "Saved Gen\n" << flush ;
    }

void Save::write_ham ()					// Write Ham record
    {
    auto&	record	{ global.data().ham } ;
    auto&	stream	{ global.info().sysstream } ;
    auto&	Hterms	{ global.info().Hterms } ;
    ushort	nterms	( Hterms.size() ) ;

    record.clear() ;
    for (ushort i(0) ; i < nterms ; ++i)
	{
	record.add (Poly {HamHdr {i}}, Hterms[i].cpoly) ;
	}
    if (record.entry().id != RecordID::Ham)
	fatal ("write_ham: bad record ID!") ;

    record.entry().nelem = nterms ;
    record.writerec (stream) ;
    record.free() ;
    if (stream.fail()) ioerror ("write_ham: I/O error!") ;
    if (record.size() && Blab::level(Blab::SAVE))
	cout << "Saved Ham\n" << flush ;
    }

void Save::write_grad ()				// Write Grad record
    {
    auto&	stream	{ global.info().sysstream } ;
    auto&	record	{ global.data().grad } ;

    if (record.entry().id != RecordID::Grad)
	fatal ("write_grad: bad record ID!") ;

    record.writerec (stream) ;
    if (stream.fail()) ioerror ("write_grad: I/O error!") ;
    if (record.size() && Blab::level(Blab::SAVE))
	cout << "Saved Grad\n" << flush ;
    }

void Save::write_curv ()				// Write Curv records
    {
    auto&	stream	{ global.info().sysstream } ;

    for (int rep(0) ; rep < global.data().curv.size() ; ++rep)
	{
	auto&	record	{ global.data().curv[rep] } ;
	auto&	repname { Rep::list[rep].name } ;

	if (record.entry().id != RecordID::Curv)
	    fatal (format("write_curv {}: bad record ID!",repname)) ;

	record.writerec (stream) ;
	if (stream.fail())
	    ioerror (format("write_curv {}: I/O error!", repname)) ;
	else if (record.size() && Blab::level(Blab::SAVE))
	    cout << "Saved Curv " << repname << "\n" << flush ;
	}
    }

void Save::write_lagr ()					// Write Lagr records
    {
    auto&	stream	{ global.info().sysstream } ;

    for (int rep(0) ; rep < global.data().lagr.size() ; ++rep)
	{
	auto&	record	{ global.data().lagr[rep] } ;
	auto&	repname { Rep::list[rep].name } ;

	if (record.entry().id != RecordID::Lagr)
	    fatal (format("write_lagr {}: bad record ID!",repname)) ;

	record.writerec (stream) ;
	if (stream.fail())
	    ioerror (format("write_lagr {}: I/O error!",repname)) ;
	else if (record.size() && Blab::level(Blab::SAVE))
	    cout << "Saved Lagr " << Rep::list[rep].name << "\n" << flush ;
	}
    }

void Save::write_geos ()					// Write Geo records
    {
    auto	nbckts	{ global.info().bckt.size() } ;
    int		bcktnum (0) ;

    for (; bcktnum < nbckts ; ++bcktnum)
	{
	auto&	stream	{ global.info().sysstream } ;
	auto&	record	{ global.data().geos[bcktnum] } ;

	if (record.entry().id != RecordID::Geos)
	    fatal ("write_geos: bad record ID!") ;

	record.writerec (stream) ;
	if (stream.fail())
	    ioerror ("write_geos: I/O error!") ;
	else if (record.size() && Blab::level(Blab::SAVE))
	    cout << "Saved Geos [" << bcktnum << "/" << nbckts << "]\n" << flush ;
	}
    for (; bcktnum < global.data().geos.size() ; ++bcktnum)
	{
	auto&	record	{ global.data().geos[bcktnum] } ;
	record.entry().reclen = 0 ;
	}
    }

void Save::write_geo_bckt (int bcktnum)			// Write single Geo bucket
    {
    auto&	stream	{ global.info().sysstream } ;
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
    if (Blab::level(Blab::SAVE))
	cout << "Saved Geos [" << bcktnum << "/" << nbckts << "]\n" << flush ;
    write_sysindex () ;
    rewrite_stat () ;
    if (global.geoswap) record.free() ;
    }

void Save::write_stat ()				// Write Stat record
    {
    auto&	stream	{ global.info().sysstream } ;
    auto&	record	{ global.data().stat } ;
    auto&	count	{ global.info().count } ;
    uint	nelem	( sizeof (count) / sizeof (Counter) ) ;

    if (record.entry().id != RecordID::Stat)
	fatal ("write_stat: bad record ID!") ;

    record.entry().nelem   = nelem ;
    record.entry().reclen  = nelem * sizeof (Counter) / sizeof (RecHdr) ;
    record.entry().filepos = static_cast<std::streamoff>(stream.tellp()) ;

    stream.write (cast_to<const char*>(&count), nelem * sizeof (Counter)) ;
    if (stream.fail()) ioerror ("write_stat: I/O error!") ;
    if (Blab::level(Blab::SAVE)) cout << "Saved Stat\n" << flush ;
    }

void Save::rewrite_stat ()				// Rewrite Stat record
    {
    auto&	stream	{ global.info().sysstream } ;
    auto&	record	{ global.data().stat } ;
    auto&	count	{ global.info().count } ;
    uint	nelem	( sizeof (count) / sizeof (Counter) ) ;

    if (record.entry().id != RecordID::Stat)
	fatal ("rewrite_stat: bad record ID!") ;

    stream.seekg (record.entry().filepos, ios_base::beg) ;
    stream.write (cast_to<const char*>(&count), nelem * sizeof (Counter)) ;
    if (stream.fail()) ioerror ("rewrite_stat: I/O error!") ;
    if (Blab::level(Blab::SAVE)) cout << "Rewrote Stat\n" ;
    }

void Save::write_sysindex ()				// Write SysIndex
    {
    auto&	stream	{ global.info().sysstream } ;
    auto&	syspath	{ global.info().syspath   } ;
    auto&	index	{ global.info().sysindex  } ;
    auto&	last	{ index.back() } ;

    last.id	 = RecordID::Indx ;
    last.nelem   = 1 ;
    last.reclen  = sizeof (index) / sizeof (RecHdr) ;
    last.filepos = static_cast<std::streamoff>(stream.tellp()) ;
    
    stream.write (cast_to<const char*>(index.data()), sizeof index) ;
    stream.flush() ;
    if (stream.fail()) ioerror ("write_sysindex: I/O error!") ;
    std::filesystem::resize_file (syspath,stream.tellp()) ;
    if (Blab::level(Blab::SAVE) > 1) cout << "Wrote SysIndex\n" ;
    }

void Save::write_coup ()				// Write Coupling's
    {
    auto&	stream	{ global.info().vevstream } ;
    auto	ptr	{ Coupling::list.data() } ;

    stream.seekp (0, ios::end) ;
    stream.write (cast_to<char*>(ptr), coupsize()) ;
    if (stream.fail()) ioerror ("write_coup: I/O error!") ;
    if (Blab::level(Blab::SAVE)) cout << "Wrote Couplings\n" << flush ;
    }

void Save::write_vev ()					// Write Vev's
    {
    uint	start	{ global.stage ? global.info(0).nobs : 0 } ;
    auto&	stream	{ global.info().vevstream } ;
    auto	ptr	{ memptr(numerics.vev) + start } ;

    stream.seekp (0, ios::end) ;
    stream.write (cast_to<const char*>(ptr), vevsize()) ;
    if (stream.fail()) ioerror ("write_vev: I/O error!") ;
    if (Blab::level(Blab::SAVE)) cout << "Wrote Vevs\n" << flush ;
    }

void Save::load_save (int set, string file)		// Load save file
    {
    const auto& savedir { global.savedir } ;
    string	path	{ global.addsubdir (savedir) + file } ;

    if (!std::filesystem::exists (path))
	path = global.addsubdir (savedir, !global.stage) + file ;

    if (!std::filesystem::exists (path))
	gripe (format("Cannot find save file {}", file)) ;

    fstream stream { path, ios::in | ios::out | ios::binary } ;
    if (stream.is_open())
	{
	global.stageinit (read_header (stream)) ;
	if (filehdr.is_sysfile())			// sys info file
	    {
	    if (set >= 0) gripe ("No set numbers in sys files") ;

	    cout << "Loading sys info from " << path << "\n" ;
	    global.info().syspath   = std::move (path) ;
	    global.info().sysstream = std::move (stream) ;

	    read_sysindex () ;
	    read_op   () ;
	    read_obs  () ;
	    read_gen  () ;
	    read_ham  () ;
	    read_grad () ;
	    read_curv () ;
	    read_lagr () ;
	    read_stat () ;
	    read_geos () ;
	    }
	else /* filehdr.is_vevfile() */		// vev data file
	    {
	    cout << "Loading vev data from " << path ;
	    if (set < 0)	cout << "\n" ;
	    else		cout << ", set " << set << "\n" ;

	    global.info().vevpath   = std::move (path) ;
	    global.info().vevstream = std::move (stream) ;
	    load_vev (set) ;
	    }
	}
    else gripe (format("Cannot open save file {}", path)) ;
    }

void Save::load_vev (int set)				// Load vev data set
    {
    if (global.info(global.stage).vevstream.is_open())
	{
	if (read_coup (set,nullptr)) read_vev (set) ;
	else gripe (format("Cannot load vev set {}", set)) ;
	}
    else gripe ("No open vev file") ;
    }

int Save::read_header (fstream& stream)		// Read save file header
    {
    const auto& hashG	{ global.info(0).obshash } ;
    const auto& hashF	{ global.info(1).obshash } ;
    const auto& nobsG	{ global.info(0).nobs } ;
    const auto& nobsF	{ global.info(1).nobs } ;
    auto&	hdr	{ filehdr } ;
    auto	badobs	{ format("Inconsistent Obs set in {}", hdr.name.data()) } ;

    stream.read (cast_to<char*>(&hdr), sizeof hdr) ;
    if (stream.fail()) ioerror ("read_header: I/O error!") ;

    if (hdr.version.incompat()) 
	gripe (format("Inconsistent save file {}: wrong real num type",
	    hdr.name.data())) ;
    if (hdr.version.newer())
	cout << format("Warning: save file {} from newer program version",
	    hdr.name.data()) ;

    if (theory.parent() == hdr.name)				// YM save file
	{
	if (!hdr.ncoup && !hdr.nvev && !hdr.hashF)
	    {
	    if (!hdr.hashF) return 0 ;				// YM sys file
	    else gripe (badobs) ;
	    }
	if (hdr.ncoup == Coupling::ncoup(0) && hdr.nvev  == nobsG)
	    {
	    if (hdr.hashG == hashG && !hdr.hashF) return 0 ;	// YM vev file
	    else gripe (badobs) ;
	    }
	}
    else if (theory.name == hdr.name && theory.nf)		// QCD save file
	{
	if (hdr.ncoup == 0 && hdr.nvev == 0)
	    {
	    if (hdr.hashG == hashG) return 1 ;			// QCD sys file
	    else gripe (badobs) ;
	    }
	if (hdr.ncoup == Coupling::ncoup(1) && hdr.nvev == nobsF)
	    {
	    if (hdr.hashG == hashG && hdr.hashF == hashF) return 1 ; // QCD vev file
	    else gripe (badobs) ;
	    }
	}
    gripe (format("Inappropriate save file from theory {}", hdr.name.data())) ;
    }

void Save::read_sysindex ()					// Load SysIndex
    {
    auto&	stream	{ global.info().sysstream } ;
    auto&	index	{ global.info().sysindex } ;
    auto	ptr	{ index.data() } ;

    stream.seekg (-sizeof index, ios_base::end) ;
    stream.read (cast_to<char*>(ptr), sizeof index) ;
    if (stream.fail()) ioerror ("read_sysindex: I/O error!") ;
    if (Blab::level(Blab::SAVE)) cout << "Read Index\n" ;
    }

void Save::read_op ()						// Read Op record
    {
    auto&	stream	{ global.info().sysstream } ;
    auto&	record	{ global.data().op } ;

    if (record.entry().id != RecordID::Op)
	gripe ("read_op: bad record ID!") ;

    record.readrec (stream) ;
    if (stream.fail()) ioerror ("read_op: I/O error!") ;

    uint	indx	{ 0 } ;
    uint	nop	{ record.entry().nelem } ;
    RecHdr*	recptr	{ record.data() } ;
    RecHdr*	recend	{ recptr + record.size() } ;
    OpList&	oplist	{ global.info().ops } ;

    oplist.clear () ;
    oplist.reserve (nop) ;
    while (recptr < recend)
	{
	OpHdr	hdr	( *recptr++ ) ;
	symb*	ptr	{ cast_to<symb*> (recptr) } ;
	Str	s	{ ptr, ptr + hdr.len } ;
	Op	o	{ s, hdr } ;
	auto	k	{ oplist.store(o) } ;

	if (k != indx++) gripe ("read_op: Inconsistent OpList") ;

	recptr += 1 + (hdr.len - 1) / sizeof (RecHdr) ;
	}
    record.free() ;
    if (indx != nop)
	gripe ("read_op: Inconsistent save record!") ;

    if (Blab::level(Blab::SAVE)) cout << "Loaded Op\n" << flush ;
    }

void Save::read_obs ()					// Read Obs record
    {
    auto&	stream	{ global.info().sysstream } ;
    auto&	record	{ global.data().obs } ;
    auto&	obslist	{ ObsList::obs } ;

    if (record.entry().id != RecordID::Obs)
	gripe ("read_obs: bad record ID!") ;
    record.readrec (stream) ;
    if (stream.fail())
	ioerror ("read_obs: I/O error!") ;

    uint	start	{ global.stage ? global.info(0).nobs : 0 } ;
    uint	indx    { start } ;
    uint	nobs	{ record.entry().nelem } ;
    RecHdr*	recptr	{ record.data() } ;
    RecHdr*	recend	{ recptr + record.size() } ;

    obslist.purge   (start) ;
    obslist.reserve (start + nobs) ;
    while (recptr < recend)
	{
	ObsHdr	hdr	( *recptr++ ) ;
	symb*	ptr	{ cast_to<symb*> (recptr) } ;
	Str	s	{ ptr, ptr + hdr.len } ;
	Obs	o	{ s, hdr } ;

	if (obslist.store (o) != indx++)
	    gripe ("read_obs: Inconsistent ObsList::obs") ;

	if (global.info().maxord < o.order() && o.corder == o.xorder)
	    global.info().maxord = o.order() ;

	recptr += 1 + (hdr.len - 1) / sizeof (RecHdr) ;
	}
    record.free() ;
    if (indx - start != nobs)
	gripe ("read_obs: Inconsistent save record!") ;
    if (global.info(0).nobs + global.info(1).nobs != obslist.size())
	abort ("read_obs: Inconsistent ObsList size") ;

    if (Blab::level(Blab::SAVE)) cout << "Loaded Obs\n" << flush ;
    global.mk_bcktlist  () ;
    numerics.initialize () ;
    Canon::cache.reload() ;
    global.info(0).obshash = filehdr.hashG ;
    global.info(1).obshash = filehdr.hashF ;
    }

void Save::read_gen ()					// Read Gen record
    {
    auto&	stream	{ global.info().sysstream } ;
    auto&	record	{ global.data().gen } ;
    int		optsiz	{ sizeof (OpTerm) / sizeof (RecHdr) } ;

    if (record.entry().id != RecordID::Gen)
	gripe ("read_gen: bad record ID!") ;
    record.readrec (stream) ;
    if (stream.fail())
	ioerror ("read_gen: I/O error!") ;

    RecHdr*	recptr	{ record.data() } ;
    RecHdr*	recend	{ recptr + record.size() } ;
    uint	nelem   ( 0 ) ;

    for (auto& gens : global.info().gens)  gens.clear() ;
    for (auto& even : global.info().neven) even = 0 ;

    while (recptr < recend)
	{
	GenHdr	hdr	( *recptr++ ) ;
	doub	coeff	{ *cast_to<doub*> (recptr++) } ;
	OpTerm*	ptr	{ cast_to<OpTerm*> (recptr) } ;
	OpList& list	{ global.info().ops } ;
	OpSum	s	( ptr, ptr + hdr.len, list ) ;
	Gen	gen	{ s, hdr, coeff } ;

	global.info().gens[hdr.rep].push_back (gen) ;

	if (!gen.T_odd) ++global.info().neven[hdr.rep] ;
	if (global.info().maxgen < gen.order)
	    global.info().maxgen = gen.order ;
	++nelem ;
	recptr += hdr.len * optsiz ;
	}
    record.free() ;

    if (nelem != record.entry().nelem)
	gripe ("read_gen: Inconsistent save record!") ;

    if (Blab::level(Blab::SAVE)) cout << "Loaded Gen\n" << flush ;
    }

void Save::read_ham ()				// Read Ham record
    {
    auto&	stream	{ global.info().sysstream } ;
    auto& 	record	{ global.data().ham } ;
    auto&	Hterms	{ global.info().Hterms } ;

    if (record.entry().id != RecordID::Ham)
	gripe ("read_ham: bad record ID!") ;
    record.readrec (stream) ;
    if (stream.fail())
	ioerror ("read_ham: I/O error!") ;

    if (Hterms.size() != record.entry().nelem)
	gripe ("read_ham: Inconsistent save record!") ;

    for (const auto& poly : record)
	{
	const HamHdr& info ( poly ) ;
	auto& Hterm { Hterms[info.Hterm] } ;
	Hterm.cpoly.clear() ;
	Hterm.cpoly.add (poly) ;
	}
    record.free() ;
    if (Blab::level(Blab::SAVE)) cout << "Loaded Ham\n" << flush ;
    }

void Save::read_grad ()				// Read Grad record
    {
    auto&	stream	{ global.info().sysstream } ;
    auto& 	record	{ global.data().grad } ;

    if (record.entry().id != RecordID::Grad)
	gripe ("read_grad: bad record ID!") ;
    record.readrec (stream) ;
    if (stream.fail())
	ioerror ("read_grad: I/O error!") ;
    else if (record.size() && Blab::level(Blab::SAVE))
	cout << "Loaded Grad\n" << flush ;
    }

void Save::read_curv ()				// Read Curv records
    {
    auto&	stream	{ global.info().sysstream } ;

    for (int rep(0) ; rep < theory.nrep ; ++rep)
	{
	auto&	record	{ global.data().curv[rep] } ;
	auto&	repname	{ Rep::list[rep].name } ;

	if (record.entry().id != RecordID::Curv)
	    gripe (format("read_curv {}: bad record ID!", repname)) ;

	record.readrec (stream) ;
	if (stream.fail())
	    ioerror (format("read_curv {}: I/O error!", repname)) ;
	else if (record.size() && Blab::level(Blab::SAVE))
	    cout << "Loaded Curv " << repname << "\n" << flush ;
	}
    }

void Save::read_lagr ()					// Read Lagr records
    {
    auto&	stream	{ global.info().sysstream } ;

    for (int rep(0) ; rep < theory.nrep ; ++rep)
	{
	auto&	record	{ global.data().lagr[rep] } ;
	auto&	repname	{ Rep::list[rep].name } ;

	if (record.entry().id != RecordID::Lagr)
	    gripe (format("read_lagr {}: bad record ID!", repname)) ;

	record.readrec (stream) ;
	if (stream.fail())
	    ioerror (format("read_lagr {}: I/O error!", repname)) ;
	else if (record.size() && Blab::level(Blab::SAVE))
	    cout << "Loaded Lagr " << repname << "\n" << flush ;
	}
    }

void Save::read_geos ()					// Read Geo records
    {
    auto&	stream	{ global.info().sysstream } ;
    auto	nbckts	{ global.info().bckt.size() } ;

    if (global.geoswap) return ;
    for (int bcktnum(0) ; bcktnum < nbckts ; ++bcktnum)
	{
	auto&	record	{ global.data().geos[bcktnum] } ;

	if (record.entry().id != RecordID::Geos)
	    gripe (format("read_geos [{}]: bad record ID!",bcktnum)) ;

	record.readrec (stream) ;
	if (stream.fail())
	    ioerror (format("read_geos [{}]: I/O error!",bcktnum)) ;
	else if (record.size() && Blab::level(Blab::SAVE))
	    cout << "Loaded Geos [" << bcktnum << "/" << nbckts << "]\n" << flush ;
	}
    }

void Save::read_geo_bckt (int bcktnum)			// Read Geo bucket
    {
    const auto&				syspath  { global.info().syspath } ;
    thread_local string			mypath   { syspath } ;
    thread_local bool			checkout { false } ;
    thread_local fstream		mystream ;
    thread_local vector<RecHdr>		myvec ;

    if (bcktnum >= 0)
	{
	auto	 	nbckts	{ global.info().bckt.size() } ;
	auto& 	 	record	{ global.data().geos[bcktnum] } ;
//	cout << "read_geo_bckt: bcktnum " << bcktnum << " stage " << global.stage << " nbckts " << nbckts << " ncol " << record.entry().ncol << " nrow " << record.entry().nrow << "\n" << flush ;

	if (!global.info().sysstream.is_open())
	    gripe ("Must write or load save file!") ;

	if (mypath != syspath) mystream.close() ;

	if (!mystream.is_open())
	    {
	    alignas (PAGESIZE) thread_local char mybuf [512*1024] ;
	    mystream.rdbuf()->pubsetbuf (mybuf, sizeof mybuf) ;
	    mystream.open (mypath = syspath, ios::in | ios::binary) ;
	    }

	if (record.entry().id != RecordID::Geos)
	    gripe (format("read_geo_bckt [{}]: bad record ID!",bcktnum)) ;
	if (record.size())
	    gripe (format("read_geo_bckt [{}]: non-empty record!",bcktnum)) ;
	if (!record.entry().reclen)
	    gripe (format("Empty geodesic equation bucket {}!",bcktnum)) ;

	if (!checkout)
	    {
	    record.swap (myvec) ;
	    checkout = true ;
	    }
	else gripe (format("read_geo_bckt [{}]: bad check out!",bcktnum)) ;

	record.readrec (mystream) ;
	if (mystream.fail())
	    ioerror (format("read_geo_bckt [{}]: I/O error!", bcktnum)) ;
	if (Blab::level(Blab::SAVE))
	    cout << "Loaded Geos [" << bcktnum << "/" << nbckts << "]\n" << flush ;
	}
    else if (checkout)
	{
	myvec.swap (global.data().geos[-bcktnum-1]) ;
	checkout = false ;
	}
    else gripe (format("read_geo_bckt [{}]: bad check in!",bcktnum)) ;
    }

void Save::read_stat ()					// Read Stat record
    {
    auto&	stream	{ global.info().sysstream } ;
    auto&	record	{ global.data().stat } ;
    uint	nelem	{ record.entry().nelem } ;
    auto	ptr	{ &global.info().count } ;

    if (record.entry().id != RecordID::Stat)
	gripe ("read_stat: bad record ID!") ;

    stream.seekg (record.entry().filepos, ios_base::beg) ;
    stream.read (cast_to<char*>(ptr), nelem * sizeof (Counter)) ;
    if (stream.fail()) ioerror ("read_stat: I/O error!") ;
    if (Blab::level(Blab::SAVE)) cout << "Loaded Stat\n" << flush ;
    }

ulong Save::vevsize ()				// Vev record bytes
    {
    return global.info().nobs * sizeof (real) ;
    }

ulong Save::coupsize ()				// Coupling record bytes
    {
    return Coupling::ncoup() * sizeof (Coupling) ;
    }

ulong Save::cvsetsize ()			// coup + vev record size
    {
    return coupsize() + vevsize() ;
    }

bool Save::read_coup (int set, Couplings* ptr)		// Read Coupling set
    {
    auto&	stream	{ global.info().vevstream } ;
    auto&	list	{ ptr ? *ptr : Coupling::list } ;
    ulong	offset	{ set * cvsetsize() } ;

    if (set < 0) stream.seekg (offset, ios_base::end) ;
    else	 stream.seekg (offset + sizeof filehdr, ios_base::beg) ;

    stream.read (cast_to<char*>(list.data()), coupsize()) ;
    if (stream.eof()) { stream.clear() ; return false ; }
    if (stream.fail()) ioerror ("read_coup: I/O error!") ;
    if (Blab::level(Blab::SAVE)) cout << "Loaded Coup\n" << flush ;
    return true ;
    }

void Save::read_vev (int set)				// Read Vev data
    {
    auto&	stream	{ global.info().vevstream } ;
    uint	start	{ global.stage ? global.info(0).nobs : 0 } ;
    ulong	offset	{ coupsize() + set * cvsetsize() } ;
    auto	ptr	{ memptr(numerics.vev) + start } ;

    if (set < 0) stream.seekg (offset, ios_base::end) ;
    else	 stream.seekg (offset + sizeof filehdr, ios_base::beg) ;

    stream.read (cast_to<char*>(ptr), filehdr.nvev * sizeof (real)) ;
    if (stream.fail()) ioerror ("read_vev: I/O error!") ;
    if (Blab::level(Blab::SAVE)) cout << "Loaded Vev\n" << flush ;
    }
