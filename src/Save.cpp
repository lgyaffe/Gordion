#include "Save.h"
#include "Rep.h"
#include "Canon.h"
#include "Numerics.h"
#include "Blab.h"
#include "Gripe.h"
#include <filesystem>

using std::ios_base ;

static constexpr int		rhdrsiz = sizeof (RecHdr) ;
alignas (rhdrsiz) static symb	symbbuf [rhdrsiz] ;	// alignment buffer

void Save::save_sys (string sysfile)			// Save sys info
    {
    if (sysfile.size())
	{
	auto pos { sysfile.find_last_of ('/') } ;
	if (pos != sysfile.npos)			// separate directory
	    {
	    global.savedir.assign (sysfile, 0, pos) ;
	    sysfile.erase (0, pos+1) ;
	    }
	}
    else sysfile = Global::dfltfilename("sys") ;

    string	path	{ global.savedir + "/" + sysfile } ;
    fstream	stream ;

    if (global.sysstream.is_open()) global.sysstream.close() ;
    stream.open (path, ios::out | ios::trunc | ios::binary) ;
    if (stream.is_open())
	{
	cout << "Writing sys info file " << path << "\n" ;
	global.sysfile   = sysfile ;
	global.sysstream = std::move (stream) ;
	syspath = path ;
	write_header (global.sysstream) ;
	for (int stage(0) ; stage < 2 ; ++stage)
	    {
	    write_op	(stage) ;
	    write_obs	(stage) ;
	    write_gen 	(stage) ;
	    write_grad	(stage) ;
	    write_curv	(stage) ;
	    write_lagr	(stage) ;
	    write_stat	(stage) ;
	    write_geos	(stage) ;
	    if (!theory.nf) break ;
	    }
	write_sysindex  () ;
	}
    else gripe (format("Cannot write sys info file {}", path)) ;
    }

void Save::save_vev (string vevfile)			// Save vev data
    {
    alignas (16*1024) char mybuf [512*1024] ;

    if (Blab::blablevel[BLAB::SAVE]) cout << "Saving vev's\n" << flush ;

    if (vevfile.size())
	{
	auto pos { vevfile.find_last_of ('/') } ;
	if (pos != vevfile.npos)			// separate directory
	    {
	    global.savedir.assign (vevfile, 0, pos) ;
	    vevfile.erase (0, pos+1) ;
	    }
	global.vevstream.close() ;
	}
    else vevfile = Global::dfltfilename("vev") ;

    if (!global.vevstream.is_open())
	{
	string	path	{ global.savedir + "/" + vevfile } ;
	bool	exists	{ std::filesystem::exists (path) } ;
	auto	mode	{ ios::in | ios::out | ios::binary } ;
	fstream	stream ;
	stream.rdbuf()->pubsetbuf (mybuf, sizeof mybuf) ;

	if (!exists)			// touch file into existence
	    {
	    stream.open (path, ios::out | ios::binary) ;
	    stream.close() ;
	    }
	else mode |= (global.vevappend ? ios::app : ios::trunc) ;

	stream.open (path, mode) ;
	if (stream.is_open())
	    {
	    if (exists && global.vevappend)
		{
		if (read_header (stream,true) && filehdr.is_vevfile())
		    cout << "Appending vev data to " << path << "\n" ;
		else
		    gripe (format("Inappropriate vev data file {}", path)) ;
		}
	    else
		{
		cout << "Writing vev data to " << path << "\n" ;
		write_header (stream, Coupling::list.size(), global.nobs()) ;
		}
	    global.vevfile   = vevfile ;
	    global.vevstream = std::move (stream) ;
	    }
	else gripe (format("Cannot write vev data file {}", path)) ;
	}
    write_coup () ;
    write_vev  () ;
    }

void Save::write_header (fstream& stream, int ncoup, int nvev)	// Write file header
    {
    filehdr.name    = theory.name ;
    filehdr.version = global.version ;
    filehdr.ncoup   = ncoup ;
    filehdr.nvev    = nvev ;
    stream.write (cast_to<char*>(&filehdr), sizeof filehdr) ;
    if (stream.fail()) ioerror ("write_header: I/O error!") ;
    }

void Save::write_op (int stage)				// Write Op record
    {
    auto&	stream	{ global.sysstream } ;
    auto&	record	{ global.data(stage).op  } ;
    uint	nelem	{ global.info(stage).nop } ;
    uint	start	{ stage ? global.nopG() : 0 } ;

    record.clear() ;
    for (int indx(start) ; indx < start + nelem ; ++indx)
	{
	const Op& op	{ Op::list[indx] } ;
	OpHdr hdr {op.order, (char) op.type, op.primary, 0, (ushort) op.size()} ;

	record.emplace_back (hdr) ;
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
    if (Blab::blablevel[BLAB::SAVE]) cout << "Saved Op\n" << flush ;
    }

void Save::write_obs (int stage)			// Write Obs record
    {
    auto&	stream	{ global.sysstream } ;
    auto&	record	{ global.data(stage).obs  } ;
    uint	nelem	{ global.info(stage).nobs } ;
    uint	start	{ stage ? global.nobsG() : 0 } ;

    record.clear() ;
    for (uint indx(start) ; indx < start + nelem ; ++indx)
	{
	const Obs& obs { ObsList::obs(indx) } ;
	ObsHdr hdr {obs.corder, obs.xorder, (char) obs.type, 0, (ushort) obs.size()} ;
	record.emplace_back (hdr) ;
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
    if (Blab::blablevel[BLAB::SAVE]) cout << "Saved Obs\n" << flush ;
    }

void Save::write_gen (int stage)			// Write Gen record
    {
    auto&	stream	{ global.sysstream } ;
    auto&	record	{ global.data(stage).gen } ;
    int		optsiz	{ sizeof (OpTerm) / sizeof (RecHdr) } ;
    uint	nelem	( 0 ) ;
    
    record.clear() ;
    for (short rep(0) ; rep < theory.nrep ; ++rep)
	{
	for (const auto& gen : global.info(stage).gens[rep])
	    {
	    GenHdr hdr {rep, gen.order, (char) gen.type, gen.T_odd, (ushort) gen.size()} ;
	    record.emplace_back (hdr) ;
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
    if (Blab::blablevel[BLAB::SAVE]) cout << "Saved Gen\n" << flush ;
    }

void Save::write_grad (int stage)			// Write Grad record
    {
    auto&	stream	{ global.sysstream } ;
    auto&	record	{ global.data(stage).grad } ;

    if (record.entry().id != RecordID::Grad)
	fatal ("write_grad: bad record ID!") ;

    record.writerec (stream) ;
    if (stream.fail()) ioerror ("write_grad: I/O error!") ;
    if (record.size() && Blab::blablevel[BLAB::SAVE])
	cout << "Saved Grad\n" << flush ;
    }

void Save::write_curv (int stage)			// Write Curv records
    {
    for (int rep(0) ; rep < theory.nrep ; ++rep)
	{
	auto&	stream	{ global.sysstream } ;
	auto&	record	{ global.data(stage).curv[rep] } ;
	auto&	repname { Rep::list[rep].name } ;

	if (record.entry().id != RecordID::Curv)
	    fatal (format("write_curv {}: bad record ID!",repname)) ;

	record.writerec (stream) ;
	if (stream.fail())
	    ioerror (format("write_curv {}: I/O error!", repname)) ;
	else if (record.size() && Blab::blablevel[BLAB::SAVE])
	    cout << "Saved Curv " << repname << "\n" << flush ;
	}
    }

void Save::write_lagr (int stage)			// Write Lagr records
    {
    for (int rep(0) ; rep < global.data(stage).lagr.size() ; ++rep)
	{
	auto&	stream	{ global.sysstream } ;
	auto&	record	{ global.data(stage).lagr[rep] } ;
	auto&	repname { Rep::list[rep].name } ;

	if (record.entry().id != RecordID::Lagr)
	    fatal (format("write_lagr {}: bad record ID!",repname)) ;

	record.writerec (stream) ;
	if (stream.fail())
	    ioerror (format("write_lagr {}: I/O error!",repname)) ;
	else if (record.size() && Blab::blablevel[BLAB::SAVE])
	    cout << "Saved Lagr " << Rep::list[rep].name << "\n" << flush ;
	}
    }

void Save::write_geos (int stage)			// Write Geo records
    {
    int	 bcktnum (0) ;
    auto nbckts	{ global.info(stage).bckt.size() } ;
    for (; bcktnum < global.info(stage).bckt.size() ; ++bcktnum)
	{
	auto&	stream	{ global.sysstream } ;
	auto&	record	{ global.data(stage).geos[bcktnum] } ;

	if (record.entry().id != RecordID::Geos)
	    fatal ("write_geos: bad record ID!") ;

	record.writerec (stream) ;
	if (stream.fail())
	    ioerror ("write_geos: I/O error!") ;
	else if (record.size() && Blab::blablevel[BLAB::SAVE])
	    cout << "Saved Geos [" << bcktnum << "/" << nbckts << "]\n" << flush ;
	}
    for (; bcktnum < global.data(stage).geos.size() ; ++bcktnum)
	{
	auto&	record	{ global.data(stage).geos[bcktnum] } ;
	record.entry().reclen = 0 ;
	}
    }

void Save::write_geo_bckt (int bcktnum)			// Write single Geo bucket
    {
    auto&	stream	{ global.sysstream } ;
    auto&	record	{ global.data().geos[bcktnum] } ;
    auto	nbckts	{ global.info().bckt.size() } ;
    int		indxsiz	{ global.sysindex.size() * sizeof (RecIndx) } ;

    if (global.interrupt) return ;
    if (record.entry().id != RecordID::Geos)
	fatal (format("write_geos [{}]: bad record ID!",bcktnum)) ;

    std::lock_guard<std::mutex> lock (savemutex) ;
    stream.seekg (-indxsiz, ios_base::end) ;
    record.writerec (stream) ;
    if (stream.fail())
	ioerror (format("write_geo_bckt [{}]: I/O error!",bcktnum)) ;
    if (Blab::blablevel[BLAB::SAVE])
	cout << "Saved Geos [" << bcktnum << "/" << nbckts << "]\n" << flush ;
    write_sysindex () ;
    rewrite_stat () ;
    if (global.geoswap) record.free() ;
    }

void Save::write_stat (int stage)			// Write Stat record
    {
    auto&	stream	{ global.sysstream } ;
    auto&	record	{ global.data(stage).stat } ;
    auto&	count	{ global.info(stage).count } ;
    uint	nelem	( sizeof (count) / sizeof (Counter) ) ;

    if (record.entry().id != RecordID::Stat)
	fatal ("write_stat: bad record ID!") ;

    record.entry().nelem   = nelem ;
    record.entry().reclen  = nelem * sizeof (Counter) / sizeof (RecHdr) ;
    record.entry().filepos = static_cast<std::streamoff>(stream.tellp()) ;

    stream.write (cast_to<const char*>(&count), nelem * sizeof (Counter)) ;
    if (stream.fail()) ioerror ("write_stat: I/O error!") ;
    if (Blab::blablevel[BLAB::SAVE]) cout << "Saved Stat\n" << flush ;
    }

void Save::rewrite_stat ()				// Rewrite Stat record
    {
    auto&	stream	{ global.sysstream } ;
    auto&	record	{ global.data().stat } ;
    auto&	count	{ global.info().count } ;
    uint	nelem	( sizeof (count) / sizeof (Counter) ) ;

    if (record.entry().id != RecordID::Stat)
	fatal ("rewrite_stat: bad record ID!") ;

    stream.seekg (record.entry().filepos, ios_base::beg) ;
    stream.write (cast_to<const char*>(&count), nelem * sizeof (Counter)) ;
    if (stream.fail()) ioerror ("rewrite_stat: I/O error!") ;
    if (Blab::blablevel[BLAB::SAVE]) cout << "Rewrote Stat\n" ;
    }

void Save::write_sysindex ()				// Write SysIndex
    {
    auto&	stream	{ global.sysstream } ;
    auto&	index	{ global.sysindex } ;
    auto&	last	{ index.back() } ;

    last.id	 = RecordID::Indx ;
    last.nelem   = 1 ;
    last.reclen  = sizeof (index) / sizeof (RecHdr) ;
    last.filepos = static_cast<std::streamoff>(stream.tellp()) ;
    
    stream.write (cast_to<const char*>(index.data()), sizeof index) ;
    stream.flush() ;
    if (stream.fail()) ioerror ("write_sysindex: I/O error!") ;
    std::filesystem::resize_file (syspath,stream.tellp()) ;
    if (Blab::blablevel[BLAB::SAVE] > 1) cout << "Wrote SysIndex\n" ;
    }

void Save::write_coup ()				// Write Coupling's
    {
    auto&	stream	{ global.vevstream } ;
    uint	ncoup	( Coupling::list.size() ) ;
    auto	ptr	{ Coupling::list.data() } ;

    stream.seekp (0, ios::end) ;
    stream.write (cast_to<const char*>(ptr), ncoup * sizeof (Coupling)) ;

    if (stream.fail()) ioerror ("write_coup: I/O error!") ;
    if (Blab::blablevel[BLAB::SAVE]) cout << "Wrote Couplings\n" << flush ;
    }

void Save::write_vev ()					// Write Vev's
    {
    auto&	stream	{ global.vevstream } ;
    uint	nvev	( numerics.vev.n_elem ) ;
    auto	ptr	{ numerics.vev.memptr() } ;

    stream.seekp (0, ios::end) ;
    stream.write (cast_to<const char*>(ptr), nvev * sizeof (real)) ;
    if (stream.fail()) ioerror ("write_vev: I/O error!") ;
    if (Blab::blablevel[BLAB::SAVE]) cout << "Wrote Vevs\n" << flush ;
    }

void Save::load_save (int set, string file)		// Load save file
    {
    auto pos { file.find_last_of ('/') } ;
    if (pos != file.npos)				// separate directory
	{
	global.savedir.assign (file, 0, pos) ;
	file.erase (0, pos+1) ;
	}
    string  path   { global.savedir + "/" + file } ;
    fstream stream { path, ios::in | ios::out | ios::binary } ;
    if (stream.is_open())
	{
	if (read_header (stream,false))
	    {
	    if (filehdr.is_sysfile())			// sys info file
		{
		cout << "Loading sys info from " << path << "\n" ;
		global.sysfile = file ;
		global.sysstream = std::move (stream) ;
		syspath = path ;

		read_sysindex () ;
		for (int stage(0) ; stage < 2 ; ++stage)
		    {
		    read_op   (stage) ;
		    read_obs  (stage) ;
		    read_gen  (stage) ;
		    read_grad (stage) ;
		    read_curv (stage) ;
		    read_lagr (stage) ;
		    read_stat (stage) ;
		    read_geos (stage) ;
		    if (!theory.nf) break ;
		    }
		}
	    else /* filehdr.is_vevfile() */		// vev data file
		{
		cout << "Loading vev data from " << path ;
		if (set < 0)	cout << "\n" ;
		else		cout << ", set " << set << "\n" ;

		global.vevfile = file ;
		global.vevstream = std::move (stream) ;
		load_vev (set) ;
		}
	    }
	else gripe (format("Inappropriate save file {}", path)) ;
	}
    else gripe (format("Cannot open save file {}", path)) ;
    }

void Save::load_vev (int set)				// Load vev data set
    {
    if (global.vevstream.is_open())
	{
	if (read_coup (set,nullptr)) read_vev (set) ;
	else gripe (format("Cannot load vev set {}", set)) ;
	}
    else gripe ("No open vev file") ;
    }

bool Save::read_header (fstream& stream, bool write)	// Read save file header
    {
    auto& hdr	{ filehdr } ;
    auto  nvev  { global.nobs() } ;
    auto  ncoup { Coupling::list.size() } ;

    stream.read (cast_to<char*>(&hdr), sizeof hdr) ;
    if (stream.fail()) ioerror ("read_header: I/O error!") ;

    if (hdr.version.incompat()) 
	gripe (format("Inconsistent save file {}: wrong real num type", hdr.name.data())) ;
    if (hdr.version.newer())
	cout << format("Warning: save file {} from newer program version", hdr.name.data()) ;

    if (theory.name == hdr.name || theory.parent() == hdr.name)
	{
	if (hdr.ncoup == 0 && hdr.nvev == 0)	    return true ; // sys file
	if (hdr.ncoup == ncoup && hdr.nvev == nvev) return true ; // vev file

	return !write && ( hdr.nvev == global.nobs()	// parent thy vev file
			|| hdr.nvev == global.nobsG() )
		      && hdr.ncoup <= ncoup 
		      && hdr.nvev  <= numerics.vev.n_elem ;
	}
    else gripe (format("Inconsistent save file theory {}", hdr.name.data())) ;
    }

void Save::read_sysindex ()				// Load SysIndex
    {
    auto&	stream	{ global.sysstream } ;
    auto&	index	{ global.sysindex } ;
    auto	ptr	{ index.data() } ;

    stream.seekg (-sizeof index, ios_base::end) ;
    stream.read (cast_to<char*>(ptr), sizeof index) ;
    if (stream.fail()) ioerror ("read_sysindex: I/O error!") ;
    if (Blab::blablevel[BLAB::SAVE]) cout << "Read Index\n" ;
    }

void Save::read_op (int stage)				// Read Op record
    {
    auto&	stream	{ global.sysstream } ;
    auto&	record	{ global.data(stage).op } ;

    if (record.entry().id != RecordID::Op)
	gripe ("read_op: bad record ID!") ;

    record.readrec (stream) ;
    if (stream.fail()) ioerror ("read_op: I/O error!") ;

    uint	start	{ stage ? global.nopG() : 0 } ;
    uint	indx	{ start } ;
    uint	nop	{ record.entry().nelem } ;
    RecHdr*	recptr	{ record.data() } ;
    RecHdr*	recend	{ recptr + record.size() } ;

    Op::purge (start) ;
    Op::list.reserve (start + nop) ;
    while (recptr < recend)
	{
	OpHdr	hdr	( *recptr++ ) ;
	symb*	ptr	{ cast_to<symb*> (recptr) } ;
	SymbStr	s	{ ptr, ptr + hdr.len } ;
	Op	o	{ s, hdr } ;
	auto	k	{ Op::store(o) } ;

	if (k != indx++) gripe ("read_op: Inconsistent Op::list") ;

	recptr += 1 + (hdr.len - 1) / sizeof (RecHdr) ;
	}
    record.free() ;
    if (indx - start != nop)
	gripe ("read_op: Inconsistent save record!") ;

    global.info(stage).nop = nop ;
    if (Blab::blablevel[BLAB::SAVE]) cout << "Loaded Op\n" << flush ;
    }

void Save::read_obs (int stage)				// Read Obs record
    {
    auto&	stream	{ global.sysstream } ;
    auto&	record	{ global.data(stage).obs } ;

    if (record.entry().id != RecordID::Obs)
	gripe ("read_obs: bad record ID!") ;
    record.readrec (stream) ;
    if (stream.fail())
	ioerror ("read_obs: I/O error!") ;

    uint	start	{ stage ? global.nobsG() : 0 } ;
    uint	indx    { start } ;
    uint	nobs	{ record.entry().nelem } ;
    RecHdr*	recptr	{ record.data() } ;
    RecHdr*	recend	{ recptr + record.size() } ;

    ObsList::obs.purge (indx) ;
    ObsList::obs.map.reserve (start + nobs) ;
    ObsList::obs.reserve     (start + nobs) ;
    while (recptr < recend)
	{
	ObsHdr	hdr	( *recptr++ ) ;
	symb*	ptr	{ cast_to<symb*> (recptr) } ;
	SymbStr	s	{ ptr, ptr + hdr.len } ;
	Obs	o	{ s, hdr } ;

	if (ObsList::obs.store (o) != indx++)
	    gripe ("read_obs: Inconsistent ObsList::obs") ;

	if (global.info(stage).maxord < o.order() && o.corder == o.xorder)
	    global.info(stage).maxord = o.order() ;

	recptr += 1 + (hdr.len - 1) / sizeof (RecHdr) ;
	}
    record.free() ;
    if (indx - start != nobs)
	gripe ("read_obs: Inconsistent save record!") ;

    global.info(stage).nobs = nobs ;
    Global::mk_bcktlist (stage) ;
    Numerics::numericsinit() ;
    Canon::cache.reload() ;
    if (Blab::blablevel[BLAB::SAVE]) cout << "Loaded Obs\n" << flush ;
    }

void Save::read_gen (int stage)				// Read Gen record
    {
    auto&	stream	{ global.sysstream } ;
    auto&	record	{ global.data(stage).gen } ;
    int		optsiz	{ sizeof (OpTerm) / sizeof (RecHdr) } ;

    if (record.entry().id != RecordID::Gen)
	gripe ("read_gen: bad record ID!") ;
    record.readrec (stream) ;
    if (stream.fail())
	ioerror ("read_gen: I/O error!") ;

    RecHdr*	recptr	{ record.data() } ;
    RecHdr*	recend	{ recptr + record.size() } ;
    uint	nelem   ( 0 ) ;
    uint	opstart	{ stage ? global.nopG() : 0 } ;

    for (auto& gens : global.info(stage).gens)  gens.clear() ;
    for (auto& even : global.info(stage).neven) even = 0 ;

    while (recptr < recend)
	{
	GenHdr	hdr	( *recptr++ ) ;
	doub	coeff	{ *cast_to<doub*> (recptr++) } ;
	OpTerm*	ptr	{ cast_to<OpTerm*> (recptr) } ;
	OpSum	s	( ptr, ptr + hdr.len ) ;
	Gen	gen	{ s, hdr, coeff } ;

	global.info(stage).gens[hdr.rep].push_back (gen) ;

	if (!gen.T_odd) ++global.info(stage).neven[hdr.rep] ;
	if (global.info(stage).maxgen < gen.order)
	    global.info(stage).maxgen = gen.order ;
	++nelem ;
	recptr += hdr.len * optsiz ;
	}
    record.free() ;

    if (nelem != record.entry().nelem)
	gripe ("read_gen: Inconsistent save record!") ;

    //Op::clean (opstart,global.info(stage).maxgen) ;
    if (Blab::blablevel[BLAB::SAVE]) cout << "Loaded Gen\n" << flush ;
    }

void Save::read_grad (int stage)			// Read Grad record
    {
    auto&	stream	{ global.sysstream } ;
    auto& 	record	{ global.data(stage).grad } ;

    if (record.entry().id != RecordID::Grad)
	gripe ("read_grad: bad record ID!") ;
    record.readrec (stream) ;
    if (stream.fail())
	ioerror ("read_grad: I/O error!") ;
    else if (record.size() && Blab::blablevel[BLAB::SAVE])
	cout << "Loaded Grad\n" << flush ;
    }

void Save::read_curv (int stage)			// Read Curv records
    {
    for (int rep(0) ; rep < theory.nrep ; ++rep)
	{
	auto&	stream	{ global.sysstream } ;
	auto&	record	{ global.data(stage).curv[rep] } ;
	auto&	repname	{ Rep::list[rep].name } ;

	if (record.entry().id != RecordID::Curv)
	    gripe (format("read_curv {}: bad record ID!", repname)) ;

	record.readrec (stream) ;
	if (stream.fail())
	    ioerror (format("read_curv {}: I/O error!", repname)) ;
	else if (record.size() && Blab::blablevel[BLAB::SAVE])
	    cout << "Loaded Curv " << repname << "\n" << flush ;
	}
    }

void Save::read_lagr (int stage)			// Read Lagr records
    {
    for (int rep(0) ; rep < theory.nrep ; ++rep)
	{
	auto&	stream	{ global.sysstream } ;
	auto&	record	{ global.data(stage).lagr[rep] } ;
	auto&	repname	{ Rep::list[rep].name } ;

	if (record.entry().id != RecordID::Lagr)
	    gripe (format("read_lagr {}: bad record ID!", repname)) ;

	record.readrec (stream) ;
	if (stream.fail())
	    ioerror (format("read_lagr {}: I/O error!", repname)) ;
	else if (record.size() && Blab::blablevel[BLAB::SAVE])
	    cout << "Loaded Lagr " << repname << "\n" << flush ;
	}
    }

void Save::read_geos (int stage)			// Read Geo records
    {
    if (global.geoswap) return ;
    for (int bcktnum(0) ; bcktnum < global.info(stage).bckt.size() ; ++bcktnum)
	{
	auto&	stream	{ global.sysstream } ;
	auto&	record	{ global.data(stage).geos[bcktnum] } ;
	auto	nbckts	{ global.info(stage).bckt.size() } ;

	if (record.entry().id != RecordID::Geos)
	    gripe (format("read_geos [{}]: bad record ID!",bcktnum)) ;

	record.readrec (stream) ;
	if (stream.fail())
	    ioerror (format("read_geos [{}]: I/O error!",bcktnum)) ;
	else if (record.size() && Blab::blablevel[BLAB::SAVE])
	    cout << "Loaded Geos [" << bcktnum << "/" << nbckts << "]\n" << flush ;
	}
    }

void Save::read_geo_bckt (int bcktnum)			// Read Geo bucket
    {
    alignas (16*1024) thread_local char	mybuf [512*1024] ;
    thread_local vector<RecHdr>		myvec ;
    thread_local fstream		mystream ;
    thread_local string			mypath { syspath } ;
    thread_local bool			checkout { false } ;

    if (bcktnum >= 0)
	{
	auto	 nbckts { global.info().bckt.size() } ;
	auto& 	 record	{ global.data().geos[bcktnum] } ;

	if (!global.sysstream.is_open())
	    gripe ("Must write or load save file!") ;

	if (mypath != syspath) mystream.close() ;

	if (!mystream.is_open())
	    {
	    mystream.rdbuf()->pubsetbuf (mybuf, sizeof mybuf) ;
	    mystream.open (mypath = syspath, ios::in | ios::binary) ;
	    }

	if (record.entry().id != RecordID::Geos)
	    gripe (format("read_geo_bckt [{}]: bad record ID!",bcktnum)) ;
	if (record.size())
	    gripe (format("read_geo_bckt [{}]: non-empty record!",bcktnum)) ;

	if (!checkout)
	    {
	    record.swap (myvec) ;
	    checkout = true ;
	    }
	else gripe (format("read_geo_bckt [{}]: bad check out!",bcktnum)) ;

	record.readrec (mystream) ;
	if (mystream.fail())
	    ioerror (format("read_geo_bckt [{}]: I/O error!", bcktnum)) ;
	if (Blab::blablevel[BLAB::SAVE])
	    cout << "Loaded Geos [" << bcktnum << "/" << nbckts << "]\n" << flush ;
	}
    else if (checkout)
	{
	myvec.swap (global.data().geos[-bcktnum-1]) ;
	checkout = false ;
	}
    else gripe (format("read_geo_bckt [{}]: bad check in!",bcktnum)) ;
    }

void Save::read_stat (int stage)			// Read Stat record
    {
    auto&	stream	{ global.sysstream } ;
    auto&	record	{ global.data(stage).stat } ;
    uint	nelem	{ record.entry().nelem } ;
    auto	ptr	{ &global.info(stage).count } ;

    if (record.entry().id != RecordID::Stat)
	gripe ("read_stat: bad record ID!") ;

    stream.seekg (record.entry().filepos, ios_base::beg) ;
    stream.read (cast_to<char*>(ptr), nelem * sizeof (Counter)) ;
    if (stream.fail()) ioerror ("read_stat: I/O error!") ;
    if (Blab::blablevel[BLAB::SAVE]) cout << "Loaded Stat\n" << flush ;
    }

bool Save::read_coup (int set, Couplings* ptr)		// Read Coupling set
    {
    auto&	stream	{ global.vevstream } ;
    ulong	offset	{ set * filehdr.cvsetsize() } ;
    auto&	list	{ ptr ? *ptr : Coupling::list } ;

    if (set < 0) stream.seekg (offset, ios_base::end) ;
    else	 stream.seekg (offset + sizeof filehdr, ios_base::beg) ;

    stream.read (cast_to<char*>(list.data()), filehdr.coupsize()) ;
    if (stream.eof()) { stream.clear() ; return false ; }
    if (stream.fail()) ioerror ("read_coup: I/O error!") ;
    if (Blab::blablevel[BLAB::SAVE]) cout << "Loaded Coup\n" << flush ;
    return true ;
    }

void Save::read_vev (int set)				// Read Vev data
    {
    auto&	stream	{ global.vevstream } ;
    ulong	offset	{ filehdr.coupsize() + set * filehdr.cvsetsize() } ;
    auto	ptr	{ numerics.vev.memptr() } ;

    if (set < 0) stream.seekg (offset, ios_base::end) ;
    else	 stream.seekg (offset + sizeof filehdr, ios_base::beg) ;

    stream.read (cast_to<char*>(ptr), filehdr.nvev * sizeof (real)) ;
    if (stream.fail()) ioerror ("read_vev: I/O error!") ;
    if (Blab::blablevel[BLAB::SAVE]) cout << "Loaded Vev\n" << flush ;
    }
