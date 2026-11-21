#ifndef DATA_H
#define DATA_H
#include "Types.h"
#include "Gripe.h"
#include <functional>

#ifdef NUM64
using usmall = uint ;
using unumb  = ulong ;
#define LENBITS	29
#define ORDBITS 16
#else
using usmall = ushort ;
using unumb  = uint ;
#define LENBITS	13
#define ORDBITS  8
#endif

class Op ;
class Obs ;
class Gen ;

#ifdef DELPHI
    struct OpHdr			// Op record item header
	{
	usmall	len ;
	usmall	type ;
	OpHdr (const Op&) ;
	} ;
    struct ObsHdr			// Obs record item header
	{
	usmall	len ;
	usmall	type ;
	ObsHdr (const Obs&) ;
	} ;
    struct GenHdr			// Gen record item header
	{
	unumb 	len   : LENBITS ;
	unumb	type  : 2 ;
	unumb	T_odd : 1 ;
	unumb	group : ORDBITS ;
	unumb	rep   : 6 ;
	GenHdr (const Gen&, int) ;
	} ;
#else
    struct OpHdr			// Op record item header
	{
	unumb	len   : LENBITS ;
	unumb	type  : 2  ;
	unumb	prim  : 1  ;
	unumb	order : ORDBITS ;
	OpHdr (const Op&) ;
	} ;
    struct ObsHdr			// Obs record item header
	{
	unumb	len    : LENBITS ;
	unumb	type   : 3 ;
	unumb	corder : ORDBITS ;
	unumb	xorder : ORDBITS ;
	ObsHdr (const Obs&) ;
	} ;
    struct GenHdr			// Gen record item header
	{
	unumb	len   : LENBITS ;
	unumb	type  : 2  ;
	unumb	T_odd : 1 ;
	unumb	order : ORDBITS ;
	unumb	rep   : 6 ;
	GenHdr (const Gen&, int) ;
	} ;				
#endif

struct PolyHdr { unumb len ; } ;	// Poly record header

struct RecHdr				// Data record header
    {
    union
	{
	OpHdr	op   ;
	ObsHdr	obs  ;
	GenHdr	gen  ;
	PolyHdr	poly ;
	} ;
    } ;

struct Element				// Data record element
    {
    union
	{
	RecHdr	hdr ;		// Record header
	real	coeff ;		// Poly coefficient
	numb	index ;		// Poly Obs index
	} ;

    Element () {}
    Element (OpHdr  h) { hdr.op  = h ; }
    Element (ObsHdr h) { hdr.obs = h ; }
    Element (GenHdr h) { hdr.gen = h ; }
    Element (real   c) { coeff   = c ; }
    Element (numb   n) { index   = n ; }
    } ;

static_assert (sizeof (Element) == sizeof (real)) ;

enum class RecordID : char		// Data record type
    {
    Null, Op,   Obs,  Gen,
    Ham,  Grad, Curv, Lagr,
    Geos, Stat, Indx
    } ;

class RecIndx				// Data index entry
    {
    public:
    RecordID	id { RecordID::Null } ;	// record ID
    uchar	nslice  { 1 } ;		// # slices if multi-dim block
    uint	ncol    { 1 } ;		// # cols if multi-dim block
    union
	{
	ulong	nrow    { 0 } ;		// # rows if multi-dim block
	ulong	nelem ;			// # items in block
	} ;
    ulong	reclen  { 0 } ;		// data length in Element units
    ulong	filepos { 0 } ;		// file offset of data

    ulong items () const		// total number of items
	{ return nrow * ncol * (int) nslice ; }

    inline static string idname[] { "Null", "Op", "Obs", "Gen", "Ham", "Grad",
				    "Curv", "Lagr", "Geos", "Stat", "Indx" } ;
    } ;

template <size_t N>
class RecIndxArr : public array<RecIndx,N>	// RecIndx array
    {
    public:
    RecIndx& next()			// return next unused entry
	{
	for (auto& d : *this)
	    if (d.id == RecordID::Null) return d ;
	abort ("No unused RecIndx entry!") ;
	}
    } ;

class SysIndex ;

class DataRec : public vector<Element>		// Data record
    {
    public:
    std::reference_wrapper<RecIndx> indexref ;	// Data index entry

    DataRec (SysIndex&, RecordID) ;		// Constructor

    RecIndx& entry() const { return indexref ;} // Return index entry

    char* recptr()				// char* ptr to data
	{
	return cast_to<char*>(data()) ;
	}
    const char* reccptr()			// const char* ptr to data
	{
	return cast_to<const char*>(data()) ;
	}
    void writerec (fstream& stream)		// Write record to stream
	{
	entry().reclen  = size() ;
	entry().filepos = static_cast<std::streamoff> (stream.tellp()) ;
	if (size()) stream.write (reccptr(), size() * sizeof (Element)) ;
	}
    void readrec (fstream& stream)		// Read record from stream
	{
	resize (entry().reclen) ;
	if (!entry().reclen) return ;
	stream.seekg (entry().filepos, std::ios_base::beg) ;
	stream.read (recptr(), size() * sizeof (Element)) ;
	if (stream.gcount() != entry().reclen * sizeof (Element))
	    ioerror ("I/O error: short read") ;
	}
    void clear()				// Clear data record
	{
	vector<Element>::clear() ;
	entry().nslice	= 1 ;
	entry().ncol	= 1 ;
	entry().nrow	= 0 ;
	}
    void free()					// Free data record
	{
	vector<Element>().swap(*this) ;
	}
    ulong filepos() const			// File position
	{
	return entry().filepos ;
	}
    } ;

#endif
