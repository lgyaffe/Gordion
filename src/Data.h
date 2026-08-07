#ifndef DATA_H
#define DATA_H
#include "Gordion.h"
#include "Gripe.h"
#include <functional>

#ifdef NUM32
    struct OpHdr				// Op record item header
	{
	ushort	order : 6 ;
	ushort	type  : 2 ;
	ushort	prim  : 1 ;
	} ;

    struct ObsHdr				// Obs record item header
	{
	ushort	corder : 6 ;
	ushort	xorder : 6 ;
	ushort	type   : 3 ;
	} ;

    struct GenHdr				// Gen record item header
	{
	ushort	rep   : 6 ;
	ushort	order : 6 ;
	ushort	type  : 2 ;
	ushort	T_odd : 1 ;
	} ;
#else
    struct OpHdr				// Op record item header
	{
	uint	order : 10 ;
	uint	type  : 2  ;
	uint	prim  : 1  ;
	} ;

    struct ObsHdr				// Obs record item header
	{
	uint	corder : 10 ;
	uint	xorder : 10 ;
	uint	type   : 4  ;
	} ;

    struct GenHdr				// Gen record item header
	{
	uint	rep   : 8 ;
	uint	order : 10 ;
	uint	type  : 2 ;
	uint	T_odd : 1 ;
	} ;
#endif

struct RecHdr				// Data record header
    {
    usmall	len ;			// Record length (w/o header)
    union
	{
	OpHdr	op  ;
	ObsHdr	obs ;
	GenHdr	gen ;
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

    const auto& len() const { return hdr.len ; }
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
    } ;

#endif
