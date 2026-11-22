#ifndef DATA_H
#define DATA_H
#include "Types.h"
#include "Gripe.h"
#include <functional>

#ifdef NUM64
using uqutr = ushort ;		// quarter size unsigned int
using uhalf = uint ;		// half size unsigned int
using unumb = ulong ;		// full size unsigned int
#else
using uqutr = uchar ;		// quarter size unsigned int
using uhalf = ushort ;		// half size unsigned int
using unumb = uint ;		// full size unsigned int
#endif

class Op ;
class Obs ;
class Gen ;

inline constexpr int ordbits { 8 * sizeof (uqutr) } ;
inline constexpr int lenbits { 8 * sizeof (uhalf) - 3 } ;
inline constexpr int genbits { 8 * sizeof (unumb) - 3 - ordbits } ;

#ifdef DELPHI
    struct OpHdr			// Op record item header
	{
	uhalf	len ;
	uhalf	type ;
	OpHdr (const Op&) ;
	} ;
    struct ObsHdr			// Obs record item header
	{
	uhalf	len ;
	uhalf	type ;
	ObsHdr (const Obs&) ;
	} ;
    struct GenHdr			// Gen record item header
	{
	unumb 	len   : genbits ;
	unumb	type  : 2 ;
	unumb	T_odd : 1 ;
	unumb	group : ordbits ;
	GenHdr (const Gen&) ;
	} ;
#else
    struct OpHdr			// Op record item header
	{
	unumb	len   : lenbits ;
	unumb	type  : 2  ;
	unumb	prim  : 1  ;
	unumb	order : ordbits ;
	OpHdr (const Op&) ;
	} ;
    struct ObsHdr			// Obs record item header
	{
	unumb	len    : lenbits ;
	unumb	type   : 3 ;
	unumb	corder : ordbits ;
	unumb	xorder : ordbits ;
	ObsHdr (const Obs&) ;
	} ;
    struct GenHdr			// Gen record item header
	{
	unumb	len   : genbits ;
	unumb	type  : 2  ;
	unumb	T_odd : 1 ;
	unumb	order : ordbits ;
	GenHdr (const Gen&) ;
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

template <size_t N>
class DataArr : public array<DataRec,N>			// Array of DataRec's
    {
    public:
    DataArr (SysIndex& indx, RecordID id)		// Public constructor
	: DataArr (std::make_index_sequence<N>{}, indx, id) {}

    private:
    template <size_t... Is> constexpr			// Private constructor
    DataArr (std::index_sequence<Is...>, SysIndex&, RecordID) ;
    } ;

#endif
