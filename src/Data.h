#ifndef DATA_H
#define DATA_H
#include "Gordion.h"
#include "Gripe.h"
#include <functional>

struct OpHdr					// Operator record item header
    {
    short	order ;
    char	type ;
    bool	prim ;
    ushort	pad ;
    ushort	len ;
    doub&	data()	{ return *cast_to<doub*>(this) ; }
    } ;

struct ObsHdr					// Observable record item header
    {
    short	corder ;
    short	xorder ;
    char	type ;
    char	pad ;
    ushort	len ;
    doub&	data()	{ return *cast_to<doub*>(this) ; }
    } ;

struct GenHdr					// Generator record item header
    {
    short	rep ;
    short	order ;
    char	type ;
    bool	T_odd ;
    ushort	len ;
    doub&	data()	{ return *cast_to<doub*>(this) ; }
    } ;

struct GradHdr					// Gradient record poly header
    {
    ushort	gen ;
    ushort	Hterm ;
    ushort	pad ;
    ushort	len ;
    doub&	data()	{ return *cast_to<doub*>(this) ; }
    } ;

struct CurvHdr					// Curvature record poly header
    {
    ushort	gen2 ;
    ushort	gen1 ;
    ushort	Hterm ;
    ushort	len ;
    doub&	data()	{ return *cast_to<doub*>(this) ; }
    } ;

struct LagrHdr					// Lagrange record poly header
    {
    ushort	gen2 ;
    ushort	gen1 ;
    ushort	pad ;
    ushort	len ;
    doub&	data()	{ return *cast_to<doub*>(this) ; }
    } ;

struct GeoHdr					// Geodesic record poly header
    {
    uint	indx ;
    ushort	gen ;
    ushort	len ;
    doub&	data()	{ return *cast_to<doub*>(this) ; }
    } ;

struct RecHdr				// Data record header
    {
    ushort	info[3] ;			// type-specific info
    ushort	len ;				// subsequent # data elements

    operator const ObsHdr&  () const { return *cast_to<const ObsHdr*>(this)  ; }
    operator const OpHdr&   () const { return *cast_to<const OpHdr*>(this)   ; }
    operator const GenHdr&  () const { return *cast_to<const GenHdr*>(this)  ; }
    operator const GradHdr& () const { return *cast_to<const GradHdr*>(this) ; }
    operator const CurvHdr& () const { return *cast_to<const CurvHdr*>(this) ; }
    operator const LagrHdr& () const { return *cast_to<const LagrHdr*>(this) ; }
    operator const GeoHdr&  () const { return *cast_to<const GeoHdr*>(this)  ; }

    doub& data() { return *cast_to<doub*>(this) ; }

    RecHdr () {}
    RecHdr (OpHdr&    h) { data() = h.data() ; }
    RecHdr (OpHdr&&   h) { data() = h.data() ; }
    RecHdr (ObsHdr&   h) { data() = h.data() ; }
    RecHdr (ObsHdr&&  h) { data() = h.data() ; }
    RecHdr (GenHdr&   h) { data() = h.data() ; }
    RecHdr (GenHdr&&  h) { data() = h.data() ; }
    RecHdr (GradHdr&  h) { data() = h.data() ; }
    RecHdr (GradHdr&& h) { data() = h.data() ; }
    RecHdr (CurvHdr&  h) { data() = h.data() ; }
    RecHdr (CurvHdr&& h) { data() = h.data() ; }
    RecHdr (LagrHdr&  h) { data() = h.data() ; }
    RecHdr (LagrHdr&& h) { data() = h.data() ; }
    RecHdr (GeoHdr&   h) { data() = h.data() ; }
    RecHdr (GeoHdr&&  h) { data() = h.data() ; }
    } ;

static_assert (sizeof (RecHdr) == sizeof (doub)) ;

enum class RecordID : char			// Data record type
    {
    Null,
    Op,
    Obs,
    Gen,
    Grad,
    Curv,
    Lagr,
    Geos,
    Stat,
    Indx
    } ;

class RecIndx				// Data index entry
    {
    public:
    RecordID	id { RecordID::Null } ;	// record ID
    uchar	nslice  { 1 } ;		// # slices if multi-dim block
    uint	ncol    { 1 } ;		// # cols if multi-dim block
    union
	{
	uint	nrow    { 0 } ;		// # rows if multi-dim block
	uint	nelem ;			// # items in block
	} ;
    ulong	reclen  { 0 } ;		// data length in RecHdr units
    ulong	filepos { 0 } ;		// file offset of data

    ulong items () const		// total number of items
	{ return nrow * ncol * (int) nslice ; }

    inline static string idname[] { "Null", "Op", "Obs", "Gen", "Grad",
				    "Curv", "Lagr", "Geos", "Stat", "Indx" } ;
    } ;

template <size_t N>
class RecIndxArr : public array<RecIndx,N>	// RecIndx array
    {
    public:
    RecIndx& next()			// return next unused entry
	{
	for (auto& d : *this)
	    {
	    if (d.id == RecordID::Null) return d ;
	    }
	gripe ("No unused RecIndx entry!") ;
	}
    } ;

class DataRec : public vector<RecHdr>		// Data record
    {
    public:
    std::reference_wrapper<RecIndx> indexref ;	// Data index entry

    DataRec (RecordID = RecordID::Null) ;	// Constructor

    RecIndx& entry() const { return indexref ;} // Data index entry

    char* recptr()				// char* ptr to data
	{
	return cast_to<char*> (data()) ;
	}
    const char* reccptr()			// const char* ptr to data
	{
	return cast_to<const char*> (data()) ;
	}
    void writerec (fstream& stream)		// Write record to stream
	{
	entry().reclen  = size() ;
	entry().filepos = static_cast<std::streamoff> (stream.tellp()) ;
	if (size()) stream.write (reccptr(), size() * sizeof (RecHdr)) ;
	}
    void readrec (fstream& stream)		// Read record from stream
	{
	resize (entry().reclen) ;
	if (!entry().reclen) return ;
	stream.seekg (entry().filepos, std::ios_base::beg) ;
	stream.read (recptr(), size() * sizeof (RecHdr)) ;
	}
    void clear()				// Clear data record
	{
	vector<RecHdr>::clear() ;
	entry().nslice	= 1 ;
	entry().ncol	= 1 ;
	entry().nrow	= 0 ;
	}
    void free()					// Free data record
	{
	vector<RecHdr>().swap(*this) ;
	}
    } ;

#endif
