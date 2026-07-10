#include "Theory.h"
#include "Print.h"
#include "Global.h"
#include "Obs.h"
#ifndef MASSSCALE
#define MASSSCALE 0
#endif

void Theory::theoryinit()			// Verity theory applicability
    {
    global.savedir.append (theory.name.data()) ;
    global.MMAdir.append  (theory.name.data()) ;
    Print::print_theory() ;
    Print::print_version() ;
    }

void Theory::theorydefn (int stage)		// Define hamiltonian or action
    {
    Coeff	unitcoeff {} ;
    char8	lambda	  {"lambda"} ;
    int		lamindx   ( Coupling::indx (lambda) ) ;
    PolyMap	map	  { ObsList::base } ;
    ObsList&	baslist   { ObsList::base } ;
    symb	link[4]   { 'x', 'y', 'z', 'w' } ;
    symb	Link[4]   { 'X', 'Y', 'Z', 'W' } ;
    symb	ferm[4]   { 'f', 'g', 'h', 'i' } ;
    symb	Ferm[4]   { 'F', 'G', 'H', 'I' } ;
    bool	isham	  { !theory.euclid } ;
    bool	iseuc	  {  theory.euclid } ;

    if (lamindx < 0)
	{
	Coupling::list.emplace_back (lambda) ;
	lamindx = Coupling::indx (lambda) ;
	Coupling::list[lamindx].value = 1000 ;
	Coupling::list[lamindx].stage = 0 ;
	}
    Coeff lamcoeff	{{lamindx,1}} ;
    Coeff laminvcoeff	{{lamindx,-1}} ;

    if (stage == 0)
	{
	global.info(0).Hterms.clear() ;

	if (isham)				// Hamiltonian gauge kinetic
	    {
	    ObsPoly	kinetic  { baslist } ;
	    ObsPoly	ckinetic { baslist } ;
	    doub	coeff    { theory.dim > 1 ? 0.25 : 1.0 } ;
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		uint	indx	{ baslist.find (Str (KinG + i)) } ;
		Obs	EE	{ baslist(indx) } ; EE.canon() ;
		uint	cindx	{ baslist.find (EE) } ;
		if (indx == UINT_MAX) fatal ("Baselist missing EE") ;
		kinetic.push_back (PolyTerm (indx,  coeff)) ;
		map.add		  (PolyTerm (cindx, coeff)) ;
		}
	    ckinetic.push_map (map) ;
	    global.info(0).Hterms.emplace_back (lamcoeff, kinetic, ckinetic) ;
	    }
	else					// Euclidean gauge entropy
	    {
	    ObsPoly	gauge_ent (baslist) ;
	    uint	indx	{ baslist.find (Str (EntrG)) } ;
	    if (indx == UINT_MAX) fatal ("Baselist missing S") ;
	    gauge_ent.push_back (PolyTerm (indx, -1.0)) ;
	    global.info(1).Hterms.emplace_back (unitcoeff, gauge_ent, gauge_ent) ;
	    }

	ObsPoly plaquette  (baslist) ;
	ObsPoly cplaquette (baslist) ;
	if (theory.dim == 1 && theory.box.comp[0])	// Polyakov loop
	    {
	    plaquette.push_back  (PolyTerm(0, 2)) ;
	    cplaquette.push_back (PolyTerm(0, 2)) ;
	    Str		lll	 { string (theory.box.comp[0], link[0]) } ;
	    Str		LLL	 { string (theory.box.comp[0], Link[0]) } ;
	    uint	indx	 { baslist.find (lll) } ;
	    uint	Indx	 { baslist.find (LLL) } ;
	    Obs		polyakov { baslist(indx) } ; polyakov.canon() ;
	    Obs		Polyakov { baslist(Indx) } ; Polyakov.canon() ;
	    uint	cindx	 { baslist.find (polyakov) } ;
	    uint	cIndx	 { baslist.find (Polyakov) } ;
	    if (indx == UINT_MAX) fatal ("Baselist missing polyakov") ;
	    if (Indx == UINT_MAX) fatal ("Baselist missing Polyakov") ;
	    plaquette.push_back (PolyTerm (indx,  -1.0)) ;
	    plaquette.push_back (PolyTerm (Indx,  -1.0)) ;
	    map.add		(PolyTerm (cindx, -1.0)) ;
	    map.add		(PolyTerm (cIndx, -1.0)) ;
	    }
	else						// Plaquette terms
	    {
	    if (theory.dim * (theory.dim-1))
		{
		plaquette.push_back  (PolyTerm(0, theory.dim * (theory.dim-1))) ;
		cplaquette.push_back (PolyTerm(0, theory.dim * (theory.dim-1))) ;
		}
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		for (int j(i) ; ++j < theory.dim ;)
		    {
		    Str		xyXY   { string{link[i],link[j],Link[i],Link[j]} } ;
		    Str		xYXy   { string{link[i],Link[j],Link[i],link[j]} } ;
		    uint	indx  { baslist.find (xyXY) } ;
		    uint	Indx  { baslist.find (xYXy) } ;
		    Obs		plaq  { baslist(indx) } ; plaq.canon() ;
		    Obs		Plaq  { baslist(Indx) } ; Plaq.canon() ;
		    uint	cindx { baslist.find (plaq) } ;
		    uint	cIndx { baslist.find (Plaq) } ;
		    if (indx == UINT_MAX) fatal ("Baselist missing plaq") ;
		    if (Indx == UINT_MAX) fatal ("Baselist missing Plaq") ;
		    plaquette.push_back (PolyTerm (indx,  -1.0)) ;
		    plaquette.push_back (PolyTerm (Indx,  -1.0)) ;
		    map.add		(PolyTerm (cindx, -1.0)) ;
		    map.add		(PolyTerm (cIndx, -1.0)) ;
		    }
		}
	    }
	cplaquette.push_map (map) ;
	if (plaquette.size())
	    global.info(0).Hterms.emplace_back (laminvcoeff, plaquette, cplaquette) ;
	}
    else if (theory.nf)
	{
	char8	mass	 {"mass"} ;
	int	massindx { Coupling::indx (mass) } ;
	if (massindx < 0)
	    {
	    Coupling::list.emplace_back (mass) ;
	    massindx  = Coupling::indx (mass) ;
	    Coupling::list[massindx].value = 1 ;
	    Coupling::list[massindx].stage = 1 ;
	    }
	Coeff masscoeff {{massindx,1},{lamindx,MASSSCALE}} ;
	global.info(1).Hterms.clear() ;
	
	if (isham)
	    {
	    ObsPoly kinetic_F  (baslist) ;
	    ObsPoly ckinetic_F (baslist) ;
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		uint    indx	{ baslist.find (Str (KinF + i)) } ;
		Obs	ee	{ baslist(indx) } ; ee.canon() ;
		uint	cindx	{ baslist.find (ee) } ;
		if (indx == UINT_MAX) fatal ("Baselist missing ee") ;
		kinetic_F.push_back (PolyTerm (indx,  0.25)) ;
		map.add		    (PolyTerm (cindx, 0.25)) ;
		}
	    ckinetic_F.push_map (map) ;
	    global.info(1).Hterms.emplace_back (lamcoeff, kinetic_F, ckinetic_F) ;
	    }

	ObsPoly hop_term   (baslist) ;
	ObsPoly chop_term  (baslist) ;
	for (int i(0) ; i < theory.nf ; i+=2)
	    {
	    for (int j(0) ; j < theory.dim ; ++j)
		{
		Str 	Fxf	{ string{Ferm[i],link[j],ferm[i]} } ;
		Str 	FXf	{ string{Ferm[i],Link[j],ferm[i]} } ;
		uint	indx	{ baslist.find (Fxf) } ;
		uint	Indx	{ baslist.find (FXf) } ;
		Obs	hop	{ baslist(indx) } ; hop.canon() ;
		Obs	Hop	{ baslist(Indx) } ; Hop.canon() ;
		uint	cindx	{ baslist.find (hop) } ;
		uint	cIndx	{ baslist.find (Hop) } ;
		if (indx == UINT_MAX) fatal ("Baselist missing Fxf") ;
		if (Indx == UINT_MAX) fatal ("Baselist missing FXf") ;
		hop_term.push_back (PolyTerm (indx,  0.5)) ;
		hop_term.push_back (PolyTerm (Indx,  0.5)) ;
		map.add		   (PolyTerm (cindx, 0.5)) ;
		map.add		   (PolyTerm (cIndx, 0.5)) ;
		}
	    }
	chop_term.push_map (map) ;
	global.info(1).Hterms.emplace_back (unitcoeff, hop_term, chop_term) ;

	ObsPoly mass_term  (baslist) ;
	ObsPoly cmass_term (baslist) ;

	for (int i(0) ; i < theory.nf ; i+=2)
	    {
	    uint indx { baslist.find (string {Ferm[i+isham],ferm[i]}) } ;
	    if (indx == UINT_MAX) fatal ("Baselist missing Ff/Gf") ;
	    mass_term.push_back (PolyTerm (indx, 1.0)) ;
	    map.add		(PolyTerm (indx, 1.0)) ;
	    }
	cmass_term.push_map (map) ;
	global.info(1).Hterms.emplace_back (masscoeff, mass_term, cmass_term, iseuc) ;

	if (iseuc)
	    {
	    ObsPoly	fermi_ent (baslist) ;
	    uint	indx	{ baslist.find (Str {EntrF}) } ;
	    if (indx == UINT_MAX) fatal ("Baselist missing s") ;
	    fermi_ent.push_back (PolyTerm (indx, -1.0)) ;
	    global.info(1).Hterms.emplace_back (unitcoeff, fermi_ent, fermi_ent) ;
	    }
	}
    global.info().MMAlimit = global.nobs() ;
    global.info().MMAlist.clear() ;
    }
