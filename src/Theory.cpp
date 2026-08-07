#include "Theory.h"
#include "Print.h"
#include "Global.h"
#include "Obs.h"
#ifndef MASSSCALE
#define MASSSCALE 0
#endif

void Theory::theoryinit()			// Report theory & version
    {
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
	Coupling::list.emplace_back (lambda, 0) ;
	lamindx = Coupling::indx (lambda) ;
	++Coupling::ncoupG ;
	}
    Coeff lamcoeff	{{lamindx, 1}} ;
    Coeff laminvcoeff	{{lamindx,-1}} ;

    if (stage == 0)
	{
	global.info(0).Hterms.clear() ;

	if (isham)				// Hamiltonian gauge kinetic
	    {
	    ObsPoly	kinetic  { baslist } ;
	    doub	coeff    { theory.dim > 1 ? 0.25 : 1.0 } ;
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		numb	indx	{ baslist.find (Str (KinG + i)) } ;
		if (indx == MAXNUM) fatal ("Baselist missing EE") ;
		kinetic.push_back (PolyTerm (indx, coeff)) ;
		}
	    global.info(0).Hterms.emplace_back (lamcoeff, kinetic) ;
	    }
	else					// Euclidean gauge entropy
	    {
	    ObsPoly	gauge_ent (baslist) ;
	    numb	indx	{ baslist.find (Str (EntrG)) } ;
	    if (indx == MAXNUM) fatal ("Baselist missing S") ;
	    gauge_ent.push_back (PolyTerm (indx, -1.0)) ;
	    global.info(1).Hterms.emplace_back (unitcoeff, gauge_ent) ;
	    }

	ObsPoly plaquette  (baslist) ;
	if (theory.dim == 1 && theory.box.comp[0])	// Polyakov loop
	    {
	    plaquette.push_back  (PolyTerm(0, 2)) ;
	    Str		lll	 { string (theory.box.comp[0], link[0]) } ;
	    Str		LLL	 { string (theory.box.comp[0], Link[0]) } ;
	    numb	indx	 { baslist.find (lll) } ;
	    numb	Indx	 { baslist.find (LLL) } ;
	    if (indx == MAXNUM) fatal ("Baselist missing polyakov") ;
	    if (Indx == MAXNUM) fatal ("Baselist missing Polyakov") ;
	    plaquette.push_back (PolyTerm (indx, -1.0)) ;
	    plaquette.push_back (PolyTerm (Indx, -1.0)) ;
	    }
	else						// Plaquette terms
	    {
	    if (theory.dim * (theory.dim-1))
		{
		plaquette.push_back  (PolyTerm(0, theory.dim * (theory.dim-1))) ;
		}
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		for (int j(i) ; ++j < theory.dim ;)
		    {
		    Str		xyXY   { string{link[i],link[j],Link[i],Link[j]} } ;
		    Str		xYXy   { string{link[i],Link[j],Link[i],link[j]} } ;
		    numb	indx  { baslist.find (xyXY) } ;
		    numb	Indx  { baslist.find (xYXy) } ;
		    if (indx == MAXNUM) fatal ("Baselist missing plaq") ;
		    if (Indx == MAXNUM) fatal ("Baselist missing Plaq") ;
		    plaquette.push_back (PolyTerm (indx, -1.0)) ;
		    plaquette.push_back (PolyTerm (Indx, -1.0)) ;
		    }
		}
	    }
	if (plaquette.size())
	    global.info(0).Hterms.emplace_back (laminvcoeff, plaquette) ;
	}
    else if (theory.nf)
	{
	char8	mass	 {"mass"} ;
	int	massindx { Coupling::indx (mass) } ;
	if (massindx < 0)
	    {						// N.B.: fermion 
	    Coupling::list.emplace_back (mass, 1) ;	// coupling must follow
	    massindx = Coupling::indx (mass) ;		// all gauge couplings
	    }
	Coeff masscoeff {{massindx,1},{lamindx,MASSSCALE}} ;
	global.info(1).Hterms.clear() ;
	
	if (isham)
	    {
	    ObsPoly kinetic_F  (baslist) ;
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		numb    indx	{ baslist.find (Str (KinF + i)) } ;
		if (indx == MAXNUM) fatal ("Baselist missing ee") ;
		kinetic_F.push_back (PolyTerm (indx,  0.25)) ;
		}
	    global.info(1).Hterms.emplace_back (lamcoeff, kinetic_F) ;
	    }

	ObsPoly hop_term   (baslist) ;
	for (int i(0) ; i < theory.nf ; i+=2)
	    {
	    for (int j(0) ; j < theory.dim ; ++j)
		{
		Str 	Fxf	{ string{Ferm[i],link[j],ferm[i]} } ;
		Str 	FXf	{ string{Ferm[i],Link[j],ferm[i]} } ;
		numb	indx	{ baslist.find (Fxf) } ;
		numb	Indx	{ baslist.find (FXf) } ;
		if (indx == MAXNUM) fatal ("Baselist missing Fxf") ;
		if (Indx == MAXNUM) fatal ("Baselist missing FXf") ;
		hop_term.push_back (PolyTerm (indx,  0.5)) ;
		hop_term.push_back (PolyTerm (Indx,  0.5)) ;
		}
	    }
	global.info(1).Hterms.emplace_back (unitcoeff, hop_term) ;

	ObsPoly mass_term  (baslist) ;
	for (int i(0) ; i < theory.nf ; i+=2)
	    {
	    numb indx { baslist.find (string {Ferm[i+isham],ferm[i]}) } ;
	    if (indx == MAXNUM) fatal ("Baselist missing Ff/Gf") ;
	    mass_term.push_back (PolyTerm (indx, 1.0)) ;
	    }
	global.info(1).Hterms.emplace_back (masscoeff, mass_term, iseuc) ;

	if (iseuc)
	    {
	    ObsPoly	fermi_ent (baslist) ;
	    numb	indx	{ baslist.find (Str {EntrF}) } ;
	    if (indx == MAXNUM) fatal ("Baselist missing s") ;
	    fermi_ent.push_back (PolyTerm (indx, -1.0)) ;
	    global.info(1).Hterms.emplace_back (unitcoeff, fermi_ent) ;
	    }
	}
    }
