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

void Theory::theorydefn ()			// Define hamiltonian or action
    {
    Coeff	unitcoeff {} ;
    char8	lambda	  {"lambda"} ;
    int		lamindx   ( Coupling::indx (lambda) ) ;
    PolyMap	map	  { ObsList::obs } ;
    ObsList&	obslist   { ObsList::obs } ;
    ObsList&	baslist   { ObsList::base } ;
    symb	link[4]   { 0x00, 0x01, 0x02, 0x03 } ;
    symb	Link[4]   { 0x04, 0x05, 0x06, 0x07 } ;
    symb	ferm[4]   { 0x28, 0x29, 0x2a, 0x2b } ;
    symb	Ferm[4]   { 0x2c, 0x2d, 0x2e, 0x2f } ;

    if (lamindx < 0)
	{
	Coupling::list.emplace_back (lambda) ;
	lamindx = Coupling::indx (lambda) ;
	Coupling::list[lamindx].value = 1000 ;
	}
    Coeff lamcoeff	{{lamindx,1}} ;
    Coeff laminvcoeff	{{lamindx,-1}} ;

    ObsList::freeze = false ;
    if (global.stage == 0)
	{
	global.info().Hterms.clear() ;

	if (!theory.euclid)
	    {
	    ObsPoly	kinetic  { baslist } ;
	    ObsPoly	ckinetic { obslist } ;
	    doub	coeff    { theory.dim > 1 ? 0.25 : 1.0 } ;
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		SymbStr	EE ( KinG + i ) ;
		uint	indx { baslist.find (EE) } ;
		if (indx == UINT_MAX) fatal ("Baselist missing EE") ;
		kinetic.push_back (PolyTerm (indx, coeff)) ;
		map.add (PolyTerm (obslist.catalog (baslist(indx)), coeff)) ;
		}
	    ckinetic.push_map (map) ;
	    global.info().Hterms.emplace_back (lamcoeff, kinetic, ckinetic) ;
	    }
	else
	    {
	    ObsPoly gauge_ent  (baslist) ;
	    ObsPoly cgauge_ent (obslist) ;
	    SymbStr S ( EntrG ) ;
	    uint    indx { baslist.find (S) } ;
	    if (indx == UINT_MAX) fatal ("Baselist missing S") ;
	    gauge_ent.push_back (PolyTerm(indx, -1.0)) ;
	    map.add (PolyTerm (obslist.catalog (baslist(indx)), -1.0)) ;
	    cgauge_ent.push_map (map) ;
	    global.info().Hterms.emplace_back (unitcoeff, gauge_ent, cgauge_ent) ;
	    }

	ObsPoly plaquette  (baslist) ;
	ObsPoly cplaquette (obslist) ;
	if (theory.dim == 1 && theory.box.comp[0])
	    {
	    plaquette.push_back  (PolyTerm(0, 2)) ;
	    cplaquette.push_back (PolyTerm(0, 2)) ;
	    SymbStr	polyakov { string (theory.box.comp[0], link[0]) } ;
	    SymbStr	Polyakov { string (theory.box.comp[0], Link[0]) } ;
	    uint	indx { baslist.find (polyakov) } ;
	    uint	Indx { baslist.find (Polyakov) } ;
	    if (indx == UINT_MAX) fatal ("Baselist missing polyakov") ;
	    if (Indx == UINT_MAX) fatal ("Baselist missing Polyakov") ;
	    plaquette.push_back (PolyTerm(indx, -1.0)) ;
	    plaquette.push_back (PolyTerm(Indx, -1.0)) ;
	    map.add (PolyTerm (obslist.catalog (baslist(indx)), -1.0)) ;
	    map.add (PolyTerm (obslist.catalog (baslist(Indx)), -1.0)) ;
	    }
	else
	    {
	    plaquette.push_back  (PolyTerm(0, theory.dim * (theory.dim-1))) ;
	    cplaquette.push_back (PolyTerm(0, theory.dim * (theory.dim-1))) ;
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		for (int j(i) ; ++j < theory.dim ;)
		    {
		    SymbStr plaq { string{link[i],link[j],Link[i],Link[j]} } ;
		    SymbStr Plaq { string{link[i],Link[j],Link[i],link[j]} } ;
		    uint    indx { baslist.find (plaq) } ;
		    uint    Indx { baslist.find (Plaq) } ;
		    if (indx == UINT_MAX) fatal ("Baselist missing plaq") ;
		    if (Indx == UINT_MAX) fatal ("Baselist missing Plaq") ;
		    plaquette.push_back (PolyTerm(indx, -1.0)) ;
		    plaquette.push_back (PolyTerm(Indx, -1.0)) ;
		    map.add (PolyTerm (obslist.catalog (baslist(indx)), -1.0)) ;
		    map.add (PolyTerm (obslist.catalog (baslist(Indx)), -1.0)) ;
		    }
		}
	    }
	cplaquette.push_map (map) ;
	global.info().Hterms.emplace_back (laminvcoeff, plaquette, cplaquette) ;
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
	    }
	Coeff masscoeff {{massindx,1},{lamindx,MASSSCALE}} ;
	global.info().Hterms.clear() ;
	
	if (!theory.euclid)
	    {
	    ObsPoly kinetic_F  (baslist) ;
	    ObsPoly ckinetic_F (obslist) ;
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		SymbStr ee ( KinF + i ) ;
		uint    indx { baslist.find (ee) } ;
		if (indx == UINT_MAX) fatal ("Baselist missing ee") ;
		kinetic_F.push_back (PolyTerm(indx, 0.25)) ;
		map.add (PolyTerm (obslist.catalog (baslist(indx)), 0.25)) ;
		}
	    ckinetic_F.push_map (map) ;
	    global.info().Hterms.emplace_back (lamcoeff, kinetic_F, ckinetic_F) ;
	    }

	ObsPoly hop_term   (baslist) ;
	ObsPoly chop_term  (obslist) ;
	for (int i(0) ; i < theory.nf ; i+=2)
	    {
	    for (int j(0) ; j < theory.dim ; ++j)
		{
		SymbStr Fxf { string{Ferm[i],link[j],ferm[i]} } ;
		SymbStr FXf { string{Ferm[i],Link[j],ferm[i]} } ;
		uint	indx { baslist.find(Fxf) } ;
		uint	Indx { baslist.find(FXf) } ;
		if (indx == UINT_MAX) fatal ("Baselist missing Fxf") ;
		if (Indx == UINT_MAX) fatal ("Baselist missing FXf") ;
		hop_term.push_back (PolyTerm(indx, 0.5)) ;
		hop_term.push_back (PolyTerm(Indx, 0.5)) ;
		map.add (PolyTerm (obslist.catalog (baslist(indx)), 0.5)) ;
		map.add (PolyTerm (obslist.catalog (baslist(Indx)), 0.5)) ;
		}
	    }
	chop_term.push_map (map) ;
	global.info().Hterms.emplace_back (unitcoeff, hop_term, chop_term) ;

	ObsPoly mass_term  (baslist) ;
	ObsPoly cmass_term (obslist) ;
	bool	iseuc	   { theory.euclid } ;

	for (int i(0) ; i < theory.nf ; i+=2)
	    {
	    uint indx ;
	    if (theory.euclid)
		{
		SymbStr Ff { string{Ferm[i],ferm[i]} } ;
		indx = baslist.find (Ff) ;
		if (indx == UINT_MAX) fatal ("Baselist missing Ff") ;
		}
	    else
		{
		SymbStr Gf { string{Ferm[i+1],ferm[i]} } ;
		indx = baslist.find (Gf) ;
		if (indx == UINT_MAX) fatal ("Baselist missing Gf") ;
		}
	    mass_term.push_back (PolyTerm(indx, 1.0)) ;
	    map.add (PolyTerm (obslist.catalog (baslist(indx)), 1.0)) ;
	    }
	cmass_term.push_map (map) ;
	global.info().Hterms.emplace_back (masscoeff, mass_term, cmass_term, iseuc) ;

	if (theory.euclid)
	    {
	    ObsPoly fermi_ent  (baslist) ;
	    ObsPoly cfermi_ent (obslist) ;
	    SymbStr s { EntrF } ;
	    uint    indx { baslist.find (s) } ;
	    if (indx == UINT_MAX) fatal ("Baselist missing s") ;
	    fermi_ent.push_back (PolyTerm(indx, -1.0)) ;
	    map.add (PolyTerm (obslist.catalog (baslist(indx)), -1.0)) ;
	    cfermi_ent.push_map (map) ;
	    global.info().Hterms.emplace_back (unitcoeff, fermi_ent, cfermi_ent) ;
	    }
	}
    ObsList::freeze = true ;
    global.info().MMAlimit = global.info().nobs ;
    global.info().MMAlist.clear() ;
    }
