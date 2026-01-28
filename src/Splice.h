#ifndef SPLICE_H
#define SPLICE_H

static constexpr char Y { -1 } ;

struct Splice					// Commute splice info
    {
    signed char	coeff  = 0 ;
    signed char	prefix = Y ;
    signed char	suffix = Y ;
    } ;

static const Splice splicetbl[10][10][2]	// Commute splice table
    //
    // 	Same tnR codes as in Symb::ligtable:
    //		0 = link,	1 = Link
    //		2 = E,		3 = LEl
    //		4 = El,		5 = LE
    //		6 = EE,		7 = LEEl
    //		8 = EEl,	9 = LEE
    //		Y = -1 := omit symb
    {{
	{{  0, Y, Y }, {  0, Y, Y }}	// [ link, link ]	-
    },{
    	{{  0, Y, Y }, {  0, Y, Y }},	// [ Link, link ]	-
    	{{  0, Y, Y }, {  0, Y, Y }}	// [ Link, Link ]	-
    },{
    	{{ +1, Y, 0 }, {  0, Y, Y }},	// [ E, link ]		ok
    	{{ -1, 1, Y }, {  0, Y, Y }},	// [ E, Link ]		ok
	{{ +1, Y, 2 }, {  0, Y, Y }}	// [ E, E ]		ok
    },{
    	{{ +1, 0, Y }, {  0, Y, Y }},	// [ LEl, link ]	ok
    	{{ -1, Y, 1 }, {  0, Y, Y }},	// [ LEl, Link ]	ok
	{{  0, Y, Y }, {  0, Y, Y }},	// [ LEl, E ]		ok
	{{ +1, 3, Y }, {  0, Y, Y }}	// [ LEl, LEl ]		ok
    },{
    	{{ +1, 0, 0 }, {  0, Y, Y }},	// [ El, link ]		ok
    	{{ -1, Y, Y }, {  0, Y, Y }},	// [ El, Link ]		ok
	{{ -1, 4, Y }, {  0, Y, Y }},	// [ El, E ]		ok
	{{ -1, Y, 4 }, {  0, Y, Y }},	// [ El, LEl ]		ok
	{{  0, Y, Y }, {  0, Y, Y }},	// [ El, El ]		ok
	{{  0, Y, Y }, {  0, Y, Y }},	// [ El, LE ]		-
	{{  0, Y, Y }, {  0, Y, Y }},	// [ El, EE ]		-
	{{  0, Y, Y }, {  0, Y, Y }},	// [ El, LEEl ]		-
	{{  0, Y, Y }, {  0, Y, Y }},	// [ El, EEl ]		-
	{{ -1, Y, 6 }, {  0, Y, Y }}	// [ El, LEE ]		ok
    },{
    	{{ +1, Y, Y }, {  0, Y, Y }},	// [ LE, link ]		ok
    	{{ -1, 1, 1 }, {  0, Y, Y }},	// [ LE, Link ]		ok
	{{ +1, Y, 5 }, {  0, Y, Y }},	// [ LE, E ]		ok
	{{ +1, 5, Y }, {  0, Y, Y }},	// [ LE, LEl ]		ok
	{{ +1, Y, 3 }, { +1, 2, Y }},	// [ LE, El ]		ok
	{{  0, Y, Y }, {  0, Y, Y }},	// [ LE, LE ]		ok
	{{  0, Y, Y }, {  0, Y, Y }},	// [ LE, EE ]		-
	{{  0, Y, Y }, {  0, Y, Y }},	// [ LE, LEEl ]		-
	{{ +1, 6, Y }, {  0, Y, Y }}	// [ LE, EEl ]		ok
    },{
    	{{ +1, 2, 0 }, { +1, Y, 4 }},	// [ EE, link ]		ok
    	{{ -1, 1, 2 }, { -1, 5, Y }},	// [ EE, Link ]		ok
	{{ +1, Y, 6 }, { -1, 6, Y }},	// [ EE, E ]		ok
	{{  0, Y, Y }, {  0, Y, Y }},	// [ EE, LEl ]		ok
	{{ +1, Y, 8 }, { +1, 2, 4 }},	// [ EE, El ]		ok
	{{ -1, 9, Y }, { -1, 5, 2 }}	// [ EE, LE ]		ok
    },{
    	{{ +1, 4, Y }, { +1, 0, 3 }},	// [ LEEl, link ]	ok
    	{{ -1, Y, 5 }, { -1, 3, 1 }},	// [ LEEl, Link ]	ok
	{{  0, Y, Y }, {  0, Y, Y }},	// [ LEEl, E ]		ok
	{{ +1, 7, Y }, { -1, Y, 7 }},	// [ LEEl, LEl ]	ok
	{{ +1, 4, 3 }, { +1, 8, Y }},	// [ LEEl, El ]		ok
	{{ -1, Y, 9 }, { -1, 3, 5 }}	// [ LEEl, LE ]		ok
    },{
    	{{ +1, 4, 0 }, { +1, 0, 4 }},	// [ EEl, link ]	ok
    	{{ -1, Y, 2 }, { -1, 3, Y }},	// [ EEl, Link ]	ok
	{{ -1, 8, Y }, {  0, Y, Y }},	// [ EEl, E ]		ok
	{{ -1, Y, 8 }, {  0, Y, Y }},	// [ EEl, LEl ]		ok
	{{ +1, 4, 4 }, {  0, Y, Y }},	// [ EEl, El ]		ok
	{{ -1, 3, 2 }, { -1, 7, Y }},	// [ EEl, LE ]		ok
    },{
    	{{ +1, 2, Y }, { +1, Y, 3 }},	// [ LEE, link ]	ok
    	{{ -1, 1, 5 }, { -1, 5, 1 }},	// [ LEE, Link ]	ok
	{{ +1, Y, 9 }, {  0, Y, Y }},	// [ LEE, E ]		ok
	{{ +1, 9, Y }, {  0, Y, Y }},	// [ LEE, LEl ]		ok
	{{ +1, 2, 3 }, { +1, Y, 7 }},	// [ LEE, El ]		ok
	{{ -1, 5, 5 }, {  0, Y, Y }}	// [ LEE, LE ]		ok
    }} ;

#endif
