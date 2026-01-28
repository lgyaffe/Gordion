#include "RK.h"

namespace RK				// Runge-Kutta method definitions
    {
    using dvec = const double[] ;

    double constexpr do_sqrt (double x, double curr, double prev)
	{ return curr == prev ? curr : do_sqrt (x, 0.5 * (curr + x/curr), curr) ; }

    double constexpr mysqrt (double x)			// constexpr sqrt function
	{ return x >= 0 && x < std::numeric_limits<double>::infinity()
	     ? do_sqrt(x,x,0) : std::numeric_limits<double>::quiet_NaN() ; }

    //
    // Bogacki & Shampine 3(2) method, "A 3(2) pair of Runge-Kutta formulas"
    // Applied Mathematics Letters, Volume 2, 1989, pg 321
    // https://www.sciencedirect.com/science/article/pii/0893965989900797
    //
    static constexpr char bs32_n[] { "bs32: Bogacki/Shampine 3(2)" } ;
    static constexpr dvec bs32_c  { 0.,     1./2,   3./4, 1. } ;
    static constexpr dvec bs32_a2 { 1./2 } ;
    static constexpr dvec bs32_a3 { 0.,     3./4 } ;
    static constexpr dvec bs32_a4 { 2./9,   1./3,   4./9 } ;
    static constexpr dvec bs32_b  { 2./9,   1./3,   4./9,  0. } ;
    static constexpr dvec bs32_bh { 7./24,  1./4,   1./3,  1./8 } ;
    static constexpr dvec bs32_db { 5./72, -1./12, -1./9,  1./8 } ;
    static constexpr const double* const bs32_a[] { bs32_a2, bs32_a3, bs32_a4 } ;

    //
    // Sofroniou & Spaletta 3(2) method
    // Construction of explicit runge-kutta pairs with stiffness detection
    // Mathematical and Computer Modelling, Volume 40, December 2004, pg 1157
    // https://www.sciencedirect.com/science/article/pii/S0895717705000117
    // https://reference.wolfram.com/language/tutorial/NDSolveExplicitRungeKutta.html
    //
    static constexpr char ss32_n[] { "ss32: Sofroniou/Spaletta 3(2)" } ;
    static constexpr dvec ss32_c  {  0.,   1./2, 1., 1. } ;
    static constexpr dvec ss32_a2 {  1./2 } ;
    static constexpr dvec ss32_a3 { -1.,   2. } ;
    static constexpr dvec ss32_a4 {  1./6, 2./3, 1./6 } ;
    static constexpr dvec ss32_b  {  1./6, 2./3, 1./6, 0. } ;
    static constexpr dvec ss32_bh { ( 22-mysqrt(82.))/72, ( 14+mysqrt(82.))/36,
				   ( -4+mysqrt(82.))/144,( 16-mysqrt(82.))/48 } ;
    static constexpr dvec ss32_db { ( 10-mysqrt(82.))/72, (-10+mysqrt(82.))/36,
				   (-28+mysqrt(82.))/144,( 16-mysqrt(82.))/48 } ;
    static constexpr const double* const ss32_a[] { ss32_a2, ss32_a3, ss32_a4 } ;

    //
    // Sofroniou & Spaletta 4(3) method
    // "Construction of explicit runge-kutta pairs with stiffness detection"
    // Mathematical and Computer Modelling, Volume 40, December 2004, pg 1157
    // https://www.sciencedirect.com/science/article/pii/S0895717705000117
    // https://reference.wolfram.com/language/tutorial/NDSolveExplicitRungeKutta.html
    //
    static constexpr char ss43_n[] { "ss43: Sofroniou/Spaletta 4(3)" } ;
    static constexpr dvec ss43_c  {  0.,      2./5,   3./5,  1.,  1. } ;
    static constexpr dvec ss43_a2 {  2./5 } ;
    static constexpr dvec ss43_a3 { -3./20,   3./4 } ;
    static constexpr dvec ss43_a4 { 19./44, -15./44, 10./11 } ;
    static constexpr dvec ss43_a5 { 11./72,  25./72, 25./72, 11./72 } ;
    static constexpr dvec ss43_b  { 11./72,  25./72, 25./72, 11./72, 0. } ;
    static constexpr dvec ss43_bh { 1251515./8970912, 3710105./8970912,
				   2519695./8970912,   61105./8970912, 119041./747576 } ;
    static constexpr dvec ss43_db { -119041./8970912,  595205./8970912,
				   -595205./8970912, -1309451./8970912, 119041./747576 } ;
    static constexpr const double* const ss43_a[] { ss43_a2, ss43_a3, ss43_a4, ss43_a5 } ;

    //
    // Dormand & Prince RK5(4)7M, "A family of embedded Runge-Kutta formulae"
    // Journal of Computational and Applied Mathematics, Volume 6, March 1980, 
    // https://www.sciencedirect.com/science/article/pii/0771050X80900133
    // Also in: J. Butcher, Numerical Methods for Ordinary Differential Equations 
    // 3rd ed, 2016, pg 224, ISBN: 9781119121503, EISBN: 9781119121527
    //
    static constexpr char dp54_n[] { "dp54: Dormand/Prince 5(4)7M" } ;
    static constexpr dvec dp54_c  {     0.,     1./5,    3./10,      4./5,        8./9,        1.,     1.  } ;
    static constexpr dvec dp54_a2 {     1./5 } ;
    static constexpr dvec dp54_a3 {     3./40,           9./40 } ;
    static constexpr dvec dp54_a4 {    44./45,         -56./15,      32./9 } ;
    static constexpr dvec dp54_a5 { 19372./6561,    -25360./2187, 64448./6561,  -212./729 } ;
    static constexpr dvec dp54_a6 {  9017./3168,      -355./33,   46732./5247,    49./176,  -5103./18656 } ;
    static constexpr dvec dp54_a7 {    35./384,   0.,  500./1113,   125./192,  -2187./6784,    11./84 } ;
    static constexpr dvec dp54_b  {    35./384,   0.,  500./1113,   125./192,  -2187./6784,    11./84,  0. };
    static constexpr dvec dp54_bh {  5179./57600, 0., 7571./16695 , 393./640, -92097./339200, 187./2100, 1./40 } ;
    static constexpr dvec dp54_db {   -71./57600, 0.,   71./16695 , -71./1920, 17253./339200, -22./525,  1./40 } ;
    static constexpr const double* const dp54_a[] { dp54_a2, dp54_a3, dp54_a4, dp54_a5, dp54_a6, dp54_a7 } ;

    //
    // J. Butcher, 5(4) method, Numerical methods for ordinary differential equations 
    // Wiley, 3rd ed, 2016, pg 225, ISBN: 9781119121503, EISBN: 9781119121527
    //
    static constexpr char jb54_n[] { "jb54: John Butcher 5(4)" } ;
    static constexpr dvec jb54_c  {   0.,     2./9,  1./3,       5./9,       2./3,      1.,      1.  } ;
    static constexpr dvec jb54_a2 {   2./9 } ;
    static constexpr dvec jb54_a3 {   1./12,    1./4 } ;
    static constexpr dvec jb54_a4 {  55./324, -25./108, 50./81 } ;
    static constexpr dvec jb54_a5 {  83./330, -13./22,  61./66,      9./110 } ;
    static constexpr dvec jb54_a6 { -19./28,    9./4,    1./7,     -27./7,      22./7 } ;
    static constexpr dvec jb54_a7 {  19./200,   0.,      3./5,    -243./400,    33./40,     7./80 } ;
    static constexpr dvec jb54_b  {  19./200,   0.,      3./5,    -243./400,    33./40,     7./80,    0. } ;
    static constexpr dvec jb54_bh { 431./5000,  0.,    333./500, -7857./10000, 957./1000, 193./2000, -1./50 } ;
    static constexpr dvec jb54_db { -11./1250,  0.,     33./500,  -891./5000,   33./250,    9./1000, -1./50 } ;
    static constexpr const double* const jb54_a[] { jb54_a2, jb54_a3, jb54_a4, jb54_a5, jb54_a6, jb54_a7 } ;

    //
    // Bogacki & Shampine 5(4) method
    // An efficient Runge-Kutta (4,5) pair, Rept. 89-20, Math. Dept.
    // Southern Methodist University, Dallas, Texas, USA, 1989
    // https://www.sciencedirect.com/science/article/pii/0898122196001411
    // https://www.netlib.org/ode/rksuite/rksuite.f (lines 2527-2629
    // http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK5/RKcoeff5p_1.pdf
    //
    static constexpr char bs54_n[] { "bs54: Bogacki/Shampine 5(4)" } ;
    static constexpr dvec bs54_c  { 0., 1./6, 2./9, 3./7, 2./3, 3./4, 1., 1. } ;
    static constexpr dvec bs54_a2 {   1./6 } ;
    static constexpr dvec bs54_a3 {   2./27,       4./27 } ;
    static constexpr dvec bs54_a4 { 183./1372,  -162./343,  1053./1372 } ;
    static constexpr dvec bs54_a5 {  68./297,     -4./11,     42./143,     1960./3861 } ;
    static constexpr dvec bs54_a6 { 597./22528,   81./352, 63099./585728, 58653./366080, 4617./20480 } ;
    static constexpr dvec bs54_a7 { 174197./959244,   -30942./79937,  8152137./19744439,
				   666106./1039181,   -29421./29068,   482048./414219 } ;
    static constexpr dvec bs54_a8 {   587./8064,  0, 4440339./15491840, 24353./124800,
				      387./44800,       2152./5985,      7267./94080 } ;
    static constexpr dvec bs54_b  {   587./8064,  0, 4440339./15491840, 24353./124800,
				      387./44800,       2152./5985,      7267./94080,      0. } ;
    static constexpr dvec bs54_bh {  2479./34992, 0.,    123./416,     612941./3411720,
				       43./1440,         2272./6561,    79937./1113912, 3293./556956 } ;
    static constexpr dvec bs54_db  { -3817./1959552, 0., 140181./15491840, -4224731./272937600,
				      8557./403200,  -57928./4363065,  -23930231./4366535040, 3293./556956  } ;
    static constexpr const double* const bs54_a[] { bs54_a2, bs54_a3, bs54_a4, bs54_a5, bs54_a6, bs54_a7, bs54_a8 } ;

    //
    // J. Verner RK6(5) IIIX-65 method, "Numerically optimal RungeKutta pairs with interpolants"
    // Numerical algorithms, 2010-03, Vol.53 (2-3), p.383-396; Boston: Springer US
    // https://www.sfu.ca/~jverner/
    //
    static constexpr char jv65_n[] { "jv65: Jim Verner 6(5)" } ;
    static constexpr dvec jv65_c
	{
	/* c[1] = */  0 ,
	/* c[2] = */  .6e-1 ,
	/* c[3] = */  .9593333333333333333333333333333333333333e-1 ,
	/* c[4] = */  .1439 ,
	/* c[5] = */  .4973 ,
	/* c[6] = */  .9725 ,
	/* c[7] = */  .9995 ,
	/* c[8] = */  1 ,
	/* c[9] = */  1
	} ;
    static constexpr dvec jv65_a2
	{
	/* a[2,1] = */  .6e-1
	} ;
    static constexpr dvec jv65_a3
	{
	/* a[3,1] = */  .1923996296296296296296296296296296296296e-1 ,
	/* a[3,2] = */  .7669337037037037037037037037037037037037e-1
	} ;
    static constexpr dvec jv65_a4
	{
	/* a[4,1] = */  .35975e-1 ,
	/* a[4,2] = */  0 ,
	/* a[4,3] = */  .107925
	} ;
    static constexpr dvec jv65_a5
	{
	/* a[5,1] = */  1.318683415233148260919747276431735612861 ,
	/* a[5,2] = */  0 ,
	/* a[5,3] = */ -5.042058063628562225427761634715637693344 ,
	/* a[5,4] = */  4.220674648395413964508014358283902080483
	} ;
    static constexpr dvec jv65_a6
	{
	/* a[6,1] = */ -41.87259166432751461803757780644346812905 ,
	/* a[6,2] = */  0 ,
	/* a[6,3] = */  159.4325621631374917700365669070346830453 ,
	/* a[6,4] = */ -122.1192135650100309202516203389242140663 ,
	/* a[6,5] = */  5.531743066200053768252631238332999150076
	} ;
    static constexpr dvec jv65_a7
	{
	/* a[7,1] = */ -54.43015693531650433250642051294142461271 ,
	/* a[7,2] = */  0 ,
	/* a[7,3] = */  207.0672513650184644273657173866509835987 ,
	/* a[7,4] = */ -158.6108137845899991828742424365058599469 ,
	/* a[7,5] = */  6.991816585950242321992597280791793907096 ,
	/* a[7,6] = */ -.1859723106220323397765171799549294623692e-1
	} ;
    static constexpr dvec jv65_a8
	{
	/* a[8,1] = */ -54.66374178728197680241215648050386959351 ,
	/* a[8,2] = */  0 ,
	/* a[8,3] = */  207.9528062553893734515824816699834244238 ,
	/* a[8,4] = */ -159.2889574744995071508959805871426654216 ,
	/* a[8,5] = */  7.018743740796944434698170760964252490817 ,
	/* a[8,6] = */ -.1833878590504572306472782005141738268361e-1 ,
	/* a[8,7] = */ -.5119484997882099077875432497245168395840e-3
	} ;
    static constexpr dvec jv65_a9
	{
	/* a[9,1] = */  .3438957868357036009278820124728322386520e-1 ,
	/* a[9,2] = */  0 ,
	/* a[9,3] = */  0 ,
	/* a[9,4] = */  .2582624555633503404659558098586120858767 ,
	/* a[9,5] = */  .4209371189673537150642551514069801967032 ,
	/* a[9,6] = */  4.405396469669310170148836816197095664891 ,
	/* a[9,7] = */ -176.4831190242986576151740942499002125029 ,
	/* a[9,8] = */  172.3641334014150730294022582711902413315
	} ;
    static constexpr dvec jv65_b
	{
	/* b[1] = */  .3438957868357036009278820124728322386520e-1 ,
	/* b[2] = */  0 ,
	/* b[3] = */  0 ,
	/* b[4] = */  .2582624555633503404659558098586120858767 ,
	/* b[5] = */  .4209371189673537150642551514069801967032 ,
	/* b[6] = */  4.405396469669310170148836816197095664891 ,
	/* b[7] = */ -176.4831190242986576151740942499002125029 ,
	/* b[8] = */  172.3641334014150730294022582711902413315 ,
	/* b[9] = */  0
	} ;
    static constexpr dvec jv65_bh
	{
	/* bh[ 1] = */  .4909967648382489730906854927971225836479e-1 ,
	/* bh[ 2] = */  0 ,
	/* bh[ 3] = */  0 ,
	/* bh[ 4] = */  .2251112229516524153401395320539875329485 ,
	/* bh[ 5] = */  .4694682253029562039431948525047387412553 ,
	/* bh[ 6] = */  .8065792249988867707634161808995217981443 ,
	/* bh[ 7] = */  0 ,
	/* bh[ 8] = */ -.6071194891777959797672951465256217122488 ,
	/* bh[ 9] = */  .5686113944047569241147603178766138153594e-1
	} ;
    static constexpr dvec jv65_db
	{
	/* bh[ 1] = */  0.01471009780025453721628034803242903449959 ,
	/* bh[ 2] = */  0 ,
	/* bh[ 3] = */  0 ,
	/* bh[ 4] = */  -0.0331512326116979251258162778046245529282 ,
	/* bh[ 5] = */  0.0485311063356024888789397010977585445521 ,
	/* bh[ 6] = */  -3.5988172446704233993854206352975738667467 ,
	/* bh[ 7] = */  176.4831190242986576151740942499002125029 ,
	/* bh[ 8] = */ -172.9712528905928690091695534177158630437488 ,
	/* bh[ 9] = */  .5686113944047569241147603178766138153594e-1
	} ;
    static constexpr const double* const jv65_a[] { jv65_a2, jv65_a3, jv65_a4, jv65_a5, jv65_a6, jv65_a7, jv65_a8, jv65_a9 } ;
    } ;

const std::array<RKdef,7> RKdef::list
    {
    RKdef { 5, 4, RK::ss43_c, RK::ss43_a, RK::ss43_b, RK::ss43_bh, RK::ss43_db, RK::ss43_n },
    RKdef { 7, 5, RK::jb54_c, RK::jb54_a, RK::jb54_b, RK::jb54_bh, RK::jb54_db, RK::jb54_n },
    RKdef { 7, 5, RK::dp54_c, RK::dp54_a, RK::dp54_b, RK::dp54_bh, RK::dp54_db, RK::dp54_n },
    RKdef { 8, 5, RK::bs54_c, RK::bs54_a, RK::bs54_b, RK::bs54_bh, RK::bs54_db, RK::bs54_n },
    RKdef { 4, 3, RK::bs32_c, RK::bs32_a, RK::bs32_b, RK::bs32_bh, RK::bs32_db, RK::bs32_n },
    RKdef { 9, 6, RK::jv65_c, RK::jv65_a, RK::jv65_b, RK::jv65_bh, RK::jv65_db, RK::jv65_n },
    RKdef { 4, 3, RK::ss32_c, RK::ss32_a, RK::ss32_b, RK::ss32_bh, RK::ss32_db, RK::ss32_n }
    } ;
