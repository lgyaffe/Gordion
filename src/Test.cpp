#include "Test.h"
#include "Commute.h"
#include "Global.h"
#include "Rep.h"
#include "Gripe.h"

void Test::irreps ()					// Test irrep projector validity
    {
    Proj* preva { nullptr } ;
    for (const auto& a : Proj::list)
	{
	auto	[ai,aj] { a.indices () } ;
	auto	repa { a.repname() } ;

	for (const auto& b : Proj::list)
	    {
	    auto	[bi,bj] { b.indices () } ;
	    auto	repb { b.repname() } ;
	    Proj	ans ;

	    if (repa == repb && aj == bi) ans = *(&a + bj-aj) ;

	    Proj prod (a(b)) ;
	    Proj diff (prod - ans) ;
	    if (!diff.allzero())
		{
		cout << a.name << " * " << b.name << " error " << diff ;
		fatal (format("Inconsistent projectors {} * {}", a.name, b.name)) ;
		}
	    }
	}
    cout << "All projectors orthonormal\n" ;

    Proj v ;
    for (const auto& proj : Proj::list) if (proj[0]) v += proj ;
    v[0] -= v.denom ;
    cout << "Sum of diagonal projectors " << (v.allzero() ? "is" : "is not") << " complete\n" ;
    }

void Test::jacobi (const string& w1, const string& w2, numb obsindx) // Test specific Jacobi identity
    {
    ObsList	list	{ "JacobiTemp" } ;
    ObsPoly	poly	{ obsindx, ObsList::obs } ;
    Obs		obs	{ ObsList::obs(obsindx) } ;
    PolyMap	ans	{ list } ;
    PolyMap	tmp1	{ list } ;
    PolyMap	tmp2	{ list } ;
    Gen		g1	{ Op {w1,-1} } ;
    Gen		g2	{ Op {w2,-1} } ;
    bool	isF	( g1.is_Fermion() ^ g2.is_Fermion() ) ;
    auto&	oplist	{ global.info(isF).ops } ;
    Gen		g21	{ oplist } ;
    uint	opnum	( oplist.size() ) ;

    Commute::commute_poly (g2, poly, tmp1) ;
    Commute::commute_poly (g1, tmp1, ans) ;
    Commute::commute_poly (Commute::commute_gen (g2, g1, g21), poly, ans) ;
    Commute::commute_poly (g1, poly.negate(), tmp2) ;
    Commute::commute_poly (g2, tmp2, ans) ;
    oplist.purge (opnum) ;

    cout << " jacobi(" << w1 << "," << w2 << ',' << obs << ")" ;
    if (ans.allzero())	cout << " OK\n" ;
    else		cout << " = " << ans << "\n" ;
    }

void Test::jacobi (numb obsindx)				// Test Jacobi identities on Obs
    {
    auto&	gens	{ global.info().gens[global.repnum] } ;
    auto&	oplist	{ global.info().ops } ;
    uint	ngens	( gens.size() ) ;
    uint	opnum	( oplist.size() ) ;
    ObsList	list	{ "JacobiTemp" } ;
    Obs		obs	{ ObsList::obs(obsindx) } ;

    for (int i(0) ; i < ngens ; ++i)
	{
	Gen& g1 { gens[i] } ;

	for (int j(i) ; j < ngens ; ++j)
	    {
	    Gen&	g2	{ gens[j] } ;
	    Gen		g21	{ oplist } ;
	    ObsPoly	poly	{ obsindx, list } ;
	    PolyMap	ans	{ list } ;
	    PolyMap	tmp1	{ list } ;
	    PolyMap	tmp2	{ list } ;

	    Commute::commute_poly (g2, poly, tmp1) ;
	    Commute::commute_poly (g1, tmp1, ans) ;
	    Commute::commute_poly (Commute::commute_gen (g2, g1, g21), poly, ans) ;
	    Commute::commute_poly (g1, poly.negate(), tmp2) ;
	    Commute::commute_poly (g2, tmp2, ans) ;

	    cout << " jacobi(" << global.fg() << i << ","
			       << global.fg() << j << ',' << obs << ")" ;
	    if (ans.allzero())	cout << " OK\n" ;
	    else		cout << " = " << ans << "\n" ;
	    }
	}
    oplist.purge (opnum) ;
    }

void Test::jacobi()					// Test all Gen Jacobi identities
    {
    auto&	gens	{ global.info().gens[global.repnum] } ;
    auto&	oplist	{ global.info().ops } ;
    uint	ngens	( gens.size() ) ;
    uint	opnum	( oplist.size() ) ;

    for (int i(0) ; i < ngens ; ++i)
	{
	Gen& g1 { gens[i] } ;

	for (int j(i) ; j < ngens ; ++j)
	    {
	    Gen& g2 { gens[j] } ;

	    for (int k(j) ; k < ngens ; ++k)
		{
		Gen& g3   { gens[k] } ;
		Gen  ans  { oplist } ;
		Gen  tmp1 { oplist } ;
		Gen  tmp2 { oplist } ;
		Gen  tmp3 { oplist } ;
		Commute::commute_gen (g1, Commute::commute_gen (g2, g3, tmp1), ans) ;
		Commute::commute_gen (g2, Commute::commute_gen (g3, g1, tmp2), ans) ;
		Commute::commute_gen (g3, Commute::commute_gen (g1, g2, tmp3), ans) ;

		cout << " jacobi(" << global.fg() << i << ","
				   << global.fg() << j << ","
				   << global.fg() << k << ")" ;
		if (ans.allzero())	cout << " OK\n" << flush ;
		else 			cout << " = " << ans << "\n" << flush ;
		}
	    }
	}
    oplist.purge (opnum) ;
    }
