#include "Test.h"
#include "Global.h"
#include "Commute.h"
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

void Test::jacobi (string w1, string w2, uint obsindx)	 // Test specfic Jacobi identity
    {
    ObsList	list  { "JacobiTemp" } ;
    auto	opnum { Op::list.size() } ;
    Obs		obs   { ObsList::obs(obsindx) } ;
    ObsPoly	poly  { obs, list } ;
    PolyMap	ans   { list } ;
    PolyMap	tmp1  { list } ;
    PolyMap	tmp2  { list } ;
    Gen		g1    { Op {w1,-1} } ;
    Gen		g2    { Op {w2,-1} } ;
    Gen		g21 ;

    Commute::commute_poly (g2, poly, tmp1) ;
    Commute::commute_poly (g1, tmp1, ans) ;
    Commute::commute_poly (Commute::commute_gen (g2, g1, g21), poly, ans) ;
    Commute::commute_poly (g1, poly.negate(), tmp2) ;
    Commute::commute_poly (g2, tmp2, ans) ;
    Op::purge (opnum) ;

    cout << " jacobi(" << w1 << "," << w2 << ',' << obs << ")" ;
    if (ans.allzero())	cout << " OK\n" ;
    else		cout << " = " << ans << "\n" ;
    }

void Test::jacobi (uint obsindx)				// Test Jacobi identities on Obs
    {
    auto&	gens  { global.info().gens[global.repnum] } ;
    auto	opnum { Op::list.size() } ;
    uint	ngens ( gens.size() ) ;
    ObsList	list  { "JacobiTemp" } ;
    Obs		obs   { ObsList::obs(obsindx) } ;

    for (int i(0) ; i < ngens ; ++i)
	{
	Gen& g1 { gens[i] } ;

	for (int j(i) ; j < ngens ; ++j)
	    {
	    Gen&	g2 { gens[j] } ;
	    ObsPoly	poly  { obs, list } ;
	    PolyMap	ans   { list } ;
	    PolyMap	tmp1  { list } ;
	    PolyMap	tmp2  { list } ;
	    Gen		g21 ;

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
    Op::purge (opnum) ;
    }

void Test::jacobi()					// Test all Gen Jacobi identities
    {
    auto&	gens  { global.info().gens[global.repnum] } ;
    auto	opnum { Op::list.size() } ;
    uint	ngens ( gens.size() ) ;

    for (int i(0) ; i < ngens ; ++i)
	{
	Gen& g1 { gens[i] } ;

	for (int j(i) ; j < ngens ; ++j)
	    {
	    Gen& g2 { gens[j] } ;

	    for (int k(j) ; k < ngens ; ++k)
		{
		Gen& g3 { gens[k] } ;
		Gen  tmp1, tmp2, tmp3, ans ;
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
    Op::purge (opnum) ;
    }
