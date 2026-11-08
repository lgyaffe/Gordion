#ifndef TYPES_H
#define TYPES_H
#include <cstdint>
#include <set>
#include <list>
#include <span>
#include <array>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <string_view>
#include <iomanip>
#include <complex>
#include <algorithm>
#include <climits>
#include <version>
#include <chrono>

using std::ios ;			// I/O classes & functions
using std::cout ;
using std::flush ;
using std::fstream ;
using std::istream ;
using std::ostream ;
using std::ifstream ;
using std::ofstream ;
using std::istringstream ;
using std::ostringstream ;
using std::setprecision ;

using std::pair ;			// STL classes, etc.
using std::list ;
using std::span ;
using std::tuple ;
using std::array ;
using std::string ;
using std::string_view ;
using std::vector ;
using std::complex ;
using std::unordered_map ;
using std::unordered_set ;

using std::exception ;			// Other standard types
using std::chrono::milliseconds ;
using std::chrono::high_resolution_clock ;

using symb   = char ;			// Basic symbol type
using doub   = double ;			// Spectral matrices, etc.
using cmplx  = complex<float> ;		// Complex matrix elements
using dcmplx = complex<doub> ;		// Complex eigenvalues
using uchar  = uint8_t ;		// Hterm #
using ushort = uint16_t ;		// Gen #, Obs/Poly length
using uint   = uint32_t ;		// various indices
using ulong  = uint64_t ;		// Statistics counters
using int2   = pair<int,int> ;		// Canonicalization trials
using uint2  = pair<uint,uint> ;	// Rep projector indices
using char8  = array<char,8> ;		// Theory/Coupling name
using sv     = string_view ;		// Constexpr strings

#ifdef NUM64
using numb   = long ;			// ObsList index
using real   = doub ;			// Poly coeffs & vev's
using usmall = uint ;			// Half size integer
inline constexpr int  MAXORD = 1023 ;	// Max Op/Obs order
inline constexpr long MAXNUM = LONG_MAX; // Max Obs index
#else
using numb   = int ;			// Observable indices
using real   = float ;			// Poly coeffs & matrix elements
using usmall = ushort ;			// Half size integer
inline constexpr int  MAXGRP = 63 ;	// Mas Op group
inline constexpr int  MAXORD = 63 ;	// Max Op/Obs order
inline constexpr long MAXNUM = INT_MAX; // Max Obs index
#endif
using numb2  = pair<numb,numb> ;	// Index range, fermion init
using numb3  = array<numb,3> ;		// Obs bucket position

static_assert (sizeof (ulong) >= 4 * sizeof (short), "Need 8 byte longs!") ;
static_assert (sizeof (real) == sizeof (numb)) ;

template <typename T, typename U>	// reinterpret_cast alias
    T cast_to(U* ptr) { return reinterpret_cast<T>(ptr) ; }

inline constexpr ulong ipow(ulong base, ulong exp, ulong ans = 1) // Integer powers
    {
    return exp < 1 ? ans : ipow(base*base, exp/2, (exp % 2) ? ans * base : ans) ;
    }

#ifdef __cpp_lib_format
#include <format>
using std::format ;
#else
#include <sstream>
template <typename T>
void apply_format_specs(ostringstream& oss, sv spec)
    {
    if (spec.empty()) return ;
    oss << std::defaultfloat;
    size_t dot_pos { spec.find('.') } ;
    if (dot_pos != sv::npos)	// decimal point
	{
	if (dot_pos > 0)
	    {
	    int width = 0;
	    for (size_t i = 0; i < dot_pos; ++i)
		{
		if (std::isdigit(spec[i])) width = width * 10 + (spec[i] - '0') ;
		}
	    if (width > 0) oss << std::setw(width) ;
	    }
	int precision { 0 } ;
	size_t i { dot_pos + 1 } ;
	for (; i < spec.length() && std::isdigit(spec[i]) ; ++i)
	    precision = precision * 10 + (spec[i] - '0') ;
	oss << std::setprecision(precision);
	if (i < spec.length())
	    {
	    if      (spec[i] == 'f') oss << std::fixed ;
	    else if (spec[i] == 'e') oss << std::scientific ;
	    }
	}
    else			// no decimal point
	{
	int width { 0 } ;
	size_t i { 0 } ;
	for (; i < spec.length() && std::isdigit(spec[i]) ; ++i)
	    width = width * 10 + (spec[i] - '0') ;
	if (width > 0) oss << std::setw(width) ;
	if (i < spec.length() && spec[i] == 'f') oss << std::fixed;
	}
    }
template <typename T>
void format_helper(ostringstream& oss, sv& str, const T& value)
    {
    std::size_t openBracket = str.find('{');
    if (openBracket == string::npos) { return; }
    std::size_t closeBracket = str.find('}', openBracket + 1);
    if (closeBracket == string::npos) { return; }
    oss << str.substr(0, openBracket) ;
    sv spec = str.substr (openBracket+1, closeBracket - openBracket -1) ;
    if (!spec.empty() && spec[0] == ':') spec.remove_prefix(1) ;
    apply_format_specs<T>(oss, spec) ;
    oss << value ;
    str = str.substr(closeBracket + 1);
    }
template <typename... Targs>
string format(sv str, Targs...args)
    {
    ostringstream oss;
    (format_helper(oss, str, args),...);
    oss << str;
    return oss.str();
    }
#endif

#endif
