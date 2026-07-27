#ifndef GORDION_H
#define GORDION_H
#include <cstdint>
#include <set>
#include <array>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <algorithm>
#include <climits>
#include <version>

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

using std::pair ;			// STL classes, etc.
using std::array ;
using std::string ;
using std::vector ;
using std::complex ;
using std::unordered_map ;
using std::unordered_set ;
using std::exception ;

using symb   = char ;			// Basic symbol type
using doub   = double ;			// Spectral matrices, etc.
using cmplx  = complex<doub> ;		// Complex eigenvalues
using uchar  = uint8_t ;		// Hterm #
using ushort = uint16_t ;		// Gen #, Obs/Poly length
using uint   = uint32_t ;		// ObsList index
using ulong  = uint64_t ;		// Statistics counters
using uint2  = pair<uint,uint> ;	// Rep projector indices

#ifdef NUM32
using real   = float ;			// Poly coeffs & vev's
using numb   = int ;			// ObsList index
using usmall = ushort ;			// Half size integer
using small  = short ;			// Half size integer
constexpr int  MAXORD  = 63 ;		// Max Op/Obs order
constexpr long MAXOBS  = INT_MAX ;	// Max Obs index
constexpr doub DFLTTOL = 1.e-6 ;	// Default ODE tolerance
#else
using real   = doub ;			// Poly coeffs & vev's
using numb   = long ;			// ObsList index
using usmall = uint ;			// Half size integer
using small  = int ;			// Half size integer
constexpr int  MAXORD  = 1023 ;		// Max Op/Obs order
constexpr long MAXOBS  = LONG_MAX ;	// Max Obs index
constexpr doub DFLTTOL = 1.e-10 ;	// Default ODE tolerance
#endif

using numb2  = pair<numb,numb> ;	// Fermion initalization map
using numb3  = array<numb,3> ;		// Obs bucket position
using char8  = array<char,8> ;		// Theory/Coupling name
using strview = std::string_view ;	// Substrings

static_assert (sizeof (ulong) >= 4 * sizeof (short), "Need 8 byte longs!") ;
static_assert (sizeof (real) == sizeof (numb)) ;

template <typename Key, typename... Value>
    using hash = unordered_map<Key, Value...> ;	// unordered_map alias

template <typename T, typename U>		// reinterpret_cast alias
    T cast_to(U* ptr) { return reinterpret_cast<T>(ptr) ; }

static string program ;
static string cmdargs { " [-f <startup_file>] [sys-info or vev-data files]\n" } ;

static constexpr ulong ipow(ulong base, ulong exp, ulong ans = 1) // Integer power function
    {
    return exp < 1 ? ans : ipow(base*base, exp/2, (exp % 2) ? ans * base : ans) ;
    }

#ifdef __cpp_lib_format
#include <format>
using std::format ;
#else
#include <sstream>
#include <string_view>

template <typename T>
void format_helper(ostringstream& oss, std::string_view& str, const T& value)
    {
    std::size_t openBracket = str.find('{');
    if (openBracket == std::string::npos) { return; }
    std::size_t closeBracket = str.find('}', openBracket + 1);
    if (closeBracket == std::string::npos) { return; }
    oss << str.substr(0, openBracket) << value;
    str = str.substr(closeBracket + 1);
    }

template <typename... Targs>
string format(std::string_view str, Targs...args)
    {
    ostringstream oss;
    (format_helper(oss, str, args),...);
    oss << str;
    return oss.str();
    }

#endif
#endif
