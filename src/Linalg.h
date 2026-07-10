#ifndef LINALG
#define LINALG
#include "Gordion.h"
#include <complex>
#include <numeric>
#include <iostream>

#ifdef ARMA				// Use Armadillo linear algrbra library
#define ARMA_WARN_LEVEL 1
#include <armadillo>

using arma::abs ;
using arma::eig_gen ;
using arma::eig_pair ;
using arma::index_max ;
using arma::pinv ;
using arma::rank ;
using arma::sort ;
using arma::sqrt ;
using arma::svd ;

using Uvec = arma::Col<arma::uword> ;
using Rvec = arma::Col<real>  ;
using Dvec = arma::Col<doub>  ;
using Cvec = arma::Col<cmplx> ;
using Dmtx = arma::Mat<doub>  ;
using Cmtx = arma::Mat<cmplx> ;

template <typename T> uint ncol      (T& x)		{ return x.n_cols ; }
template <typename T> uint nrow      (T& x)		{ return x.n_rows ; }
template <typename T> uint nelem     (T& x)		{ return x.n_elem ; }
template <typename T> auto memptr    (T& x)		{ return x.memptr() ; }
template <typename T> void set_zero  (T&  x)		{ x.zeros() ; }
template <typename T> void set_zero  (T&& x)		{ x.zeros() ; }
template <typename T> void set_zero  (T&  x,int n)	{ x.zeros(n) ; }
template <typename T> void set_zero  (T&& x,int n)	{ x.zeros(n) ; }
template <typename T> void set_zero  (T& x,int n,int m)	{ x.zeros(n,m) ; }
template <typename T> void set_size  (T& x,int n)	{ x.set_size(n) ; }
template <typename T> void resize    (T& x,int n)	{ x.resize(n) ; }
template <typename T> void raw_print (T&  x,string s)	{ x.raw_print (s) ; }
template <typename T> void raw_print (T&& x,string s)	{ x.raw_print (s) ; }

inline auto opts	{ arma::solve_opts::refine + arma::solve_opts::equilibrate } ;
inline Dvec linsolve	(const Dmtx& a, const Dvec& x) { return arma::solve(a,x) ; }
inline doub infnorm	(const Dvec& v) { return arma::norm(v,"inf") ; }
inline doub l2norm	(const Dvec& v) { return arma::norm(v,2) ; }
inline Dvec realpart	(const Cvec& v) { return arma::real(v) ; }
inline Dvec imagpart	(const Cvec& v) { return arma::imag(v) ; }
inline bool has_nan	(const Cvec& v) { return v.has_nan() ; }
inline Dmtx inv		(const Dmtx& m) { return m.i() ; }
inline Dmtx transpose	(const Dmtx& m) { return m.t() ; }
inline Rvec aliasvec	(real* mem, uint n) { return Rvec (mem, n, false, true) ; }
inline Uvec sort_index	(const Dvec& v, bool up=true)
    {
    if (up)	return arma::sort_index (v,"ascend") ;
    else	return arma::sort_index (v,"descend") ;
    }

#elif EIGEN				// Use Eigen linear algebra library
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/SVD>

using Uvec = Eigen::Matrix<uint,  Eigen::Dynamic, 1> ;
using Rvec = Eigen::Matrix<real,  Eigen::Dynamic, 1> ;
using Dvec = Eigen::Matrix<doub,  Eigen::Dynamic, 1> ;
using Cvec = Eigen::Matrix<cmplx, Eigen::Dynamic, 1> ;
using Dmtx = Eigen::Matrix<doub,  Eigen::Dynamic, Eigen::Dynamic> ;
using Cmtx = Eigen::Matrix<cmplx, Eigen::Dynamic, Eigen::Dynamic> ;

template <typename T> auto ncol      (T& x)		{ return x.cols() ; }
template <typename T> auto nrow      (T& x)		{ return x.rows() ; }
template <typename T> auto nelem     (T& x)		{ return x.size() ; }
template <typename T> auto memptr    (T& x)		{ return x.data() ; }
template <typename T> void set_zero  (T&  x)		{ x.setZero() ; }
template <typename T> void set_zero  (T&& x)		{ x.setZero() ; }
template <typename T> void set_zero  (T&  x,int n)	{ x.setZero(n) ; }
template <typename T> void set_zero  (T&& x,int n)	{ x.setZero(n) ; }
template <typename T> void set_zero  (T& x,int n,int m)	{ x.setZero(n,m) ; }
template <typename T> void set_size  (T& x,int n)	{ x.resize(n) ; }
template <typename T> void resize    (T& x,int n)	{ x.conservativeResize(n) ; }
template <typename T> void raw_print (T&  x,string s)	{ std::cout<< s<< "\n"<< x<< "\n";}
template <typename T> void raw_print (T&& x,string s)	{ std::cout<< s<< "\n"<< x<< "\n";}

inline doub infnorm	(const Dvec& v) { return v.lpNorm<Eigen::Infinity>(); }
inline doub l2norm	(const Dvec& v) { return v.norm() ; }
inline Dvec realpart	(const Cvec& v) { return v.real() ; }
inline Dvec imagpart	(const Cvec& v) { return v.imag() ; }
inline Dvec sqrt	(const Dvec& v) { return v.cwiseSqrt() ; }
inline Cvec sqrt	(const Cvec& v) { return v.cwiseSqrt() ; }
inline Dvec abs		(const Dvec& v) { return v.cwiseAbs() ; }
inline bool has_nan	(const Cvec& v) { return v.hasNaN() ; }
inline Dmtx inv		(const Dmtx& m) { return m.inverse() ; }
inline Dmtx transpose	(const Dmtx& m) { return m.transpose() ; }
inline Rvec aliasvec	(doub* mem, uint n) { return Rvec {Eigen::Map<Rvec> (mem,n)} ; }

inline Cvec sort	(const Cvec& v)
    {
    Cvec w { v } ;
    std::sort (w.begin(),w.end(), [](cmplx a,cmplx b){ return std::abs(a) < std::abs(b);});
    return w ;
    }
inline Cvec sort	(const Cvec&& v)
    {
    return sort (v) ;
    }
inline Uvec sort_index	(const Dvec& v, bool up=true)
    {
    Uvec i	(v.size()) ; std::iota (i.begin(), i.end(), 0) ;
    auto lt	{ [&](doub j,doub k) { return v[j] < v[k] ;} } ;
    auto gt	{ [&](doub j,doub k) { return v[j] > v[k] ;} } ;
    if (up)	std::sort (i.begin(), i.end(), lt) ;
    else	std::sort (i.begin(), i.end(), gt) ;
    return i ;
    }
inline uint index_max	(const Dvec& v)
    {
    uint i ; v.maxCoeff (&i) ; return i ;
    }
inline Dvec svd		(const Dmtx& m)
    {
    Eigen::BDCSVD<Dmtx, Eigen::ComputeThinU | Eigen::ComputeThinV> dcomp (m) ;
    return dcomp.singularValues() ;
    }
inline uint rank	(const Dmtx& m)
    {
    Eigen::BDCSVD<Dmtx, Eigen::ComputeThinU | Eigen::ComputeThinV> dcomp (m) ;
    return dcomp.rank() ;
    }
inline Dmtx pinv	(const Dmtx& m, doub cut)
    {
    Eigen::BDCSVD<Dmtx, Eigen::ComputeFullU | Eigen::ComputeFullV> dcomp(m) ;
    Dvec sv = dcomp.singularValues() ;
    Dvec inv (sv.size()) ;
    for (int i = 0; i < sv.size(); ++i) inv(i) = (sv(i) > cut ? 1.0 / sv(i) : 0.0) ;
    return dcomp.matrixV() * inv.asDiagonal() * dcomp.matrixU().transpose() ;
    }
inline Cvec eig_gen	(const Dmtx& m)
    {
    Eigen::EigenSolver<Dmtx> es(m);
    return es.eigenvalues() ;
    }
inline bool eig_pair	(Cvec& val, Cmtx& vec, const Dmtx& a, const Dmtx& b)
    {
    Eigen::GeneralizedEigenSolver<Dmtx> ges(a,b) ;
    vec = ges.eigenvectors() ;
    val = ges.alphas().cwiseQuotient (ges.betas()) ;
    return true ;
    }
inline Dvec linsolve (const Dmtx& a, const Dvec& x)
    {
    return a.colPivHouseholderQr().solve(x) ;
    }

#endif
#endif
