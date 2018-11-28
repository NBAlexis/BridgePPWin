/*!
        @file    $Id: bridge_complex.h #$

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2017-10-05 17:07:53 #$

        @version $LastChangedRevision: 1668 $
*/

#ifndef BRIDGE_COMPLEX_INCLUDED
#define BRIDGE_COMPLEX_INCLUDED

//- chose the complex type
// std::complex (c++) -> USE_STD_COMPLEX
// _Complex (C99)     -> USE_C99_COMPLEX
//
// When you make original complex class, you have to implement such functions:
// cmplx(T real, T imag), return your complex type, real + i imag,
// real(complex& z), imag(complex& z), return real(imaginary) part of the complex number z.
// abs(complex& z), return absolute value |z|,
// norm(complex& z), return squared norm |z|^2,
// arg(complex& z), return phase angle of the complex number z,
// conj(complex& z), return complex conjugate of the complex number z.

#if defined USE_STD_COMPLEX
// std::complex
#include <complex>
using std::abs;

typedef std::complex<int>      icomplex;
typedef std::complex<float>    fcomplex;
typedef std::complex<double>   dcomplex;
// typedef std::complex<long double> lcomplex;

inline icomplex cmplx(const int r, const int i) { return icomplex(r, i); }
inline fcomplex cmplx(const float r, const float i) { return fcomplex(r, i); }
inline dcomplex cmplx(const double r, const double i) { return dcomplex(r, i); }
// inline lcomplex cmplx (const long double r, const long double i) {
//   return lcomplex(r, i);
// }

#ifdef LIB_CPP11
// do nothing
#else
inline int real(const icomplex& z) { return z.real(); }
inline int imag(const icomplex& z) { return z.imag(); }

inline float real(const fcomplex& z) { return z.real(); }
inline float imag(const fcomplex& z) { return z.imag(); }

inline double real(const dcomplex& z) { return z.real(); }
inline double imag(const dcomplex& z) { return z.imag(); }

// inline long double real (const lcomplex& z) { return z.real(); }
// inline long double imag (const lcomplex& z) { return z.imag(); }
#endif

#elif defined USE_C99_COMPLEX
// C99 complex
#include <complex.h>

class icomplex {
 public:
  icomplex() : m_r(0), m_i(0) {}
  icomplex(const int r, const int i = 0) : m_r(r), m_i(i) {}
  icomplex(const icomplex& rhs) : m_r(rhs.m_r), m_i(rhs.m_i) {}
  ~icomplex() {}

  icomplex& operator=(const icomplex& rhs)
  {
    m_r = rhs.m_r;
    m_i = rhs.m_i;
    return *this;
  }

  icomplex& operator*=(const icomplex& rhs)
  {
    double re = m_r;
    double im = m_i;

    m_r = re * rhs.m_r - im * rhs.m_i;
    m_i = re * rhs.m_i + im * rhs.m_r;
    return *this;
  }

  friend icomplex operator*(const icomplex& lhs, const icomplex& rhs)
  {
    return icomplex(lhs.m_r * rhs.m_r - lhs.m_i * rhs.m_i,
                    lhs.m_r * rhs.m_i + lhs.m_i * rhs.m_r);
  }

  friend int real(const icomplex& z) { return z.m_r; }
  friend int imag(const icomplex& z) { return z.m_i; }

 private:
  int m_r;
  int m_i;
};

inline icomplex cmplx(const int real, const int imag)
{
  return icomplex(real, imag);
}


typedef float _Complex    fcomplex;
typedef double _Complex   dcomplex;
//typedef long double _Complex lcomplex;

inline float real(const fcomplex& z) { return crealf(z); }
inline float imag(const fcomplex& z) { return cimagf(z); }
inline float abs(const fcomplex& z) { return cabsf(z); }
inline fcomplex cmplx(const float real, const float imag)
{
  return real + imag * _Complex_I;
}


inline double real(const dcomplex& z) { return creal(z); }
inline double imag(const dcomplex& z) { return cimag(z); }
inline double abs(const dcomplex& z) { return cabs(z); }
inline dcomplex cmplx(const double real, const double imag)
{
  return real + imag * _Complex_I;
}


/*
inline long double real(const lcomplex& z) { return creall(z); }
inline long double imag(const lcomplex& z) { return cimagl(z); }
inline long double abs(const lcomplex& z) { return cabsl(z); }
inline lcomplex cmplx(const long double real, const long double imag) {
  return real + imag * _Complex_I;
}
*/

/*
dcomplex operator-(const dcomplex& rhs)
  {
    return cmplx(-real(rhs),-imag(rhs));
  }
*/
#endif
#endif // BRIDGE_COMPLEX_INCLUDED
