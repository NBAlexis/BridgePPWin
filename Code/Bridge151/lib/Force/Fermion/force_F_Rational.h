/*!
        @file    force_F_Rational.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FORCE_F_RATIONAL_INCLUDED
#define FORCE_F_RATIONAL_INCLUDED

#include "force_F.h"

#include "Field/field_F.h"
#include "Fopr/fopr_Rational.h"
#include "Tools/math_Rational.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for smeared fermion operators.

/*!
    This class BAPI determines the force of a rational approximation
    for a given fermion operator.
                                    [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
 */


class BAPI Force_F_Rational : public Force
{
 public:
  static const std::string class_name;

 private:
  int m_Np;                // number of poles in rational approx.
  int m_n_exp, m_d_exp;    // numerator and denominator of the exponent
  double m_x_min, m_x_max; // valid range of approximate sign function
  int m_Niter;             // max iteration of shiftsolver
  double m_Stop_cond;      // stopping condition of shift solver

  Field_G *m_U;
  Fopr *m_fopr;
  Force *m_force;

  // rational approx. coefficients
  double m_a0;
  std::vector<double> m_bl;
  std::vector<double> m_cl;

 public:
  Force_F_Rational(Fopr *fopr, Force *force)
    : Force(), m_fopr(fopr), m_force(force) {}

  Force_F_Rational(unique_ptr<Fopr>& fopr, unique_ptr<Force>& force)
    : Force(), m_fopr(fopr.get()), m_force(force.get()) {}

  ~Force_F_Rational() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Np, const int n_exp, const int d_exp, const double x_min, const double x_max,
                      const int Niter, const double Stop_cond);

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr->set_config(U);
    m_force->set_config(U);
  }

  void force_udiv(Field&, const Field&);

  void force_core1(Field&, const Field&, const Field&);  // dummy entry
  void force_udiv1(Field&, const Field&, const Field&);  // dummy entry

 private:
  void force_udiv_impl(Field_G&, const Field_F&);
};
#endif
