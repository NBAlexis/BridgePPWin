/*!
        @file    corr2pt_Wilson_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef CORR2PT_WILSON_SF_INCLUDED
#define CORR2PT_WILSON_SF_INCLUDED

#include "Field/field_F.h"
#include "Field/index_lex.h"

#include "Tools/gammaMatrixSet.h"

#include "bridge_complex.h"
using Bridge::vout;

//! Two-point correlator for Wilson-element_type fermions with SF BC.

/*!
  <ul>
  <li>Evaluate the axial current - pseudo scalar two point funciton fA, fA'.
  <li>Evaluate the pseudo scalar - pseudo scalar two point funciton fP, fP'
  <li>These two point functions are used to evaluate the PCAC mass and \f$\Delta M\f$ needed for the cSW determination.
  <li>[19 Apr 2012 Y.Taniguchi]
  </ul>
 */

class BAPI Corr2pt_Wilson_SF
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  Index_lex m_index;
  GammaMatrixSet *m_gmset;

 public:
  Corr2pt_Wilson_SF(GammaMatrixSet *gmset)
    : m_vl(CommonParameters::Vlevel()), m_gmset(gmset) {}

  Corr2pt_Wilson_SF(unique_ptr<GammaMatrixSet>& gmset)
    : m_vl(CommonParameters::Vlevel()), m_gmset(gmset.get()) {}

 private:
  // non-copyable
  Corr2pt_Wilson_SF(const Corr2pt_Wilson_SF&);
  Corr2pt_Wilson_SF& operator=(const Corr2pt_Wilson_SF&);

 public:
  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

/*!
  The axial current:
  \f[
  A_\mu^a(x)={\overline{\psi}}(x)\gamma_\mu\gamma_5\frac{\tau^a}{2}\psi(x)
  \f]
  The pseudo scalar density:
  \f[
  P^a(x)={\overline{\psi}}(x)\gamma_5\frac{\tau^a}{2}\psi(x)
  \f]
  The boundary pseudo scalar density at t=0
  \f[
  {\cal O}^a=a^6\sum_{\vec{x},\vec{y}}\overline{\zeta}\left(\vec{x}\right)
  \gamma_5\frac{\tau^a}{2}\zeta\left(\vec{y}\right)
  \f]
  The boundary pseudo scalar density at t=T
  \f[
  {\cal O'}^a=a^6\sum_{\vec{x},\vec{y}}\overline{\zeta}'\left(\vec{x}\right)
  \gamma_5\frac{\tau^a}{2}\zeta'\left(\vec{y}\right)
  \f]
  The two point functions are
  \f[
  f_A(x_0)=-\frac{1}{N_f^2-1}\left\langle{A_0^a(x_0){\cal O}^a}\right\rangle
  =\frac{1}{2}\sum_{\vec{v},\vec{y}}{\rm tr}\left(
  \left\langle{\zeta(\vec{v}){\overline{\psi}}(x)}\right\rangle\gamma_0\gamma_5
  \left\langle{\psi(x)\overline{\zeta}(\vec{y})}\right\rangle\gamma_5\right)
  \f]
  \f[
  f_P(x_0)=-\frac{1}{N_f^2-1}\left\langle{P^a(x_0){\cal O}^a}\right\rangle
  =\frac{1}{2}\sum_{\vec{v},\vec{y}}{\rm tr}\left(
  \left\langle{\zeta(\vec{v}){\overline{\psi}}(x)}\right\rangle\gamma_5
  \left\langle{\psi(x)\overline{\zeta}(\vec{y})}\right\rangle\gamma_5\right)
  \f]
  \f[
  f_A'(x_0)=+\frac{1}{N_f^2-1}\left\langle{A_0^a(T-x_0){\cal O'}^a}\right\rangle
  =-\frac{1}{2}\sum_{\vec{v},\vec{y}}{\rm tr}\left(
  \left\langle{\zeta'(\vec{v}){\overline{\psi}}(T-x_0)}\right\rangle\gamma_0\gamma_5
  \left\langle{\psi(T-x_0)\overline{\zeta}'(\vec{y})}\right\rangle\gamma_5\right)
  \f]
  \f[
  f_P'(x_0)=-\frac{1}{N_f^2-1}\left\langle{P^a(T-x_0){\cal O'}^a}\right\rangle
  =\frac{1}{2}\sum_{\vec{v},\vec{y}}{\rm tr}\left(
  \left\langle{\zeta'(\vec{v}){\overline{\psi}}(T-x_0)}\right\rangle\gamma_5
  \left\langle{\psi(T-x_0)\overline{\zeta}'(\vec{y})}\right\rangle\gamma_5\right)
  \f]
*/

  double fAfP(
    const std::vector<Field_F>& sq1,
    const std::vector<Field_F>& sq2);

  double meson_corr(
    std::vector<dcomplex>& meson,
    const GammaMatrix& gm_sink,
    const GammaMatrix& gm_src,
    const std::vector<Field_F>& sq1,
    const std::vector<Field_F>& sq2);
};
#endif
