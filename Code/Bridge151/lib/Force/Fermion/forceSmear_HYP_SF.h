/*!
        @file    forceSmear_HYP_SF.h

        @brief

        @author  Yusuke Taniguchi  (taniguchi)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FORCESMEAR_HYP_SF_INCLUDED
#define FORCESMEAR_HYP_SF_INCLUDED

#include "forceSmear.h"

#include "Field/field_G_SF.h"
#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Recursive calculation of HYP smeared fermion force with SF BC.

/*!
                                [26 May 2012 Y.Taniguchi]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
 */


class BAPI ForceSmear_HYP_SF : public ForceSmear
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim, m_Nvol;
  double m_alpha1, m_alpha2, m_alpha3;                 // HYP smearing parameters
  Projection *m_proj;
  std::vector<Field_G> m_U;
  std::vector<Field_G> m_v1, m_v2;
  std::vector<Field_G> m_Sigma3, m_Sigma2;
  std::vector<Field_G> m_iTheta3, m_iTheta2, m_iTheta1;
  ShiftField_lex m_shift;

  //! SF boundary condition at t=0
  double m_phi[3];
  //! SF boundary condition at t=Nt
  double m_phipr[3];

  Field_G_SF set_wk;

 public:
  ForceSmear_HYP_SF(Projection *proj)
    : ForceSmear(), m_proj(proj)
  {
    init();
  }

  ForceSmear_HYP_SF(unique_ptr<Projection>& proj)
    : ForceSmear(), m_proj(proj.get())
  {
    init();
  }

  //  ~ForceSmear_HYP_SF(){ };

  void set_parameters(const Parameters& params);
  void set_parameters(const double alpha1, const double alpha2, const double alpha3,
                      double *phi, double *phipr);

  void force_udiv(Field_G& Sigma, const Field_G& Sigma_p, const Field_G& U);

  // old interface
  //Field force_udiv(const Field_G& Sigma, const Field_G& U);

 private:
  void init();

  void force_step3(Field_G&, const Field_G&);
  void force_step2(Field_G&);
  void force_step1(Field_G&);

  void force_each(Field_G&, const Field_G&, const Field_G&,
                  const Field_G&, const Field_G&, const int mu, const int nu);

  void smear_step1();
  void smear_step2();

  void staple(Field_G&, const Field_G&, const Field_G&,
              const int mu, const int nu);

  int idx1(const int mu, const int nu, const int rho)
  {
    int sig = 6 - mu - nu - rho;

    if (sig > mu) --sig;
    return mu + m_Ndim * sig;
  }

  int idx1b(const int mu, int nu, int rho)
  {
    if (nu > mu) --nu;
    if (rho > mu) --rho;
    if (rho > nu) --rho;
    return mu + m_Ndim * (nu + (m_Ndim - 1) * rho);
  }

  int idx2(const int mu, int nu)
  {
    if (nu > mu) --nu;
    return mu + m_Ndim * nu;
  }

  int size1()
  {
    return m_Ndim * (m_Ndim - 1);
  }

  int size1b()
  {
    return m_Ndim * (m_Ndim - 1) * (m_Ndim - 2);
  }

  int size2()
  {
    return m_Ndim * (m_Ndim - 1);
  }

#ifdef USE_FACTORY
 private:
  static ForceSmear *create_object(Projection *proj)
  {
    return new ForceSmear_HYP_SF(proj);
  }

 public:
  static bool register_factory()
  {
    return ForceSmear::Factory::Register("HYP_SF", create_object);
  }
#endif
};
#endif
