/*!
        @file    smear_HYP_SF.h

        @brief

        @author  Yusuke Tanigchi  (tanigchi)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef SMEAR_HYP_SF_INCLUDED
#define SMEAR_HYP_SF_INCLUDED

#include "smear.h"

#include "Field/field_G_SF.h"
#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HYP smearing of link variables with SF BC

/*!
                                    [26 May 2012 Y.Taniguchi]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
 */

class BAPI Smear_HYP_SF : public Smear
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim, m_Nvol;
  double m_alpha1, m_alpha2, m_alpha3;                 // HYP smearing parameters
  Projection *m_proj;
  std::vector<Field_G> m_U;
  std::vector<Field_G> m_v1;
  std::vector<Field_G> m_v2;
  ShiftField_lex m_shift;

  Field_G_SF set_wk;

 public:
  Smear_HYP_SF(Projection *proj)
    : Smear(), m_proj(proj)
  {
    init();
  }

  Smear_HYP_SF(unique_ptr<Projection>& proj)
    : Smear(), m_proj(proj.get())
  {
    init();
  }

  ~Smear_HYP_SF() {}

  void init();

  void set_parameters(const Parameters& params);
  void set_parameters(const double alpha1, const double alpha2, const double alpha3,
                      double *phi, double *phipr);

  void smear(Field_G& Usmear, const Field_G& U);

 private:
  void staple(Field_G&, const Field_G&, const Field_G&,
              const int mu, const int nu);

  void step1();
  void step2();
  void step3(Field_G&);

  int index_v1(const int mu, const int nu, const int rho)
  {
    int sig = 6 - mu - nu - rho;

    if (sig > mu) --sig;
    return mu + m_Ndim * sig;
  }

  int index_v2(const int mu, int nu)
  {
    if (nu > mu) --nu;
    return mu + m_Ndim * nu;
  }

  int size_v1()
  {
    return m_Ndim * (m_Ndim - 1);
  }

  int size_v2()
  {
    return m_Ndim * (m_Ndim - 1);
  }

#ifdef USE_FACTORY
 private:
  static Smear *create_object(Projection *proj)
  {
    return new Smear_HYP_SF(proj);
  }

 public:
  static bool register_factory()
  {
    return Smear::Factory::Register("HYP_SF", create_object);
  }
#endif
};
#endif
