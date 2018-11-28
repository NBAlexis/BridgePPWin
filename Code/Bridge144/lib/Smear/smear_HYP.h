/*!
        @file    $Id:: smear_HYP.h #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/
#ifndef SMEAR_HYP_INCLUDED
#define SMEAR_HYP_INCLUDED

#include "smear.h"
#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HYP smearing of link variables.

/*!
                                    [18 Mar 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
 */


class Smear_HYP : public Smear
{
 public:
  static const std::string class_name;

 private:
  int                  m_Ndim, m_Nvol;
  double               m_alpha1, m_alpha2, m_alpha3;   // HYP smearing parameters
  Projection           *m_proj;
  std::vector<Field_G> m_U;
  std::vector<Field_G> m_v1;
  std::vector<Field_G> m_v2;
  ShiftField_lex       m_shift;

 public:
  Smear_HYP(Projection *proj)
    : Smear(), m_proj(proj)
  {
    init();
  }

  Smear_HYP(unique_ptr<Projection>& proj)
    : Smear(), m_proj(proj.get())
  {
    init();
  }

  ~Smear_HYP() {}

  void init();

  void set_parameters(const Parameters& params);
  void set_parameters(double alpha1, double alpha2, double alpha3);

  void smear(Field_G& Usmear, const Field_G& U);

 private:
  void staple(Field_G&, const Field_G&, const Field_G&,
              int mu, int nu);

  void step1();
  void step2();
  void step3(Field_G&);

  int index_v1(int mu, int nu, int rho)
  {
    int sig = 6 - mu - nu - rho;

    if (sig > mu) --sig;
    return mu + m_Ndim * sig;
  }

  int index_v2(int mu, int nu)
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
};
#endif