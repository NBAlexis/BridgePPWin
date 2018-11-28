/*!
        @file    $Id:: forceSmear_APE.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FORCESMEAR_APE_INCLUDED
#define FORCESMEAR_APE_INCLUDED

#include "forceSmear.h"
#include "Smear/smear_APE.h"

#include "bridge_complex.h"
#include "Field/shiftField_lex.h"

#include "Smear/projection.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Recursive calculation for APE smeared fermion force.

/*!
                                [08 Apr 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
 */


class ForceSmear_APE : public ForceSmear
{
 public:
  static const std::string class_name;

 private:
  int                  m_Ndim, m_Nvol;
  std::vector<double>  m_rho;
  Projection           *m_proj;
  ShiftField_lex       m_shift;
  std::vector<Field_G> m_U;
  std::vector<Field_G> m_iTheta;

 public:
  ForceSmear_APE(Projection *proj)
    : ForceSmear(), m_proj(proj)
  {
    init();
  }

  ForceSmear_APE(unique_ptr<Projection>& proj)
    : ForceSmear(), m_proj(proj.get())
  {
    init();
  }

  //  ~ForceSmear_APE(){
  //  };

  // Setting parameters with Parameters object.
  void set_parameters(const Parameters& params);

  // Setting parameters with uniform smearing parameter.
  void set_parameters(const double rho1);

  // Setting parameters with anisotropic smearing parameter.
  void set_parameters(const std::vector<double>& rho);

  // Force computation.
  void force_udiv(Field_G& Sigma, const Field_G& Sigma_p, const Field_G& U);

 private:

  void init();

  double rho(int mu, int nu)
  {
    return m_rho[mu + nu * m_Ndim];
  }

  void force_each(Field_G&, const Field_G&, const Field_G&,
                  const Field_G&, const Field_G&, int mu, int nu);

  void staple(Field_G&, const Field_G&, const Field_G&,
              int mu, int nu);
};
#endif
