/*!
        @file    forceSmear_APE_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FORCESMEAR_APE_ALT_INCLUDED
#define FORCESMEAR_APE_ALT_INCLUDED

#include "forceSmear.h"

#include "Field/field_G_SF.h"
#include "Field/shiftField_lex.h"

#include "Smear/smear_APE_SF.h"
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


class BAPI ForceSmear_APE_SF : public ForceSmear
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim, m_Nvol;
  std::vector<double> m_rho;
  Projection *m_proj;
  ShiftField_lex m_shift;
  std::vector<Field_G> m_U;
  std::vector<Field_G> m_iTheta;

  //! SF boundary condition at t=0
  double m_phi[3];
  //! SF boundary condition at t=Nt
  double m_phipr[3];

  Field_G_SF set_wk;

 public:
  ForceSmear_APE_SF(Projection *proj)
    : ForceSmear(), m_proj(proj)
  {
    init();
  }

  ForceSmear_APE_SF(unique_ptr<Projection>& proj)
    : ForceSmear(), m_proj(proj.get())
  {
    init();
  }

  //  ~ForceSmear_APE_SF(){
  //  };

  void set_parameters(const Parameters& params);

  void set_parameters(const double rho1, double *phi, double *phipr);
  void set_parameters(const std::vector<double>& rho, double *phi, double *phipr);

  void force_udiv(Field_G& Sigma, const Field_G& Sigma_p, const Field_G& U);

  // old interface
  //Field force_udiv(const Field_G& Sigma, const Field_G& U);

 private:
  void init();

  double rho(const int mu, const int nu)
  {
    return m_rho[mu + nu * m_Ndim];
  }

  void force_each(Field_G&,
                  const Field_G&, const Field_G&,
                  const Field_G&, const Field_G&, const int mu, const int nu);

  void staple(Field_G&,
              const Field_G&, const Field_G&,
              const int mu, const int nu);

#ifdef USE_FACTORY
 private:
  static ForceSmear *create_object(Projection *proj)
  {
    return new ForceSmear_APE_SF(proj);
  }

 public:
  static bool register_factory()
  {
    return ForceSmear::Factory::Register("APE_SF", create_object);
  }
#endif
};
#endif
