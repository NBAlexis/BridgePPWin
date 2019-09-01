/*!
        @file    smear_APE_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef SMEAR_APE_SF_INCLUDED
#define SMEAR_APE_SF_INCLUDED

#include "smear.h"

#include "Measurements/Gauge/staple_SF.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! APE element_type smearing of link variables.

/*!
    This class BAPI is alternative to the Smear_APE class BAPI.
                            [08 Apr 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.    [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                            [21 Mar 2015 Y.Namekawa]
 */


class BAPI Smear_APE_SF : public Smear
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim;
  std::valarray<double> m_rho;
  Projection *m_proj;

  //! SF boundary condition at t=0
  double m_phi[3];
  //! SF boundary condition at t=Nt
  double m_phipr[3];

 public:
  Smear_APE_SF(Projection *proj)
    : Smear(),
    m_Ndim(CommonParameters::Ndim()), m_rho(0.0, m_Ndim * m_Ndim),
    m_proj(proj) {}

  Smear_APE_SF(unique_ptr<Projection>& proj)
    : Smear(),
    m_Ndim(CommonParameters::Ndim()), m_rho(0.0, m_Ndim * m_Ndim),
    m_proj(proj.get()) {}

  ~Smear_APE_SF() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const double rho1, double *phi, double *phipr);
  void set_parameters(const std::vector<double>& rho, double *phi, double *phipr);

  void smear(Field_G& Usmear, const Field_G& U);

#ifdef USE_FACTORY
 private:
  static Smear *create_object(Projection *proj)
  {
    return new Smear_APE_SF(proj);
  }

 public:
  static bool register_factory()
  {
    return Smear::Factory::Register("APE_SF", create_object);
  }
#endif
};
#endif
