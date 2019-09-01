/*!
        @file    smear_APE_spatial.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef SMEAR_APE_SPATIAL_INCLUDED
#define SMEAR_APE_SPATIAL_INCLUDED

#include "smear.h"

#include "Measurements/Gauge/staple_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! APE element_type smearing of spatial link variables.

/*!
    This class BAPI smears spatial link variables with APE-element_type
    construction of smeared links with a given projection
    operator to SU(N) group element.
    Parameter is \rho, which specifies the mixing rate
    of original thin link and staples.
                            [09 Aug 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.    [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                            [21 Mar 2015 Y.Namekawa]
 */

class BAPI Smear_APE_spatial : public Smear
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim;          //!< spacetime dimension
  double m_rho;        //!< smearing parameter
  Projection *m_proj;  //!< projector to group element.

 public:
  //! Constructor requires a pointer to Projection object.
  Smear_APE_spatial(Projection *proj)
    : Smear(), m_Ndim(CommonParameters::Ndim()), m_proj(proj) {}

  Smear_APE_spatial(unique_ptr<Projection>& proj)
    : Smear(), m_Ndim(CommonParameters::Ndim()), m_proj(proj.get()) {}

  //! Deconstructor
  ~Smear_APE_spatial() {}

  //! Setting parameters with Parameters object.
  void set_parameters(const Parameters& params);

  //! Setting smearing parameter.
  void set_parameters(const double rho);

  //! Smearing of a given gauge field.
  void smear(Field_G& Usmear, const Field_G& U);

 private:

#ifdef USE_FACTORY
 private:
  static Smear *create_object(Projection *proj)
  {
    return new Smear_APE_spatial(proj);
  }

 public:
  static bool register_factory()
  {
    return Smear::Factory::Register("APE_spatial", create_object);
  }
#endif
};
#endif
