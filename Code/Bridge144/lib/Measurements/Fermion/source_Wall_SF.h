/*!
        @file    $Id:: source_Wall_SF.h #$

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef SOURCE_WALL_SF_INCLUDED
#define SOURCE_WALL_SF_INCLUDED

#include "Field/field_F.h"
#include "Smear/director_Smear.h"
#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wall source for SF boundary propagator.

/*!
    This class sets a wall source vector for the 4-spinor
    (Wilson-type) fermion.
                                 [12 Apr 2012 Y.Taniguchi]
    (Coding history will be recovered from trac.)
    YAML is implemented.         [14 Nov 2012 Y.Namekawa]
 */


class Source_Wall_SF
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  //! Pointer to gauge field
  Field_G *m_U;
  //! Boundary O(a) improvement factor \f$\widetilde{c}_t\f$ for the Wilson (clover) fermion.
  double m_ct_tilde;

  Index_lex m_index;

 public:

  Source_Wall_SF()
    : m_vl(CommonParameters::Vlevel()) {}

 private:
  Source_Wall_SF(const Source_Wall_SF&);
  Source_Wall_SF& operator=(const Source_Wall_SF&);

 public:
  void set_parameters(const Parameters& params);

  void set_parameters(double ct_tilde);
  void set_parameters(Field_G *U, double ct_tilde);
  void set_parameters(Field_G *U, Director_Smear *dr_smear, double ct_tilde);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  //! setting pointer to configuration
  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
  }

  //! Set the 3D wall source at t=1.
  void set_t0(Field_F& src, int ic, int id);

  //! Set the 3D wall source at t=T-1.
  void set_tT(Field_F& src, int ic, int id);
};
#endif
