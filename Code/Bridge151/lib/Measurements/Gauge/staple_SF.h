/*!
        @file    staple_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef STAPLE_SF_INCLUDED
#define STAPLE_SF_INCLUDED

#include "Field/field_G_SF.h"
#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Staple construction.

/*!
  \brief Evaluate staple with SF BC.
  <ul>
  <li>The BC parameters phi, phipr are set into SU(3) matrix wk, wkpr for boundary spatial link variable.
  <li> [26 Jan. 2012 Y.Taniguchi]
  </ul>
    (Coding history will be recovered from trac.)
    YAML is implemented.         [14 Nov 2012 Y.Namekawa]
 */


class BAPI Staple_SF
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  Index_lex m_index;
  ShiftField_lex m_shift;

  Mat_SU_N m_wk, m_wkpr, m_i_omega0;
  int m_initialized;

 public:
  Staple_SF()
    : m_vl(CommonParameters::Vlevel()),
    m_wk(CommonParameters::Nc()),
    m_wkpr(CommonParameters::Nc()),
    m_i_omega0(CommonParameters::Nc()),
    m_initialized(0) {}

 private:
  // non-copyable
  Staple_SF(const Staple_SF&);
  Staple_SF& operator=(const Staple_SF&);

 public:
  void set_parameters(const Parameters& params);

  void set_parameters(const double *phi, const double *phipr);
  void set_parameters(const double *phi, const double *phipr, const double *pomega);
  void set_parameters(const std::vector<double>& phi, const std::vector<double>& phipr, const std::vector<double>& p_omega);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  void upper(Field_G_SF&, const Field_G&, const int, const int);
  void lower(Field_G_SF&, const Field_G&, const int, const int);
  double plaq_s(const Field_G&);
  double plaq_t(const Field_G&);
  double plaq_t_ct(const Field_G&, const double ct);
  double plaquette(const Field_G&);
  double plaquette_ct(const Field_G&, const double ct);

  double sf_coupling_plaq(const Field_G&, const double ct);
  double sf_coupling_rect(const Field_G&, const double ctr);

  void staple(Field_G&, const Field_G&, const int);
  void staple_ct(Field_G&, const Field_G&, const int, const double ct);

  void print_plaquette(const Field_G&);
};
#endif
