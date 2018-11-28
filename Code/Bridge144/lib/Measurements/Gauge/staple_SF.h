/*!
        @file    $Id:: staple_SF.h #$

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/


#ifndef STAPLE_SF_INCLUDED
#define STAPLE_SF_INCLUDED

#include "Parameters/parameters.h"
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


class Staple_SF
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  int            Nc;
  int            Ndim;
  int            Nvol, Lvol;
  Field_G_SF     Umu, Unu, v, w;
  Index_lex      index;
  ShiftField_lex shift;

  int      Nx, Ny, Nz, Nt;
  int      Lx, Ly, Lz, Lt;
  int      NPEt;
  Mat_SU_N wk, wkpr, iomega0;
  int      initialized;

 public:

  Staple_SF()
    : m_vl(CommonParameters::Vlevel()),
      Nc(CommonParameters::Nc()),
      Ndim(CommonParameters::Ndim()),
      Nvol(CommonParameters::Nvol()),
      Lvol(CommonParameters::Lvol()),
      Nx(CommonParameters::Nx()),
      Ny(CommonParameters::Ny()),
      Nz(CommonParameters::Nz()),
      Nt(CommonParameters::Nt()),
      Lx(CommonParameters::Lx()),
      Ly(CommonParameters::Ly()),
      Lz(CommonParameters::Lz()),
      Lt(CommonParameters::Lt()),
      NPEt(CommonParameters::NPEt()),
      wk(Nc), wkpr(Nc), iomega0(Nc),
      initialized(0) {}

 private:
  // non-copyable
  Staple_SF(const Staple_SF&);
  Staple_SF& operator=(const Staple_SF&);

 public:
  void set_parameters(const Parameters& params);

  void set_parameters(double *phi, double *phipr);
  void set_parameters(const double *phi, const double *phipr, const double *pomega);
  void set_parameters(std::vector<double>& phi, std::vector<double>& phipr, std::vector<double>& pomega);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  void upper(Field_G_SF&, const Field_G&, const int, const int);
  void lower(Field_G_SF&, const Field_G&, const int, const int);
  double plaq_s(const Field_G&);
  double plaq_t(const Field_G&);
  double plaq_t_ct(const Field_G&, double ct);
  double plaquette(const Field_G&);
  double plaquette_ct(const Field_G&, double ct);

  double sf_coupling_plaq(const Field_G&, double ct);
  double sf_coupling_rect(const Field_G&, double ctr);

  void staple(Field_G&, const Field_G&, const int);
  void staple_ct(Field_G&, const Field_G&, const int, double ct);

  void print_plaquette(const Field_G&);
};
#endif
