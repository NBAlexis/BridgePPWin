/*!
        @file    fopr_Clover_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FOPR_CLOVER_SF_INCLUDED
#define FOPR_CLOVER_SF_INCLUDED

#include "fopr_Wilson_SF.h"
#include "Measurements/Gauge/staple_SF.h"

#include "Field/shiftField_lex.h"
#include "Tools/gammaMatrixSet.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Clover fermion operator.

/*!
    This class BAPI implements the clover (improved Wilson) fermion
    operator with SF BC.
    <ul>
    <li>The field strength is calculate when the function set_config() is called.
    <li>Dirac representation only!
    <li>[10 Apr 2012 Y.Taniguchi]
    </ul>
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
 */


class BAPI Fopr_Clover_SF : public Fopr
{
 public:
  static const std::string class_name;

 private:
  int m_Nvol, m_Ndim, m_Nc, m_Nd, m_NinF;
  double m_kappa, m_cSW;
  std::vector<int> m_boundary;
  std::string m_repr;
  std::string m_mode;

  void (Fopr_Clover_SF::*m_csw)(Field_F&, const Field_F&);

  //  Index_lex m_idx;
  Fopr_Wilson_SF *m_fopr_w;
  const Field_G *m_U;
  ShiftField_lex m_shift;

  Field_G m_Bx, m_By, m_Bz, m_Ex, m_Ey, m_Ez;
  // Bx = -iF(1,2), By = -iF(2,1), -iBz = F(0,1)
  // Ex = -iF(4,0), Ey = -iF(4,1), Ez = -iF(4,2)

  std::vector<GammaMatrix> m_GM, m_SG;

  //! SF boundary condition at t=0
  double m_phi[3];
  //! SF boundary condition at t=Nt
  double m_phipr[3];

  //! In order to set the boundary field to zero.
  Field_F_SF setzero;

 public:
  Fopr_Clover_SF()
  {
    init("Dirac");
  }

  ~Fopr_Clover_SF()
  {
    tidyup();
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const double kappa, const double cSW, const std::vector<int> bc,
                      double *phi, double *phipr);

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_w->set_config(U);

    set_csw();
  }

  void set_config(unique_ptr<Field_G>& U)
  {
    m_U = U.get();
    m_fopr_w->set_config(U.get());

    set_csw();
  }

  void set_mode(const std::string mode)
  {
    m_mode = mode;
  }

  std::string get_mode() const
  {
    return m_mode;
  }

  void mult(Field& v, const Field& f)
  {
    if (m_mode == "D") {
      D(v, f);
    } else if (m_mode == "DdagD") {
      DdagD(v, f);
    } else if (m_mode == "Ddag") {
      Ddag(v, f);
    } else if (m_mode == "H") {
      H(v, f);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n", class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }

  void mult_dag(Field& v, const Field& f)
  {
    if (m_mode == "D") {
      Ddag(v, f);
    } else if (m_mode == "DdagD") {
      DdagD(v, f);
    } else if (m_mode == "Ddag") {
      D(v, f);
    } else if (m_mode == "H") {
      H(v, f);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n", class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }

  void DdagD(Field&, const Field&);
  void D(Field&, const Field&);
  void Ddag(Field&, const Field&);
  void H(Field&, const Field&);

  void mult_gm5(Field& v, const Field& w)
  {
    m_fopr_w->mult_gm5(v, w);
  }

  void mult_isigma(Field_F&, const Field_F&,
                   const int mu, const int nu);

  int field_nvol() { return m_Nvol; }
  int field_nin() { return 2 * m_Nc * m_Nd; }
  int field_nex() { return 1; }

  //! this returns the number of floating point number operations.
  double flop_count();

 private:
  void init(const std::string repr);
  void tidyup();

  void set_csw();
  void mult_csw(Field_F&, const Field_F&);
  void set_fieldstrength(Field_G&, const int, const int);

  void mult_csw_dirac(Field_F&, const Field_F&);
  void mult_csw_chiral(Field_F&, const Field_F&);

  int sg_index(const int mu, const int nu) { return mu * m_Ndim + nu; }

#ifdef USE_FACTORY
 private:
  static Fopr *create_object()
  {
    return new Fopr_Clover_SF();
  }

 public:
  static bool register_factory()
  {
    return Fopr::Factory_noarg::Register("Clover_SF", create_object);
  }
#endif
};
#endif
