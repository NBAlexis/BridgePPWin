/*!
        @file    fopr_Wilson_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/


#ifndef FOPR_WILSON_SF_INCLUDED
#define FOPR_WILSON_SF_INCLUDED

#include "fopr_Wilson.h"

#include "Field/field_F_SF.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wilson fermion operator with SF BC.

/*!
    This class BAPI implements the Wilson fermion operator with the SF BC.
    <ul>
    <li>The SF BC for the fermion field is:
    <ul>
    <li>The fermion field at t=0 is always set to zero (Dirichlet BC).
    <li>The field at t=1,...,Lt-1 is active.
    </ul>
    <li>Implemented by delegation of the Fopr_Wilson class BAPI like Fopr_Clover.
    <li>The modification is only in Fopr_Wilson_SF::D to set the boundary fermion field to zero before and after multiplication of the standard Wilson Dirac operator.
    <ul>
    <li>By this manipulation the SF BC is introduced in the fermion field.
    </ul>
    <li>A private function set_boundary_zero(Field&) is introduced for this bounadry manipulation.
    <li>A few private members are added for this function.
    <li> [04 Apr 2012 Y.Taniguchi]
    </ul>
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
 */

class BAPI Fopr_Wilson_SF : public Fopr
{
 public:
  static const std::string class_name;

 private:
  int m_Nvol, m_Ndim, m_Nc, m_Nd, m_NinF;
  double m_kappa;
  std::vector<int> m_boundary;
  std::string m_mode;

  Fopr_Wilson *m_fopr_w;
  const Field_G *m_U;

  /*
  //! Needed to know a node at the temporal boundary.
  Communicator* comm;
  //! A spatial volume in a node.
  int Svol;
  //! num of the double color elements
  int m_Nc2;
  */

  //! In order to set the boundary field to zero.
  Field_F_SF set_zero;

 public:
  Fopr_Wilson_SF()
  {
    m_Nvol = CommonParameters::Nvol();
    m_Ndim = CommonParameters::Ndim();
    m_Nc   = CommonParameters::Nc();
    m_Nd   = CommonParameters::Nd();
    m_NinF = 2 * m_Nc * m_Nd;

    m_boundary.resize(m_Ndim);
    m_U      = 0;
    m_fopr_w = new Fopr_Wilson;

    /*
    comm = Communicator::init();
    m_Nc2 = 2*m_Nc;
    int Nt = CommonParameters::Nt();
    Svol=m_Nvol/Nt;
    */
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const double kappa, const std::vector<int> bc);

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_w->set_config(U);
  }

  void set_config(unique_ptr<Field_G>& U)
  {
    m_U = U.get();
    m_fopr_w->set_config(U.get());
  }

  ~Fopr_Wilson_SF()
  {
    delete m_fopr_w;
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
      vout.crucial(m_vl, "Error at %s: mode undefined.\n", class_name.c_str());
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
      vout.crucial(m_vl, "Error at %s: mode undefined.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }

  void set_mode(const std::string mode)
  {
    m_mode = mode;
  }

  std::string get_mode() const
  {
    return m_mode;
  }

  void DdagD(Field&, const Field&);
  void D(Field&, const Field&);
  void Ddag(Field&, const Field&);
  void H(Field&, const Field&);

  void mult_gm5(Field& v, const Field& w)
  {
    m_fopr_w->mult_gm5(v, w);
  }

  void mult_gm5p(const int mu, Field_F& v, const Field_F& w);

  int field_nvol() { return m_Nvol; }
  int field_nin() { return 2 * m_Nc * m_Nd; }
  int field_nex() { return 1; }

  //! this returns the number of floating point operations.
  double flop_count();

 private:
  // A function to set the fermion field to zero at the t=0 boundary.
  //  void set_boundary_zero(Field&);

#ifdef USE_FACTORY
 private:
  static Fopr *create_object()
  {
    return new Fopr_Wilson_SF();
  }

 public:
  static bool register_factory()
  {
    return Fopr::Factory_noarg::Register("Wilson_SF", create_object);
  }
#endif
};
#endif
