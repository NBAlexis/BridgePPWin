/*!
        @file    fopr_Wilson_Isochemical.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FOPR_WILSON_ISOCHEMICAL_INCLUDED
#define FOPR_WILSON_ISOCHEMICAL_INCLUDED

#include "fopr_Wilson.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wilson fermion operator with isospin chemical potential.

/*!
                                    [22 Aug 2012 H,Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    Multi-threaded.                 [12 Sep 2014 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
 */

class BAPI Fopr_Wilson_Isochemical : public Fopr
{
 public:
  static const std::string class_name;

 private:
  int m_Nvol, m_Ndim;
  double m_kappa;                //!< hopping parameter
  double m_mu;                   //!< isospin chemical potential
  std::vector<int> m_boundary;
  std::string m_mode;
  Fopr_Wilson *m_fopr_w;

  double m_exp_mu;               //!< exp(mu)

  std::string m_repr;
  void (Fopr_Wilson_Isochemical::*m_mult)(Field&, const Field&);
  void (Fopr_Wilson_Isochemical::*m_mult_dag)(Field&, const Field&);
  void (Fopr_Wilson_Isochemical::*m_D)(Field&, const Field&);
  void (Fopr_Wilson_Isochemical::*m_gm5)(Field&, const Field&);

  const Field_G *m_U;

  Field m_vt;  //!< working field.
  Field m_w2;  //!< working field.

 public:

  Fopr_Wilson_Isochemical() : Fopr()
  {
    init("Dirac");
  }

  Fopr_Wilson_Isochemical(const std::string repr) : Fopr()
  {
    init(repr);
  }

  ~Fopr_Wilson_Isochemical()
  {
    tidyup();
  }

  // this method is temporary. After all calls specify bc, to be removed.
  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa, const double mu,
                      const std::vector<int> bc);

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

  void mult(Field& v, const Field& f)
  {
    (this->*m_mult)(v, f);
  }

  void mult_dag(Field& v, const Field& f)
  {
    (this->*m_mult_dag)(v, f);
  }

  void set_mode(const std::string mode)
  {
    m_mode = mode;

    if (m_mode == "D") {
      m_mult     = &Fopr_Wilson_Isochemical::D;
      m_mult_dag = &Fopr_Wilson_Isochemical::Ddag;
    } else if (m_mode == "Ddag") {
      m_mult     = &Fopr_Wilson_Isochemical::Ddag;
      m_mult_dag = &Fopr_Wilson_Isochemical::D;
    } else if (m_mode == "DdagD") {
      m_mult     = &Fopr_Wilson_Isochemical::DdagD;
      m_mult_dag = &Fopr_Wilson_Isochemical::DdagD;
    } else if (m_mode == "H") {
      m_mult     = &Fopr_Wilson_Isochemical::H;
      m_mult_dag = &Fopr_Wilson_Isochemical::Hdag;
    } else {
      vout.crucial(m_vl, "Error at %s: input mode is undefined.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }

  std::string get_mode() const
  {
    return m_mode;
  }

  void mult_gm5p(const int mu, Field_F& v, const Field_F& w);

  void mult_gm5(Field&, const Field&);
  void D(Field&, const Field&);
  void Dminmu(Field&, const Field&);
  void Dspc(Field&, const Field&);

  void DdagD(Field&, const Field&);
  void Ddag(Field&, const Field&);
  void H(Field&, const Field&);
  void Hdag(Field&, const Field&);

  void mult_undef(Field&, const Field&)
  {
    vout.crucial(m_vl, "Error at %s: mode undefined.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  void mult_up(const int mu, Field& v, const Field& w)
  {
    m_fopr_w->mult_up(mu, v, w);
  }

  void mult_dn(const int mu, Field& v, const Field& w)
  {
    m_fopr_w->mult_dn(mu, v, w);
  }

  int field_nvol() { return CommonParameters::Nvol(); }
  int field_nin() { return 2 * CommonParameters::Nc() * CommonParameters::Nd(); }
  int field_nex() { return 1; }

  //! this returns the number of floating point operations.
  double flop_count();

 private:
  void init(const std::string repr);
  void tidyup();

#ifdef USE_FACTORY
 private:
  static Fopr *create_object()
  {
    return new Fopr_Wilson_Isochemical();
  }

  static Fopr *create_object_with_repr(const std::string& repr)
  {
    return new Fopr_Wilson_Isochemical(repr);
  }

 public:
  static bool register_factory()
  {
    bool init1 = Fopr::Factory_noarg::Register("Wilson_Isochemical", create_object);
    bool init2 = Fopr::Factory_string::Register("Wilson_Isochemical", create_object_with_repr);

    return init1 && init2;
  }
#endif
};
#endif
