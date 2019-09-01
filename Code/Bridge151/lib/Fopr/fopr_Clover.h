/*!
        @file    fopr_Clover.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FOPR_CLOVER_INCLUDED
#define FOPR_CLOVER_INCLUDED

#include "fopr_Wilson.h"
#include "fopr_CloverTerm.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Clover fermion operator.

/*!
    This class BAPI implements the clover (improved Wilson) fermion
    operator.
    Wilson kernel and clover term are implemented in other classes,
    and this class BAPI holds them as objects.
    (The implementation was modified after revision 645: before
     that, clover term was being implemented inside this class BAPI.)
    The `mode' controls which of D, Ddag, H, DdagD are multiplied
    when mult or mult_dag is called.
       [first ver. 24 Dec 2011/ modified 28 Aug 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
 */


class BAPI Fopr_Clover : public Fopr
{
 public:
  static const std::string class_name;

 private:
  double m_kappa;                         //!< hopping parameter
  double m_cSW;                           //!< clover coefficient
  std::vector<int> m_boundary;            //!< boundary conditions
  std::string m_repr;                     //!<  gamma matrix representation
  std::string m_mode;                     //!<  mode of multiplication

  int m_Nvol, m_Ndim, m_Nc, m_Nd, m_NinF; //!< internal parameters

  Fopr_Wilson *m_fopr_w;                  //!< Wilson fermion kernel
  Fopr_CloverTerm *m_fopr_csw;            //!< Clover term operator
  const Field_G *m_U;                     //!< gauge configuration (pointer)

  //  Fopr_Clover_imp* m_imp;    //!< pimple prescription
  Field m_v1, m_v2;              //!< working field.

 public:
  Fopr_Clover()
    : Fopr()
  {
    init("Dirac");
  }

  Fopr_Clover(const std::string repr)
    : Fopr()
  {
    init(repr);
  }

  ~Fopr_Clover()
  {
    tidyup();
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const double kappa, const double cSW, const std::vector<int> bc);

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_w->set_config(U);
    m_fopr_csw->set_config(U);
  }

  void set_config(unique_ptr<Field_G>& U)
  {
    m_U = U.get();
    m_fopr_w->set_config(U.get());
    m_fopr_csw->set_config(U.get());
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
    } else if (m_mode == "Ddag") {
      Ddag(v, f);
    } else if (m_mode == "DdagD") {
      DdagD(v, f);
    } else if (m_mode == "DDdag") {
      DDdag(v, f);
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

  void D(Field&, const Field&);
  void Ddag(Field&, const Field&);
  void DdagD(Field&, const Field&);
  void DDdag(Field&, const Field&);
  void H(Field&, const Field&);

  void mult_gm5(Field& v, const Field& w)
  {
    m_fopr_w->mult_gm5(v, w);
  }

  void mult_up(const int mu, Field& v, const Field& w)
  {
    m_fopr_w->mult_up(mu, v, w);
  }

  void mult_dn(const int mu, Field& v, const Field& w)
  {
    m_fopr_w->mult_dn(mu, v, w);
  }

  void mult_isigma(Field_F&, const Field_F&,
                   const int mu, const int nu);

  int field_nvol() { return m_Nvol; }
  int field_nin() { return 2 * m_Nc * m_Nd; }
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
    return new Fopr_Clover();
  }

  static Fopr *create_object_with_arg(const std::string& repr)
  {
    return new Fopr_Clover(repr);
  }

 public:
  static bool register_factory()
  {
    bool init1 = Fopr::Factory_noarg::Register("Clover", create_object);
    bool init2 = Fopr::Factory_string::Register("Clover", create_object_with_arg);

    return init1 && init2;
  }
#endif
};
#endif
