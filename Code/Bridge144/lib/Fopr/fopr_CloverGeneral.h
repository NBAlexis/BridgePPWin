/*!
        @file    $Id:: fopr_CloverGeneral.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FOPR_CLOVERGENERAL_INCLUDED
#define FOPR_CLOVERGENERAL_INCLUDED

#include "fopr_WilsonGeneral.h"
#include "fopr_CloverTerm_General.h"


#include "IO/bridgeIO.h"
using Bridge::vout;

//! Clover general fermion operator.

/*!
    This class implements Clover General fermion operator,
    including Clover fermion on anisotropic lattice,
      T. R. Klassen, Nucl. Phys. B509 391 (1998)
      QCD-TARO Collab., Nucl. Phys. Proc. Suppl. 83 411 (2000)
    and relativistic heavy quarks,
      A. X. El-Khadra et al., Phys. Rev. D 55, 3933 (1997)
      S. Aoki et al., Prog. Theor. Phys. 109, 383 (2003)
      N. H. Christ et al., Phys. Rev. D 76, 074505 (2007)
                                     [21 Mar 2015 Y.Namekawa]
 */


class Fopr_CloverGeneral : public Fopr
{
 public:
  static const std::string class_name;

 private:
  double           m_kappa_s;             //!< spatial hopping parameter
  double           m_kappa_t;             //!< temporal hopping parameter
  double           m_nu_s;                //!< spatial anisotopy parameter
  double           m_r_s;                 //!< spatial Wilson parameter
  double           m_cSW_s;               //!< spatial clover coefficient
  double           m_cSW_t;               //!< temporal clover coefficient
  std::vector<int> m_boundary;            //!< boundary conditions
  std::string      m_repr;                //!<  gamma matrix representation
  std::string      m_mode;                //!<  mode of multiplication

  int m_Nvol, m_Ndim, m_Nc, m_Nd, m_NinF; //!< internal parameters

  Fopr_WilsonGeneral      *m_fopr_w;      //!< Wilson fermion kernel
  Fopr_CloverTerm_General *m_fopr_csw;    //!< Clover term operator
  const Field_G           *m_U;           //!< gauge configuration (pointer)

  //  Fopr_Clover_imp* m_imp;    //!< pimple prescription
  Field m_v1, m_v2;              //!< working field.

 public:
  Fopr_CloverGeneral()
    : Fopr()
  {
    init("Dirac");
  }

  Fopr_CloverGeneral(std::string repr)
    : Fopr()
  {
    init(repr);
  }

  ~Fopr_CloverGeneral()
  {
    tidyup();
  }

  void set_parameters(const Parameters& params);
  void set_parameters(double kappa_s, double kappa_t,
                      double nu_s, double r_s,
                      double cSW_s, double cSW_t,
                      std::vector<int> bc);

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

  void set_mode(std::string mode)
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

  void mult_up(int mu, Field& v, const Field& w)
  {
    m_fopr_w->mult_up(mu, v, w);
  }

  void mult_dn(int mu, Field& v, const Field& w)
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

  void init(std::string repr);
  void tidyup();
};
#endif
