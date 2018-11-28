/*!
        @file    $Id:: fopr_Clover_eo.h #$

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FOPR_CLOVER_EO_INCLUDED
#define FOPR_CLOVER_EO_INCLUDED

#include <vector>
#include <string>

#include "fopr_Wilson_eo.h"
//#include "fopr_Wilson_eo_impl.h"
#include "fopr_CloverTerm_eo.h"

//#include "Solver/solver_CG.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Even-odd Clover fermion operator.

/*!
    This class is an even-odd version of Clover fermion operator.
    At present this is rough implementation, while correctly
    works, and to be updated by supplying complete functionality.
    Only the functions needed for even-odd preconditioned solver
    is ready.
                                        [20 June 2012 S.UEDA]
    (Coding history will be recovered from trac.)
    Modify this code to work.           [03 Mar 2013 Y.Namekawa]
    Multi-threaded.                     [12 Jul 2014 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                        [21 Mar 2015 Y.Namekawa]
 */

class Fopr_Clover_eo : public Fopr_eo
{
 public:
  static const std::string class_name;

 private:
  int m_Nc, m_Nd, m_NinF, m_Ndim;
  int m_Nvol, m_Nvol2;

  double           m_kappa;    //!< hopping parameter.
  double           m_cSW;      //!< clover coefficient.
  std::vector<int> m_boundary; //!< boundary condition.

  std::string m_mode;

  Fopr_Wilson_eo     *m_fopr_w;
  Fopr_CloverTerm_eo *m_fopr_csw;

  Index_eo m_idx;
  Field_G  *m_Ueo;

  Field   m_w1;                  //!< working field.
  Field_F m_vF1, m_vF2, m_vF3;   //!< working field.

  void (Fopr_Clover_eo::*m_mult)(Field&, const Field&);
  void (Fopr_Clover_eo::*m_mult_dag)(Field&, const Field&);
  void (Fopr_Clover_eo::*m_preProp)(Field&, Field&, const Field&);
  void (Fopr_Clover_eo::*m_postProp)(Field&, const Field&, const Field&);

 public:
  Fopr_Clover_eo(std::string repr) { init(repr); }

  ~Fopr_Clover_eo() { tidyup(); }

  void set_parameters(const Parameters& params);
  void set_parameters(const double kappa, const double cSW,
                      const std::vector<int> bc);

  void set_config(Field *U);

  void set_config(unique_ptr<Field_G>& U)
  {
    set_config(U.get());
  }

  void set_mode(std::string mode)
  {
    m_mode = mode;

    if (m_mode == "D") {
      m_mult     = &Fopr_Clover_eo::D;
      m_mult_dag = &Fopr_Clover_eo::Ddag;
      m_preProp  = &Fopr_Clover_eo::prePropD;
      m_postProp = &Fopr_Clover_eo::postPropD;
    } else if (m_mode == "Ddag") {
      m_mult     = &Fopr_Clover_eo::Ddag;
      m_mult_dag = &Fopr_Clover_eo::D;
      m_preProp  = &Fopr_Clover_eo::prePropDag;
      m_postProp = &Fopr_Clover_eo::postPropDag;
    } else if (m_mode == "DdagD") {
      m_mult     = &Fopr_Clover_eo::DdagD;
      m_mult_dag = &Fopr_Clover_eo::DdagD;
    } else if (m_mode == "DDdag") {
      m_mult     = &Fopr_Clover_eo::DDdag;
      m_mult_dag = &Fopr_Clover_eo::DDdag;
    } else if (m_mode == "H") {
      m_mult     = &Fopr_Clover_eo::H;
      m_mult_dag = &Fopr_Clover_eo::H;
    } else {
      vout.crucial("Error at %s: undefined mode = %s\n", class_name.c_str(), mode.c_str());
      exit(EXIT_FAILURE);
    }
  }

  std::string get_mode() const
  { return m_mode; }

  void mult(Field& v, const Field& f)
  { (this->*m_mult)(v, f); }

  void mult_dag(Field& v, const Field& f)
  { (this->*m_mult_dag)(v, f); }

  //- method for even odd fermion operator
  void preProp(Field& Be, Field& bo, const Field& b)
  { (this->*m_preProp)(Be, bo, b); }

  void postProp(Field& x, const Field& xe, const Field& bo)
  { (this->*m_postProp)(x, xe, bo); }

  void prePropD(Field&, Field&, const Field&);
  void postPropD(Field&, const Field&, const Field&);
  void prePropDag(Field&, Field&, const Field&);
  void postPropDag(Field&, const Field&, const Field&);

  const Field_F mult_csw_inv(const Field_F&, const int ieo);

  void D(Field& v, const Field& f);
  void Ddag(Field& v, const Field& f);
  void DdagD(Field& v, const Field& f);
  void DDdag(Field& v, const Field& f);
  void H(Field& v, const Field& f);
  void mult_gm5(Field& v, const Field& f);
  void MeoMoe(Field& v, const Field& f);

  // ieo=0: even <-- odd
  // ieo=1: odd  <-- even
  void Meo(Field&, const Field&, const int ieo);
  void Meo_gm5(Field_F&, const Field_F&, const int ieo);
  void Mdageo(Field_F&, const Field_F&, const int ieo);

  void mult_isigma(Field_F& w, const Field_F& f,
                   const int mu, const int nu);

  inline std::vector<double> csmatrix(const int& site)
  { return m_fopr_csw->csmatrix(site); }

  int field_nvol() { return m_Nvol; }
  int field_nin()  { return 2 * m_Nc * m_Nd; }
  int field_nex()  { return 1; }

  //! this returns the number of floating point operations.
  double flop_count();

 private:
  void init(const std::string repr);
  void tidyup();
};
#endif
