/*!
        @file    fopr_Wilson_eo_impl.h

        @brief

        @author  UEDA, Satoru
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2019-01-25 20:04:20 #$

        @version $LastChangedRevision: 1934 $
*/

#ifndef FOPR_WILSON_EO_IMPL_ORG_INCLUDED
#define FOPR_WILSON_EO_IMPL_ORG_INCLUDED

#include "Fopr/fopr_eo.h"

#include "Field/shiftField_eo.h"
#include "Tools/gammaMatrixSet.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Implementation of even-odd Wilson fermion operator.

/*!
    This class BAPI is a subclass of Fopr_Wilson_eo and implements
    an even-odd version of the Wilson fermion operator.
    This is rather straightforward and readable while slower
    version that was coded by S.Ueda [20 Jun 2012 S.UEDA].
    The implementation class BAPI was separated.
                                     [07 Jul 2014 H.Matsufuru]
 */

namespace Org {
  class BAPI Fopr_Wilson_eo : public Fopr_eo
  {
   public:
    static const std::string class_name;

   private:
    int m_Nvol, m_Nvol2, m_Ndim;
    int m_Nc, m_Nd;
    double m_kappa;
    std::vector<int> m_boundary;
    Index_eo m_index;
    Field_G *m_Ueo;

    ShiftField_eo shift;
    Field_F trf, trf2;
    Vec_SU_N v1, v2;

    std::vector<GammaMatrix> m_GM;

    std::string m_mode;
    std::string m_repr;

   public:
    Fopr_Wilson_eo()
      : m_Nvol(CommonParameters::Nvol()),
      m_Nvol2(CommonParameters::Nvol() / 2),
      m_Ndim(CommonParameters::Ndim()),
      m_Nc(CommonParameters::Nc()),
      m_Nd(CommonParameters::Nd()),
      trf(m_Nvol2, 1), trf2(m_Nvol2, 1)
    {
      init("Dirac");
    }

    Fopr_Wilson_eo(const std::string repr)
      : m_Nvol(CommonParameters::Nvol()),
      m_Nvol2(CommonParameters::Nvol() / 2),
      m_Ndim(CommonParameters::Ndim()),
      m_Nc(CommonParameters::Nc()),
      m_Nd(CommonParameters::Nd()),
      trf(m_Nvol2, 1), trf2(m_Nvol2, 1)
    {
      init(repr);
    }

    ~Fopr_Wilson_eo()
    {
      delete m_Ueo;
    }

    void set_parameters(const Parameters& params);
    void set_parameters(const double kappa, const std::vector<int> bc);

    void set_config(Field *U);

    void set_config(unique_ptr<Field_G>& U)
    { return set_config(U.get()); }

    void set_mode(std::string mode) override;
    std::string get_mode() const;

    // method for even odd fermion operator
    void preProp(Field& Be, Field& bo, const Field& b)
    { (this->*m_preProp)(Be, bo, b); }

    void postProp(Field& x, const Field& xe, const Field& bo)
    { (this->*m_postProp)(x, xe, bo); }

    void prePropD(Field&, Field&, const Field&);
    void postPropD(Field&, const Field&, const Field&);
    void prePropDag(Field&, Field&, const Field&);
    void postPropDag(Field&, const Field&, const Field&);

    void mult(Field& v, const Field& f)
    { (this->*m_mult)(v, f); }

    void mult_dag(Field& v, const Field& f)
    { (this->*m_mult_dag)(v, f); }

    void D(Field& v, const Field& f);
    void Ddag(Field& v, const Field& f);
    void DdagD(Field& v, const Field& f);
    void DDdag(Field& v, const Field& f);
    void H(Field& v, const Field& f);

    // ieo=0: even <-- odd
    // ieo=1: odd  <-- even

    void Meo(Field&, const Field&, const int ieo);
    void Mdageo(Field&, const Field&, const int ieo);
    void MeoMoe(Field& v, const Field& f);
    void Meo_gm5(Field&, const Field&, const int ieo);

    void mult_gm5(Field&, const Field&);

    //! gamma_5 (1 - gamma_mu) v(x + mu)
    void gm5p(const int mu, Field&, const Field& v);

    int field_nvol() { return CommonParameters::Nvol() / 2; }
    int field_nin() { return 2 * CommonParameters::Nc() * CommonParameters::Nd(); }
    int field_nex() { return 1; }

    //! this returns the number of floating point operations of Meo.
    double flop_count();

   private:
    void init(const std::string);

    void (Fopr_Wilson_eo::*m_mult)(Field&, const Field&);
    void (Fopr_Wilson_eo::*m_mult_dag)(Field&, const Field&);
    void (Fopr_Wilson_eo::*m_D)(Field&, const Field&);
    void (Fopr_Wilson_eo::*m_gm5)(Field&, const Field&);
    void (Fopr_Wilson_eo::*m_preProp)(Field&, Field&, const Field&);
    void (Fopr_Wilson_eo::*m_postProp)(Field&, const Field&, const Field&);

    void mult_p(const int mu, Field_F&, const Field_F&, const int ieo);
    void mult_m(const int mu, Field_F&, const Field_F&, const int ieo);

    /*
    void mult_xp(Field_F&, const Field_F&, const int ieo);
    void mult_xm(Field_F&, const Field_F&, const int ieo);
    void mult_yp(Field_F&, const Field_F&, const int ieo);
    void mult_ym(Field_F&, const Field_F&, const int ieo);
    void mult_zp(Field_F&, const Field_F&, const int ieo);
    void mult_zm(Field_F&, const Field_F&, const int ieo);
    void mult_tp(Field_F&, const Field_F&, const int ieo);
    void mult_tm(Field_F&, const Field_F&, const int ieo);
    */

#ifdef USE_FACTORY
   private:
    static Fopr *create_object()
    {
      return new Fopr_Wilson_eo();
    }

    static Fopr *create_object_with_repr(const std::string& repr)
    {
      return new Fopr_Wilson_eo(repr);
    }

   public:
    static bool register_factory()
    {
      bool init1 = Fopr::Factory_noarg::Register("Wilson_eo/Org", create_object);
      bool init2 = Fopr::Factory_string::Register("Wilson_eo/Org", create_object_with_repr);

      return init1 && init2;
    }
#endif
  };
}
#endif /* FOPR_WILSON_EO_IMPL_ORG_INCLUDED */
