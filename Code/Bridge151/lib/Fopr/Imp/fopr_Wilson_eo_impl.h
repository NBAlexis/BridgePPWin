/*!
        @file    fopr_Wilson_eo_impl.h

        @brief

        @author  UEDA, Satoru (sueda)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2019-01-25 20:04:20 #$

        @version $LastChangedRevision: 1934 $
*/

#ifndef FOPR_WILSON_EO_IMPL_IMP_INCLUDED
#define FOPR_WILSON_EO_IMPL_IMP_INCLUDED

#include "Fopr/fopr_eo.h"
#include "Field/index_eo.h"

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

namespace Imp {
  class BAPI Fopr_Wilson_eo : public Fopr_eo
  {
   public:
    static const std::string class_name;

   private:
    int m_Nc, m_Nd, m_Nvc, m_Ndf;
    int m_Nvol, m_Nvol2, m_Ndim;
    int m_Nx, m_Ny, m_Nz, m_Nt, m_Nx2;

    double m_kappa;                           //!< hopping parameter.
    std::vector<int> m_boundary;              //!< boundary condition.
    std::vector<double> m_boundary_each_node; //!< b.c. for each node.

    std::string m_mode;                       //!< mult mode.
    std::string m_repr;                       //!< Dirac matrix representation.

    Index_eo m_index;
    Field_G *m_Ueo;
    Field_G *m_U;   //!< dummy: pointing m_Ueo.

    std::vector<int> m_yzt_eo;
    Field m_v1, m_v2, m_w1, m_w2;            //!< working field.

    //! arrays for data transfer.
    double *vcp1_xp, *vcp2_xp, *vcp1_xm, *vcp2_xm;
    double *vcp1_yp, *vcp2_yp, *vcp1_ym, *vcp2_ym;
    double *vcp1_zp, *vcp2_zp, *vcp1_zm, *vcp2_zm;
    double *vcp1_tp, *vcp2_tp, *vcp1_tm, *vcp2_tm;

    void (Fopr_Wilson_eo::*m_gm5)(Field&, const Field&);
    void (Fopr_Wilson_eo::*m_gm5_self)(Field&);
    void (Fopr_Wilson_eo::*m_mult_tp)(Field&, const Field&, const int ieo);
    void (Fopr_Wilson_eo::*m_mult_tm)(Field&, const Field&, const int ieo);

    void (Fopr_Wilson_eo::*m_mult)(Field&, const Field&);
    void (Fopr_Wilson_eo::*m_mult_dag)(Field&, const Field&);
    void (Fopr_Wilson_eo::*m_D)(Field&, const Field&);
    void (Fopr_Wilson_eo::*m_preProp)(Field&, Field&, const Field&);
    void (Fopr_Wilson_eo::*m_postProp)(Field&, const Field&, const Field&);

   public:
    Fopr_Wilson_eo()
    { init("Dirac"); }

    Fopr_Wilson_eo(const std::string repr)
    { init(repr); }

    ~Fopr_Wilson_eo()
    { tidyup(); }

    void set_parameters(const Parameters& params);
    void set_parameters(const double kappa, const std::vector<int> bc);

    void set_config(Field *U);

    void set_config(unique_ptr<Field_G>& U)
    { return set_config(U.get()); }

    void set_mode(const std::string mode);
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
    void MeoMoe(Field&, const Field&);
    void Meo_gm5(Field&, const Field&, const int ieo);

    void mult_gm5(Field&, const Field&);
    void mult_gm5(Field&);
    void gm5_dirac(Field&, const Field&);
    void gm5_chiral(Field&, const Field&);
    void gm5_self_dirac(Field&);
    void gm5_self_chiral(Field&);

    //! gamma_5 (1 - gamma_mu) v(x + mu) used in force calculation.
    void gm5p(const int mu, Field&, const Field& v);

    int field_nvol() { return m_Nvol2; }
    int field_nin()  { return m_Nvc * m_Nd; }
    int field_nex()  { return 1; }

    //! this returns the number of floating point operations of Meo.
    double flop_count();

   private:
    void init(const std::string);
    void tidyup();

    void mult_xp(Field&, const Field&, const int ieo);
    void mult_xm(Field&, const Field&, const int ieo);
    void mult_yp(Field&, const Field&, const int ieo);
    void mult_ym(Field&, const Field&, const int ieo);
    void mult_zp(Field&, const Field&, const int ieo);
    void mult_zm(Field&, const Field&, const int ieo);
    void mult_tp_dirac(Field&, const Field&, const int ieo);
    void mult_tm_dirac(Field&, const Field&, const int ieo);
    void mult_tp_chiral(Field&, const Field&, const int ieo);
    void mult_tm_chiral(Field&, const Field&, const int ieo);

    void clear_impl(Field&);
    void scal_impl(Field&, const double);

    // member data for threading
    int m_Mz, m_Mt;
    int m_Nthread, m_Ntask;
    int m_Ntask_z, m_Ntask_t;

    struct mult_arg
    {
      int isite;
      int isite_cp_x, isite_cp_y, isite_cp_z, isite_cp_t;
      int kz0, kz1, kt0, kt1;
    };
    std::vector<mult_arg> m_arg;

    // member functions for threading
    void setup_thread();

    void mult_xp1_thread(const int, double *, const double *, const int);
    void mult_xp2_thread(const int, double *, const double *, const int);
    void mult_xpb_thread(const int, double *, const double *, const int);

    void mult_xm1_thread(const int, double *, const double *, const int);
    void mult_xm2_thread(const int, double *, const double *, const int);
    void mult_xmb_thread(const int, double *, const double *, const int);

    void mult_yp1_thread(const int, double *, const double *, const int);
    void mult_yp2_thread(const int, double *, const double *, const int);
    void mult_ypb_thread(const int, double *, const double *, const int);

    void mult_ym1_thread(const int, double *, const double *, const int);
    void mult_ym2_thread(const int, double *, const double *, const int);
    void mult_ymb_thread(const int, double *, const double *, const int);

    void mult_zp1_thread(const int, double *, const double *, const int);
    void mult_zp2_thread(const int, double *, const double *, const int);
    void mult_zpb_thread(const int, double *, const double *, const int);

    void mult_zm1_thread(const int, double *, const double *, const int);
    void mult_zm2_thread(const int, double *, const double *, const int);
    void mult_zmb_thread(const int, double *, const double *, const int);

    void mult_tp1_dirac_thread(const int, double *, const double *, const int);
    void mult_tp2_dirac_thread(const int, double *, const double *, const int);
    void mult_tpb_dirac_thread(const int, double *, const double *, const int);

    void mult_tm1_dirac_thread(const int, double *, const double *, const int);
    void mult_tm2_dirac_thread(const int, double *, const double *, const int);
    void mult_tmb_dirac_thread(const int, double *, const double *, const int);

    void mult_tp1_chiral_thread(const int, double *, const double *, const int);
    void mult_tp2_chiral_thread(const int, double *, const double *, const int);
    void mult_tpb_chiral_thread(const int, double *, const double *, const int);

    void mult_tm1_chiral_thread(const int, double *, const double *, const int);
    void mult_tm2_chiral_thread(const int, double *, const double *, const int);
    void mult_tmb_chiral_thread(const int, double *, const double *, const int);

    void scal_thread(const int, double *, const double);
    void clear_thread(const int, double *);

    void gm5_dirac_thread(const int, double *, const double *);
    void gm5_chiral_thread(const int, double *, const double *);
    void gm5_dirac_thread(const int, double *);
    void gm5_chiral_thread(const int, double *);

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
      bool init1 = Fopr::Factory_noarg::Register("Wilson_eo/Imp", create_object);
      bool init2 = Fopr::Factory_string::Register("Wilson_eo/Imp", create_object_with_repr);

      return init1 && init2;
    }
#endif
  };
}
#endif /* FOPR_WILSON_EO_IMPL_IMP_INCLUDED */
