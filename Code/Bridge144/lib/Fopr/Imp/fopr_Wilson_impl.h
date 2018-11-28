/*!
        @file    $Id:: fopr_Wilson_impl.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FOPR_WILSON_IMPL_INCLUDED
#define FOPR_WILSON_IMPL_INCLUDED

//#include "Fopr/fopr_Wilson.h"
#include "Fopr/fopr.h"

#include "Tools/gammaMatrixSet.h"
#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wilson fermion operator.

/*!
    This fermion operator defines the standard Wilson fermion.
    The gamma matrix representation is given as control string
    "Dirac"(default) or "Chiral" at the construction, which is
    used to construct the Fopr_Wilson instance.
    The `mode', which of D, Ddag, H, DdagD are multiplied, is
    controlled by setting the pointers to these functions,
    m_mult and m_mult_dag.
    At the beginning, they are set to point mult_undef() which
    just represent the mode has not been set.
    set_mode(string) must be called before mult() is called.
                                    [24 Dec 2011 H,Matsufuru]
 */

namespace Imp {
  class Fopr_Wilson : public Fopr
  {
   public:
    static const std::string class_name;

   private:
    // lattice parameters
    int m_Nc, m_Nd, m_Nvc, m_Ndf;
    int m_Nx, m_Ny, m_Nz, m_Nt;
    int m_Nvol, m_Ndim;

    // physical parameters
    double              m_kappa;
    std::vector<int>    m_boundary;  //!< boundary condition.
    std::vector<double> m_boundary2; //!< b.c. for each node.

    const Field_G            *m_U;   //!< gauge configuration.
    std::vector<GammaMatrix> m_GM;   //!< gamma matrices.

    Bridge::VerboseLevel m_vl;

    std::string m_mode;
    std::string m_repr;

    Field m_w1, m_w2; //!< temporary fields.

    //! arrays for data transfer.
    double *vcp1_xp, *vcp2_xp, *vcp1_xm, *vcp2_xm;
    double *vcp1_yp, *vcp2_yp, *vcp1_ym, *vcp2_ym;
    double *vcp1_zp, *vcp2_zp, *vcp1_zm, *vcp2_zm;
    double *vcp1_tp, *vcp2_tp, *vcp1_tm, *vcp2_tm;

   public:

    Fopr_Wilson() { init("Dirac"); }

    Fopr_Wilson(std::string repr) { init(repr); }

    ~Fopr_Wilson() { tidyup(); }

    void init(std::string repr);

    void tidyup();

    void set_parameters(const Parameters& params);
    void set_parameters(const double kappa, const std::vector<int> bc);

    void set_config(Field *U)
    { m_U = (Field_G *)U; }

    void set_config(unique_ptr<Field_G>& U)
    { return set_config(U.get()); }

    void set_mode(std::string mode);

    std::string get_mode() const;

    inline void mult(Field& v, const Field& f)
    { (this->*m_mult)(v, f); }

    inline void mult_dag(Field& v, const Field& f)
    { (this->*m_mult_dag)(v, f); }

    inline void mult_gm5(Field& v, const Field& f)
    { (this->*m_gm5)(v, f); }

    inline void D(Field& v, const Field& f)
    { (this->*m_D)(v, f); }

    inline void Ddag(Field& w, const Field& f)
    {
      mult_gm5(w, f);
      D(m_w1, w);
      mult_gm5(w, m_w1);
    }

    inline void DdagD(Field& w, const Field& f)
    {
      D(m_w1, f);

      // Ddag(w, m_w1);
      mult_gm5(w, m_w1);
      D(m_w1, w);
      mult_gm5(w, m_w1);
    }

    inline void DDdag(Field& w, const Field& f)
    {
      // Ddag(m_w1, f);
      mult_gm5(m_w1, f);
      D(w, m_w1);
      mult_gm5(m_w1, w);

      D(w, m_w1);
    }

    inline void H(Field& w, const Field& f)
    {
      D(m_w1, f);
      mult_gm5(w, m_w1);
    }

    inline void mult_undef(Field&, const Field& f)
    {
      vout.crucial(m_vl, "Error at Fopr_Wilson: mode undefined.\n");
      exit(EXIT_FAILURE);
    }

    inline void D_ex(Field& v, const int ex1,
                     const Field& f, const int ex2)
    { (this->*m_D_ex)(v, ex1, f, ex2); }

    const Field_F mult_gm5p(int mu, const Field_F& w)
    {
      mult_gm5p(mu, m_w1, w);
      return (Field_F)m_w1;
    }

    void mult_gm5p(int mu, Field&, const Field&);

    // void gm5p(Field_F& v, const int mu, const Field_F& w);

    void proj_chiral(Field& w, const int ex1,
                     const Field& v, const int ex2, const int ipm);

    void mult_up(int mu, Field&, const Field&);
    void mult_dn(int mu, Field&, const Field&);

    int field_nvol()
    { return CommonParameters::Nvol(); }
    int field_nin()
    { return 2 * CommonParameters::Nc() * CommonParameters::Nd(); }
    int field_nex()
    { return 1; }

    double flop_count();

    inline void fprop_normalize(Field& v)
    { scal(v, 2.0 * m_kappa); } //v *= (2.0 * m_kappa);

    inline void fopr_normalize(Field& v)
    { scal(v, 1.0 / (2.0 * m_kappa)); } //v *= 1.0 / (2.0 * m_kappa);

    // const double get_fprop_normfactor(){ return 2.0 * m_kappa; }
    // const double get_fopr_normfactor(){  return  1.0/(2.0 * m_kappa); }

   private:

    // prohibit copy
    Fopr_Wilson(const Fopr_Wilson&) {}
    Fopr_Wilson& operator=(const Fopr_Wilson&);

    void (Fopr_Wilson::*m_mult)(Field&, const Field&);
    void (Fopr_Wilson::*m_mult_dag)(Field&, const Field&);
    void (Fopr_Wilson::*m_D)(Field&, const Field&);
    void (Fopr_Wilson::*m_gm5)(Field&, const Field&);
    void (Fopr_Wilson::*m_mult_tp)(Field&, const Field&);
    void (Fopr_Wilson::*m_mult_tm)(Field&, const Field&);
    void (Fopr_Wilson::*m_D_ex)(Field&, const int,
                                const Field&, const int);

    void D_chiral(Field&, const Field&);
    void D_dirac(Field&, const Field&);
    void gm5_chiral(Field&, const Field&);
    void gm5_dirac(Field&, const Field&);

    void D_ex_chiral(Field&, const int ex1, const Field&, const int ex2);
    void D_ex_dirac(Field&, const int ex1, const Field&, const int ex2);

    void mult_p(int mu, Field_F&, const Field_F&);
    void mult_m(int mu, Field_F&, const Field_F&);

    void mult_xp(Field&, const Field&);
    void mult_xm(Field&, const Field&);
    void mult_yp(Field&, const Field&);
    void mult_ym(Field&, const Field&);
    void mult_zp(Field&, const Field&);
    void mult_zm(Field&, const Field&);

    void mult_tp_dirac(Field&, const Field&);
    void mult_tm_dirac(Field&, const Field&);
    void mult_tp_chiral(Field&, const Field&);
    void mult_tm_chiral(Field&, const Field&);

    void daypx(Field&, double, const Field&);
    void clear(Field&);

    // member data for threading
    int m_Mz, m_Mt;
    int m_Nthread, m_Ntask;
    int m_Ntask_z, m_Ntask_t;

    struct mult_arg
    {
      int isite;
      int isite_cpx, isite_cpy, isite_cpz, isite_cpt;
      int kz0, kz1, kt0, kt1;
    };
    std::vector<mult_arg> m_arg;

    // member functions for threading
    void setup_thread();
    void mult_xp1_thread(int, double *, const double *);
    void mult_xp2_thread(int, double *, const double *);
    void mult_xpb_thread(int, double *, const double *);
    void mult_xm1_thread(int, double *, const double *);
    void mult_xm2_thread(int, double *, const double *);
    void mult_xmb_thread(int, double *, const double *);
    void mult_yp1_thread(int, double *, const double *);
    void mult_yp2_thread(int, double *, const double *);
    void mult_ypb_thread(int, double *, const double *);
    void mult_ym1_thread(int, double *, const double *);
    void mult_ym2_thread(int, double *, const double *);
    void mult_ymb_thread(int, double *, const double *);
    void mult_zp1_thread(int, double *, const double *);
    void mult_zp2_thread(int, double *, const double *);
    void mult_zpb_thread(int, double *, const double *);
    void mult_zm1_thread(int, double *, const double *);
    void mult_zm2_thread(int, double *, const double *);
    void mult_zmb_thread(int, double *, const double *);
    void mult_tp1_dirac_thread(int, double *, const double *);
    void mult_tp2_dirac_thread(int, double *, const double *);
    void mult_tpb_dirac_thread(int, double *, const double *);
    void mult_tm1_dirac_thread(int, double *, const double *);
    void mult_tm2_dirac_thread(int, double *, const double *);
    void mult_tmb_dirac_thread(int, double *, const double *);
    void mult_tp1_chiral_thread(int, double *, const double *);
    void mult_tp2_chiral_thread(int, double *, const double *);
    void mult_tpb_chiral_thread(int, double *, const double *);
    void mult_tm1_chiral_thread(int, double *, const double *);
    void mult_tm2_chiral_thread(int, double *, const double *);
    void mult_tmb_chiral_thread(int, double *, const double *);
    void daypx_thread(int, double *, double, const double *);
    void clear_thread(int, double *);
    void gm5_dirac_thread(int, double *, const double *);
    void gm5_chiral_thread(int, double *, const double *);

    // for derivatives
    void mult_xpu(Field&, const Field&);
    void mult_ypu(Field&, const Field&);
    void mult_zpu(Field&, const Field&);
    void mult_tpu_dirac(Field&, const Field&);
    void mult_tpu_chiral(Field&, const Field&);
  };
}
#endif
