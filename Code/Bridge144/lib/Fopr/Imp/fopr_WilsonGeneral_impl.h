/*!
        @file    $Id:: fopr_WilsonGeneral_impl.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FOPR_WILSON_GENERAL_IMP_INCLUDED
#define FOPR_WILSON_GENERAL_IMP_INCLUDED

#include "Fopr/fopr.h"

#include "Tools/gammaMatrixSet.h"
#include "Field/shiftField_lex.h"

#include "ResourceManager/threadManager_OpenMP.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wilson_General fermion operator.

/*!
    This class implements the Wilson_General fermion operator,
    including Wilson fermion on anisotropic lattice, and
    relativistic heavy quarks.
                                    [21 Mar 2015 Y.Namekawa]
 */

namespace Imp {
  class Fopr_WilsonGeneral : public Fopr
  {
   public:
    static const std::string class_name;

   private:
    // lattice parameters
    int m_Nc, m_Nd, m_Nvc, m_Ndf;
    int m_Nx, m_Ny, m_Nz, m_Nt;
    int m_Nvol, m_Ndim;

    // physical parameters
    double              m_kappa_s, m_kappa_t;
    double              m_nu_s, m_r_s;
    std::vector<int>    m_boundary;  //!< boundary condition.
    std::vector<double> m_boundary2; //!< b.c. for each node.

    const Field_G            *m_U;   //!< gauge configuration.
    std::vector<GammaMatrix> m_GM;   //!< gamma matrices.

    Bridge::VerboseLevel m_vl;

    std::string m_mode;
    std::string m_repr;

    Field m_w1, m_w2; //!< temporary fields.

    //! arrays for data transfer.
    double *vcp1_x_plus, *vcp2_x_plus, *vcp1_x_minus, *vcp2_x_minus;
    double *vcp1_y_plus, *vcp2_y_plus, *vcp1_y_minus, *vcp2_y_minus;
    double *vcp1_z_plus, *vcp2_z_plus, *vcp1_z_minus, *vcp2_z_minus;
    double *vcp1_t_plus, *vcp2_t_plus, *vcp1_t_minus, *vcp2_t_minus;

   public:
    Fopr_WilsonGeneral() { init("Dirac"); }
    Fopr_WilsonGeneral(std::string repr) { init(repr); }
    ~Fopr_WilsonGeneral() { tidyup(); }

    void init(std::string repr);

    void tidyup();

    void set_parameters(const Parameters& params);
    void set_parameters(const double kappa_s, const double kappa_t,
                        const double nu_s, const double r_s,
                        const std::vector<int> bc);

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
      vout.crucial(m_vl, "Error at %s: mode undefined.\n", class_name.c_str());
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

    // inline void fprop_normalize(Field& v)
    // { v *= (2.0 * m_kappa); }
    //
    // inline void fopr_normalize(Field& v)
    // { v *= 1.0 / (2.0 * m_kappa); }

    // const double get_fprop_normfactor(){ return 2.0 * m_kappa; }
    // const double get_fopr_normfactor(){  return  1.0/(2.0 * m_kappa); }

   private:
    // prohibit copy
    Fopr_WilsonGeneral(const Fopr_WilsonGeneral&) {}
    Fopr_WilsonGeneral& operator=(const Fopr_WilsonGeneral&);

    void (Fopr_WilsonGeneral::*m_mult)(Field&, const Field&);
    void (Fopr_WilsonGeneral::*m_mult_dag)(Field&, const Field&);
    void (Fopr_WilsonGeneral::*m_D)(Field&, const Field&);
    void (Fopr_WilsonGeneral::*m_gm5)(Field&, const Field&);
    void (Fopr_WilsonGeneral::*m_mult_t_plus)(Field&, const Field&);
    void (Fopr_WilsonGeneral::*m_mult_t_minus)(Field&, const Field&);
    void (Fopr_WilsonGeneral::*m_D_ex)(Field&, const int,
                                       const Field&, const int);

    void D_chiral(Field&, const Field&);
    void D_dirac(Field&, const Field&);
    void gm5_chiral(Field&, const Field&);
    void gm5_dirac(Field&, const Field&);

    void D_ex_chiral(Field&, const int ex1, const Field&, const int ex2);
    void D_ex_dirac(Field&, const int ex1, const Field&, const int ex2);

    void mult_p(int mu, Field_F&, const Field_F&);
    void mult_m(int mu, Field_F&, const Field_F&);

    void mult_x_plus(Field&, const Field&);
    void mult_x_minus(Field&, const Field&);
    void mult_y_plus(Field&, const Field&);
    void mult_y_minus(Field&, const Field&);
    void mult_z_plus(Field&, const Field&);
    void mult_z_minus(Field&, const Field&);

    void mult_t_plus_dirac(Field&, const Field&);
    void mult_t_minus_dirac(Field&, const Field&);
    void mult_t_plus_chiral(Field&, const Field&);
    void mult_t_minus_chiral(Field&, const Field&);

    void daxpy(Field&, double, const Field&);
    void daypx(Field&, double, const Field&);
    void scal(Field&, double);
    void clear(Field&);

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
    void mult_x_plus1_thread(int, double *, const double *);
    void mult_x_plus2_thread(int, double *, const double *);
    void mult_x_plus_bulk_thread(int, double *, const double *);
    void mult_x_minus1_thread(int, double *, const double *);
    void mult_x_minus2_thread(int, double *, const double *);
    void mult_x_minus_bulk_thread(int, double *, const double *);
    void mult_y_plus1_thread(int, double *, const double *);
    void mult_y_plus2_thread(int, double *, const double *);
    void mult_y_plus_bulk_thread(int, double *, const double *);
    void mult_y_minus1_thread(int, double *, const double *);
    void mult_y_minus2_thread(int, double *, const double *);
    void mult_y_minus_bulk_thread(int, double *, const double *);
    void mult_z_plus1_thread(int, double *, const double *);
    void mult_z_plus2_thread(int, double *, const double *);
    void mult_z_plus_bulk_thread(int, double *, const double *);
    void mult_z_minus1_thread(int, double *, const double *);
    void mult_z_minus2_thread(int, double *, const double *);
    void mult_z_minus_bulk_thread(int, double *, const double *);
    void mult_t_plus1_dirac_thread(int, double *, const double *);
    void mult_t_plus2_dirac_thread(int, double *, const double *);
    void mult_t_plus_bulk_dirac_thread(int, double *, const double *);
    void mult_t_minus1_dirac_thread(int, double *, const double *);
    void mult_t_minus2_dirac_thread(int, double *, const double *);
    void mult_t_minus_bulk_dirac_thread(int, double *, const double *);
    void mult_t_plus1_chiral_thread(int, double *, const double *);
    void mult_t_plus2_chiral_thread(int, double *, const double *);
    void mult_t_plus_bulk_chiral_thread(int, double *, const double *);
    void mult_t_minus1_chiral_thread(int, double *, const double *);
    void mult_t_minus2_chiral_thread(int, double *, const double *);
    void mult_t_minus_bulk_chiral_thread(int, double *, const double *);

    void daxpy_thread(int, double *, double, const double *);
    void daypx_thread(int, double *, double, const double *);
    void scal_thread(int, double *, double);
    void clear_thread(int, double *);

    void gm5_dirac_thread(int, double *, const double *);
    void gm5_chiral_thread(int, double *, const double *);

    //- for derivatives
    // void mult_xpu(Field&, const Field&);
    // void mult_ypu(Field&, const Field&);
    // void mult_zpu(Field&, const Field&);
    // void mult_tpu_dirac(Field&, const Field&);
    // void mult_tpu_chiral(Field&, const Field&);
  };
}
#endif
