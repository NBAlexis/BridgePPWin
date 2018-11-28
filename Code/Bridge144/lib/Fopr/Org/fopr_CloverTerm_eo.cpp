#include "BridgeLib_Private.h"
#if USE_ORG

/*!
        @file    $Id:: fopr_CloverTerm_eo.cpp #$

        @brief

        @author  UEDA, Satoru (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "Fopr/fopr_CloverTerm_eo.h"
#include "Measurements/Gauge/staple_eo.h"

#include "ResourceManager/threadManager_OpenMP.h"
#include "Solver/solver_CG.h"

namespace Org {
//====================================================================

  const std::string Fopr_CloverTerm_eo::class_name = "Org::Fopr_CloverTerm_eo";

//====================================================================
  void Fopr_CloverTerm_eo::init(std::string repr)
  {
    m_repr = repr;

    m_Nc    = CommonParameters::Nc();
    m_Nd    = CommonParameters::Nd();
    m_Ndim  = CommonParameters::Ndim();
    m_NinF  = 2 * m_Nc * m_Nd;
    m_Nvol  = CommonParameters::Nvol();
    m_Nvol2 = m_Nvol / 2;

    m_boundary.resize(m_Ndim);

    m_Ueo = 0;

    m_GM.resize(m_Ndim + 1);
    m_SG.resize(m_Ndim * m_Ndim);

    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(m_repr));

    m_GM[0] = gmset->get_GM(gmset->GAMMA1);
    m_GM[1] = gmset->get_GM(gmset->GAMMA2);
    m_GM[2] = gmset->get_GM(gmset->GAMMA3);
    m_GM[3] = gmset->get_GM(gmset->GAMMA4);
    m_GM[4] = gmset->get_GM(gmset->GAMMA5);

    m_SG[sg_index(0, 1)] = gmset->get_GM(gmset->SIGMA12);
    m_SG[sg_index(1, 2)] = gmset->get_GM(gmset->SIGMA23);
    m_SG[sg_index(2, 0)] = gmset->get_GM(gmset->SIGMA31);
    m_SG[sg_index(3, 0)] = gmset->get_GM(gmset->SIGMA41);
    m_SG[sg_index(3, 1)] = gmset->get_GM(gmset->SIGMA42);
    m_SG[sg_index(3, 2)] = gmset->get_GM(gmset->SIGMA43);

    m_SG[sg_index(1, 0)] = m_SG[sg_index(0, 1)].mult(-1);
    m_SG[sg_index(2, 1)] = m_SG[sg_index(1, 2)].mult(-1);
    m_SG[sg_index(0, 2)] = m_SG[sg_index(2, 0)].mult(-1);
    m_SG[sg_index(0, 3)] = m_SG[sg_index(3, 0)].mult(-1);
    m_SG[sg_index(1, 3)] = m_SG[sg_index(3, 1)].mult(-1);
    m_SG[sg_index(2, 3)] = m_SG[sg_index(3, 2)].mult(-1);

    m_SG[sg_index(0, 0)] = gmset->get_GM(gmset->UNITY);
    m_SG[sg_index(1, 1)] = gmset->get_GM(gmset->UNITY);
    m_SG[sg_index(2, 2)] = gmset->get_GM(gmset->UNITY);
    m_SG[sg_index(3, 3)] = gmset->get_GM(gmset->UNITY);
    // these 4 gamma matrices are actually not used.

    m_fee_inv = new Field_F(m_Nvol2, m_Nc * m_Nd);
    m_foo_inv = new Field_F(m_Nvol2, m_Nc * m_Nd);

    m_vf.reset(m_Nvol2, 1);
    m_ff.reset(m_Nvol2, 1);

    int Nfst = 6;
    m_T2.resize(2); // 0: even, 1: odd.
    m_T2[0].reset(m_Nvol2, Nfst);
    m_T2[1].reset(m_Nvol2, Nfst);
  }


//====================================================================
  void Fopr_CloverTerm_eo::tidyup()
  {
    delete m_foo_inv;
    delete m_fee_inv;
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_parameters(const Parameters& params)
  {
    const string str_vlevel = params.get_string("verbose_level");

    m_vl = vout.set_verbose_level(str_vlevel);

    //- fetch and check input parameters
    double           kappa, cSW;
    std::vector<int> bc;

    int err = 0;
    err += params.fetch_double("hopping_parameter", kappa);
    err += params.fetch_double("clover_coefficient", cSW);
    err += params.fetch_int_vector("boundary_condition", bc);

    if (err) {
      vout.crucial(m_vl, "Error at %s: input parameter not found\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    set_parameters(kappa, cSW, bc);
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_parameters(const double kappa, const double cSW,
                                          const std::vector<int> bc)
  {
    //- print input parameters
    vout.general(m_vl, "%s:\n", class_name.c_str());
    vout.general(m_vl, "  kappa = %12.8f\n", kappa);
    vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
    }

    //- range check
    // NB. kappa,cSW == 0 is allowed.
    assert(bc.size() == m_Ndim);

    //- store values
    m_kappa = kappa;
    m_cSW   = cSW;
    assert(bc.size() == m_Ndim);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_config(Field *Ueo)
  {
    m_Ueo = (Field_G *)Ueo;

    set_csw();
    solve_csw_inv();
  }


//====================================================================
  void Fopr_CloverTerm_eo::solve_csw_inv()
  {
    double eps2 = CommonParameters::epsilon_criterion2();

#if 1
    Parameters params_solver;

    params_solver.set_string("solver_type", "CG");
    params_solver.set_int("maximum_number_of_iteration", 100);
    params_solver.set_int("maximum_number_of_restart", 40);
    params_solver.set_double("convergence_criterion_squared", 1.0e-30);
    //- NB. set VerboseLevel to CRUCIAL to suppress frequent messages.
    params_solver.set_string("verbose_level", "Crucial");
#else
    //
#endif

    int    Nconv;
    double diff;

    unique_ptr<Solver> solver(new Solver_CG(this));

#if 1
    solver->set_parameters(params_solver);
#else
    const int    Niter              = 100;
    const int    Nrestart           = 40;
    const double Stopping_condition = 1.0e-30;

    solver->set_parameters(Niter, Nrestart, Stopping_condition);
    solver->set_parameter_verboselevel(Bridge::CRUCIAL);
#endif

    Field_F w(m_Nvol2);
    Field_F w2(m_Nvol2);

    for (int ispin = 0; ispin < m_Nd; ++ispin) {
      for (int icolor = 0; icolor < m_Nc; ++icolor) {
        int spin_color = icolor + m_Nc * ispin;
        w.set(0.0);
        for (int isite = 0; isite < m_Nvol2; ++isite) {
          w.set_ri(icolor, ispin, isite, 0, 1, 0);
        }

        if (m_cSW * m_cSW < eps2) {
          m_fee_inv->setpart_ex(spin_color, w, 0);
          m_foo_inv->setpart_ex(spin_color, w, 0);
        } else {
          //        m_fopr_csw->set_mode("even");
          set_mode("even");
          solver->solve(w2, w, Nconv, diff);
          m_fee_inv->setpart_ex(spin_color, w2, 0);

          vout.detailed(m_vl, "  Nconv,diff = %d %12.4e\n", Nconv, diff);

          set_mode("odd");
          solver->solve(w2, w, Nconv, diff);
          m_foo_inv->setpart_ex(spin_color, w2, 0);

          vout.detailed(m_vl, "  Nconv,diff = %d %12.4e\n", Nconv, diff);
        }
      }
    }

    // redefine the inverse matrix with its dagger.
    double re, im;
    for (int ics = 0; ics < m_Nc * m_Nd; ++ics) {
      for (int site = 0; site < m_Nvol2; ++site) {
        for (int id = 0; id < m_Nd; ++id) {
          for (int ic = 0; ic < m_Nc; ++ic) {
            re = m_foo_inv->cmp_r(ic, id, site, ics);
            im = m_foo_inv->cmp_i(ic, id, site, ics);
            m_foo_inv->set_ri(ic, id, site, ics, re, -im);
            re = m_fee_inv->cmp_r(ic, id, site, ics);
            im = m_fee_inv->cmp_i(ic, id, site, ics);
            m_fee_inv->set_ri(ic, id, site, ics, re, -im);
          }
        }
      }
    }
  }


//====================================================================

/*
const Field_F Fopr_CloverTerm_eo::mult_csw_inv(const Field_F& f, const int ieo)
{
  int     nex = f.nex();
  Field_F w(m_Nvol2, nex);

  mult_csw_inv(w, f, ieo);

  return w;
}
*/

//====================================================================
  void Fopr_CloverTerm_eo::mult_csw_inv(Field& v,
                                        const Field& f, const int ieo)
  {
    int     nex = f.nex();
    Field_F v2(m_Nvol2, nex), f2(m_Nvol2, nex);
    Field_F *csw_inv;

    copy(f2, f);
    v2.set(0.0);

    if (ieo == 0) {
      csw_inv = m_fee_inv;
    } else if (ieo == 1) {
      csw_inv = m_foo_inv;
    } else {
      vout.crucial(m_vl, "Error at %s: wrong parameter, ieo=%d.\n",
                   class_name.c_str(), ieo);
      exit(EXIT_FAILURE);
    }

    for (int iex = 0; iex < nex; ++iex) {
      for (int isite = 0; isite < m_Nvol2; ++isite) {
        for (int ispin = 0; ispin < m_Nd; ++ispin) {
          for (int icolor = 0; icolor < m_Nc; ++icolor) {
            double re = 0.0;
            double im = 0.0;

            for (int jspin = 0; jspin < m_Nd; ++jspin) {
              for (int jcolor = 0; jcolor < m_Nc; ++jcolor) {
                int spin_color = jcolor + m_Nc * jspin;
                // Hermiteness of clover term is used here.
                re += csw_inv->cmp_r(icolor, ispin, isite, spin_color) *
                      f2.cmp_r(jcolor, jspin, isite, iex);
                re += csw_inv->cmp_i(icolor, ispin, isite, spin_color) *
                      f2.cmp_i(jcolor, jspin, isite, iex);

                im += csw_inv->cmp_r(icolor, ispin, isite, spin_color) *
                      f2.cmp_i(jcolor, jspin, isite, iex);
                im -= csw_inv->cmp_i(icolor, ispin, isite, spin_color) *
                      f2.cmp_r(jcolor, jspin, isite, iex);
              }
            }

            v2.set_ri(icolor, ispin, isite, iex, re, im);
          }
        }
      }
    }
    copy(v, v2);
  }


//====================================================================
  std::vector<double> Fopr_CloverTerm_eo::csmatrix(const int& site)
  {
    std::vector<double> matrix(m_Nc * m_Nc * m_Nd * m_Nd * 2);

    for (int ispin = 0; ispin < m_Nd / 2; ++ispin) {
      for (int icolor = 0; icolor < m_Nc; ++icolor) {
        int ics = icolor + ispin * m_Nc;
        for (int jspin = 0; jspin < m_Nd; ++jspin) {
          int js2 = (jspin + m_Nd / 2) % m_Nd;
          for (int jcolor = 0; jcolor < m_Nc; ++jcolor) {
            int cs1 = jcolor + m_Nc * (jspin + m_Nd * ics);
            int cs2 = jcolor + m_Nc * (jspin + m_Nd * (ics + m_Nc * m_Nd / 2));
            int cc  = jcolor + icolor * m_Nc;
            int ss1 = jspin + ispin * m_Nd;
            int ss2 = js2 + ispin * m_Nd;

            matrix[2 * cs1]     = m_T.cmp_r(cc, site, ss1);
            matrix[2 * cs1 + 1] = m_T.cmp_i(cc, site, ss1);

            matrix[2 * cs2]     = m_T.cmp_r(cc, site, ss2);
            matrix[2 * cs2 + 1] = m_T.cmp_i(cc, site, ss2);
          }
        }
      }
    }

    return matrix;
  }


//====================================================================
  void Fopr_CloverTerm_eo::D(Field& v, const Field& f, const int ieo)
  {
    Field_F vf(v.nvol(), v.nex());
    Field_F ff(f.nvol(), f.nex());
    Field_F wt(f.nvol(), f.nex());

    int     Nfst = 6;
    Field_G T2(m_Nvol2, Nfst);
    int     NinG = T2.nin();

    /*
    for(int ex = 0; ex < Nfst; ++ex){
      for(int isite = 0; isite < m_Nvol2; ++isite){
        for(int in = 0; in < NinG; ++in){
        T2.set(in, isite, ex, m_T.cmp(in, idx.site(isite, ieo), ex));
      }
     }
    }
    */

    copy(ff, f); //  ff = (Field_F)f;
    copy(vf, ff);
    double coeff = -m_kappa * m_cSW;

    // i s_23 F_23
    mult_iGM(wt, m_SG[sg_index(1, 2)], ff);
    multadd_Field_Gn(vf, 0, m_T2[ieo], 0, wt, 0, coeff);

    // i s_31 F_31
    mult_iGM(wt, m_SG[sg_index(2, 0)], ff);
    multadd_Field_Gn(vf, 0, m_T2[ieo], 1, wt, 0, coeff);

    // i s_12 F_12
    mult_iGM(wt, m_SG[sg_index(0, 1)], ff);
    multadd_Field_Gn(vf, 0, m_T2[ieo], 2, wt, 0, coeff);

    // i s_41 F_41
    mult_iGM(wt, m_SG[sg_index(3, 0)], ff);
    multadd_Field_Gn(vf, 0, m_T2[ieo], 3, wt, 0, coeff);

    // i s_42 F_42
    mult_iGM(wt, m_SG[sg_index(3, 1)], ff);
    multadd_Field_Gn(vf, 0, m_T2[ieo], 4, wt, 0, coeff);

    // i s_43 F_43
    mult_iGM(wt, m_SG[sg_index(3, 2)], ff);
    multadd_Field_Gn(vf, 0, m_T2[ieo], 5, wt, 0, coeff);

    copy(v, vf);

    /*
    Field_F vf(v.nvol(), v.nex());
    Field_F ff(f.nvol(), f.nex());

    ff = (Field_F)f;

    int n_ex = f.nex();

    Vec_SU_N u_vec, d_vec;
    for (int iex = 0; iex < n_ex; ++iex) {
      for (int isite = 0; isite < m_Nvol2; ++isite) {
        for (int ispin = 0; ispin < m_Nd / 2; ++ispin) {
          u_vec.zero();
          d_vec.zero();
          for (int jspin = 0; jspin < m_Nd; ++jspin) {
            int u_spin = jspin + ispin * m_Nd;
            u_vec += m_T.mat(idx.site(isite, ieo), u_spin)
                     * ff.vec(jspin, isite, iex);
            int d_spin = (jspin + m_Nd / 2) % m_Nd + ispin * m_Nd;
            d_vec += m_T.mat(idx.site(isite, ieo), d_spin)
                     * ff.vec(jspin, isite, iex);
          }
          vf.set_vec(ispin, isite, iex, u_vec);
          vf.set_vec(ispin + m_Nd / 2, isite, iex, d_vec);
        }
      }
    }

    v = (Field)vf;
    */
  }


//====================================================================
  void Fopr_CloverTerm_eo::mult_isigma(Field_F& v, const Field_F& w,
                                       const int mu, const int nu)
  {
    assert(mu != nu);
    mult_iGM(v, m_SG[sg_index(mu, nu)], w);
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_csw()
  {
    m_T.set(0.0);

    Field_G F;

    //- F_23
    set_fieldstrength(F, 1, 2);
    copy(m_T, 0, F, 0);

    //- F_31
    set_fieldstrength(F, 2, 0);
    copy(m_T, 1, F, 0);

    //- F_12
    set_fieldstrength(F, 0, 1);
    copy(m_T, 2, F, 0);

    //- F_41
    set_fieldstrength(F, 3, 0);
    copy(m_T, 3, F, 0);

    //- F_42
    set_fieldstrength(F, 3, 1);
    copy(m_T, 4, F, 0);

    //- F_43
    set_fieldstrength(F, 3, 2);
    copy(m_T, 5, F, 0);

    int Nfst = 6;
    int NinG = m_T2[0].nin();
    for (int ieo = 0; ieo < 2; ++ieo) {
      for (int ex = 0; ex < Nfst; ++ex) {
        for (int isite = 0; isite < m_Nvol2; ++isite) {
          for (int in = 0; in < NinG; ++in) {
            m_T2[ieo].set(in, isite, ex, m_T.cmp(in, idx.site(isite, ieo), ex));
          }
        }
      }
    }


    /*
    //- sigma23
    Field_G F;
    set_fieldstrength(F, 1, 2);
    F.xI();
    m_T.addpart_ex(1, F, 0);
    m_T.addpart_ex(4, F, 0);

    //- sigma31
    set_fieldstrength(F, 2, 0);
    m_T.addpart_ex(1, F, 0);
    scal(F, -1.0);
    m_T.addpart_ex(4, F, 0);

    //- sigma12
    set_fieldstrength(F, 0, 1);
    F.xI();
    m_T.addpart_ex(0, F, 0);
    scal(F, -1.0);
    m_T.addpart_ex(5, F, 0);

    //- sigma41
    set_fieldstrength(F, 3, 0);
    F.xI();
    scal(F, -1.0);  //  F *= -1.0;
    m_T.addpart_ex(3, F, 0);
    m_T.addpart_ex(6, F, 0);

    //- sigma42
    set_fieldstrength(F, 3, 1);
    m_T.addpart_ex(6, F, 0);
    scal(F, -1.0);
    m_T.addpart_ex(3, F, 0);

    //- sigma43
    set_fieldstrength(F, 3, 2);
    F.xI();
    m_T.addpart_ex(7, F, 0);
    scal(F, -1.0);
    m_T.addpart_ex(2, F, 0);

    scal(m_T, -m_kappa*m_cSW);

    //- add unit color matrix for diag part of dirac matrix;
    for (int ispin = 0; ispin < m_Nd / 2; ++ispin) {
      int spin_diag = ispin + m_Nd * ispin;
      for (int isite = 0; isite < m_Nvol; ++isite) {
        for (int icolor = 0; icolor < m_Nc; ++icolor) {
          int cc2 = 2 * (icolor * m_Nc + icolor);
          m_T.add(cc2, isite, spin_diag, 1.0);
        }
      }
    }
    */
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_fieldstrength(Field_G& Fst,
                                             const int mu, const int nu)
  {
    // Staple_eo staple;
    unique_ptr<Staple> staple(Staple::New("EvenOdd"));

    Field_G Cup(m_Nvol, 1), Cdn(m_Nvol, 1);
    Field_G Umu(m_Nvol, 1);
    Field_G w(m_Nvol, 1), v(m_Nvol, 1), v2(m_Nvol, 1);

    staple->upper(Cup, *m_Ueo, mu, nu);
    staple->lower(Cdn, *m_Ueo, mu, nu);
    Umu.setpart_ex(0, *m_Ueo, mu);

    for (int site = 0; site < m_Nvol; ++site) {
      w.set_mat(site, 0, Umu.mat(site) * Cup.mat_dag(site));
    }

    for (int site = 0; site < m_Nvol; ++site) {
      v2.set_mat(site, 0, Umu.mat(site) * Cdn.mat_dag(site));
    }

    axpy(w, -1.0, v2);

    for (int site = 0; site < m_Nvol; ++site) {
      v.set_mat(site, 0, Cup.mat_dag(site) * Umu.mat(site));
    }

    for (int site = 0; site < m_Nvol; ++site) {
      v2.set_mat(site, 0, Cdn.mat_dag(site) * Umu.mat(site));
    }

    axpy(v, -1.0, v2);

    m_shift_eo.forward(v2, v, mu);

    axpy(w, 1.0, v2);

    for (int site = 0; site < m_Nvol; ++site) {
      Fst.set_mat(site, 0, w.mat(site).ah());
    }

    scal(Fst, 0.25);
  }


//====================================================================
  void Fopr_CloverTerm_eo::trSigmaInv(Field_G& tr_sigma_inv, const int mu, const int nu)
  {
    int      nex_finv = m_fee_inv->nex();
    Vec_SU_N v;
    Field_F  sigma_inv(m_Nvol, nex_finv);

    assert(tr_sigma_inv.nvol() == m_Nvol);
    assert(tr_sigma_inv.nex() == 1);

    {
      Field_F sigma_eo_inv(m_Nvol2, nex_finv);
      mult_isigma(sigma_eo_inv, *m_fee_inv, mu, nu);
      m_idx.reverseField(sigma_inv, sigma_eo_inv, 0);
      mult_isigma(sigma_eo_inv, *m_foo_inv, mu, nu);
      m_idx.reverseField(sigma_inv, sigma_eo_inv, 1);
    }

    for (int isite = 0; isite < m_Nvol; ++isite) {
      for (int ispin = 0; ispin < m_Nd; ++ispin) {
        for (int icolor = 0; icolor < m_Nc; ++icolor) {
          v = sigma_inv.vec(ispin, isite, icolor + m_Nc * ispin);
          for (int jcolor = 0; jcolor < m_Nc; ++jcolor) {
            int cc = icolor + m_Nc * jcolor;
            tr_sigma_inv.set_r(cc, isite, 0, v.r(jcolor));
            tr_sigma_inv.set_i(cc, isite, 0, v.i(jcolor));
          }
        }
      }
    }
  }


//====================================================================
  double Fopr_CloverTerm_eo::flop_count()
  {
    // The following counting explicitly depends on the implementation
    // and to be recalculated when the code is modified.
    // Present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

    int    Lvol      = CommonParameters::Lvol();
    double flop_site =
      static_cast<double>(8 * m_Nc * m_Nc * m_Nd * m_Nd);

    double flop = flop_site * static_cast<double>(Lvol / 2);

    return flop;
  }


//====================================================================
}
//============================================================END=====

#endif
