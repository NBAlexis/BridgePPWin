#include "BridgeLib_Private.h"
#if USE_IMP

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

namespace Imp {
#if defined USE_GROUP_SU3
#include "fopr_Wilson_impl_SU3.inc"
#elif defined USE_GROUP_SU2
#include "fopr_Wilson_impl_SU2.inc"
#elif defined USE_GROUP_SU_N
#include "fopr_Wilson_impl_SU_N.inc"
#endif

    //====================================================================

    const std::string Fopr_CloverTerm_eo::class_name = "Imp::Fopr_CloverTerm_eo";

    //====================================================================
    void Fopr_CloverTerm_eo::init(std::string repr)
    {
        m_repr = repr;

        m_Nc = CommonParameters::Nc();
        m_Nd = CommonParameters::Nd();
        m_Ndim = CommonParameters::Ndim();
        m_NinF = 2 * m_Nc * m_Nd;
        m_Nvol = CommonParameters::Nvol();
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
            vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
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
        m_cSW = cSW;
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

        unique_ptr<Solver_CG> solver(new Solver_CG(this));

#if 1
        solver->set_parameters(params_solver);
#else
        const int    Niter = 100;
        const int    Nrestart = 40;
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
                }
                else {
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
    void Fopr_CloverTerm_eo::mult_csw_inv(Field& v,
        const Field& f, const int ieo)
    {
        if (m_repr == "Dirac") {
            mult_csw_inv_dirac(v, f, ieo);
        }
        else if (m_repr == "Chiral") {
            mult_csw_inv_chiral(v, f, ieo);
        }
    }


    //====================================================================
    void Fopr_CloverTerm_eo::mult_csw_inv_dirac(Field& v,
        const Field& f, const int ieo)
    {
        int Nvc = 2 * m_Nc;

        const double *v1_f = f.ptr(0);
        double       *v2_f = v.ptr(0);

        double *csw_inv = 0;

        if (ieo == 0) {
            csw_inv = m_fee_inv->ptr(0);
        }
        else if (ieo == 1) {
            csw_inv = m_foo_inv->ptr(0);

            /*
            } else {
            vout.crucial(m_vl, "Error at %s: wrong parameter, ieo = %d.\n",
            class_name.c_str(), ieo);
            exit(EXIT_FAILURE);
            */
        }

        // threadding applied.
        int Nthread = ThreadManager_OpenMP::get_num_threads();
        int i_thread = ThreadManager_OpenMP::get_thread_id();
        int is = m_Nvol2 * i_thread / Nthread;
        int ns = m_Nvol2 * (i_thread + 1) / Nthread;

        int Nd2 = m_Nd / 2;
        for (int site = is; site < ns; ++site) {
            for (int icd = 0; icd < m_Nc * Nd2; ++icd) {
                int iv2 = 2 * icd + m_NinF * site;
                v2_f[iv2] = 0.0;
                v2_f[iv2 + 1] = 0.0;
                for (int jd = 0; jd < m_Nd; ++jd) {
                    int jcd = Nvc * jd;
                    int iv = jcd + m_NinF * site;
                    int ig = jcd + m_NinF * (site + m_Nvol2 * icd);
                    v2_f[iv2] += mult_uv_r(&csw_inv[ig], &v1_f[iv], m_Nc);
                    v2_f[iv2 + 1] += mult_uv_i(&csw_inv[ig], &v1_f[iv], m_Nc);
                }
            }

            for (int icd = 0; icd < m_Nc * Nd2; ++icd) {
                int iv2 = 2 * (icd + m_Nc * Nd2) + m_NinF * site;
                v2_f[iv2] = 0.0;
                v2_f[iv2 + 1] = 0.0;
                for (int jd = 0; jd < m_Nd; ++jd) {
                    int jd2 = (jd + Nd2) % m_Nd;
                    int iv = Nvc * jd + m_NinF * site;
                    int ig = Nvc * jd2 + m_NinF * (site + m_Nvol2 * icd);
                    v2_f[iv2] += mult_uv_r(&csw_inv[ig], &v1_f[iv], m_Nc);
                    v2_f[iv2 + 1] += mult_uv_i(&csw_inv[ig], &v1_f[iv], m_Nc);
                }
            }
        }
#pragma omp barrier
    }


    //====================================================================
    void Fopr_CloverTerm_eo::mult_csw_inv_chiral(Field& v,
        const Field& f, const int ieo)
    {
        int Nvc = 2 * m_Nc;

        const double *v1_f = f.ptr(0);
        double       *v2_f = v.ptr(0);

        double *csw_inv = 0;

        if (ieo == 0) {
            csw_inv = m_fee_inv->ptr(0);
        }
        else if (ieo == 1) {
            csw_inv = m_foo_inv->ptr(0);

            /*
            } else {
            vout.crucial(m_vl, "Error at %s: wrong parameter, ieo = %d.\n",
            class_name.c_str(), ieo);
            exit(EXIT_FAILURE);
            */
        }

        // threadding applied.
        int Nthread = ThreadManager_OpenMP::get_num_threads();
        int i_thread = ThreadManager_OpenMP::get_thread_id();
        int is = m_Nvol2 * i_thread / Nthread;
        int ns = m_Nvol2 * (i_thread + 1) / Nthread;

        for (int site = is; site < ns; ++site) {
            for (int icd = 0; icd < m_Nc * m_Nd / 2; ++icd) {
                int iv2 = 2 * icd + m_NinF * site;
                v2_f[iv2] = 0.0;
                v2_f[iv2 + 1] = 0.0;

                for (int jd = 0; jd < m_Nd / 2; ++jd) {
                    int jcd = Nvc * jd;
                    int iv = jcd + m_NinF * site;
                    int ig = jcd + m_NinF * (site + m_Nvol2 * icd);
                    v2_f[iv2] += mult_uv_r(&csw_inv[ig], &v1_f[iv], m_Nc);
                    v2_f[iv2 + 1] += mult_uv_i(&csw_inv[ig], &v1_f[iv], m_Nc);
                }
            }

            for (int icd = m_Nc * m_Nd / 2; icd < m_Nc * m_Nd; ++icd) {
                int iv2 = 2 * icd + m_NinF * site;
                v2_f[iv2] = 0.0;
                v2_f[iv2 + 1] = 0.0;

                for (int jd = m_Nd / 2; jd < m_Nd; ++jd) {
                    int jcd = Nvc * jd;
                    int iv = jcd + m_NinF * site;
                    int ig = jcd + m_NinF * (site + m_Nvol2 * icd);
                    v2_f[iv2] += mult_uv_r(&csw_inv[ig], &v1_f[iv], m_Nc);
                    v2_f[iv2 + 1] += mult_uv_i(&csw_inv[ig], &v1_f[iv], m_Nc);
                }
            }
        }
#pragma omp barrier
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
                        int cc = jcolor + icolor * m_Nc;
                        int ss1 = jspin + ispin * m_Nd;
                        int ss2 = js2 + ispin * m_Nd;

                        matrix[2 * cs1] = m_T.cmp_r(cc, site, ss1);
                        matrix[2 * cs1 + 1] = m_T.cmp_i(cc, site, ss1);

                        matrix[2 * cs2] = m_T.cmp_r(cc, site, ss2);
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
        if (m_repr == "Dirac") {
            D_dirac(v, f, ieo);
        }
        else if (m_repr == "Chiral") {
            D_chiral(v, f, ieo);
        }
    }


    //====================================================================
    void Fopr_CloverTerm_eo::D_dirac(Field& v, const Field& f, const int ieo)
    {
        //  assert(f.nvol() == m_Nvol2);
        //  assert(f.nex()  == 1);
        //  assert(v.nvol() == m_Nvol2);
        //  assert(v.nex()  == 1);

        const double *fp = f.ptr(0);
        double       *vp = v.ptr(0);
        double       *tp = m_T.ptr(0, m_Nvol2 * ieo, 0);

        int Nthread = ThreadManager_OpenMP::get_num_threads();
        int i_thread = ThreadManager_OpenMP::get_thread_id();
        int is = m_Nvol2 * i_thread / Nthread;
        int ns = m_Nvol2 * (i_thread + 1) / Nthread;

        int Nvc = 2 * m_Nc;
        int Nd2 = m_Nd / 2;
        int NinF = 2 * m_Nc * m_Nd;
        int NinG = 2 * m_Nc * m_Nc;

        for (int site = is; site < ns; ++site) {
            for (int id = 0; id < Nd2; ++id) {
                for (int ic = 0; ic < m_Nc; ++ic) {
                    int icd = ic + m_Nc * id;

                    int iv2 = 2 * icd + NinF * site;
                    vp[iv2] = 0.0;
                    vp[iv2 + 1] = 0.0;
                    for (int jd = 0; jd < m_Nd; ++jd) {
                        int iv = Nvc * jd + NinF * site;
                        int ig = Nvc * ic + NinG * (site + m_Nvol * (id * m_Nd + jd));
                        vp[iv2] += mult_uv_r(&tp[ig], &fp[iv], m_Nc);
                        vp[iv2 + 1] += mult_uv_i(&tp[ig], &fp[iv], m_Nc);
                    }

                    iv2 += Nvc * Nd2;
                    vp[iv2] = 0.0;
                    vp[iv2 + 1] = 0.0;
                    for (int jd = 0; jd < m_Nd; ++jd) {
                        int jd2 = (2 + jd) % m_Nd;
                        int iv = Nvc * jd + NinF * site;
                        int ig = Nvc * ic + NinG * (site + m_Nvol * (id * m_Nd + jd2));
                        vp[iv2] += mult_uv_r(&tp[ig], &fp[iv], m_Nc);
                        vp[iv2 + 1] += mult_uv_i(&tp[ig], &fp[iv], m_Nc);
                    }
                }
            }
        }
#pragma omp barrier
    }


    //====================================================================
    void Fopr_CloverTerm_eo::D_chiral(Field& v, const Field& f, const int ieo)
    {
        const double *fp = f.ptr(0);
        double       *vp = v.ptr(0);
        double       *tp = m_T.ptr(0, m_Nvol2 * ieo, 0);

        int Nthread = ThreadManager_OpenMP::get_num_threads();
        int i_thread = ThreadManager_OpenMP::get_thread_id();
        int is = m_Nvol2 * i_thread / Nthread;
        int ns = m_Nvol2 * (i_thread + 1) / Nthread;

        int Nvc = 2 * m_Nc;
        int Nd2 = m_Nd / 2;
        int NinF = 2 * m_Nc * m_Nd;
        int NinG = 2 * m_Nc * m_Nc;

        for (int site = is; site < ns; ++site) {
            for (int id = 0; id < Nd2; ++id) {
                for (int ic = 0; ic < m_Nc; ++ic) {
                    int icd = ic + m_Nc * id;

                    int iv2 = 2 * icd + NinF * site;
                    vp[iv2] = 0.0;
                    vp[iv2 + 1] = 0.0;
                    for (int jd = 0; jd < Nd2; ++jd) {
                        int iv = Nvc * jd + NinF * site;
                        int ig = Nvc * ic + NinG * (site + m_Nvol * (id * Nd2 + jd));
                        vp[iv2] += mult_uv_r(&tp[ig], &fp[iv], m_Nc);
                        vp[iv2 + 1] += mult_uv_i(&tp[ig], &fp[iv], m_Nc);
                    }

                    iv2 += Nvc * Nd2;
                    vp[iv2] = 0.0;
                    vp[iv2 + 1] = 0.0;
                    for (int jd = 0; jd < Nd2; ++jd) {
                        int iv = Nvc * (Nd2 + jd) + NinF * site;
                        int ig = Nvc * ic + NinG * (site + m_Nvol * (m_Nd + id * Nd2 + jd));
                        vp[iv2] += mult_uv_r(&tp[ig], &fp[iv], m_Nc);
                        vp[iv2 + 1] += mult_uv_i(&tp[ig], &fp[iv], m_Nc);
                    }
                }
            }
        }
#pragma omp barrier
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
        if (m_repr == "Dirac") {
            set_csw_dirac();
        }
        else if (m_repr == "Chiral") {
            set_csw_chiral();
        }
        else {
            vout.crucial(m_vl, "Error at %s: unsupported gamma matrix repr. %s.\n",
                class_name.c_str(), m_repr.c_str());
            exit(EXIT_FAILURE);
        }
    }


    //====================================================================
    void Fopr_CloverTerm_eo::set_csw_dirac()
    {
        // The clover term in the Dirac representation is as spin-space
        // matrix
        //   [ P Q ]
        //   [ Q P ],
        // where P and Q are 2x2 block matrices as
        //   P =  [          iF(1,2)   F(3,1) + iF(2,3) ]
        //        [-F(3,1) + iF(2,3)          - iF(1,2) ]
        // and
        //   Q =  [        - iF(4,3)  -F(4,2) - iF(4,1) ]
        //        [ F(4,2) - iF(4,1)            iF(4,3) ]
        // up to the coefficient.
        // in the following what defined is
        // [ P Q ] = [ T(0) T(1)  T(2) T(3) ]
        //           [ T(4) T(5)  T(6) T(7) ].

        m_T.set(0.0);

        //- sigma23
        Field_G F;
        set_fieldstrength(F, 1, 2);
        F.xI();
        axpy(m_T, 1, 1.0, F, 0);
        axpy(m_T, 4, 1.0, F, 0);

        //- sigma31
        set_fieldstrength(F, 2, 0);
        axpy(m_T, 1, 1.0, F, 0);
        axpy(m_T, 4, -1.0, F, 0);

        //- sigma12
        set_fieldstrength(F, 0, 1);
        F.xI();
        axpy(m_T, 0, 1.0, F, 0);
        axpy(m_T, 5, -1.0, F, 0);

        //- sigma41
        set_fieldstrength(F, 3, 0);
        F.xI();
        axpy(m_T, 3, -1.0, F, 0);
        axpy(m_T, 6, -1.0, F, 0);

        //- sigma42
        set_fieldstrength(F, 3, 1);
        axpy(m_T, 3, -1.0, F, 0);
        axpy(m_T, 6, 1.0, F, 0);

        //- sigma43
        set_fieldstrength(F, 3, 2);
        F.xI();
        axpy(m_T, 2, -1.0, F, 0);
        axpy(m_T, 7, 1.0, F, 0);

        scal(m_T, -m_kappa * m_cSW);

        Field_G Unity(m_Nvol, 1);
        Unity.set_unit();
        axpy(m_T, 0, 1.0, Unity, 0);
        axpy(m_T, 5, 1.0, Unity, 0);
    }


    //====================================================================
    void Fopr_CloverTerm_eo::set_csw_chiral()
    {
        // The clover term in the Dirac representation is
        // as spin-space matrix
        //  [ P+Q  0  ]
        //  [ 0   P-Q ],
        // where P and Q are 2x2 block matrices as
        //        [          iF(1,2) |  F(3,1) + iF(2,3) ]
        //   P =  [ -----------------+------------------ ]
        //        [-F(3,1) + iF(2,3) |         - iF(1,2) ]
        // and
        //        [        - iF(4,3) | -F(4,2) - iF(4,1) ]
        //   Q =  [ -----------------+------------------ ]
        //        [ F(4,2) - iF(4,1) |           iF(4,3) ]
        // up to the coefficient.
        // in the following what defined is
        //        [ T(0) | T(1) ]          [ T(4) | T(5) ]
        //  P+Q = [ -----+----- ]  P - Q = [ -----+----- ]
        //        [ T(2) | T(3) ]          [ T(6) | T(7) ]

        m_T.set(0.0);

        Field_G F;

        //- sigma23
        set_fieldstrength(F, 1, 2);
        F.xI();
        axpy(m_T, 1, 1.0, F, 0);
        axpy(m_T, 2, 1.0, F, 0);
        axpy(m_T, 5, 1.0, F, 0);
        axpy(m_T, 6, 1.0, F, 0);

        //- sigma31
        set_fieldstrength(F, 2, 0);
        axpy(m_T, 1, 1.0, F, 0);
        axpy(m_T, 2, -1.0, F, 0);
        axpy(m_T, 5, 1.0, F, 0);
        axpy(m_T, 6, -1.0, F, 0);

        //- sigma12
        set_fieldstrength(F, 0, 1);
        F.xI();
        axpy(m_T, 0, 1.0, F, 0);
        axpy(m_T, 3, -1.0, F, 0);
        axpy(m_T, 4, 1.0, F, 0);
        axpy(m_T, 7, -1.0, F, 0);

        //- sigma41
        set_fieldstrength(F, 3, 0);
        F.xI();
        axpy(m_T, 1, -1.0, F, 0);
        axpy(m_T, 2, -1.0, F, 0);
        axpy(m_T, 5, 1.0, F, 0);
        axpy(m_T, 6, 1.0, F, 0);

        //- sigma42
        set_fieldstrength(F, 3, 1);
        axpy(m_T, 1, -1.0, F, 0);
        axpy(m_T, 2, 1.0, F, 0);
        axpy(m_T, 5, 1.0, F, 0);
        axpy(m_T, 6, -1.0, F, 0);

        //- sigma43
        set_fieldstrength(F, 3, 2);
        F.xI();
        axpy(m_T, 0, -1.0, F, 0);
        axpy(m_T, 3, 1.0, F, 0);
        axpy(m_T, 4, 1.0, F, 0);
        axpy(m_T, 7, -1.0, F, 0);

        scal(m_T, -m_kappa * m_cSW);

        Field_G Unity(m_Nvol, 1);
        Unity.set_unit();
        axpy(m_T, 0, 1.0, Unity, 0);
        axpy(m_T, 3, 1.0, Unity, 0);
        axpy(m_T, 4, 1.0, Unity, 0);
        axpy(m_T, 7, 1.0, Unity, 0);
    }


    //====================================================================
    void Fopr_CloverTerm_eo::set_fieldstrength(Field_G& Fst,
        const int mu, const int nu)
    {
        // Staple_eo staple;
        unique_ptr<Staple> staple(Staple::New("EvenOdd"));

        Field_G Cup(m_Nvol, 1), Cdn(m_Nvol, 1);
        Field_G Umu(m_Nvol, 1);
        Field_G w(m_Nvol, 1), v(m_Nvol, 1), v2_f(m_Nvol, 1);

        staple->upper(Cup, *m_Ueo, mu, nu);
        staple->lower(Cdn, *m_Ueo, mu, nu);
        Umu.setpart_ex(0, *m_Ueo, mu);

        for (int site = 0; site < m_Nvol; ++site) {
            w.set_mat(site, 0, Umu.mat(site) * Cup.mat_dag(site));
        }

        for (int site = 0; site < m_Nvol; ++site) {
            v2_f.set_mat(site, 0, Umu.mat(site) * Cdn.mat_dag(site));
        }

        axpy(w, -1.0, v2_f);

        for (int site = 0; site < m_Nvol; ++site) {
            v.set_mat(site, 0, Cup.mat_dag(site) * Umu.mat(site));
        }

        for (int site = 0; site < m_Nvol; ++site) {
            v2_f.set_mat(site, 0, Cdn.mat_dag(site) * Umu.mat(site));
        }

        axpy(v, -1.0, v2_f);

        m_shift_eo.forward(v2_f, v, mu);

        axpy(w, 1.0, v2_f);

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

        int    Lvol = CommonParameters::Lvol();
        double flop_site = 0.0;

        if (m_repr == "Dirac") {
            flop_site = static_cast<double>(8 * m_Nc * m_Nc * m_Nd * m_Nd);
        }
        else if (m_repr == "Chiral") {
            flop_site = static_cast<double>(4 * m_Nc * m_Nc * m_Nd * m_Nd);
        }

        double flop = flop_site * static_cast<double>(Lvol / 2);

        return flop;
    }


    //====================================================================
}
//============================================================END=====

#endif
