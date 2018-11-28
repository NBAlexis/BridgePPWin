#include "BridgeLib_Private.h"
#if USE_ORG

/*!
        @file    $Id:: fopr_CloverTerm.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "Fopr/fopr_CloverTerm.h"

namespace Org {
//====================================================================

  const std::string Fopr_CloverTerm::class_name = "Org::Fopr_CloverTerm";

//====================================================================
  void Fopr_CloverTerm::set_parameters(const Parameters& params)
  {
    const std::string str_vlevel = params.get_string("verbose_level");

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
  void Fopr_CloverTerm::set_parameters(double kappa, double cSW,
                                       std::vector<int> bc)
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

    // m_boundary.resize(m_Ndim);  // already resized in init.
    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
    }
  }


//====================================================================
  void Fopr_CloverTerm::set_config(Field *U)
  {
    m_U = (Field_G *)U;
    set_csw();
  }


//====================================================================
  void Fopr_CloverTerm::init(std::string repr)
  {
    m_Nvol = CommonParameters::Nvol();
    m_Ndim = CommonParameters::Ndim();
    m_Nc   = CommonParameters::Nc();
    m_Nd   = CommonParameters::Nd();
    m_NinF = 2 * m_Nc * m_Nd;

    m_U = 0;

    m_repr = repr;

    m_boundary.resize(m_Ndim);
    m_SG.resize(m_Ndim * m_Ndim);

    GammaMatrixSet *gmset = GammaMatrixSet::New(m_repr);

    m_GM5 = gmset->get_GM(gmset->GAMMA5);

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

    delete gmset;
  }


//====================================================================
  void Fopr_CloverTerm::tidyup()
  {
    // nothing to do.
  }


//====================================================================
  void Fopr_CloverTerm::mult_gm5(Field& v, const Field& f)
  {
    assert(v.nvol() == f.nvol());
    assert(v.nex() == f.nex());
    assert(v.nin() == f.nin());

    Field_F vt(f.nvol(), f.nex());

    mult_GM(vt, m_GM5, (Field_F)f);
    v = (Field)vt;
  }


//====================================================================
  void Fopr_CloverTerm::mult_isigma(Field_F& v, const Field_F& w,
                                    const int mu, const int nu)
  {
    assert(mu != nu);
    mult_iGM(v, m_SG[sg_index(mu, nu)], w);
  }


//====================================================================
  void Fopr_CloverTerm::mult_sigmaF(Field& v, const Field& f)
  {
    mult_csw(v, f);
  }


//====================================================================
  void Fopr_CloverTerm::mult_csw(Field& v, const Field& w)
  {
    int Nc   = CommonParameters::Nc();
    int Nd   = CommonParameters::Nd();
    int Nvol = w.nvol();

    Field_F wt(Nvol, 1), vt(Nvol, 1);

    vt.set(0.0);

    mult_iGM(wt, m_SG[sg_index(1, 2)], (Field_F)w);
    multadd_Field_Gn(vt, 0, m_Bx, 0, wt, 0, 1.0);

    mult_iGM(wt, m_SG[sg_index(2, 0)], (Field_F)w);
    multadd_Field_Gn(vt, 0, m_By, 0, wt, 0, 1.0);

    mult_iGM(wt, m_SG[sg_index(0, 1)], (Field_F)w);
    multadd_Field_Gn(vt, 0, m_Bz, 0, wt, 0, 1.0);

    mult_iGM(wt, m_SG[sg_index(3, 0)], (Field_F)w);
    multadd_Field_Gn(vt, 0, m_Ex, 0, wt, 0, 1.0);

    mult_iGM(wt, m_SG[sg_index(3, 1)], (Field_F)w);
    multadd_Field_Gn(vt, 0, m_Ey, 0, wt, 0, 1.0);

    mult_iGM(wt, m_SG[sg_index(3, 2)], (Field_F)w);
    multadd_Field_Gn(vt, 0, m_Ez, 0, wt, 0, 1.0);

    scal(vt, m_kappa * m_cSW);

    v = (Field)vt;
  }


//====================================================================
  void Fopr_CloverTerm::set_csw()
  {
    set_fieldstrength(m_Bx, 1, 2);
    set_fieldstrength(m_By, 2, 0);
    set_fieldstrength(m_Bz, 0, 1);
    set_fieldstrength(m_Ex, 3, 0);
    set_fieldstrength(m_Ey, 3, 1);
    set_fieldstrength(m_Ez, 3, 2);
  }


//====================================================================
  void Fopr_CloverTerm::set_fieldstrength(Field_G& Fst,
                                          const int mu, const int nu)
  {
    int Nvol = CommonParameters::Nvol();

    Staple_lex staple;

    Field_G Cup(Nvol, 1), Cdn(Nvol, 1);
    Field_G Umu(Nvol, 1);
    Field_G v(Nvol, 1), v2(Nvol, 1);

    staple.upper(Cup, *m_U, mu, nu);
    staple.lower(Cdn, *m_U, mu, nu);
    Umu.setpart_ex(0, *m_U, mu);

    mult_Field_Gnd(Fst, 0, Umu, 0, Cup, 0);
    multadd_Field_Gnd(Fst, 0, Umu, 0, Cdn, 0, -1.0);

    mult_Field_Gdn(v, 0, Cup, 0, Umu, 0);
    multadd_Field_Gdn(v, 0, Cdn, 0, Umu, 0, -1.0);

    m_shift.forward(v2, v, mu);

    axpy(Fst, 1, v2);

    ah_Field_G(Fst, 0);
    scal(Fst, 0.25);
  }


//====================================================================
  double Fopr_CloverTerm::flop_count()
  {
    // The following counting explicitly depends on the implementation
    // and to be recalculated when the code is modified.
    // Present counting is based on rev.1131. [10 Sep 2014 H.Matsufuru]

    int Lvol = CommonParameters::Lvol();

    double flop_site
      = static_cast<double>(m_Nc * m_Nd * (2 + 12 + 48 * m_Nc));
    double flop = flop_site * static_cast<double>(Lvol);

    return flop;
  }


//====================================================================
}
//============================================================END=====

#endif
