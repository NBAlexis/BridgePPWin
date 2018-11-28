#include "BridgeLib_Private.h"

/*!
@file    $Id:: action_F_Ratio_eo.cpp #$

@brief

@author  Yusuke Namekawa (namekawa)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 00:49:46 #$

@version $LastChangedRevision: 1561 $
*/

#include "action_F_Ratio_eo.h"

const std::string Action_F_Ratio_eo::class_name = "Action_F_Ratio_eo";

//====================================================================
void Action_F_Ratio_eo::set_parameters(const Parameters& params)
{
    Bridge::VerboseLevel vl = params.get_VerboseLevel();

    set_parameter_verboselevel(vl);
}


//====================================================================
void Action_F_Ratio_eo::set_parameters()
{
    //int Nc   = CommonParameters::Nc();
    //int Nvol = CommonParameters::Nvol();
    //int Ndim = CommonParameters::Ndim();
    //int NinG = 2 * Nc * Nc;

    vout.general(m_vl, "%s:\n", class_name.c_str());
}


//====================================================================
void Action_F_Ratio_eo::set_config(Field *U)
{
    m_U = U;

    //- NB. only solver part is even-odd preconditioned.
    m_fopr_prec->set_config(U);
    m_fopr_prec_force->set_config(U);
    m_fopr->set_config(U);
    m_fopr_force->set_config(U);
}


//====================================================================
double Action_F_Ratio_eo::langevin(RandomNumbers *rand)
{
    int Nvol = CommonParameters::Nvol();
    //int Ndim = CommonParameters::Ndim();

    int NinF = m_fopr_prec->field_nin();
    int NvolF = m_fopr_prec->field_nvol();
    int NexF = m_fopr_prec->field_nex();
    int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

    assert(NvolF == Nvol);
    m_psf.reset(NinF, NvolF, NexF);

    vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

    Field xi(NinF, NvolF, NexF);
    rand->gauss_lex_global(xi);

    m_fopr_prec->set_config(m_U);
    m_fopr->set_config(m_U);

    Field v1(NinF, NvolF, NexF), v2(NinF, NvolF, NexF);

    m_fopr->set_mode("H");
    m_fopr->mult_dag(v2, xi);

    m_fopr_prec->set_mode("H");
    m_fopr_prec->mult_dag(v1, v2);

    int    Nconv;
    double diff;

    m_fprop_H_prec->set_config(m_U);
    m_fprop_H_prec->invert_DdagD(m_psf, v1, Nconv, diff);
    vout.general(m_vl, "    Nconv = %d  diff  = %.8e\n", Nconv, diff);

    double xi2 = xi.norm();
    double H_psf = xi2 * xi2;

    vout.general(m_vl, "    H_Fratio     = %18.8f\n", H_psf);
    vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

    return H_psf;
}


//====================================================================
double Action_F_Ratio_eo::calcH()
{
    //int Nvol = CommonParameters::Nvol();
    //int Ndim = CommonParameters::Ndim();

    int NinF = m_fopr_prec->field_nin();
    int NvolF = m_fopr_prec->field_nvol();
    int NexF = m_fopr_prec->field_nex();
    int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

    Field v1(NinF, NvolF, NexF), v2(NinF, NvolF, NexF);

    vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

    m_fopr_prec->set_config(m_U);
    m_fopr->set_config(m_U);

    m_fopr_prec->set_mode("H");
    m_fopr_prec->mult(v1, m_psf);

    int    Nconv;
    double diff;

    m_fprop_H->set_config(m_U);
    m_fprop_H->invert_DdagD(v2, v1, Nconv, diff);
    vout.general(m_vl, "    Nconv = %d  diff  = %.8e\n", Nconv, diff);

    double H_psf = dot(v1, v2);

    vout.general(m_vl, "    H_Fratio     = %18.8f\n", H_psf);
    vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

    return H_psf;
}


//====================================================================
void Action_F_Ratio_eo::force(Field& force)
{
    int Nin = m_U->nin();
    int Nvol = m_U->nvol();
    int Nex = m_U->nex();
    //int Nc = CommonParameters::Nc();
    //int Ndim = CommonParameters::Ndim();

    assert(force.nin() == Nin);
    assert(force.nvol() == Nvol);
    assert(force.nex() == Nex);

    int   NinF = m_fopr_prec->field_nin();
    int   NvolF = m_fopr_prec->field_nvol();
    int   NexF = m_fopr_prec->field_nex();
    Field eta(NinF, NvolF, NexF);

    vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

    m_fopr_prec->set_config(m_U);
    m_fopr->set_config(m_U);
    m_fopr_prec_force->set_config(m_U);
    m_fopr_force->set_config(m_U);

    Field v1(NinF, NvolF, NexF), v2(NinF, NvolF, NexF);
    Field force_tmp(Nin, Nvol, Nex);

    m_fopr_prec->set_mode("H");
    m_fopr_prec->mult(v1, m_psf);

    int    Nconv;
    double diff;

    m_fprop_MD->set_config(m_U);
    m_fprop_MD->invert_DdagD(v2, v1, Nconv, diff);
    vout.general(m_vl, "    Solver: Nconv = %6d  diff  = %12.6e\n", Nconv, diff);

    m_fopr_force->force_core(force, v2);

    m_fopr_prec_force->set_mode("H");
    m_fopr_prec_force->force_core1(force_tmp, v2, m_psf);
    axpy(force, -1.0, force_tmp);

    m_fopr_prec_force->set_mode("Hdag");
    m_fopr_prec_force->force_core1(force_tmp, m_psf, v2);
    axpy(force, -1.0, force_tmp);

    double Fave, Fmax, Fdev;
    force.stat(Fave, Fmax, Fdev);
    vout.general(m_vl, "    Fratio_ave = %12.6f  Fratio_max = %12.6f  Fratio_dev = %12.6f\n",
        Fave, Fmax, Fdev);
}


//====================================================================
//============================================================END=====
