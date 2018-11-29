#include "BridgeLib_Private.h"

/*!
@file    $Id:: action_F_Rational.cpp #$

@brief

@author  Hideo Matsufuru (matsufuru)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 00:49:46 #$

@version $LastChangedRevision: 1561 $
*/

#include "action_F_Rational.h"



const std::string Action_F_Rational::class_name = "Action_F_Rational";

//====================================================================
void Action_F_Rational::set_parameters(const Parameters& params)
{
    const string str_vlevel = params.get_string("verbose_level");

    m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
void Action_F_Rational::setup()
{
    //what are you doing?
    //int Nc   = CommonParameters::Nc();
    //int Nvol = CommonParameters::Nvol();
    //int Ndim = CommonParameters::Ndim();
    //int NinG = 2 * Nc * Nc;
}


//====================================================================
double Action_F_Rational::langevin(RandomNumbers *rand)
{
    int NinF = m_fopr_langev->field_nin();
    int NvolF = m_fopr_langev->field_nvol();
    int NexF = m_fopr_langev->field_nex();
    int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

    m_psf.reset(NinF, NvolF, NexF);

    vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

    Field xi(NinF, NvolF, NexF);
    rand->gauss_lex_global(xi);

    m_fopr_langev->set_config(m_U);
    m_fopr_langev->mult(m_psf, xi);

    double xi2 = xi.norm();
    double H_psf = xi2 * xi2;

    vout.general(m_vl, "    H_Frational  = %18.8f\n", H_psf);
    vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

    return H_psf;
}


//====================================================================
double Action_F_Rational::calcH()
{
    int NinF = m_fopr_H->field_nin();
    int NvolF = m_fopr_H->field_nvol();
    int NexF = m_fopr_H->field_nex();
    int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

    Field v1(NinF, NvolF, NexF);

    vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

    m_fopr_H->set_config(m_U);
    m_fopr_H->mult(v1, m_psf);

    double H_psf = dot(v1, m_psf);

    vout.general(m_vl, "    H_Frational  = %18.8f\n", H_psf);
    vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

    return H_psf;
}


//====================================================================
void Action_F_Rational::force(Field& force)
{
    //int Nin = m_U->nin();
    //int Nvol = m_U->nvol();
    //int Nex = m_U->nex();

    assert(force.nin() == m_U->nin());
    assert(force.nvol() == m_U->nvol());
    assert(force.nex() == m_U->nex());

    vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

    m_fopr_force_MD->set_config(m_U);
    m_fopr_force_MD->force_core(force, m_psf);

    double Fave, Fmax, Fdev;
    force.stat(Fave, Fmax, Fdev);
    vout.general(m_vl, "    Frational_ave = %12.6f  Frational_max = %12.6f  Frational_dev = %12.6f\n",
        Fave, Fmax, Fdev);
}


//====================================================================
//============================================================END=====
