#include "BridgeLib_Private.h"

/*!
        @file    $Id:: action_F_Standard_eo.cpp #$

        @brief

        @author  UEDA, Satoru  (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "action_F_Standard_eo.h"

const std::string Action_F_Standard_eo::class_name = "Action_F_Standard_eo";

//====================================================================
void Action_F_Standard_eo::set_parameters(const Parameters& params)
{
  Bridge::VerboseLevel vl = params.get_VerboseLevel();

  set_parameter_verboselevel(vl);
}


//====================================================================
void Action_F_Standard_eo::set_parameters()
{
  int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  int NinG = 2 * Nc * Nc;

  vout.detailed(m_vl, "%s:\n", class_name.c_str());
}


//====================================================================
void Action_F_Standard_eo::set_config(Field *U)
{
  m_U = U;

  //- NB. only solver part is even-odd preconditioned.
  m_fopr->set_config(U);
  m_fopr_force->set_config(U);
}


//====================================================================
double Action_F_Standard_eo::langevin(RandomNumbers *rand)
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF     = m_fopr->field_nin();
  int NvolF    = m_fopr->field_nvol();
  int NexF     = m_fopr->field_nex();
  int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  assert(NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  m_fopr->set_config(m_U);
  m_fopr->set_mode("Ddag");

  m_fopr->mult(m_psf, xi);

  double xi2   = xi.norm();
  double H_psf = xi2 * xi2;

  vout.general(m_vl, "    H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
double Action_F_Standard_eo::calcH()
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF     = m_fopr->field_nin();
  int NvolF    = m_fopr->field_nvol();
  int NexF     = m_fopr->field_nex();
  int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  Field v1(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  int    Nconv;
  double diff;

  m_fprop_H->set_config(m_U);
  m_fprop_H->invert_DdagD(v1, m_psf, Nconv, diff);

  double H_psf = dot(v1, m_psf);

  vout.general(m_vl, "    H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
void Action_F_Standard_eo::force(Field& force)
{
  int Nin  = m_U->nin();
  int Nvol = m_U->nvol();
  int Nex  = m_U->nex();
  int Nc   = CommonParameters::Nc();
  int Ndim = CommonParameters::Ndim();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  int   NinF  = m_fopr->field_nin();
  int   NvolF = m_fopr->field_nvol();
  int   NexF  = m_fopr->field_nex();
  Field eta(NinF, NvolF, NexF);

  vout.detailed(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  int    Nconv;
  double diff;

  m_fprop_MD->set_config(m_U);
  m_fprop_MD->invert_DdagD(eta, m_psf, Nconv, diff);

  // force of smeared fermion operator
  m_fopr_force->set_config(m_U);
  m_fopr_force->force_core(force, eta);

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl,
               "    Fstandard_ave = %12.6f  Fstandard_max = %12.6f  Fstandard_dev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//====================================================================
//============================================================END=====
