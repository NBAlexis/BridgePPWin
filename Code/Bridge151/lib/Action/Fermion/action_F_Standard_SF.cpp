/*!
        @file    action_F_Standard_SF.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "action_F_Standard_SF.h"

const std::string Action_F_Standard_SF::class_name = "Action_F_Standard_SF";

//====================================================================
void Action_F_Standard_SF::set_parameters(const Parameters& params)
{
  const Bridge::VerboseLevel vl = params.get_VerboseLevel();

  set_parameter_verboselevel(vl);
}


//====================================================================
void Action_F_Standard_SF::set_parameters()
{
  vout.general(m_vl, "%s:\n", class_name.c_str());

  const int         Niter           = 100;
  const int         Nrestart        = 40;
  const double      Stop_cond       = 1.0e-24;
  const std::string str_solver_type = "CG";

  m_solver = Solver::New(str_solver_type, m_fopr);
  m_solver->set_parameters(Niter, Nrestart, Stop_cond);
}


//====================================================================
double Action_F_Standard_SF::langevin(RandomNumbers *rand)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr->field_nin();
  const int NvolF    = m_fopr->field_nvol();
  const int NexF     = m_fopr->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  assert(NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  m_fopr->set_config(m_U);
  m_fopr->set_mode("Ddag");
  m_fopr->mult(m_psf, xi);

  //  set_boundary_zero(xi);
  Field_F_SF setzero;
  setzero.set_boundary_zero(xi);

  const double xi2   = xi.norm();
  const double H_psf = xi2 * xi2;

  vout.general(m_vl, "    H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
double Action_F_Standard_SF::calcH()
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr->field_nin();
  const int NvolF    = m_fopr->field_nvol();
  const int NexF     = m_fopr->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  m_fopr->set_config(m_U);
  m_fopr->set_mode("DdagD");

  Field  v1(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_solver->solve(v1, m_psf, Nconv, diff);

  vout.general(m_vl, "    Nconv = %d  diff  = %.8e\n", Nconv, diff);

  Field_F_SF setzero;
  setzero.set_boundary_zero(v1);
  setzero.set_boundary_zero(m_psf);

  const double H_psf = dot(v1, m_psf);

  vout.general(m_vl, "    H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
void Action_F_Standard_SF::force(Field& force)
{
  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();
  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  const int NinF  = m_fopr->field_nin();
  const int NvolF = m_fopr->field_nvol();
  const int NexF  = m_fopr->field_nex();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  //- fermion inversion for smeared gauge field
  m_fopr->set_config(m_U);
  m_fopr->set_mode("DdagD");

  Field  eta(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_solver->solve(eta, m_psf, Nconv, diff);

  vout.general(m_vl, "    Solver: Nconv = %6d  diff  = %12.6e\n", Nconv, diff);

  //- force of smeared fermion operator
  m_fopr_force->set_config(m_U);

  Field force_org(Nin, Nvol, Nex);
  m_fopr_force->force_core(force_org, eta);

  Field_G_SF Fboundary(force_org);
  Fboundary.set_boundary_spatial_link_zero();
  force = (Field)Fboundary;

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "    Fstandard_ave = %12.6f  Fstandard_max = %12.6f  Fstandard_dev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//====================================================================

/*!
  Set the boundary field to zero: \f$xi(t=0,\vec{x})=0\f$
 */

/*
void Action_F_Standard_SF::set_boundary_zero(Field& f){
  if(comm->ipe(3)==0){
    for(int site = 0; site < Svol; ++site){
      for(int s = 0; s < m_Nd; ++s){
        for(int cc = 0; cc < m_Nc2; ++cc){
          f.set(cc+m_Nc2*s, site, 0, 0.0);
        }
      }
    }
  }
}
*/

//====================================================================
//============================================================END=====
