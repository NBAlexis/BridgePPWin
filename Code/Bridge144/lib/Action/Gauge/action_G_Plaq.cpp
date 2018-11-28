#include "BridgeLib_Private.h"

/*!
        @file    $Id:: action_G_Plaq.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "action_G_Plaq.h"


#ifdef USE_FACTORY
namespace {
  Action *create_object()
  {
    return new Action_G_Plaq();
  }


  bool init = Action::Factory::Register("Action_G_Plaq", create_object);
}
#endif



const std::string Action_G_Plaq::class_name = "Action_G_Plaq";

//====================================================================
void Action_G_Plaq::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double beta;

  int err = 0;
  err += params.fetch_double("beta", beta);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(beta);

  //- post-process
  m_force_G->set_parameters(params);
}


//====================================================================
void Action_G_Plaq::set_parameters(double beta)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  beta = %8.4f\n", beta);

  //- range check
  // NB. beta == 0 is allowed.

  //- store values
  m_beta = beta;
}


//====================================================================
double Action_G_Plaq::langevin(RandomNumbers *rand)
{
  double H_U = calcH(); // calculate action H_U=beta*(1-Plaq)*Lvol*6 (SA)

  return H_U;
}


//====================================================================
double Action_G_Plaq::calcH()
{
  int Lvol = CommonParameters::Lvol();
  int Ndim = CommonParameters::Ndim();

  int Ndim2  = Ndim * (Ndim - 1) / 2;
  int size_U = Lvol * Ndim2;

  double plaq = m_staple.plaquette(*m_U);                    // calculate plaquette (SA)
  double H_U  = m_beta * (1.0 - plaq) * Lvol * Ndim2;        // action (SA)

  vout.general(m_vl, "H_Gplaq    = %18.8f\n", H_U);          // total action (SA)
  vout.general(m_vl, "H_G/dof    = %18.8f\n", H_U / size_U); // action per dof (SA)

  return H_U;
}


//====================================================================
void Action_G_Plaq::force(Field& force)
{
  //- check of argument type
  assert(force.nin() == m_U->nin());
  assert(force.nvol() == m_U->nvol());
  assert(force.nex() == m_U->nex());

  vout.general(m_vl, "  %s:  %s\n", class_name.c_str(), m_label.c_str());

  force.set(0.0);

  m_force_G->force_core(force, m_U);
}


//====================================================================
//============================================================END=====
