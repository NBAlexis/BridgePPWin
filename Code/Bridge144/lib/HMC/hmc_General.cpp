#include "BridgeLib_Private.h"

/*!
        @file    $Id:: hmc_General.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "hmc_General.h"



const std::string HMC_General::class_name = "HMC_General";

//====================================================================
//! constructor with action_list, directors, and random number generator
HMC_General::HMC_General(const ActionList& action_list,
                         std::vector<Director *> director,
                         Integrator *integrator,
                         RandomNumbers *rand)
  : m_vl(CommonParameters::Vlevel())
{
  ActionSet action_set = action_list.get_actions();

  m_action.resize(action_set.size());
  for (int i = 0; i < action_set.size(); ++i) {
    m_action[i] = action_set[i];
  }
  m_director.resize(director.size());
  for (int i = 0; i < director.size(); ++i) {
    m_director[i] = director[i];
  }
  m_integrator        = integrator;
  m_rand              = rand;
  m_staple            = new Staple_lex;
  m_Metropolis_test   = 0;
  m_Langevin_P        = new Langevin_Momentum(m_rand);
  m_trajectory_length = 0.0;
}


//====================================================================
HMC_General::HMC_General(const ActionList& action_list,
                         std::vector<Director *> director,
                         Integrator *integrator,
                         unique_ptr<RandomNumbers>& rand)
  : m_vl(CommonParameters::Vlevel())
{
  ActionSet action_set = action_list.get_actions();

  m_action.resize(action_set.size());
  for (int i = 0; i < action_set.size(); ++i) {
    m_action[i] = action_set[i];
  }
  m_director.resize(director.size());
  for (int i = 0; i < director.size(); ++i) {
    m_director[i] = director[i];
  }
  m_integrator        = integrator;
  m_rand              = rand.get();
  m_staple            = new Staple_lex;
  m_Metropolis_test   = 0;
  m_Langevin_P        = new Langevin_Momentum(m_rand);
  m_trajectory_length = 0.0;
}


//====================================================================
//! constructor with action_list when no director is necessary
HMC_General::HMC_General(const ActionList& action_list,
                         Integrator *integrator,
                         RandomNumbers *rand)
  : m_vl(CommonParameters::Vlevel())
{
  ActionSet action_set = action_list.get_actions();

  m_action.resize(action_set.size());
  for (int i = 0; i < action_set.size(); ++i) {
    m_action[i] = action_set[i];
  }
  m_integrator        = integrator;
  m_rand              = rand;
  m_staple            = new Staple_lex;
  m_Metropolis_test   = 0;
  m_Langevin_P        = new Langevin_Momentum(m_rand);
  m_trajectory_length = 0.0;
}


//====================================================================
HMC_General::HMC_General(const ActionList& action_list,
                         Integrator *integrator,
                         unique_ptr<RandomNumbers>& rand)
  : m_vl(CommonParameters::Vlevel())
{
  ActionSet action_set = action_list.get_actions();

  m_action.resize(action_set.size());
  for (int i = 0; i < action_set.size(); ++i) {
    m_action[i] = action_set[i];
  }
  m_integrator        = integrator;
  m_rand              = rand.get();
  m_staple            = new Staple_lex;
  m_Metropolis_test   = 0;
  m_Langevin_P        = new Langevin_Momentum(m_rand);
  m_trajectory_length = 0.0;
}


//====================================================================
//! destructor
HMC_General::~HMC_General()
{
  delete m_Langevin_P;
  delete m_staple;
}


//====================================================================
void HMC_General::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double traj_length;
  int    Metropolis_test;

  int err = 0;
  err += params.fetch_double("trajectory_length", traj_length);
  err += params.fetch_int("Metropolis_test", Metropolis_test);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(traj_length, Metropolis_test);
}


//====================================================================
void HMC_General::set_parameters(double trajectory_length, int Metropolis_test)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Number of actions: %4d\n", m_action.size());
  vout.general(m_vl, "  traj_length     = %8.6f\n", trajectory_length);
  vout.general(m_vl, "  Metropolis_test = %4d\n", Metropolis_test);

  //- range check
  // NB. Metropolis_test == 0 is allowed.

  //- store values
  m_trajectory_length = trajectory_length;
  m_Metropolis_test   = Metropolis_test;
}


//====================================================================
double HMC_General::update(Field_G& Uorg)
{
  int Nc   = CommonParameters::Nc();
  int Lvol = CommonParameters::Lvol();

  Field_G U(Uorg);

  for (int i = 0; i < m_action.size(); ++i) {
    m_action[i]->set_config(&U);
  }

  int Nin  = U.nin();
  int Nvol = U.nvol();
  int Nex  = U.nex();

  Field_G iP(Nvol, Nex);


  vout.general(m_vl, "HMC (general) start.\n");

  // Langevin step
  double H_total0 = langevin(iP, U);
  vout.general(m_vl, "  H_total(init)  = %20.10f\n", H_total0);

  double plaq0 = m_staple->plaquette(U);
  vout.general(m_vl, "  plaq(init)     = %20.12f\n", plaq0);

  // initial Hamiltonian
  //    H_total0 = calc_Hamiltonian(iP,U);

  // molecular dynamical integration
  m_integrator->evolve(m_trajectory_length, iP, U);

  // trial Hamiltonian
  double H_total1 = calc_Hamiltonian(iP, U);
  vout.general(m_vl, "  H_total(trial) = %20.10f\n", H_total1);

  double plaq1 = m_staple->plaquette(U);
  vout.general(m_vl, "  plaq(trial)    = %20.12f\n", plaq1);


  // Metropolis test
  vout.general(m_vl, "Metropolis test.\n");

  double diff_H           = H_total1 - H_total0;
  double exp_minus_diff_H = exp(-diff_H);
  double rand             = -1.0;
  if (m_Metropolis_test != 0) {
    rand = m_rand->get();
  }
  vout.general(m_vl, "  H(diff)       = %14.8f\n", diff_H);
  vout.general(m_vl, "  exp(-diff_H)  = %14.8f\n", exp_minus_diff_H);
  vout.general(m_vl, "  Random number = %14.8f\n", rand);

  if (rand <= exp_minus_diff_H) {  // accepted
    vout.general(m_vl, "  Accepted\n");
    Uorg = U;
  } else {  // rejected
    vout.general(m_vl, "  Rejected\n");
  }


  double plaq_final = m_staple->plaquette(Uorg);
  vout.general(m_vl, "  plaq(final)    = %20.12f\n", plaq_final);

  return plaq_final;
}


//====================================================================
double HMC_General::langevin(Field_G& iP, Field_G& U)
{
  int Nc   = CommonParameters::Nc();
  int Lvol = CommonParameters::Lvol();
  int Ndim = CommonParameters::Ndim();
  int NcA  = Nc * Nc - 1;

  int size_iP = NcA * Lvol * Ndim;

  vout.general(m_vl, "Langevin step.\n");


  // discard caches
  m_integrator->invalidate_cache();

  for (int i = 0; i < m_director.size(); ++i) {
    m_director[i]->notify_linkv();
  }

  // kinetic term
  double H_iP = m_Langevin_P->set_iP(iP);
  vout.general(m_vl, "  Kinetic term:\n");
  vout.general(m_vl, "    H_kin        = %18.8f\n", H_iP);
  vout.general(m_vl, "    H_kin/dof    = %18.8f\n", H_iP / size_iP);


  double H_actions = 0.0;
  for (int i = 0; i < m_action.size(); ++i) {
    H_actions += m_action[i]->langevin(m_rand);
  }

  double H_total = H_iP + H_actions;
  return H_total;
}


//====================================================================
double HMC_General::calc_Hamiltonian(Field_G& iP, Field_G& U)
{
  int Nin  = U.nin();
  int Nvol = U.nvol();
  int Nex  = U.nex();

  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Lvol = CommonParameters::Lvol();
  int NcA  = Nc * Nc - 1;

  int size_iP = NcA * Lvol * Nex;

  // int size_U   = Lvol * 6.0;
  // int size_psf = Lvol * Nc * Nd;


  vout.general(m_vl, "Hamiltonian calculation.\n");

  // kinetic term
  double H_iP = calcH_P(iP);
  vout.general(m_vl, "  Kinetic term:\n");
  vout.general(m_vl, "    H_kin        = %18.8f\n", H_iP);
  vout.general(m_vl, "    H_kin/dof    = %18.8f\n", H_iP / size_iP);

  double H_actions = 0.0;
  for (int i = 0; i < m_action.size(); ++i) {
    H_actions += m_action[i]->calcH();
  }

  double H_total = H_iP + H_actions;
  return H_total;
}


//====================================================================
double HMC_General::calcH_P(Field_G& iP)
{
  double hn   = iP.norm();
  double H_iP = 0.5 * hn * hn;

  return H_iP;
}


//====================================================================
//============================================================END=====
