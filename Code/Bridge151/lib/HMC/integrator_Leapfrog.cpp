/*!
        @file    integrator_Leapfrog.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "integrator_Leapfrog.h"

const std::string Integrator_Leapfrog::class_name = "Integrator_Leapfrog";

//====================================================================
void Integrator_Leapfrog::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int level;
  int Nstep;

  int err = 0;
  err += params.fetch_int("level", level);
  err += params.fetch_int("number_of_steps", Nstep);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(level, Nstep);
}


//====================================================================
void Integrator_Leapfrog::set_parameters(const int level, const int Nstep)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Level: %4d\n", level);
  vout.general(m_vl, "  Nstep    = %4d\n", Nstep);

  //- range check
  // NB. level,Estep,Nstep == 0 is allowed.

  //- store values
  m_Nstep = Nstep;
  m_level = level;
}


//====================================================================
void Integrator_Leapfrog::set_parameter_level(const int level)
{
  m_level = level;
}


//====================================================================
void Integrator_Leapfrog::set_parameter_Nstep(const int Nstep)
{
  m_Nstep = Nstep;
}


//====================================================================
void Integrator_Leapfrog::set_parameter_Nsteps(const std::vector<int>& Nsteps)
{
  if (Nsteps.size() > 0) {
    set_parameter_Nstep(Nsteps[0]);

    // transfer to lower levels
    if (m_update_U && (Nsteps.size() > 1)) {
      std::vector<int> next_steps(Nsteps.begin() + 1, Nsteps.end());
      m_update_U->set_parameter_Nsteps(next_steps);
    }
  }
}


//====================================================================
void Integrator_Leapfrog::evolve(const double step_size, Field_G& iP, Field_G& U)
{
  const int Nin  = U.nin();
  const int Nvol = U.nvol();
  const int Nex  = U.nex();
  const int Nc   = CommonParameters::Nc();

  vout.general(m_vl, "Integration level-%d start.\n", m_level);

  if (m_Nstep > 0) {
    double estep = step_size / m_Nstep;

    // initial half step
    m_update_p->evolve(estep * 0.5, iP, U);

    for (int istep = 1; istep <= m_Nstep; ++istep) {
      vout.general(m_vl, "istep = %d\n", istep);

      m_update_U->evolve(estep, iP, U);

      if (istep < m_Nstep) {
        m_update_p->evolve(estep, iP, U);
      }
    }

    // last half step
    m_update_p->evolve(estep * 0.5, iP, U);
  } else {
    vout.general(m_vl, "Nstep is zero. skip.\n");
  }

  vout.general(m_vl, "Integration level-%d finished.\n", m_level);
}


//====================================================================
//============================================================END=====
