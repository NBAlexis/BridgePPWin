/*!
        @file    $Id:: integrator_UpdateP.cpp #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "integrator_UpdateP.h"



const std::string Integrator_UpdateP::class_name = "Integrator_UpdateP";

//====================================================================
void Integrator_UpdateP::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int err = 0;

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters();
}


//====================================================================
void Integrator_UpdateP::set_parameters()
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());

  //- range check
  int err = 0;

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
}


//====================================================================
void Integrator_UpdateP::evolve(const double step_size, Field_G& iP, Field_G& U)
{
  if (not is_cache_valid()) {  // recalc force
    m_force.set(0.0);

    Field force_tmp(iP.nin(), iP.nvol(), iP.nex());
    for (unsigned int i = 0, n = m_action.size(); i < n; ++i) {
      m_action[i]->force(force_tmp);
      vout.detailed(m_vl, "updated p by action %d finished.\n", i);
      axpy(m_force, 1.0, force_tmp);
    }

    cache_validated();

    vout.detailed(m_vl, "%s: force updated.\n", class_name.c_str());
  } else {
    vout.general(m_vl, "%s: returns previous force.\n", class_name.c_str());
  }

  axpy(iP, step_size, m_force);
}


//====================================================================
//============================================================END=====
