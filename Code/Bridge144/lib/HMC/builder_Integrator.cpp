/*!
        @file    $Id:: builder_Integrator.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "builder_Integrator.h"

#include "integrator_UpdateU.h"
#include "integrator_UpdateP.h"
#include "integrator_Leapfrog.h"
#include "integrator_Omelyan.h"



const std::string Builder_Integrator::class_name = "Builder_Integrator";

//====================================================================
Builder_Integrator::Builder_Integrator(const ActionList& action_list, std::vector<Director *> director)
  : m_vl(CommonParameters::Vlevel()),
    m_Nprec(0),
    m_lambda_Omelyan(0.0),
    m_str_integrator_type(""),
    m_actions(action_list),
    m_director(director)
{
  m_str_integrator_type = action_list.get_integrator_type(0);

  m_Nprec          = 0;
  m_lambda_Omelyan = 0.0;
}


//====================================================================
void Builder_Integrator::set_parameters(const Parameters& params)
{
  const std::string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  std::string str_integrator_type;

  std::vector<int> Nstep;

  int    Nprec;
  double lambda_Omelyan;

  int err = 0;

  if (m_str_integrator_type == "") {
    err += params.fetch_string("integrator", str_integrator_type);
  } else {
    str_integrator_type = m_str_integrator_type;
  }

  err += params.fetch_int_vector("number_of_steps", Nstep);

  err += params.fetch_int("order_of_exp_iP", Nprec);

  if ((m_str_integrator_type == "") || (m_str_integrator_type == "Omelyan")) {
    err += params.fetch_double("lambda_Omelyan", lambda_Omelyan);
  } else {
    lambda_Omelyan = 0.0;
  }

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(str_integrator_type, Nstep, Nprec, lambda_Omelyan);
}


//====================================================================
void Builder_Integrator::set_parameters(std::string str_integrator_type,
                                        const std::vector<int>& Nstep,
                                        const int Nprec,
                                        const double lambda_Omelyan)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Integrator = %s\n", str_integrator_type.c_str());
  for (int lv = 0; lv < Nstep.size(); ++lv) {
    vout.general(m_vl, "   level     = %2d:  Nstep = %4d\n", lv, Nstep[lv]);
  }
  vout.general(m_vl, "  Nprec      = %4d\n", Nprec);

  //- range check
  int err = 0;
  err += ParameterCheck::non_NULL(str_integrator_type);
  // NB. Nstep == 0 is allowed.
  err += ParameterCheck::non_zero(Nprec);
  // NB. lambda_Omelyan == 0 is allowed.

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_str_integrator_type = str_integrator_type;

  m_Nstep = Nstep;

  m_Nprec          = Nprec;
  m_lambda_Omelyan = lambda_Omelyan;
}


//====================================================================
Integrator *Builder_Integrator::build()
{
  //- select leapfrog or omelyan
  Integrator *integrator;

  if (m_str_integrator_type == "Leapfrog") {
    integrator = build_leapfrog();
  } else if (m_str_integrator_type == "Omelyan") {
    integrator = build_omelyan();
  } else {
    vout.crucial("Error at %s::build : unsupported integrator type \"%s\".\n", class_name.c_str(), m_str_integrator_type.c_str());
    exit(EXIT_FAILURE);
  }

  return integrator;
}


//====================================================================
Integrator *Builder_Integrator::build_leapfrog()
{
  // #####  Following setup is for PUP order.

  // lowest level: update_U
  Integrator_UpdateU *update_U = new Integrator_UpdateU(m_director);

  m_integs.push_back(update_U);

  update_U->set_parameters(m_Nprec);

  Integrator *integrator = update_U;

  for (int lv = m_Nstep.size() - 1; lv >= 0; --lv) {
    Integrator_UpdateP *update_p = new Integrator_UpdateP(m_actions.get_actions(lv));
    m_integs.push_back(update_p);

    // append notify client
    update_U->append_notify(update_p);

    Integrator_Leapfrog *integ = new Integrator_Leapfrog(update_p, integrator);
    m_integs.push_back(integ);

    integ->set_parameters(lv, m_Nstep[lv]);

    // forward to upper level integrator
    integrator = integ;
  }

  return integrator;
}


//====================================================================
Integrator *Builder_Integrator::build_omelyan()
{
  // #####  Following setup is for PUP order.

  // lowest level: update_U
  Integrator_UpdateU *update_U = new Integrator_UpdateU(m_director);

  m_integs.push_back(update_U);

  update_U->set_parameters(m_Nprec);

  Integrator *integrator = update_U;

  for (int lv = m_Nstep.size() - 1; lv >= 0; --lv) {
    Integrator_UpdateP *update_p = new Integrator_UpdateP(m_actions.get_actions(lv));
    m_integs.push_back(update_p);

    // append notify client
    update_U->append_notify(update_p);

    Integrator_Omelyan *integ = new Integrator_Omelyan(update_p, integrator);
    m_integs.push_back(integ);

    integ->set_parameters(lv, m_Nstep[lv], m_lambda_Omelyan);

    // forward to upper level integrator
    integrator = integ;
  }

  return integrator;
}


//====================================================================
void Builder_Integrator::tidyup()
{
  for (size_t i = 0, n = m_integs.size(); i < n; ++i) {
    if (m_integs[i]) delete m_integs[i];
  }

  m_integs.resize(0);
}


//====================================================================
//============================================================END=====
