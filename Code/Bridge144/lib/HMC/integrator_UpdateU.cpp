#include "BridgeLib_Private.h"

/*!
        @file    $Id:: integrator_UpdateU.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "integrator_UpdateU.h"



const std::string Integrator_UpdateU::class_name = "Integrator_UpdateU";

//====================================================================
void Integrator_UpdateU::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int Nprec;

  int err = 0;
  err += params.fetch_int("order_of_exp_iP", Nprec);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Nprec);
}


//====================================================================
void Integrator_UpdateU::set_parameters(int Nprec)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Nprec    = %d\n", Nprec);

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(Nprec);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Nprec = Nprec;
}


//====================================================================
void Integrator_UpdateU::set_parameter_Nprec(const int Nprec)
{
  m_Nprec = Nprec;
}


//====================================================================
void Integrator_UpdateU::evolve(const double step_size, Field_G& iP, Field_G& U)
{
  //- alias
  Field_G& W = U;

  mult_exp_Field_G(W, step_size, iP, U, m_Nprec);

  notify_update();
}


//====================================================================
void Integrator_UpdateU::notify_update()
{
  for (size_t i = 0, n = m_integs.size(); i < n; ++i) {
    m_integs[i]->invalidate_cache();
  }
  for (size_t i = 0, n = m_director.size(); i < n; ++i) {
    m_director[i]->notify_linkv();
  }
}


//====================================================================
//============================================================END=====
