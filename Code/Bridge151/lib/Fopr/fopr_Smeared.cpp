/*!
        @file    fopr_Smeared.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "fopr_Smeared.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Smeared::register_factory();
}
#endif

const std::string Fopr_Smeared::class_name = "Fopr_Smeared";

//====================================================================
void Fopr_Smeared::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
void Fopr_Smeared::set_config(Field *U)
{
  m_dr_smear->set_config(U);

  const int Nsmear = m_dr_smear->get_Nsmear();
  Field     *Uptr  = m_dr_smear->getptr_smearedConfig(Nsmear);

  m_fopr->set_config(Uptr);
}


//====================================================================
//============================================================END=====
