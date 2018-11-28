/*!
        @file    $Id:: fopr_Smeared_eo.cpp #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_Smeared_eo.h"

const std::string Fopr_Smeared_eo::class_name = "Fopr_Smeared_eo";

#ifdef USE_FACTORY
namespace {
  Fopr *create_object(Fopr *fopr, Director *director)
  {
    return new Fopr_Smeared_eo(dynamic_cast<Fopr_eo *>(fopr), dynamic_cast<Director_Smear *>(director));
  }


  bool init = Fopr::Factory_fopr_director::Register("Smeared_eo", create_object);
}
#endif

//====================================================================
void Fopr_Smeared_eo::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
void Fopr_Smeared_eo::set_config(Field *U)
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  m_dr_smear->set_config(U);

  int   Nsmear = m_dr_smear->get_Nsmear();
  Field *Uptr  = m_dr_smear->getptr_smearedConfig(Nsmear);

  m_fopr_eo->set_config(Uptr);
}


//====================================================================
//============================================================END=====
