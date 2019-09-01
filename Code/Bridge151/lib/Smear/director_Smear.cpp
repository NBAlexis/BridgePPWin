/*!
        @file    director_Smear.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "director_Smear.h"

const std::string Director_Smear::class_name = "Director_Smear";

//====================================================================
void Director_Smear::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int Nsmear;

  int err = 0;
  err += params.fetch_int("number_of_smearing", Nsmear);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Nsmear);
}


//====================================================================
void Director_Smear::set_parameters(const int Nsmear)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Nsmear = %d\n", Nsmear);

  //- range check
  // NB. Nsmear == 0 is allowed.

  //- store values
  m_Nsmear = Nsmear;

  //- post-process
  m_Usmear.resize(m_Nsmear);
  if (m_Nsmear > 0) {
    for (int i_smear = 0; i_smear < m_Nsmear; ++i_smear) {
      m_Usmear[i_smear].reset(Nvol, Ndim);
    }
    vout.detailed(m_vl, " size of Usmear[i_smear] was set.\n");
  }
}


//====================================================================
void Director_Smear::set_config(Field *U)
{
  m_U = (Field_G *)U;

  if (m_status_linkv == 0) smear();
}


//====================================================================
void Director_Smear::set_config(unique_ptr<Field_G>& U)
{
  set_config(U.get());
}


//====================================================================
Field *Director_Smear::getptr_smearedConfig(const int i_smear)
{
  assert(m_U != 0);

  if (i_smear == 0) {
    return m_U;
  } else {
    return &m_Usmear[i_smear - 1];
  }
}


//====================================================================
Field_G *Director_Smear::get_config()
{
  return m_Nsmear == 0 ? m_U : &m_Usmear[m_Nsmear - 1];
}


//====================================================================
Field_G *Director_Smear::get_config(const int i_smear)
{
  if (i_smear == 0) return m_U;

  return &m_Usmear[i_smear - 1];
}


//====================================================================
void Director_Smear::smear()
{
  if (m_Nsmear > 0) {
    const int Nvol = m_U->nvol();
    const int Nex  = m_U->nex();

    for (int i_smear = 0; i_smear < m_Nsmear; ++i_smear) {
      Field_G Uprev(Nvol, Nex);
      if (i_smear == 0) {
        Uprev = *m_U;
      } else {
        Uprev = m_Usmear[i_smear - 1];
      }

      Field_G Usmear(Nvol, Nex);
      m_smear->smear(Usmear, Uprev);
      m_Usmear[i_smear] = (Field)Usmear;
    }
  }

  ++m_status_linkv;
}


//====================================================================
//============================================================END=====
