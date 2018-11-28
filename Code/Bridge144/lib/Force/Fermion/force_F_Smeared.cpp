#include "BridgeLib_Private.h"

/*!
        @file    $Id:: force_F_Smeared.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "force_F_Smeared.h"

const std::string Force_F_Smeared::class_name = "Force_F_Smeared";

//====================================================================
void Force_F_Smeared::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
void Force_F_Smeared::force_udiv(Field& force_, const Field& eta)
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);

  int Nsmear = m_director_smear->get_Nsmear();

  if (Nsmear == 0) {
    m_force->force_udiv(force, eta);
  } else {
    Field_G *Uptr = m_director_smear->get_config();

    m_force->set_config(Uptr);
    m_force->force_udiv(force, eta);

    mult_jacobian(force);
  }

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_Smeared::force_udiv1(Field& force_, const Field& zeta, const Field& eta)
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);

  int Nsmear = m_director_smear->get_Nsmear();

  if (Nsmear == 0) {
    m_force->force_udiv1(force, zeta, eta);
  } else {
    Field_G *Uptr = m_director_smear->get_config();

    m_force->set_config(Uptr);
    m_force->force_udiv1(force, zeta, eta);

    mult_jacobian(force);
  }

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_Smeared::mult_jacobian(Field_G& force)
{
  Field_G ftmp(force);  // copy to temporal field.
  int     Nsmear = m_director_smear->get_Nsmear();

  for (int ismear = Nsmear - 1; ismear >= 0; --ismear) {
    Field_G *Uptr = m_director_smear->get_config(ismear);

    m_force_smear->force_udiv(force, ftmp, *Uptr);

    if (ismear > 0) copy(ftmp, force);  // ftmp = force;
  }
}


//====================================================================
//============================================================END=====
