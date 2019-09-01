/*!
        @file    force_F_Smeared.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FORCE_F_SMEARED_INCLUDED
#define FORCE_F_SMEARED_INCLUDED

#include "force_F.h"
#include "forceSmear.h"

#include "Fopr/fopr.h"
#include "Smear/director_Smear.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for smeared fermion operators.

/*!
    This class BAPI determines the force of smeared fermion operator
    using smearing director (MultiSmear instance) and base
    fermion force instance.
                                      [28 Dec 2011 H.Matsufuru]
    Modified: set_mode() is added to incorporate non-hermitian H
                                      [21 Jan 2012 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                      [21 Mar 2015 Y.Namekawa]
*/

class BAPI Force_F_Smeared : public Force
{
 public:
  static const std::string class_name;

 private:
  Force *m_force;
  ForceSmear *m_force_smear;
  Director_Smear *m_director_smear;

 public:
  Force_F_Smeared(
    Force *force, ForceSmear *force_smear, Director_Smear *director_smear)
    : Force(), m_force(force), m_force_smear(force_smear), m_director_smear(director_smear) {}

  Force_F_Smeared(
    unique_ptr<Force>& force, unique_ptr<ForceSmear>& force_smear, unique_ptr<Director>& director_smear)
    : Force(), m_force(force.get()), m_force_smear(force_smear.get()), m_director_smear((Director_Smear *)director_smear.get()) {}

  void set_parameters(const Parameters&);

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_director_smear->set_config(U);
    m_force->set_config(m_director_smear->get_config());
  }

  void set_mode(const std::string& mode)
  {
    m_force->set_mode(mode);
  }

  void force_udiv(Field& force, const Field& eta);

  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:
  void mult_jacobian(Field_G& force);
};
#endif
