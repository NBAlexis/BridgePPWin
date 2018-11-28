/*!
        @file    $Id:: fopr_Smeared.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/


#ifndef FOPR_SMEARED_INCLUDED
#define FOPR_SMEARED_INCLUDED

#include "fopr.h"
#include "Smear/director_Smear.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! smeared fermion operator.

/*!
    This class construct a smeared fermon operator for a given
    base fermion operator together with smearing director.
    Both of them must be constructed beforehand outside this
    class and given to the constructor.
    Smearing of link configuration is triggered by call of
    set_config(), which calls set_config() of the smearing
    director and then gets the pointer to the smeared
    config. to set it as the config. of base fermion operator.
    When mult() is called, this class just call the mult()
    of base class.
                                  [24 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                  [21 Mar 2015 Y.Namekawa]
 */

class Fopr_Smeared : public Fopr
{
 public:
  static const std::string class_name;

 private:
  Fopr           *m_fopr;
  Director_Smear *m_dr_smear;

 public:
  //! constructor requires Fopr and Director_Smear objects
  Fopr_Smeared(Fopr *fopr, Director_Smear *dr_smear)
    : Fopr(), m_fopr(fopr), m_dr_smear(dr_smear) {}

  Fopr_Smeared(unique_ptr<Fopr>& fopr, unique_ptr<Director>& dr_smear)
    : Fopr(), m_fopr(fopr.get()), m_dr_smear((Director_Smear *)dr_smear.get()) {}


  void set_parameters(const Parameters&);

  //! set pointer to original thin link variable
  void set_config(Field *U);

  void set_config(unique_ptr<Field_G>& U)
  {
    set_config(U.get());
  }

  //! multiply smeared fermion operator
  void mult(Field& v, const Field& f)
  {
    m_fopr->mult(v, f);
  }

  //! multiply smeared fermion operator
  void mult_dag(Field& v, const Field& f)
  {
    m_fopr->mult_dag(v, f);
  }

  //! set the mode of fermion operator
  void set_mode(std::string mode)
  {
    m_fopr->set_mode(mode);
  }

  std::string get_mode() const
  {
    return m_fopr->get_mode();
  }

  void mult_up(int mu, Field& v, const Field& w)
  {
    m_fopr->mult_up(mu, v, w);
  }

  void mult_dn(int mu, Field& v, const Field& w)
  {
    m_fopr->mult_dn(mu, v, w);
  }

  int field_nvol() { return m_fopr->field_nvol(); }
  int field_nin() { return m_fopr->field_nin(); }
  int field_nex() { return m_fopr->field_nex(); }
};
#endif
