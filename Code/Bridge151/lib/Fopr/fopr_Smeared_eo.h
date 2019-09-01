/*!
        @file    fopr_Smeared_eo.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FOPR_SMEARED_EO_INCLUDED
#define FOPR_SMEARED_EO_INCLUDED

#include "fopr_eo.h"
#include "Smear/director_Smear.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! smeared fermion operator with even-odd preconditioning.

/*!
    This class BAPI constructs a smeared fermon operator for a given
    base fermion operator together with smearing director.
    Both of them must be constructed beforehand outside this
    class BAPI and given to the constructor.
    Smearing of link configuration is triggered by call of
    set_config(), which calls set_config() of the smearing
    director and then gets the pointer to the smeared
    config. to set it as the config. of base fermion operator.
    When mult() is called, this class BAPI just call the mult()
    of base class BAPI.
                                  [03 Mar 2013 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                  [21 Mar 2015 Y.Namekawa]
 */

class BAPI Fopr_Smeared_eo : public Fopr_eo
{
 public:
  static const std::string class_name;

 private:
  Fopr_eo *m_fopr_eo;
  Director_Smear *m_dr_smear;

 public:
  //! constructor requires Fopr and Director_Smear objects
  Fopr_Smeared_eo(Fopr_eo *fopr_eo, Director_Smear *dr_smear) : Fopr_eo()
  {
    m_fopr_eo  = fopr_eo;
    m_dr_smear = dr_smear;
  }

  Fopr_Smeared_eo(unique_ptr<Fopr_eo>& fopr_eo, unique_ptr<Director_Smear>& dr_smear) : Fopr_eo()
  {
    m_fopr_eo  = fopr_eo.get();
    m_dr_smear = dr_smear.get();
  }

  void set_parameters(const Parameters&);

  void preProp(Field& Be, Field& bo, const Field& b)
  {
    m_fopr_eo->preProp(Be, bo, b);
  }

  void postProp(Field& x, const Field& xe, const Field& bo)
  {
    m_fopr_eo->postProp(x, xe, bo);
  }

  //! set pointer to original thin link variable
  void set_config(Field *U);

  void set_config(unique_ptr<Field_G>& U)
  {
    set_config(U.get());
  }

  //! multiply smeared fermion operator
  void mult(Field& v, const Field& f)
  {
    m_fopr_eo->mult(v, f);
  }

  //! multiply smeared fermion operator
  void mult_dag(Field& v, const Field& f)
  {
    m_fopr_eo->mult_dag(v, f);
  }

  //! set the mode of fermion operator
  void set_mode(const std::string mode)
  {
    m_fopr_eo->set_mode(mode);
  }

  std::string get_mode() const
  {
    return m_fopr_eo->get_mode();
  }

  void mult_up(const int mu, Field& v, const Field& w)
  {
    m_fopr_eo->mult_up(mu, v, w);
  }

  void mult_dn(const int mu, Field& v, const Field& w)
  {
    m_fopr_eo->mult_dn(mu, v, w);
  }

  int field_nvol() { return m_fopr_eo->field_nvol(); }
  int field_nin() { return m_fopr_eo->field_nin(); }
  int field_nex() { return m_fopr_eo->field_nex(); }

#ifdef USE_FACTORY
 private:
  static Fopr *create_object(Fopr *fopr, Director *director)
  {
    return new Fopr_Smeared_eo(dynamic_cast<Fopr_eo *>(fopr), dynamic_cast<Director_Smear *>(director));
  }

 public:
  static bool register_factory()
  {
    return Fopr::Factory_fopr_director::Register("Smeared_eo", create_object);
  }
#endif
};
#endif
