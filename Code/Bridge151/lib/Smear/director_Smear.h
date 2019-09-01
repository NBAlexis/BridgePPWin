/*!
        @file    director_Smear.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef DIRECTOR_SMEAR_INCLUDED
#define DIRECTOR_SMEAR_INCLUDED

#include <cassert>

#include "smear.h"
#include "Tools/director.h"

#include "Field/field_G.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Manager of smeared configurations.

/*!
    This director class BAPI handles smeared configurations.
                            [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.    [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                            [21 Mar 2015 Y.Namekawa]
 */

class BAPI Director_Smear : public Director
{
 public:
  static const std::string class_name;

 private:
  int m_Nsmear;                        //!< number of smearing to be applied
  Smear *m_smear;                      //!< smearing operator
  Field_G *m_U;                        //!< original thin link var.
  std::vector<Field_G> m_Usmear;       //!< smeared configs.
  int m_status_linkv;                  //!< set to zero when link var. is updated

 public:
  //! constructor requires pointer to Smear object
  Director_Smear(Smear *smear)
    : Director()
  {
    m_Nsmear       = 0;
    m_smear        = smear;
    m_U            = 0;
    m_status_linkv = 0;
  }

  Director_Smear(unique_ptr<Smear>& smear)
    : Director()
  {
    m_Nsmear       = 0;
    m_smear        = smear.get();
    m_U            = 0;
    m_status_linkv = 0;
  }

  //! set parameters, must be called before set_config
  void set_parameters(const Parameters& params);
  void set_parameters(const int Nsmear);

  //! get number of applied smearing operation
  int get_Nsmear() { return m_Nsmear; }

  //! get pointer to i-th smeared config (0th is original thin link)
  Field *getptr_smearedConfig(const int i_smear);

  Field_G *get_config();                  // smeared config
  Field_G *get_config(const int i_smear); // intermediate config

  //! set pointer to original thin link variable
  void set_config(Field *U);
  void set_config(unique_ptr<Field_G>& U);

  //! to be called when configuration is updated
  void notify_linkv()
  {
    m_status_linkv = 0;
  }

 private:
  //! smearing is performed by calling a function of Smear object
  void smear();
};
#endif
