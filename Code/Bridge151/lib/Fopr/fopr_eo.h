/*!
        @file    fopr_eo.h

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FOPR_EO_INCLUDED
#define FOPR_EO_INCLUDED

#include "fopr.h"

//! Base class BAPI of fermion operator family.

/*!
    In Bridge-OF, the fermion operator implies an operator which
    transforms a field to other field, irrespective of physical
    formulation.
    This class BAPI defines the interface of the fermion operators.
    At present, void functions mult(v,w) and mult_dag(v,w) is
    not purely virtual, because some of subclass have not
    implemented them yet.
                          [20 Dec 2011  H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                          [21 Mar 2015 Y.Namekawa]
*/

class BAPI Fopr_eo : public Fopr
{
 public:
  virtual ~Fopr_eo() {}

  // methods for even odd fermion operator
  virtual void preProp(Field&, Field&, const Field&)        = 0;
  virtual void postProp(Field&, const Field&, const Field&) = 0;

  //! setting pointer to the gauge configuration.
  virtual void set_config(Field *) = 0;
  virtual void set_config(unique_ptr<Field_G>&) = 0;

  //! \brief multiplies fermion operator to a given field (2nd argument)
  //   and set the resultant field to the 1st argument.
  virtual void mult(Field&, const Field&) {}

  //! hermitian conjugate of mult(Field&, const Field&).
  virtual void mult_dag(Field&, const Field&) {}

  //! \brief setting the mode of multiplication if necessary.
  //!  Default implementation here is just to avoid irrelevant call.
  virtual void set_mode(const std::string mode)
  {
    vout.general(m_vl, "Fopr_eo: set_mode not implemented.\n");
  }

  std::string get_mode() const
  {
    vout.general(m_vl, "Fopr_eo: get_mode not implemented.\n");
    return std::string();
  }

  //! returns the volume for which the fermion operator is defined.
  virtual int field_nvol() = 0;

  //! returns the on-site d.o.f. for which the fermion operator is defined.
  virtual int field_nin() = 0;

  //! returns the external d.o.f. for which the fermion operator is defined.
  virtual int field_nex() = 0;
};
#endif
