/*!
        @file    fprop.h

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FPROP_INCLUDED
#define FPROP_INCLUDED

#include "Fopr/fopr.h"

#include "IO/bridgeIO.h"

//! Base class BAPI for fermion propagator class BAPI family.

/*!
                                        [28 Dec 2011 H.Matsufuru]
    Introduce unique_ptr to avoid memory leaks.
                                        [21 Mar 2015 Y.Namekawa]
    Add flop_count.                     [ 8 Aug 2016 Y.Namekawa]
 */

class BAPI Fprop
{
 protected:
  Bridge::VerboseLevel m_vl;

 public:
  Fprop()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~Fprop() {}

 private:
  // non-copyable
  Fprop(const Fprop&);
  Fprop& operator=(const Fprop&);

 public:
  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual void invert_D(Field&, const Field&, int&, double&)     = 0;
  virtual void invert_DdagD(Field&, const Field&, int&, double&) = 0;

  virtual void set_config(Field *) = 0;

  virtual double flop_count() = 0;
};
#endif
