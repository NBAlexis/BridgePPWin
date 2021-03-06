/*!
        @file    fopr_Chebyshev.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/


#ifndef FOPR_CHEBYSHEV_INCLUDED
#define FOPR_CHEBYSHEV_INCLUDED

#include "fopr.h"
#include "Field/field_F.h"
#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Chebyshev polynomial of fermion operator.

/*!
    This class BAPI implements Chebyshev polynomial of a given
    fermion operator by making use of Clenshow's reccurence
    formula.
    The present version is just a rough implementation which
    assumes to be used accelerate eigenvalue solver.
                                    [24 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
 */



class BAPI Fopr_Chebyshev : public Fopr
{
 public:
  static const std::string class_name;

 private:
  int m_Npcb;
  double m_Fcb1, m_Fcb2;
  Fopr *m_fopr;

 public:
  Fopr_Chebyshev(Fopr *fopr)
    : Fopr(), m_fopr(fopr) {}

  Fopr_Chebyshev(unique_ptr<Fopr>& fopr)
    : Fopr(), m_fopr(fopr.get()) {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Np, const double v_thrs, const double v_max);

  void set_config(Field *U)
  {
    m_fopr->set_config(U);
  }

  void set_config(unique_ptr<Field_G>& U)
  {
    m_fopr->set_config(U.get());
  }

  void mult(Field& v, const Field& f);

  void mult_dag(Field& v, const Field& f)
  {
    mult(v, f);
  }

  //!  evaluate for a number
  double mult(const double);

  int field_nvol() { return m_fopr->field_nvol(); }
  int field_nin() { return m_fopr->field_nin(); }
  int field_nex() { return m_fopr->field_nex(); }

#ifdef USE_FACTORY
 private:
  static Fopr *create_object(Fopr *fopr)
  {
    return new Fopr_Chebyshev(fopr);
  }

 public:
  static bool register_factory()
  {
    return Fopr::Factory_fopr::Register("Chevyshev", create_object);
  }
#endif
};
#endif
