/*!
        @file    $Id: solver_CGNR.h #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2015-03-23 23:00:52 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef SOLVER_CGNR_INCLUDED
#define SOLVER_CGNR_INCLUDED

#include "solver_CG.h"


//! CGNR solver.

/*!
    solve D x = b
    by
    (Ddag D) x = b'  and  b' = Ddag b.

    delegate to base CG solver.

    Introduce unique_ptr to avoid memory leaks.
                                   [21 Mar 2015 Y.Namekawa]
    Add flop_count.                [ 8 Aug 2016 Y.Namekawa]
 */

class Solver_CGNR : public Solver_CG
{
 private:
  Field m_b2;

 public:
  static const std::string class_name;

  Solver_CGNR(Fopr *fopr)
    : Solver_CG(fopr) {}

  Solver_CGNR(unique_ptr<Fopr>& fopr)
    : Solver_CG(fopr.get()) {}

  ~Solver_CGNR() {}

  void set_parameters(const Parameters& params);

  void solve(Field& solution, const Field& source, int& Nconv, double& diff);

  Fopr *get_fopr() { return this->Solver_CG::get_fopr(); }

  double flop_count();
};
#endif
