/*!
        @file    $Id: solver_CGNE.h #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-03-21 15:28:34 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef SOLVER_CGNE_INCLUDED
#define SOLVER_CGNE_INCLUDED

#include "solver_CG.h"


//! CGNE solver.

/*!
    solve D x = b
    by
    (D Ddag) y = b  and  x = Ddag y.

    delegate to base CG solver.

    Introduce unique_ptr to avoid memory leaks.
                               [21 Mar 2015 Y.Namekawa]
    Add flop_count.            [ 8 Aug 2016 Y.Namekawa]
 */

class Solver_CGNE : public Solver_CG
{
 private:
  Field m_y;

 public:
  static const std::string class_name;

  Solver_CGNE(Fopr *fopr)
    : Solver_CG(fopr) {}

  Solver_CGNE(unique_ptr<Fopr>& fopr)
    : Solver_CG(fopr.get()) {}

  ~Solver_CGNE() {}

  void set_parameters(const Parameters& params);

  void solve(Field& solution, const Field& source, int& Nconv, double& diff);

  Fopr *get_fopr() { return this->Solver_CG::get_fopr(); }

  double flop_count();
};
#endif
