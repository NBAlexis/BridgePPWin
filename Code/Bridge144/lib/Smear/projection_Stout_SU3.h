/*!
        @file    $Id:: projection_Stout_SU3.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/


#ifndef PROJECTION_STOUT_SU3_INCLUDED
#define PROJECTION_STOUT_SU3_INCLUDED

#include <cassert>
#include "projection.h"
#include "bridge_complex.h"
#include "Tools/mat_SU_N.h"
//#include "libkek.h"  //- only for BG/Q

#include "IO/bridgeIO.h"
using Bridge::vout;

// The following implementation only valid for Nc = 3 case.
#define  NC    3


//! Stout(exponential)-type projection to SU(N) gauge group.

/*!
    Present implpementation applies to SU(3) case only, since
    the SU(3) properties are explicitly used.
                                    [08 Apr 2012 H.Matsufuru]
 */

class Projection_Stout_SU3 : public Projection
{
 public:
  static const std::string class_name;

 private:
  unsigned long int m_flop;
  double            m_time;

 public:

  Projection_Stout_SU3()
  {
    setup();
  }

  ~Projection_Stout_SU3() {}

  void set_parameters(const Parameters& param);

  //! projection U = P[alpha, C, Uorg]
  void project(Field_G& U,
               double alpha,
               const Field_G& C, const Field_G& Uorg);

  //! determination of fields for force calculation
  void force_recursive(Field_G& Xi, Field_G& iTheta,
                       double alpha, const Field_G& Sigmap,
                       const Field_G& C, const Field_G& U);

  void print_stat();

 private:
  void exp_iQ(Field_G& e_iQ, const Field_G& iQ);
  void exp_iQ_bf(Field_G& e_iQ, const Field_G& iQ);

  void set_uw(double& u, double& w,
              const Mat_SU_N& iQ1, const Mat_SU_N& iQ2);

  void set_fj(dcomplex& f0, dcomplex& f1, dcomplex& f2,
              const double& u, const double& w);

  double func_xi0(double w);
  double func_xi1(double w);

  void setup();
};
#endif
