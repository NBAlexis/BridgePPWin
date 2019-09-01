/*!
        @file    $Id:: math_Sign_Zolotarev.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef MATH_SIGN_ZOLOTAREV_INCLUDED
#define MATH_SIGN_ZOLOTAREV_INCLUDED

#include <cmath>
#include <cassert>

//#include "defs.h"
#include "Parameters/commonParameters.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Determination of Zolotarev coefficients.

/*!
    This class determines the Zolotarev's optimal coefficients
    of partial fractional approximation to 1/sqrt(x).
    Present implementation makes use of the code in Numerical
    Recipes, and thus cannot be put public.
    To be replaced with public implementation.
                                     [28 Dec 2011 H.Matsufuru]
 */

class Math_Sign_Zolotarev
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  int                 m_Np;
  double              m_bmax;
  std::vector<double> m_cl;
  std::vector<double> m_bl;

 public:
  Math_Sign_Zolotarev(int Np, double bmax)
    : m_vl(CommonParameters::Vlevel()),
      m_Np(Np), m_bmax(bmax)
  {
    set_sign_parameters();
  }

  void get_sign_parameters(std::vector<double>& cl, std::vector<double>& bl);

  double sign(double x);

 private:

  void set_sign_parameters();

  void poly_Zolotarev(double bmax, double& UK);

  void Jacobi_elliptic(double uu, double emmc,
                       double& sn, double& cn, double& dn);
};
#endif
