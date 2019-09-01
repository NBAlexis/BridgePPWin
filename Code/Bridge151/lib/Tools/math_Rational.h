/*!
        @file    math_Rational.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2019-02-14 17:34:46 #$

        @version $LastChangedRevision: 1946 $
*/

#ifndef MATH_RATIONAL_INCLUDED
#define MATH_RATIONAL_INCLUDED

#include <fstream>
#include <cmath>
#include <cassert>

#include "Parameters/parameters.h"
#include "Parameters/commonParameters.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Determionation of coefficients of rational approximation.

/*!
    At present, the coefficients are determined outside this
    code system and saved to a file. This class BAPI reads the file
    to set the coefficients.
    In future, self-calculation should be implemented.
                                    [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
 */


class BAPI Math_Rational
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  int m_Np, m_n_exp, m_d_exp;
  double m_x_min, m_x_max;
  double m_norm, m_error;
  std::vector<double> m_res;
  std::vector<double> m_pole;

 public:
  Math_Rational()
    : m_vl(CommonParameters::Vlevel()) {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Np, const int n_exp, const int d_exp,
                      const double x_min, const double x_max);

  void get_parameters(double& norm, std::vector<double>& res,
                      std::vector<double>& pole);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  double func(const double x);

 private:
  void set_coeff();
  void read_file();
};
#endif
