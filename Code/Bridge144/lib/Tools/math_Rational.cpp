#include "BridgeLib_Private.h"

/*!
        @file    $Id:: math_Rational.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "math_Rational.h"



const std::string Math_Rational::class_name = "Math_Rational";

//====================================================================
void Math_Rational::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Np, n_exp, d_exp;
  double x_min, x_max;

  int err = 0;
  err += params.fetch_int("number_of_poles", Np);
  err += params.fetch_int("exponent_numerator", n_exp);
  err += params.fetch_int("exponent_denominator", d_exp);
  err += params.fetch_double("lower_bound", x_min);
  err += params.fetch_double("upper_bound", x_max);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Np, n_exp, d_exp, x_min, x_max);
}


//====================================================================
void Math_Rational::set_parameters(const int Np, const int n_exp, const int d_exp,
                                   const double x_min, const double x_max)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Np    = %d\n", Np);
  vout.general(m_vl, "  n_exp = %d\n", n_exp);
  vout.general(m_vl, "  d_exp = %d\n", d_exp);
  vout.general(m_vl, "  x_min = %10.6f\n", x_min);
  vout.general(m_vl, "  x_max = %10.6f\n", x_max);

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(Np);
  err += ParameterCheck::non_zero(n_exp);
  err += ParameterCheck::non_zero(d_exp);
  // NB. x_min,x_max == 0 is allowed.

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Np    = Np;
  m_n_exp = n_exp;
  m_d_exp = d_exp;
  m_x_min = x_min;
  m_x_max = x_max;

  //- post-process
  m_res.resize(m_Np);
  m_pole.resize(m_Np);

  set_coeff();
}


//====================================================================
void Math_Rational::get_parameters(double& norm, std::vector<double>& res,
                                   std::vector<double>& pole)
{
  if (res.size() != m_Np) {
    vout.crucial(m_vl, "Error at %s: size of cl is not correct\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  if (pole.size() != m_Np) {
    vout.crucial(m_vl, "Error at %s: size of bl is not correct\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  norm = m_norm;
  for (int i = 0; i < m_Np; ++i) {
    res[i]  = m_res[i];
    pole[i] = m_pole[i];
  }
}


//====================================================================
void Math_Rational::set_coeff()
{
  read_file();
}


//====================================================================
void Math_Rational::read_file()
{
  // setting input file
  char filename[FILENAME_MAX];

  snprintf(filename, FILENAME_MAX,
           "parameter_rational.%1d_%1d%s_%02d_%010.8f_%07.3f",
           abs(m_n_exp), m_d_exp,
           ((m_n_exp < 0) ? "inv" : ""),
           m_Np,
           m_x_min, m_x_max);

  vout.general(m_vl, "%s: expected filename: %s\n", class_name.c_str(), filename);

  int    Np, n_exp, d_exp;
  double x_min, x_max;

  // read parameters from file
  std::fstream parameterfile;
  parameterfile.open(filename, std::ios::in);
  if (!parameterfile.is_open()) {
    vout.crucial(m_vl, "Error at %s: failed to open parameter file. %s(%d)\n",
                 class_name.c_str(), __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  parameterfile >> Np >> n_exp >> d_exp;
  parameterfile >> x_min >> x_max;
  parameterfile >> m_error;
  parameterfile >> m_norm;
  for (int i = 0; i < Np; i++) {
    parameterfile >> m_res[i] >> m_pole[i];
  }
  parameterfile.close();

  vout.general(m_vl, "%s: read parameter file successful.\n", class_name.c_str());

  /*
  vout.general(m_vl, " Rational approximation (read from file): %s\n",
          filename.c_str());
  vout.general(m_vl, " Np = %d\n", m_Np);
  vout.general(m_vl, " n_exp = %d,  d_exp = %d\n", m_n_exp,m_d_exp);
  vout.general(m_vl, " x_min = %12.8f\n",  m_x_min);
  vout.general(m_vl, " x_max = %12.8f\n",  m_x_max);
  vout.general(m_vl, " error = %18.16e\n", m_error);
  vout.general(m_vl, " RA_a0 = %18.16e\n", m_norm );
  for(int i = 0; i < n; i++){
    vout.general(m_vl, " RA_b[%d] = %18.16e, RA_c[%d] = %18.16e\n",
             i, m_res[i], i, m_pole[i]);
  }
  */

  assert(m_Np == Np);
  assert(m_n_exp == n_exp);
  assert(m_d_exp == d_exp);
  assert(fabs(m_x_min - x_min) < 1.e-12);
  //  assert(fabs(m_x_min - x_min) < 1.e-12);
  //  assert(fabs(m_x_max - x_max) < 1.e-12);
  assert(fabs(m_x_min - x_min) < 1.e-8);
  assert(fabs(m_x_max - x_max) < 1.e-8);
}


//====================================================================
double Math_Rational::func(double x)
{
  double y = m_norm;

  for (int k = 0; k < m_Np; ++k) {
    y += m_res[k] / (x + m_pole[k]);
  }

  return y;
}


//====================================================================
//============================================================END=====
