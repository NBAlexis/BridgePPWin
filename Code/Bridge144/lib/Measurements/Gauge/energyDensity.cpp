#include "BridgeLib_Private.h"

/*!
        @file    $Id:: energyDensity.cpp #$

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2017-05-31 23:40:23 #$

        @version $LastChangedRevision: 1639 $
*/

#include "energyDensity.h"

const std::string EnergyDensity::class_name = "EnergyDensity";

//====================================================================
void EnergyDensity::set_parameters(const Parameters& params)
{
  m_filename_output = params.get_string("filename_output");
  if (m_filename_output.empty()) {
    m_filename_output = "stdout";
  }

  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double c_plaq, c_rect;

  int err = 0;
  err += params.fetch_double("c_plaq", c_plaq);
  err += params.fetch_double("c_rect", c_rect);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(c_plaq, c_rect);
}


//====================================================================
void EnergyDensity::set_parameters(double c_plaq, double c_rect)
{
  //- print input parameters
  vout.general(m_vl, "Energy density measurement:\n");
  vout.general(m_vl, "  c_plaq = %12.6f\n", c_plaq);
  vout.general(m_vl, "  c_rect = %12.6f\n", c_rect);

  //- range check
  // NB. beta,c_plaq,c_rect == 0 is allowed.

  //- store values
  m_c_plaq = c_plaq;
  m_c_rect = c_rect;
}


//====================================================================
double EnergyDensity::E_plaq(const Field_G& U)
{
  double plaq   = m_staple.plaquette(U);
  double E_plaq = 36.0 * (1.0 - plaq);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  vout.general(m_vl, "  E_plaq       = %20.16e\n", E_plaq);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  return E_plaq;
}


//====================================================================
double EnergyDensity::E_clover(const Field_G& U)
{
  int Lvol = CommonParameters::Lvol();
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  static const double eps = CommonParameters::epsilon_criterion();

  double E_clover_1x1 = 0.0;
  double E_clover_1x2 = 0.0;


  //--- 1x1 part ---
  // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
  std::vector<Field_G> Fmunu_1x1(6);

  int i_munu = 0;

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = mu + 1; nu < Ndim; ++nu) {
      m_field_strength.construct_Fmunu_1x1(Fmunu_1x1[i_munu], mu, nu, U);

      ++i_munu;
    }
  }

  double F2_1x1 = 0.0;
  for (int site = 0; site < Nvol; ++site) {
    for (int i_munu = 0; i_munu < 6; i_munu++) {
      F2_1x1 += ReTr(Fmunu_1x1[i_munu].mat(site) * Fmunu_1x1[i_munu].mat(site));
    }
  }
  E_clover_1x1 = Communicator::reduce_sum(F2_1x1) / Lvol;
  //----------------


  //--- 1x2 part ---
  // NB. skip this part, if m_c_rect = 0.0
  if (fabs(m_c_rect) > eps) {
    // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
    std::vector<Field_G> Fmunu_1x2(6);

    int i_munu = 0;
    for (int mu = 0; mu < Ndim; ++mu) {
      for (int nu = mu + 1; nu < Ndim; ++nu) {
        m_field_strength.construct_Fmunu_1x2(Fmunu_1x2[i_munu], mu, nu, U);

        ++i_munu;
      }
    }

    double F2_1x2 = 0.0;
    for (int site = 0; site < Nvol; ++site) {
      for (int i_munu = 0; i_munu < 6; i_munu++) {
        F2_1x2 += ReTr(Fmunu_1x2[i_munu].mat(site) * Fmunu_1x2[i_munu].mat(site));
      }
    }
    E_clover_1x2 = Communicator::reduce_sum(F2_1x2) / Lvol;

    // extra factor "2" for 1x2
    E_clover_1x2 *= 2.0;
  }
  //----------------


  double E_clover = (m_c_plaq * E_clover_1x1 + m_c_rect * E_clover_1x2);


  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  vout.general(m_vl, "  E_clover_1x1 = %20.16e\n", E_clover_1x1);
  if (fabs(m_c_rect) > eps) {
    vout.general(m_vl, "  E_clover_1x2 = %20.16e\n", E_clover_1x2);
  }
  vout.general(m_vl, "  E_clover     = %20.16e\n", E_clover);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }


  return E_clover;
}


//====================================================================
//============================================================END=====
