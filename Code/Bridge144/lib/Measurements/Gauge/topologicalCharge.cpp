#include "BridgeLib_Private.h"

/*!
        @file    $Id:: topologicalCharge.cpp #$

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "topologicalCharge.h"

const std::string TopologicalCharge::class_name = "TopologicalCharge";

//====================================================================
void TopologicalCharge::set_parameters(const Parameters& params)
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
void TopologicalCharge::set_parameters(double c_plaq, double c_rect)
{
  //- print input parameters
  vout.general(m_vl, "Topological Charge measurement:\n");
  vout.general(m_vl, "  c_plaq = %12.6f\n", c_plaq);
  vout.general(m_vl, "  c_rect = %12.6f\n", c_rect);

  //- range check
  // NB. beta,c_plaq,c_rect == 0 is allowed.

  //- store values
  m_c_plaq = c_plaq;
  m_c_rect = c_rect;
}


//====================================================================
double TopologicalCharge::measure(Field_G& U)
{
  int                 Ndim = CommonParameters::Ndim();
  static const double eps  = CommonParameters::epsilon_criterion();
  static const double PI   = 4.0 * atan(1.0);
  static const double PI2  = PI * PI;

  double Q_1x1 = 0.0;
  double Q_1x2 = 0.0;


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

  // (mu,nu, rho,sigma) = (1,2, 3,4)
  Q_1x1 += contract_epsilon_tensor(Fmunu_1x1[0], Fmunu_1x1[5]);

  // (mu,nu, rho,sigma) = (1,3, 2,4)
  Q_1x1 -= contract_epsilon_tensor(Fmunu_1x1[1], Fmunu_1x1[4]);

  // (mu,nu, rho,sigma) = (1,4, 2,3)
  Q_1x1 += contract_epsilon_tensor(Fmunu_1x1[2], Fmunu_1x1[3]);

  // #degeneracy of (mu,nu, rho,sigma) is 8
  Q_1x1 *= 8.0;

  // overall factor
  Q_1x1 /= (32.0 * PI2);
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

    Q_1x2 += contract_epsilon_tensor(Fmunu_1x2[0], Fmunu_1x2[5]);
    Q_1x2 -= contract_epsilon_tensor(Fmunu_1x2[1], Fmunu_1x2[4]);
    Q_1x2 += contract_epsilon_tensor(Fmunu_1x2[2], Fmunu_1x2[3]);
    Q_1x2 *= 8.0;
    // extra factor "2" for 1x2
    Q_1x2 *= 2.0;
    Q_1x2 /= (32.0 * PI2);
  }
  //----------------


  double Q_topo = m_c_plaq * Q_1x1 + m_c_rect * Q_1x2;


  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  vout.general(m_vl, "  Q_1x1  = %20.16e\n", Q_1x1);
  if (fabs(m_c_rect) > eps) {
    vout.general(m_vl, "  Q_1x2  = %20.16e\n", Q_1x2);
  }
  vout.general(m_vl, "  Q_topo = %20.16e\n", Q_topo);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }


  return Q_topo;
}


//====================================================================
double TopologicalCharge::contract_epsilon_tensor(Field_G& Fmunu_1, Field_G& Fmunu_2)
{
  int Nvol = CommonParameters::Nvol();

  double Q_topo = 0.0;

  for (int site = 0; site < Nvol; ++site) {
    Q_topo += ReTr(Fmunu_1.mat(site) * Fmunu_2.mat(site));
  }
  double result = Communicator::reduce_sum(Q_topo);

  return result;
}


//====================================================================
//============================================================END=====
