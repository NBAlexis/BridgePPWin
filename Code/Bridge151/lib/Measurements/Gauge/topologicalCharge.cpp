/*!
        @file    topologicalCharge.cpp

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
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
  int    max_mom;

  int err = 0;
  err += params.fetch_double("c_plaq", c_plaq);
  err += params.fetch_double("c_rect", c_rect);
  err += params.fetch_int("max_momentum", max_mom);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(c_plaq, c_rect, max_mom);
}


//====================================================================
void TopologicalCharge::set_parameters(const double c_plaq, const double c_rect, const int max_mom)
{
  //- print input parameters
  vout.general(m_vl, "Topological Charge measurement:\n");
  vout.general(m_vl, "  c_plaq = %12.6f\n", c_plaq);
  vout.general(m_vl, "  c_rect = %12.6f\n", c_rect);
  vout.general(m_vl, "  max_momentum = %d\n", max_mom);

  //- range check
  // NB. beta,c_plaq,c_rect == 0 is allowed.

  //- store values
  m_c_plaq  = c_plaq;
  m_c_rect  = c_rect;
  m_max_mom = max_mom;
}


//====================================================================
double TopologicalCharge::measure(const Field_G& U)
{
  const int Ndim = CommonParameters::Ndim();

  static const double eps = CommonParameters::epsilon_criterion();
  static const double PI  = 4.0 * atan(1.0);
  static const double PI2 = PI * PI;


  //--- 1x1 part ---
  double Q_1x1 = 0.0;

  {
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
    Q_1x1 += m_field_strength.contract(Fmunu_1x1[0], Fmunu_1x1[5]);

    // (mu,nu, rho,sigma) = (1,3, 2,4)
    Q_1x1 -= m_field_strength.contract(Fmunu_1x1[1], Fmunu_1x1[4]);

    // (mu,nu, rho,sigma) = (1,4, 2,3)
    Q_1x1 += m_field_strength.contract(Fmunu_1x1[2], Fmunu_1x1[3]);

    // #degeneracy of (mu,nu, rho,sigma) is 8
    Q_1x1 *= 8.0;

    // overall factor
    Q_1x1 /= (32.0 * PI2);
  }
  //----------------


  //--- 1x2 part ---
  double Q_1x2 = 0.0;

  {
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

      Q_1x2 += m_field_strength.contract(Fmunu_1x2[0], Fmunu_1x2[5]);
      Q_1x2 -= m_field_strength.contract(Fmunu_1x2[1], Fmunu_1x2[4]);
      Q_1x2 += m_field_strength.contract(Fmunu_1x2[2], Fmunu_1x2[3]);
      Q_1x2 *= 8.0;
      // extra factor "2" for 1x2
      Q_1x2 *= 2.0;
      Q_1x2 /= (32.0 * PI2);
    }
  }
  //----------------


  const double Q_topo = m_c_plaq * Q_1x1 + m_c_rect * Q_1x2;


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
void TopologicalCharge::measure_topological_charge(const double tt)
{
  assert(m_flag_field_set == true);

  static const double PI  = 4.0 * atan(1.0);
  static const double PI2 = PI * PI;
  // #degeneracy of (mu,nu, rho,sigma) is 8
  static const double factor1 = 8.0 / (32.0 * PI2);
  // extra factor "2" for 1x2
  static const double factor2 = 16.0 / (32.0 * PI2);

  //--- Field strength by one plaquette ---
  double Q_plaq = 0.0;

  Q_plaq += m_field_strength.contract(m_Fmunu_plaq[0], m_Fmunu_plaq[5]);
  Q_plaq -= m_field_strength.contract(m_Fmunu_plaq[1], m_Fmunu_plaq[4]);
  Q_plaq += m_field_strength.contract(m_Fmunu_plaq[2], m_Fmunu_plaq[3]);
  Q_plaq *= factor1;

  //--- 1x1 part ---
  double Q_1x1 = 0.0;

  // (mu,nu, rho,sigma) = (1,2, 3,4)
  Q_1x1 += m_field_strength.contract(m_Fmunu_1x1[0], m_Fmunu_1x1[5]);
  // (mu,nu, rho,sigma) = (1,3, 2,4)
  Q_1x1 -= m_field_strength.contract(m_Fmunu_1x1[1], m_Fmunu_1x1[4]);
  // (mu,nu, rho,sigma) = (1,4, 2,3)
  Q_1x1 += m_field_strength.contract(m_Fmunu_1x1[2], m_Fmunu_1x1[3]);
  Q_1x1 *= factor1;
  //----------------

  //--- 1x2 part ---
  double Q_1x2 = 0.0;

  Q_1x2 += m_field_strength.contract(m_Fmunu_1x2[0], m_Fmunu_1x2[5]);
  Q_1x2 -= m_field_strength.contract(m_Fmunu_1x2[1], m_Fmunu_1x2[4]);
  Q_1x2 += m_field_strength.contract(m_Fmunu_1x2[2], m_Fmunu_1x2[3]);
  Q_1x2 *= factor2;
  //----------------

  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;
  const double Q_imp = m_c_plaq * Q_1x1 + m_c_rect * Q_1x2;

  m_c_rect = 1.0 / 8.0;
  const double Q_rect = m_c_rect * Q_1x2;

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  vout.general(m_vl, "  Q_clover_plaq = %.8f  %.16e\n", tt, Q_1x1);
  vout.general(m_vl, "  Q_clover_rect = %.8f  %.16e\n", tt, Q_rect);
  vout.general(m_vl, "  Q_clover_imp = %.8f  %.16e\n", tt, Q_imp);
  vout.general(m_vl, "  Q_plaq = %.8f  %.16e\n", tt, Q_plaq);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
}


//====================================================================
void TopologicalCharge::measure_topological_density_t(const double tt)
{
  assert(m_flag_field_set == true);

  static const double PI      = 4.0 * atan(1.0);
  static const double PI2     = PI * PI;
  static const double factor1 = 8.0 / (32.0 * PI2);
  static const double factor2 = 16.0 / (32.0 * PI2);

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  const int Lt = CommonParameters::Lt();

  std::vector<double> corr_scr(Lt);

  //--- Field strength by one plaquette ---
  std::vector<double> corr_plaq(Lt, 0.0);

  m_field_strength.contract_at_t(corr_scr, m_Fmunu_plaq[0], m_Fmunu_plaq[5]);
  for (int t = 0; t < Lt; ++t) {
    corr_plaq[t] += corr_scr[t];
  }
  m_field_strength.contract_at_t(corr_scr, m_Fmunu_plaq[1], m_Fmunu_plaq[4]);
  for (int t = 0; t < Lt; ++t) {
    corr_plaq[t] -= corr_scr[t];
  }
  m_field_strength.contract_at_t(corr_scr, m_Fmunu_plaq[2], m_Fmunu_plaq[3]);
  for (int t = 0; t < Lt; ++t) {
    corr_plaq[t] += corr_scr[t];
  }
  //  double Q_plaq = 0.0;
  for (int t = 0; t < Lt; ++t) {
    corr_plaq[t] *= factor1;
    //    Q_plaq += corr_plaq[t];
  }

  //--- 1x1 part ---
  std::vector<double> corr_1x1(Lt, 0.0);

  // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
  // (mu,nu, rho,sigma) = (1,2, 3,4)
  m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x1[0], m_Fmunu_1x1[5]);
  for (int t = 0; t < Lt; ++t) {
    corr_1x1[t] += corr_scr[t];
  }
  // (mu,nu, rho,sigma) = (1,3, 2,4)
  m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x1[1], m_Fmunu_1x1[4]);
  for (int t = 0; t < Lt; ++t) {
    corr_1x1[t] -= corr_scr[t];
  }
  // (mu,nu, rho,sigma) = (1,4, 2,3)
  m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x1[2], m_Fmunu_1x1[3]);
  for (int t = 0; t < Lt; ++t) {
    corr_1x1[t] += corr_scr[t];
  }
  //  double Q_1x1 = 0.0;
  for (int t = 0; t < Lt; ++t) {
    corr_1x1[t] *= factor1;
    //    Q_1x1 += corr_1x1[t];
  }

  //--- 1x2 part ---
  std::vector<double> corr_1x2(Lt, 0.0);

  // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
  m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x2[0], m_Fmunu_1x2[5]);
  for (int t = 0; t < Lt; ++t) {
    corr_1x2[t] += corr_scr[t];
  }
  m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x2[1], m_Fmunu_1x2[4]);
  for (int t = 0; t < Lt; ++t) {
    corr_1x2[t] -= corr_scr[t];
  }
  m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x2[2], m_Fmunu_1x2[3]);
  for (int t = 0; t < Lt; ++t) {
    corr_1x2[t] += corr_scr[t];
  }
  //  double Q_1x2 = 0.0;
  for (int t = 0; t < Lt; ++t) {
    corr_1x2[t] *= factor2;
    //    Q_1x2 += corr_1x2[t];
  }
  //----------------

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int t = 0; t < Lt; ++t) {
    vout.general(m_vl, "  Q_clover_plaq_t = %.8f %d %.16e\n", tt, t, corr_1x1[t]);

    double scr = m_c_plaq * corr_1x1[t] + m_c_rect * corr_1x2[t];
    vout.general(m_vl, "  Q_clover_imp_t = %.8f %d %.16e\n", tt, t, scr);

    scr = l_c_rect * corr_1x2[t];
    vout.general(m_vl, "  Q_clover_rect_t = %.8f %d %.16e\n", tt, t, scr);
    vout.general(m_vl, "  Q_plaq_t = %.8f %d %.16e\n", tt, t, corr_plaq[t]);
  }

  /*
  vout.general(m_vl, "  Q_clover_plaq_sumt = %.8f  %.16e\n", tt, Q_1x1);
  const double Q_rect = l_c_rect * Q_1x2;
  vout.general(m_vl, "  Q_clover_rect_sumt = %.8f  %.16e\n", tt, Q_rect);
  const double Q_imp = m_c_plaq * Q_1x1 + m_c_rect * Q_1x2;
  vout.general(m_vl, "  Q_clover_imp_sumt = %.8f  %.16e\n", tt, Q_imp);
  vout.general(m_vl, "  Q_plaq_sumt = %.8f  %.16e\n", tt, Q_plaq);
  */
  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
}


//====================================================================
void TopologicalCharge::measure_topological_density_t_FT(const double tt)
{
  assert(m_flag_field_set == 1);
  static const double PI       = 4.0 * atan(1.0);
  static const double PI2      = PI * PI;
  static const double factor1  = 8.0 / (32.0 * PI2);
  static const double factor2  = 16.0 / (32.0 * PI2);
  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  const int Lt = CommonParameters::Lt();

  const int   Np = (2 * m_max_mom + 1);
  vector<int> source_position(4, 0);
  vector<int> momentum_sink(3);

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;

        std::vector<double> corr_plaq(Lt, 0);
        std::vector<double> corr_1x1(Lt, 0);
        std::vector<double> corr_1x2(Lt, 0);
        std::vector<double> corr_scr(Lt);

        //--- Field strength by one plaquette ---
        m_field_strength.contract_at_t(corr_scr, m_Fmunu_plaq[0], m_Fmunu_plaq[5], momentum_sink, source_position);
        for (int t = 0; t < Lt; ++t) {
          corr_plaq[t] += corr_scr[t];
        }
        m_field_strength.contract_at_t(corr_scr, m_Fmunu_plaq[1], m_Fmunu_plaq[4], momentum_sink, source_position);
        for (int t = 0; t < Lt; ++t) {
          corr_plaq[t] -= corr_scr[t];
        }
        m_field_strength.contract_at_t(corr_scr, m_Fmunu_plaq[2], m_Fmunu_plaq[3], momentum_sink, source_position);
        for (int t = 0; t < Lt; ++t) {
          corr_plaq[t] += corr_scr[t];
        }
        //	double Q_plaq = 0.0;
        for (int t = 0; t < Lt; ++t) {
          corr_plaq[t] *= factor1;
          //	  Q_plaq += corr_plaq[t];
        }

        //--- 1x1 part ---
        // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
        // (mu,nu, rho,sigma) = (1,2, 3,4)
        m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x1[0], m_Fmunu_1x1[5], momentum_sink, source_position);
        for (int t = 0; t < Lt; ++t) {
          corr_1x1[t] += corr_scr[t];
        }
        // (mu,nu, rho,sigma) = (1,3, 2,4)
        m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x1[1], m_Fmunu_1x1[4], momentum_sink, source_position);
        for (int t = 0; t < Lt; ++t) {
          corr_1x1[t] -= corr_scr[t];
        }
        // (mu,nu, rho,sigma) = (1,4, 2,3)
        m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x1[2], m_Fmunu_1x1[3], momentum_sink, source_position);
        for (int t = 0; t < Lt; ++t) {
          corr_1x1[t] += corr_scr[t];
        }
        //	double Q_1x1 = 0.0;
        for (int t = 0; t < Lt; ++t) {
          corr_1x1[t] *= factor1;
          //	  Q_1x1 += corr_1x1[t];
        }

        //--- 1x2 part ---
        // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
        m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x2[0], m_Fmunu_1x2[5], momentum_sink, source_position);
        for (int t = 0; t < Lt; ++t) {
          corr_1x2[t] += corr_scr[t];
        }
        m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x2[1], m_Fmunu_1x2[4], momentum_sink, source_position);
        for (int t = 0; t < Lt; ++t) {
          corr_1x2[t] -= corr_scr[t];
        }
        m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x2[2], m_Fmunu_1x2[3], momentum_sink, source_position);
        for (int t = 0; t < Lt; ++t) {
          corr_1x2[t] += corr_scr[t];
        }
        //	double Q_1x2 = 0.0;
        for (int t = 0; t < Lt; ++t) {
          corr_1x2[t] *= factor2;
          //	  Q_1x2 += corr_1x2[t];
        }
        //----------------

        //- output
        std::ostream& log_file_previous = vout.getStream();
        std::ofstream log_file;

        if (m_filename_output != "stdout") {
          log_file.open(m_filename_output.c_str(), std::ios::app);
          vout.init(log_file);
        }

        for (int t = 0; t < Lt; ++t) {
          vout.general(m_vl, "  Q_clover_plaq_t_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, corr_1x1[t]);
          double scr = m_c_plaq * corr_1x1[t] + m_c_rect * corr_1x2[t];
          vout.general(m_vl, "  Q_clover_imp_t_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, scr);
          scr = l_c_rect * corr_1x2[t];
          vout.general(m_vl, "  Q_clover_rect_t_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, scr);
          vout.general(m_vl, "  Q_plaq_t_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, corr_plaq[t]);
        }

        /*
        vout.general(m_vl, "  Q_clover_plaq_sumt_FT = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_1x1);
        double Q_rect = l_c_rect * Q_1x2;
        vout.general(m_vl, "  Q_clover_rect_sumt_FT = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_rect);
        double Q_imp = m_c_plaq * Q_1x1 + m_c_rect * Q_1x2;
        vout.general(m_vl, "  Q_clover_imp_sumt_FT = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_imp);
        vout.general(m_vl, "  Q_plaq_sumt = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_plaq);
        */

        if (m_filename_output != "stdout") {
          log_file.close();
          vout.init(log_file_previous);
        }
      }
    }
  }
}


//====================================================================
void TopologicalCharge::measure_topological_density_x(const double tt)
{
  assert(m_flag_field_set == true);

  static const double PI      = 4.0 * atan(1.0);
  static const double PI2     = PI * PI;
  static const double factor1 = 8.0 / (32.0 * PI2);
  static const double factor2 = 16.0 / (32.0 * PI2);

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  const int Lx = CommonParameters::Lx();

  std::vector<double> corr_scr(Lx);

  //--- Field strength by one plaquette ---
  std::vector<double> corr_plaq(Lx, 0.0);

  m_field_strength.contract_at_x(corr_scr, m_Fmunu_plaq[0], m_Fmunu_plaq[5]);
  for (int x = 0; x < Lx; ++x) {
    corr_plaq[x] += corr_scr[x];
  }
  m_field_strength.contract_at_x(corr_scr, m_Fmunu_plaq[1], m_Fmunu_plaq[4]);
  for (int x = 0; x < Lx; ++x) {
    corr_plaq[x] -= corr_scr[x];
  }
  m_field_strength.contract_at_x(corr_scr, m_Fmunu_plaq[2], m_Fmunu_plaq[3]);
  for (int x = 0; x < Lx; ++x) {
    corr_plaq[x] += corr_scr[x];
  }
  //  double Q_plaq = 0.0;
  for (int x = 0; x < Lx; ++x) {
    corr_plaq[x] *= factor1;
    //    Q_plaq += corr_plaq[x];
  }

  //--- 1x1 part ---
  std::vector<double> corr_1x1(Lx, 0.0);

  // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
  // (mu,nu, rho,sigma) = (1,2, 3,4)
  m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x1[0], m_Fmunu_1x1[5]);
  for (int x = 0; x < Lx; ++x) {
    corr_1x1[x] += corr_scr[x];
  }
  // (mu,nu, rho,sigma) = (1,3, 2,4)
  m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x1[1], m_Fmunu_1x1[4]);
  for (int x = 0; x < Lx; ++x) {
    corr_1x1[x] -= corr_scr[x];
  }
  // (mu,nu, rho,sigma) = (1,4, 2,3)
  m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x1[2], m_Fmunu_1x1[3]);
  for (int x = 0; x < Lx; ++x) {
    corr_1x1[x] += corr_scr[x];
  }
  //  double Q_1x1 = 0.0;
  for (int x = 0; x < Lx; ++x) {
    corr_1x1[x] *= factor1;
    //    Q_1x1 += corr_1x1[x];
  }

  //--- 1x2 part ---
  std::vector<double> corr_1x2(Lx, 0.0);

  // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
  m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x2[0], m_Fmunu_1x2[5]);
  for (int x = 0; x < Lx; ++x) {
    corr_1x2[x] += corr_scr[x];
  }
  m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x2[1], m_Fmunu_1x2[4]);
  for (int x = 0; x < Lx; ++x) {
    corr_1x2[x] -= corr_scr[x];
  }
  m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x2[2], m_Fmunu_1x2[3]);
  for (int x = 0; x < Lx; ++x) {
    corr_1x2[x] += corr_scr[x];
  }
  //  double Q_1x2 = 0.0;
  for (int x = 0; x < Lx; ++x) {
    corr_1x2[x] *= factor2;
    //    Q_1x2 += corr_1x2[x];
  }
  //----------------

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int x = 0; x < Lx; ++x) {
    vout.general(m_vl, "  Q_clover_plaq_x = %.8f %d %.16e\n", tt, x, corr_1x1[x]);

    double scr = m_c_plaq * corr_1x1[x] + m_c_rect * corr_1x2[x];
    vout.general(m_vl, "  Q_clover_imp_x = %.8f %d %.16e\n", tt, x, scr);

    scr = l_c_rect * corr_1x2[x];
    vout.general(m_vl, "  Q_clover_rect_x = %.8f %d %.16e\n", tt, x, scr);
    vout.general(m_vl, "  Q_plaq_x = %.8f %d %.16e\n", tt, x, corr_plaq[x]);
  }

  /*
  vout.general(m_vl, "  Q_clover_plaq_sumx = %.8f  %.16e\n", tt, Q_1x1);
  const double Q_rect = l_c_rect * Q_1x2;
  vout.general(m_vl, "  Q_clover_rect_sumx = %.8f  %.16e\n", tt, Q_rect);
  const double Q_imp = m_c_plaq * Q_1x1 + m_c_rect * Q_1x2;
  vout.general(m_vl, "  Q_clover_imp_sumx = %.8f  %.16e\n", tt, Q_imp);
  vout.general(m_vl, "  Q_plaq_sumx = %.8f  %.16e\n", tt, Q_plaq);
  */

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
}


//====================================================================
void TopologicalCharge::measure_topological_density_x_FT(const double tt)
{
  assert(m_flag_field_set == 1);
  static const double PI       = 4.0 * atan(1.0);
  static const double PI2      = PI * PI;
  static const double factor1  = 8.0 / (32.0 * PI2);
  static const double factor2  = 16.0 / (32.0 * PI2);
  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  const int Lx = CommonParameters::Lx();

  const int   Np = (2 * m_max_mom + 1);
  vector<int> source_position(4, 0);
  vector<int> momentum_sink(3);

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;

        std::vector<double> corr_plaq(Lx, 0);
        std::vector<double> corr_1x1(Lx, 0);
        std::vector<double> corr_1x2(Lx, 0);
        std::vector<double> corr_scr(Lx);

        //--- Field strength by one plaquette ---
        m_field_strength.contract_at_x(corr_scr, m_Fmunu_plaq[0], m_Fmunu_plaq[5], momentum_sink, source_position);
        for (int t = 0; t < Lx; ++t) {
          corr_plaq[t] += corr_scr[t];
        }
        m_field_strength.contract_at_x(corr_scr, m_Fmunu_plaq[1], m_Fmunu_plaq[4], momentum_sink, source_position);
        for (int t = 0; t < Lx; ++t) {
          corr_plaq[t] -= corr_scr[t];
        }
        m_field_strength.contract_at_x(corr_scr, m_Fmunu_plaq[2], m_Fmunu_plaq[3], momentum_sink, source_position);
        for (int t = 0; t < Lx; ++t) {
          corr_plaq[t] += corr_scr[t];
        }
        //	double Q_plaq = 0.0;
        for (int t = 0; t < Lx; ++t) {
          corr_plaq[t] *= factor1;
          //	  Q_plaq += corr_plaq[t];
        }

        //--- 1x1 part ---
        // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
        // (mu,nu, rho,sigma) = (1,2, 3,4)
        m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x1[0], m_Fmunu_1x1[5], momentum_sink, source_position);
        for (int t = 0; t < Lx; ++t) {
          corr_1x1[t] += corr_scr[t];
        }
        // (mu,nu, rho,sigma) = (1,3, 2,4)
        m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x1[1], m_Fmunu_1x1[4], momentum_sink, source_position);
        for (int t = 0; t < Lx; ++t) {
          corr_1x1[t] -= corr_scr[t];
        }
        // (mu,nu, rho,sigma) = (1,4, 2,3)
        m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x1[2], m_Fmunu_1x1[3], momentum_sink, source_position);
        for (int t = 0; t < Lx; ++t) {
          corr_1x1[t] += corr_scr[t];
        }
        //	double Q_1x1 = 0.0;
        for (int t = 0; t < Lx; ++t) {
          corr_1x1[t] *= factor1;
          //	  Q_1x1 += corr_1x1[t];
        }

        //--- 1x2 part ---
        // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
        m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x2[0], m_Fmunu_1x2[5], momentum_sink, source_position);
        for (int t = 0; t < Lx; ++t) {
          corr_1x2[t] += corr_scr[t];
        }
        m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x2[1], m_Fmunu_1x2[4], momentum_sink, source_position);
        for (int t = 0; t < Lx; ++t) {
          corr_1x2[t] -= corr_scr[t];
        }
        m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x2[2], m_Fmunu_1x2[3], momentum_sink, source_position);
        for (int t = 0; t < Lx; ++t) {
          corr_1x2[t] += corr_scr[t];
        }
        //	double Q_1x2 = 0.0;
        for (int t = 0; t < Lx; ++t) {
          corr_1x2[t] *= factor2;
          //	  Q_1x2 += corr_1x2[t];
        }
        //----------------

        //- output
        std::ostream& log_file_previous = vout.getStream();
        std::ofstream log_file;

        if (m_filename_output != "stdout") {
          log_file.open(m_filename_output.c_str(), std::ios::app);
          vout.init(log_file);
        }

        for (int t = 0; t < Lx; ++t) {
          vout.general(m_vl, "  Q_clover_plaq_x_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, corr_1x1[t]);
          double scr = m_c_plaq * corr_1x1[t] + m_c_rect * corr_1x2[t];
          vout.general(m_vl, "  Q_clover_imp_x_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, scr);
          scr = l_c_rect * corr_1x2[t];
          vout.general(m_vl, "  Q_clover_rect_x_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, scr);
          vout.general(m_vl, "  Q_plaq_x_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, corr_plaq[t]);
        }

        /*
        vout.general(m_vl, "  Q_clover_plaq_sumx_FT = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_1x1);
        double Q_rect = l_c_rect * Q_1x2;
        vout.general(m_vl, "  Q_clover_rect_sumx_FT = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_rect);
        double Q_imp = m_c_plaq * Q_1x1 + m_c_rect * Q_1x2;
        vout.general(m_vl, "  Q_clover_imp_sumx_FT = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_imp);
        vout.general(m_vl, "  Q_plaq_sumx = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_plaq);
        */

        if (m_filename_output != "stdout") {
          log_file.close();
          vout.init(log_file_previous);
        }
      }
    }
  }
}


//====================================================================
void TopologicalCharge::measure_topological_density_y(const double tt)
{
  assert(m_flag_field_set == true);

  static const double PI      = 4.0 * atan(1.0);
  static const double PI2     = PI * PI;
  static const double factor1 = 8.0 / (32.0 * PI2);
  static const double factor2 = 16.0 / (32.0 * PI2);

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  const int Ly = CommonParameters::Ly();

  std::vector<double> corr_scr(Ly);

  //--- Field strength by one plaquette ---
  std::vector<double> corr_plaq(Ly, 0.0);

  m_field_strength.contract_at_y(corr_scr, m_Fmunu_plaq[0], m_Fmunu_plaq[5]);
  for (int y = 0; y < Ly; ++y) {
    corr_plaq[y] += corr_scr[y];
  }
  m_field_strength.contract_at_y(corr_scr, m_Fmunu_plaq[1], m_Fmunu_plaq[4]);
  for (int y = 0; y < Ly; ++y) {
    corr_plaq[y] -= corr_scr[y];
  }
  m_field_strength.contract_at_y(corr_scr, m_Fmunu_plaq[2], m_Fmunu_plaq[3]);
  for (int y = 0; y < Ly; ++y) {
    corr_plaq[y] += corr_scr[y];
  }
  //  double Q_plaq = 0.0;
  for (int y = 0; y < Ly; ++y) {
    corr_plaq[y] *= factor1;
    //    Q_plaq += corr_plaq[y];
  }

  //--- 1x1 part ---
  std::vector<double> corr_1x1(Ly, 0.0);

  // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
  // (mu,nu, rho,sigma) = (1,2, 3,4)
  m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x1[0], m_Fmunu_1x1[5]);
  for (int y = 0; y < Ly; ++y) {
    corr_1x1[y] += corr_scr[y];
  }
  // (mu,nu, rho,sigma) = (1,3, 2,4)
  m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x1[1], m_Fmunu_1x1[4]);
  for (int y = 0; y < Ly; ++y) {
    corr_1x1[y] -= corr_scr[y];
  }
  // (mu,nu, rho,sigma) = (1,4, 2,3)
  m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x1[2], m_Fmunu_1x1[3]);
  for (int y = 0; y < Ly; ++y) {
    corr_1x1[y] += corr_scr[y];
  }
  //  double Q_1x1 = 0.0;
  for (int y = 0; y < Ly; ++y) {
    corr_1x1[y] *= factor1;
    //    Q_1x1 += corr_1x1[y];
  }

  //--- 1x2 part ---
  std::vector<double> corr_1x2(Ly, 0.0);

  // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
  m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x2[0], m_Fmunu_1x2[5]);
  for (int y = 0; y < Ly; ++y) {
    corr_1x2[y] += corr_scr[y];
  }
  m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x2[1], m_Fmunu_1x2[4]);
  for (int y = 0; y < Ly; ++y) {
    corr_1x2[y] -= corr_scr[y];
  }
  m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x2[2], m_Fmunu_1x2[3]);
  for (int y = 0; y < Ly; ++y) {
    corr_1x2[y] += corr_scr[y];
  }
  //  double Q_1x2 = 0.0;
  for (int y = 0; y < Ly; ++y) {
    corr_1x2[y] *= factor2;
    //    Q_1x2 += corr_1x2[y];
  }
  //----------------

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int y = 0; y < Ly; ++y) {
    vout.general(m_vl, "  Q_clover_plaq_y = %.8f %d %.16e\n", tt, y, corr_1x1[y]);

    double scr = m_c_plaq * corr_1x1[y] + m_c_rect * corr_1x2[y];
    vout.general(m_vl, "  Q_clover_imp_y = %.8f %d %.16e\n", tt, y, scr);

    scr = l_c_rect * corr_1x2[y];
    vout.general(m_vl, "  Q_clover_rect_y = %.8f %d %.16e\n", tt, y, scr);
    vout.general(m_vl, "  Q_plaq_y = %.8f %d %.16e\n", tt, y, corr_plaq[y]);
  }

  /*
  vout.general(m_vl, "  Q_clover_plaq_sumy = %.8f  %.16e\n", tt, Q_1x1);
  const double Q_rect = l_c_rect * Q_1x2;
  vout.general(m_vl, "  Q_clover_rect_sumy = %.8f  %.16e\n", tt, Q_rect);
  const double Q_imp = m_c_plaq * Q_1x1 + m_c_rect * Q_1x2;
  vout.general(m_vl, "  Q_clover_imp_sumy = %.8f  %.16e\n", tt, Q_imp);
  vout.general(m_vl, "  Q_plaq_sumy = %.8f  %.16e\n", tt, Q_plaq);
  */

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
}


//====================================================================
void TopologicalCharge::measure_topological_density_y_FT(const double tt)
{
  assert(m_flag_field_set == 1);
  static const double PI       = 4.0 * atan(1.0);
  static const double PI2      = PI * PI;
  static const double factor1  = 8.0 / (32.0 * PI2);
  static const double factor2  = 16.0 / (32.0 * PI2);
  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  const int Ly = CommonParameters::Ly();

  const int   Np = (2 * m_max_mom + 1);
  vector<int> source_position(4, 0);
  vector<int> momentum_sink(3);

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;

        std::vector<double> corr_plaq(Ly, 0);
        std::vector<double> corr_1x1(Ly, 0);
        std::vector<double> corr_1x2(Ly, 0);
        std::vector<double> corr_scr(Ly);

        //--- Field strength by one plaquette ---
        m_field_strength.contract_at_y(corr_scr, m_Fmunu_plaq[0], m_Fmunu_plaq[5], momentum_sink, source_position);
        for (int t = 0; t < Ly; ++t) {
          corr_plaq[t] += corr_scr[t];
        }
        m_field_strength.contract_at_y(corr_scr, m_Fmunu_plaq[1], m_Fmunu_plaq[4], momentum_sink, source_position);
        for (int t = 0; t < Ly; ++t) {
          corr_plaq[t] -= corr_scr[t];
        }
        m_field_strength.contract_at_y(corr_scr, m_Fmunu_plaq[2], m_Fmunu_plaq[3], momentum_sink, source_position);
        for (int t = 0; t < Ly; ++t) {
          corr_plaq[t] += corr_scr[t];
        }
        //	double Q_plaq = 0.0;
        for (int t = 0; t < Ly; ++t) {
          corr_plaq[t] *= factor1;
          //	  Q_plaq += corr_plaq[t];
        }

        //--- 1x1 part ---
        // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
        // (mu,nu, rho,sigma) = (1,2, 3,4)
        m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x1[0], m_Fmunu_1x1[5], momentum_sink, source_position);
        for (int t = 0; t < Ly; ++t) {
          corr_1x1[t] += corr_scr[t];
        }
        // (mu,nu, rho,sigma) = (1,3, 2,4)
        m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x1[1], m_Fmunu_1x1[4], momentum_sink, source_position);
        for (int t = 0; t < Ly; ++t) {
          corr_1x1[t] -= corr_scr[t];
        }
        // (mu,nu, rho,sigma) = (1,4, 2,3)
        m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x1[2], m_Fmunu_1x1[3], momentum_sink, source_position);
        for (int t = 0; t < Ly; ++t) {
          corr_1x1[t] += corr_scr[t];
        }
        //	double Q_1x1 = 0.0;
        for (int t = 0; t < Ly; ++t) {
          corr_1x1[t] *= factor1;
          //	  Q_1x1 += corr_1x1[t];
        }

        //--- 1x2 part ---
        // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
        m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x2[0], m_Fmunu_1x2[5], momentum_sink, source_position);
        for (int t = 0; t < Ly; ++t) {
          corr_1x2[t] += corr_scr[t];
        }
        m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x2[1], m_Fmunu_1x2[4], momentum_sink, source_position);
        for (int t = 0; t < Ly; ++t) {
          corr_1x2[t] -= corr_scr[t];
        }
        m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x2[2], m_Fmunu_1x2[3], momentum_sink, source_position);
        for (int t = 0; t < Ly; ++t) {
          corr_1x2[t] += corr_scr[t];
        }
        //	double Q_1x2 = 0.0;
        for (int t = 0; t < Ly; ++t) {
          corr_1x2[t] *= factor2;
          //	  Q_1x2 += corr_1x2[t];
        }
        //----------------

        //- output
        std::ostream& log_file_previous = vout.getStream();
        std::ofstream log_file;

        if (m_filename_output != "stdout") {
          log_file.open(m_filename_output.c_str(), std::ios::app);
          vout.init(log_file);
        }

        for (int t = 0; t < Ly; ++t) {
          vout.general(m_vl, "  Q_clover_plaq_y_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, corr_1x1[t]);
          double scr = m_c_plaq * corr_1x1[t] + m_c_rect * corr_1x2[t];
          vout.general(m_vl, "  Q_clover_imp_y_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, scr);
          scr = l_c_rect * corr_1x2[t];
          vout.general(m_vl, "  Q_clover_rect_y_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, scr);
          vout.general(m_vl, "  Q_plaq_y_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, corr_plaq[t]);
        }

        /*
        vout.general(m_vl, "  Q_clover_plaq_sumy_FT = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_1x1);
        double Q_rect = l_c_rect * Q_1x2;
        vout.general(m_vl, "  Q_clover_rect_sumy_FT = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_rect);
        double Q_imp = m_c_plaq * Q_1x1 + m_c_rect * Q_1x2;
        vout.general(m_vl, "  Q_clover_imp_sumy_FT = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_imp);
        vout.general(m_vl, "  Q_plaq_sumy = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_plaq);
        */

        if (m_filename_output != "stdout") {
          log_file.close();
          vout.init(log_file_previous);
        }
      }
    }
  }
}


//====================================================================
void TopologicalCharge::measure_topological_density_z(const double tt)
{
  assert(m_flag_field_set == true);

  static const double PI      = 4.0 * atan(1.0);
  static const double PI2     = PI * PI;
  static const double factor1 = 8.0 / (32.0 * PI2);
  static const double factor2 = 16.0 / (32.0 * PI2);

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  const int Lz = CommonParameters::Lz();

  std::vector<double> corr_scr(Lz);

  //--- Field strength by one plaquette ---
  std::vector<double> corr_plaq(Lz, 0.0);

  m_field_strength.contract_at_z(corr_scr, m_Fmunu_plaq[0], m_Fmunu_plaq[5]);
  for (int z = 0; z < Lz; ++z) {
    corr_plaq[z] += corr_scr[z];
  }
  m_field_strength.contract_at_z(corr_scr, m_Fmunu_plaq[1], m_Fmunu_plaq[4]);
  for (int z = 0; z < Lz; ++z) {
    corr_plaq[z] -= corr_scr[z];
  }
  m_field_strength.contract_at_z(corr_scr, m_Fmunu_plaq[2], m_Fmunu_plaq[3]);
  for (int z = 0; z < Lz; ++z) {
    corr_plaq[z] += corr_scr[z];
  }
  //  double Q_plaq = 0.0;
  for (int z = 0; z < Lz; ++z) {
    corr_plaq[z] *= factor1;
    //    Q_plaq += corr_plaq[z];
  }

  //--- 1x1 part ---
  std::vector<double> corr_1x1(Lz, 0.0);

  // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
  // (mu,nu, rho,sigma) = (1,2, 3,4)
  m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x1[0], m_Fmunu_1x1[5]);
  for (int z = 0; z < Lz; ++z) {
    corr_1x1[z] += corr_scr[z];
  }
  // (mu,nu, rho,sigma) = (1,3, 2,4)
  m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x1[1], m_Fmunu_1x1[4]);
  for (int z = 0; z < Lz; ++z) {
    corr_1x1[z] -= corr_scr[z];
  }
  // (mu,nu, rho,sigma) = (1,4, 2,3)
  m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x1[2], m_Fmunu_1x1[3]);
  for (int z = 0; z < Lz; ++z) {
    corr_1x1[z] += corr_scr[z];
  }
  //  double Q_1x1 = 0.0;
  for (int z = 0; z < Lz; ++z) {
    corr_1x1[z] *= factor1;
    //    Q_1x1 += corr_1x1[z];
  }

  //--- 1x2 part ---
  std::vector<double> corr_1x2(Lz, 0.0);

  // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
  m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x2[0], m_Fmunu_1x2[5]);
  for (int z = 0; z < Lz; ++z) {
    corr_1x2[z] += corr_scr[z];
  }
  m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x2[1], m_Fmunu_1x2[4]);
  for (int z = 0; z < Lz; ++z) {
    corr_1x2[z] -= corr_scr[z];
  }
  m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x2[2], m_Fmunu_1x2[3]);
  for (int z = 0; z < Lz; ++z) {
    corr_1x2[z] += corr_scr[z];
  }
  //  double Q_1x2 = 0.0;
  for (int z = 0; z < Lz; ++z) {
    corr_1x2[z] *= factor2;
    //    Q_1x2 += corr_1x2[z];
  }
  //----------------

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int z = 0; z < Lz; ++z) {
    vout.general(m_vl, "  Q_clover_plaq_z = %.8f %d %.16e\n", tt, z, corr_1x1[z]);

    double scr = m_c_plaq * corr_1x1[z] + m_c_rect * corr_1x2[z];
    vout.general(m_vl, "  Q_clover_imp_z = %.8f %d %.16e\n", tt, z, scr);

    scr = l_c_rect * corr_1x2[z];
    vout.general(m_vl, "  Q_clover_rect_z = %.8f %d %.16e\n", tt, z, scr);
    vout.general(m_vl, "  Q_plaq_z = %.8f %d %.16e\n", tt, z, corr_plaq[z]);
  }

  /*
  vout.general(m_vl, "  Q_clover_plaq_sumz = %.8f  %.16e\n", tt, Q_1x1);
  const double Q_rect = l_c_rect * Q_1x2;
  vout.general(m_vl, "  Q_clover_rect_sumz = %.8f  %.16e\n", tt, Q_rect);
  const double Q_imp = m_c_plaq * Q_1x1 + m_c_rect * Q_1x2;
  vout.general(m_vl, "  Q_clover_imp_sumz = %.8f  %.16e\n", tt, Q_imp);
  vout.general(m_vl, "  Q_plaq_sumz = %.8f  %.16e\n", tt, Q_plaq);
  */

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
}


//====================================================================
void TopologicalCharge::measure_topological_density_z_FT(const double tt)
{
  assert(m_flag_field_set == 1);
  static const double PI       = 4.0 * atan(1.0);
  static const double PI2      = PI * PI;
  static const double factor1  = 8.0 / (32.0 * PI2);
  static const double factor2  = 16.0 / (32.0 * PI2);
  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  const int Lz = CommonParameters::Lz();

  const int   Np = (2 * m_max_mom + 1);
  vector<int> source_position(4, 0);
  vector<int> momentum_sink(3);

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;

        std::vector<double> corr_plaq(Lz, 0);
        std::vector<double> corr_1x1(Lz, 0);
        std::vector<double> corr_1x2(Lz, 0);
        std::vector<double> corr_scr(Lz);

        //--- Field strength by one plaquette ---
        m_field_strength.contract_at_z(corr_scr, m_Fmunu_plaq[0], m_Fmunu_plaq[5], momentum_sink, source_position);
        for (int t = 0; t < Lz; ++t) {
          corr_plaq[t] += corr_scr[t];
        }
        m_field_strength.contract_at_z(corr_scr, m_Fmunu_plaq[1], m_Fmunu_plaq[4], momentum_sink, source_position);
        for (int t = 0; t < Lz; ++t) {
          corr_plaq[t] -= corr_scr[t];
        }
        m_field_strength.contract_at_z(corr_scr, m_Fmunu_plaq[2], m_Fmunu_plaq[3], momentum_sink, source_position);
        for (int t = 0; t < Lz; ++t) {
          corr_plaq[t] += corr_scr[t];
        }
        //	double Q_plaq = 0.0;
        for (int t = 0; t < Lz; ++t) {
          corr_plaq[t] *= factor1;
          //	  Q_plaq += corr_plaq[t];
        }

        //--- 1x1 part ---
        // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
        // (mu,nu, rho,sigma) = (1,2, 3,4)
        m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x1[0], m_Fmunu_1x1[5], momentum_sink, source_position);
        for (int t = 0; t < Lz; ++t) {
          corr_1x1[t] += corr_scr[t];
        }
        // (mu,nu, rho,sigma) = (1,3, 2,4)
        m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x1[1], m_Fmunu_1x1[4], momentum_sink, source_position);
        for (int t = 0; t < Lz; ++t) {
          corr_1x1[t] -= corr_scr[t];
        }
        // (mu,nu, rho,sigma) = (1,4, 2,3)
        m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x1[2], m_Fmunu_1x1[3], momentum_sink, source_position);
        for (int t = 0; t < Lz; ++t) {
          corr_1x1[t] += corr_scr[t];
        }
        //	double Q_1x1 = 0.0;
        for (int t = 0; t < Lz; ++t) {
          corr_1x1[t] *= factor1;
          //	  Q_1x1 += corr_1x1[t];
        }

        //--- 1x2 part ---
        // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
        m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x2[0], m_Fmunu_1x2[5], momentum_sink, source_position);
        for (int t = 0; t < Lz; ++t) {
          corr_1x2[t] += corr_scr[t];
        }
        m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x2[1], m_Fmunu_1x2[4], momentum_sink, source_position);
        for (int t = 0; t < Lz; ++t) {
          corr_1x2[t] -= corr_scr[t];
        }
        m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x2[2], m_Fmunu_1x2[3], momentum_sink, source_position);
        for (int t = 0; t < Lz; ++t) {
          corr_1x2[t] += corr_scr[t];
        }
        //	double Q_1x2 = 0.0;
        for (int t = 0; t < Lz; ++t) {
          corr_1x2[t] *= factor2;
          //	  Q_1x2 += corr_1x2[t];
        }
        //----------------

        //- output
        std::ostream& log_file_previous = vout.getStream();
        std::ofstream log_file;

        if (m_filename_output != "stdout") {
          log_file.open(m_filename_output.c_str(), std::ios::app);
          vout.init(log_file);
        }

        for (int t = 0; t < Lz; ++t) {
          vout.general(m_vl, "  Q_clover_plaq_z_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, corr_1x1[t]);
          double scr = m_c_plaq * corr_1x1[t] + m_c_rect * corr_1x2[t];
          vout.general(m_vl, "  Q_clover_imp_z_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, scr);
          scr = l_c_rect * corr_1x2[t];
          vout.general(m_vl, "  Q_clover_rect_z_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, scr);
          vout.general(m_vl, "  Q_plaq_z_FT = %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, t, corr_plaq[t]);
        }

        /*
        vout.general(m_vl, "  Q_clover_plaq_sumz_FT = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_1x1);
        double Q_rect = l_c_rect * Q_1x2;
        vout.general(m_vl, "  Q_clover_rect_sumz_FT = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_rect);
        double Q_imp = m_c_plaq * Q_1x1 + m_c_rect * Q_1x2;
        vout.general(m_vl, "  Q_clover_imp_sumz_FT = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_imp);
        vout.general(m_vl, "  Q_plaq_sumz = %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], tt, Q_plaq);
        */

        if (m_filename_output != "stdout") {
          log_file.close();
          vout.init(log_file_previous);
        }
      }
    }
  }
}


//====================================================================

/*
double TopologicalCharge::contract_epsilon_tensor(const Field_G& Fmunu_1, const Field_G& Fmunu_2)
{
  const int Nvol = CommonParameters::Nvol();

  double Q_topo = 0.0;
  for (int site = 0; site < Nvol; ++site) {
    Q_topo += ReTr(Fmunu_1.mat(site) * Fmunu_2.mat(site));
  }
  const double result = Communicator::reduce_sum(Q_topo);

  return result;
}
*/


//====================================================================
void TopologicalCharge::set_field_strength(const Field_G& U)
{
  const int Ndim = CommonParameters::Ndim();

  // NB. #(mu,nu)=6 i.e. (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
  m_Fmunu_plaq.resize(6);
  m_Fmunu_1x1.resize(6);
  m_Fmunu_1x2.resize(6);

  int i_munu = 0;
  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = mu + 1; nu < Ndim; ++nu) {
      m_field_strength.construct_Fmunu_plaq_traceless(m_Fmunu_plaq[i_munu], mu, nu, U);
      m_field_strength.construct_Fmunu_1x1_traceless(m_Fmunu_1x1[i_munu], mu, nu, U);
      m_field_strength.construct_Fmunu_1x2_traceless(m_Fmunu_1x2[i_munu], mu, nu, U);
      ++i_munu;
    }
  }
  m_flag_field_set = true;
}


//====================================================================
//============================================================END=====
