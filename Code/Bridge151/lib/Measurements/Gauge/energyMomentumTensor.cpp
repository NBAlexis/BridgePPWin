/*!
        @file    energyMomentumTensor.cpp

        @brief

        @author  Yusuke Taniguchi  (tanigchi)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "energyMomentumTensor.h"

const std::string EnergyMomentumTensor::class_name = "EnergyMomentumTensor";

//====================================================================
void EnergyMomentumTensor::set_parameters(const Parameters& params)
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
void EnergyMomentumTensor::set_parameters(const double c_plaq, const double c_rect, const int max_mom)
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
double EnergyMomentumTensor::measure_EMT(const double t_flow)
{
  assert(m_flag_field_set == 1);

  const int Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  double O2_1x1  = 0.0;
  double O2_1x2  = 0.0;
  double O2_plaq = 0.0;
  for (int mu = 0; mu < 6; ++mu) {
    O2_1x1  += m_field_strength.contract(m_Fmunu_1x1[mu], m_Fmunu_1x1[mu]);
    O2_1x2  += m_field_strength.contract(m_Fmunu_1x2[mu], m_Fmunu_1x2[mu]);
    O2_plaq += m_field_strength.contract(m_Fmunu_plaq[mu], m_Fmunu_plaq[mu]);
  }
  O2_1x1 *= 4.0 / Nvol / NPE;
  vout.general(m_vl, "  O2_clover_plaq = %.8f %.16e\n", t_flow, O2_1x1);
  O2_1x2 *= 8.0 / Nvol / NPE;
  const double O2_imp = (m_c_plaq * O2_1x1 + m_c_rect * O2_1x2);
  vout.general(m_vl, "  O2_clover_imp = %.8f %.16e\n", t_flow, O2_imp);
  const double O2_rect = l_c_rect * O2_1x2;
  vout.general(m_vl, "  O2_clover_rect = %.8f %.16e\n", t_flow, O2_rect);
  O2_plaq *= 4.0 / Nvol / NPE;
  vout.general(m_vl, "  O2_plaq = %.8f %.16e\n", t_flow, O2_plaq);

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = 0; nu < Ndim; ++nu) {
      double O1_1x1  = 0.0;
      double O1_1x2  = 0.0;
      double O1_plaq = 0.0;
      for (int rho = 0; rho < Ndim; rho++) {
        if ((mu != rho) && (nu != rho)) {
          int ind1 = index_munu2i(mu, rho);
          int ind2 = index_munu2i(nu, rho);
          int fac1 = factor(mu, rho);
          int fac2 = factor(nu, rho);

          double scr = m_field_strength.contract(m_Fmunu_1x1[ind1], m_Fmunu_1x1[ind2]);
          scr    *= fac1 * fac2;
          O1_1x1 += scr;

          scr     = m_field_strength.contract(m_Fmunu_1x2[ind1], m_Fmunu_1x2[ind2]);
          scr    *= fac1 * fac2;
          O1_1x2 += scr;

          scr      = m_field_strength.contract(m_Fmunu_plaq[ind1], m_Fmunu_plaq[ind2]);
          scr     *= fac1 * fac2;
          O1_plaq += scr;
        }
      }
      O1_1x1  *= 2.0 / Nvol / NPE;
      O1_1x2  *= 4.0 / Nvol / NPE;
      O1_plaq *= 2.0 / Nvol / NPE;
      vout.general(m_vl, "  O1_clover_plaq = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_1x1);
      double O1_rect = l_c_rect * O1_1x2;
      vout.general(m_vl, "  O1_clover_1x2 = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_1x2);
      vout.general(m_vl, "  O1_clover_rect = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_rect);
      double O1_imp = (m_c_plaq * O1_1x1 + m_c_rect * O1_1x2);
      vout.general(m_vl, "  O1_clover_imp = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_imp);
      vout.general(m_vl, "  O1_plaq = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_plaq);
    }
  }

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  return O2_imp;
}


//====================================================================
double EnergyMomentumTensor::measure_EMT_at_t(const double t_flow)
{
  assert(m_flag_field_set == 1);

  const int Lt = CommonParameters::Lt();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int Ndim = CommonParameters::Ndim();

  const double normalization = double(Lt) / Nvol / NPE;
  double       result        = 0.0; // dummy initialization

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = 0; nu < Ndim; ++nu) {
      std::vector<double> corr_plaq(Lt, 0);
      std::vector<double> corr_1x1(Lt, 0);
      std::vector<double> corr_1x2(Lt, 0);
      std::vector<double> corr_scr(Lt);
      for (int rho = 0; rho < Ndim; rho++) {
        if ((mu != rho) && (nu != rho)) {
          int ind1 = index_munu2i(mu, rho);
          int ind2 = index_munu2i(nu, rho);
          int fac1 = factor(mu, rho);
          int fac2 = factor(nu, rho);
          m_field_strength.contract_at_t(corr_scr, m_Fmunu_plaq[ind1], m_Fmunu_plaq[ind2]);
          for (int t = 0; t < Lt; ++t) {
            corr_plaq[t] += (fac1 * fac2 * corr_scr[t]);
          }
          m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x1[ind1], m_Fmunu_1x1[ind2]);
          for (int t = 0; t < Lt; ++t) {
            corr_1x1[t] += (fac1 * fac2 * corr_scr[t]);
          }
          m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x2[ind1], m_Fmunu_1x2[ind2]);
          for (int t = 0; t < Lt; ++t) {
            corr_1x2[t] += (fac1 * fac2 * corr_scr[t]);
          }
        }
      }
      //double O1_plaq = 0.0;
      //double O1_1x1 = 0.0;
      //double O1_1x2 = 0.0;
      for (int t = 0; t < Lt; ++t) {
        corr_plaq[t] *= 2 * normalization;
        //O1_plaq += corr_plaq[t];
        corr_1x1[t] *= 2 * normalization;
        //O1_1x1 += corr_1x1[t];
        corr_1x2[t] *= 4 * normalization;
        //O1_1x2 += corr_1x2[t];

        vout.general(m_vl, "  O1_clover_plaq_t = %d %d %.8f %d %.16e\n", mu, nu, t_flow, t, corr_1x1[t]);
        double O1_rect = l_c_rect * corr_1x2[t];
        vout.general(m_vl, "  O1_clover_rect_t = %d %d %.8f %d %.16e\n", mu, nu, t_flow, t, O1_rect);
        double O1_imp = (m_c_plaq * corr_1x1[t] + m_c_rect * corr_1x2[t]);
        vout.general(m_vl, "  O1_clover_imp_t = %d %d %.8f %d %.16e\n", mu, nu, t_flow, t, O1_imp);
        result = O1_imp;
        vout.general(m_vl, "  O1_plaq_t = %d %d %.8f %d %.16e\n", mu, nu, t_flow, t, corr_plaq[t]);
      }

      /*
      double O1_imp = m_c_plaq * O1_1x1 + m_c_rect * O1_1x2;
      double O1_rect = l_c_rect * O1_1x2;
      vout.general(m_vl, "  O1_clover_plaq_sumt = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_1x1/Lt);
      //vout.general(m_vl, "  O1_clover_1x2_sumt = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_1x2/Lt);
      vout.general(m_vl, "  O1_clover_rect_sumt = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_rect/Lt);
      vout.general(m_vl, "  O1_clover_imp_sumt = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_imp/Lt);
      vout.general(m_vl, "  O1_plaq_sumt = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_plaq/Lt);
      */
    }
  }

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  return result;
}


//====================================================================
double EnergyMomentumTensor::measure_EMT_at_t_FT(const double t_flow)
{
  assert(m_flag_field_set == 1);

  const int Lt = CommonParameters::Lt();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int Ndim = CommonParameters::Ndim();

  const double normalization = double(Lt) / Nvol / NPE;
  double       result        = 0.0; // dummy initialization

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  const int   Np = (2 * m_max_mom + 1);
  vector<int> source_position(4, 0);
  vector<int> momentum_sink(3);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        for (int mu = 0; mu < Ndim; ++mu) {
          for (int nu = 0; nu < Ndim; ++nu) {
            std::vector<double> corr_plaq(Lt, 0);
            std::vector<double> corr_1x1(Lt, 0);
            std::vector<double> corr_1x2(Lt, 0);
            std::vector<double> corr_scr(Lt);
            for (int rho = 0; rho < Ndim; rho++) {
              if ((mu != rho) && (nu != rho)) {
                int ind1 = index_munu2i(mu, rho);
                int ind2 = index_munu2i(nu, rho);
                int fac1 = factor(mu, rho);
                int fac2 = factor(nu, rho);
                m_field_strength.contract_at_t(corr_scr, m_Fmunu_plaq[ind1], m_Fmunu_plaq[ind2], momentum_sink, source_position);
                for (int t = 0; t < Lt; ++t) {
                  corr_plaq[t] += (fac1 * fac2 * corr_scr[t]);
                }
                m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x1[ind1], m_Fmunu_1x1[ind2], momentum_sink, source_position);
                for (int t = 0; t < Lt; ++t) {
                  corr_1x1[t] += (fac1 * fac2 * corr_scr[t]);
                }
                m_field_strength.contract_at_t(corr_scr, m_Fmunu_1x2[ind1], m_Fmunu_1x2[ind2], momentum_sink, source_position);
                for (int t = 0; t < Lt; ++t) {
                  corr_1x2[t] += (fac1 * fac2 * corr_scr[t]);
                }
              }
            }
            //double O1_plaq = 0.0;
            //	    double O1_1x1 = 0.0;
            //	    double O1_1x2 = 0.0;
            for (int t = 0; t < Lt; ++t) {
              corr_plaq[t] *= 2 * normalization;
              //O1_plaq += corr_plaq[t];
              corr_1x1[t] *= 2 * normalization;
              //O1_1x1 += corr_1x1[t];
              corr_1x2[t] *= 4 * normalization;
              //O1_1x2 += corr_1x2[t];
              vout.general(m_vl, "  O1_clover_plaq_t_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, corr_1x1[t]);
              double O1_rect = l_c_rect * corr_1x2[t];
              vout.general(m_vl, "  O1_clover_rect_t_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, O1_rect);
              double O1_imp = (m_c_plaq * corr_1x1[t] + m_c_rect * corr_1x2[t]);
              vout.general(m_vl, "  O1_clover_imp_t_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, O1_imp);
              result = O1_imp;
              vout.general(m_vl, "  O1_plaq_t_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, corr_plaq[t]);
            }

            /*
            double O1_imp = m_c_plaq * O1_1x1 + m_c_rect * O1_1x2;
            double O1_rect = l_c_rect * O1_1x2;
            vout.general(m_vl, "  O1_clover_plaq_sumt_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_1x1/Lt);
            //vout.general(m_vl, "  O1_clover_1x2_sumt_FT = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_1x2/Lt);
            vout.general(m_vl, "  O1_clover_rect_sumt_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_rect/Lt);
            vout.general(m_vl, "  O1_clover_imp_sumt_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_imp/Lt);
            vout.general(m_vl, "  O1_plaq_sumt_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_plaq/Lt);
            */
          }
        }
      }
    }
  }

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  return result;
}


//====================================================================
double EnergyMomentumTensor::measure_EMT_at_x(const double t_flow)
{
  assert(m_flag_field_set == 1);

  const int Lx = CommonParameters::Lx();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int Ndim = CommonParameters::Ndim();

  const double normalization = double(Lx) / Nvol / NPE;
  double       result        = 0.0; // dummy initialization

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = 0; nu < Ndim; ++nu) {
      std::vector<double> corr_plaq(Lx, 0);
      std::vector<double> corr_1x1(Lx, 0);
      std::vector<double> corr_1x2(Lx, 0);
      std::vector<double> corr_scr(Lx);
      for (int rho = 0; rho < Ndim; rho++) {
        if ((mu != rho) && (nu != rho)) {
          int ind1 = index_munu2i(mu, rho);
          int ind2 = index_munu2i(nu, rho);
          int fac1 = factor(mu, rho);
          int fac2 = factor(nu, rho);
          m_field_strength.contract_at_x(corr_scr, m_Fmunu_plaq[ind1], m_Fmunu_plaq[ind2]);
          for (int x = 0; x < Lx; ++x) {
            corr_plaq[x] += (fac1 * fac2 * corr_scr[x]);
          }
          m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x1[ind1], m_Fmunu_1x1[ind2]);
          for (int x = 0; x < Lx; ++x) {
            corr_1x1[x] += (fac1 * fac2 * corr_scr[x]);
          }
          m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x2[ind1], m_Fmunu_1x2[ind2]);
          for (int x = 0; x < Lx; ++x) {
            corr_1x2[x] += (fac1 * fac2 * corr_scr[x]);
          }
        }
      }
      //double O1_plaq = 0.0;
      //double O1_1x1 = 0.0;
      //double O1_1x2 = 0.0;
      for (int x = 0; x < Lx; ++x) {
        corr_plaq[x] *= 2 * normalization;
        //O1_plaq += corr_plaq[x];
        corr_1x1[x] *= 2 * normalization;
        //O1_1x1 += corr_1x1[x];
        corr_1x2[x] *= 4 * normalization;
        //O1_1x2 += corr_1x2[x];

        vout.general(m_vl, "  O1_clover_plaq_x = %d %d %.8f %d %.16e\n", mu, nu, t_flow, x, corr_1x1[x]);
        double O1_rect = l_c_rect * corr_1x2[x];
        vout.general(m_vl, "  O1_clover_rect_x = %d %d %.8f %d %.16e\n", mu, nu, t_flow, x, O1_rect);
        double O1_imp = (m_c_plaq * corr_1x1[x] + m_c_rect * corr_1x2[x]);
        vout.general(m_vl, "  O1_clover_imp_x = %d %d %.8f %d %.16e\n", mu, nu, t_flow, x, O1_imp);
        result = O1_imp;
        vout.general(m_vl, "  O1_plaq_x = %d %d %.8f %d %.16e\n", mu, nu, t_flow, x, corr_plaq[x]);
      }

      /*
      double O1_imp = m_c_plaq * O1_1x1 + m_c_rect * O1_1x2;
      double O1_rect = l_c_rect * O1_1x2;
      vout.general(m_vl, "  O1_clover_plaq_sumx = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_1x1/Lx);
      vout.general(m_vl, "  O1_clover_rect_sumx = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_rect/Lx);
      vout.general(m_vl, "  O1_clover_imp_sumx = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_imp/Lx);
      vout.general(m_vl, "  O1_plaq_sumx = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_plaq/Lx);
      */
    }
  }

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  return result;
}


//====================================================================
double EnergyMomentumTensor::measure_EMT_at_x_FT(const double t_flow)
{
  assert(m_flag_field_set == 1);

  const int Lx = CommonParameters::Lx();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int Ndim = CommonParameters::Ndim();

  const double normalization = double(Lx) / Nvol / NPE;
  double       result        = 0.0; // dummy initialization

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  const int   Np = (2 * m_max_mom + 1);
  vector<int> source_position(4, 0);
  vector<int> momentum_sink(3);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        for (int mu = 0; mu < Ndim; ++mu) {
          for (int nu = 0; nu < Ndim; ++nu) {
            std::vector<double> corr_plaq(Lx, 0);
            std::vector<double> corr_1x1(Lx, 0);
            std::vector<double> corr_1x2(Lx, 0);
            std::vector<double> corr_scr(Lx);
            for (int rho = 0; rho < Ndim; rho++) {
              if ((mu != rho) && (nu != rho)) {
                int ind1 = index_munu2i(mu, rho);
                int ind2 = index_munu2i(nu, rho);
                int fac1 = factor(mu, rho);
                int fac2 = factor(nu, rho);
                m_field_strength.contract_at_x(corr_scr, m_Fmunu_plaq[ind1], m_Fmunu_plaq[ind2], momentum_sink, source_position);
                for (int t = 0; t < Lx; ++t) {
                  corr_plaq[t] += (fac1 * fac2 * corr_scr[t]);
                }
                m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x1[ind1], m_Fmunu_1x1[ind2], momentum_sink, source_position);
                for (int t = 0; t < Lx; ++t) {
                  corr_1x1[t] += (fac1 * fac2 * corr_scr[t]);
                }
                m_field_strength.contract_at_x(corr_scr, m_Fmunu_1x2[ind1], m_Fmunu_1x2[ind2], momentum_sink, source_position);
                for (int t = 0; t < Lx; ++t) {
                  corr_1x2[t] += (fac1 * fac2 * corr_scr[t]);
                }
              }
            }
            //double O1_plaq = 0.0;
            //double O1_1x1 = 0.0;
            //double O1_1x2 = 0.0;
            for (int t = 0; t < Lx; ++t) {
              corr_plaq[t] *= 2 * normalization;
              //O1_plaq += corr_plaq[t];
              corr_1x1[t] *= 2 * normalization;
              //O1_1x1 += corr_1x1[t];
              corr_1x2[t] *= 4 * normalization;
              //O1_1x2 += corr_1x2[t];
              vout.general(m_vl, "  O1_clover_plaq_x_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, corr_1x1[t]);
              double O1_rect = l_c_rect * corr_1x2[t];
              vout.general(m_vl, "  O1_clover_rect_x_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, O1_rect);
              double O1_imp = (m_c_plaq * corr_1x1[t] + m_c_rect * corr_1x2[t]);
              vout.general(m_vl, "  O1_clover_imp_x_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, O1_imp);
              result = O1_imp;
              vout.general(m_vl, "  O1_plaq_x_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, corr_plaq[t]);
            }

            /*
            double O1_imp = m_c_plaq * O1_1x1 + m_c_rect * O1_1x2;
            double O1_rect = l_c_rect * O1_1x2;
            vout.general(m_vl, "  O1_clover_plaq_sumx_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_1x1/Lx);
            //vout.general(m_vl, "  O1_clover_1x2_sumx_FT = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_1x2/Lx);
            vout.general(m_vl, "  O1_clover_rect_sumx_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_rect/Lx);
            vout.general(m_vl, "  O1_clover_imp_sumx_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_imp/Lx);
            vout.general(m_vl, "  O1_plaq_sumx_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_plaq/Lx);
            */
          }
        }
      }
    }
  }

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  return result;
}


//====================================================================
double EnergyMomentumTensor::measure_EMT_at_y(const double t_flow)
{
  assert(m_flag_field_set == 1);

  const int Ly = CommonParameters::Ly();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int Ndim = CommonParameters::Ndim();

  const double normalization = double(Ly) / Nvol / NPE;
  double       result        = 0.0; // dummy initialization

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = 0; nu < Ndim; ++nu) {
      std::vector<double> corr_plaq(Ly, 0);
      std::vector<double> corr_1x1(Ly, 0);
      std::vector<double> corr_1x2(Ly, 0);
      std::vector<double> corr_scr(Ly);
      for (int rho = 0; rho < Ndim; rho++) {
        if ((mu != rho) && (nu != rho)) {
          int ind1 = index_munu2i(mu, rho);
          int ind2 = index_munu2i(nu, rho);
          int fac1 = factor(mu, rho);
          int fac2 = factor(nu, rho);
          m_field_strength.contract_at_y(corr_scr, m_Fmunu_plaq[ind1], m_Fmunu_plaq[ind2]);
          for (int y = 0; y < Ly; ++y) {
            corr_plaq[y] += (fac1 * fac2 * corr_scr[y]);
          }
          m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x1[ind1], m_Fmunu_1x1[ind2]);
          for (int y = 0; y < Ly; ++y) {
            corr_1x1[y] += (fac1 * fac2 * corr_scr[y]);
          }
          m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x2[ind1], m_Fmunu_1x2[ind2]);
          for (int y = 0; y < Ly; ++y) {
            corr_1x2[y] += (fac1 * fac2 * corr_scr[y]);
          }
        }
      }
      //double O1_plaq = 0.0;
      //double O1_1x1 = 0.0;
      //double O1_1x2 = 0.0;
      for (int y = 0; y < Ly; ++y) {
        corr_plaq[y] *= 2 * normalization;
        //O1_plaq += corr_plaq[y];
        corr_1x1[y] *= 2 * normalization;
        //O1_1x1 += corr_1x1[y];
        corr_1x2[y] *= 4 * normalization;
        //O1_1x2 += corr_1x2[y];

        vout.general(m_vl, "  O1_clover_plaq_y = %d %d %.8f %d %.16e\n", mu, nu, t_flow, y, corr_1x1[y]);
        double O1_rect = l_c_rect * corr_1x2[y];
        vout.general(m_vl, "  O1_clover_rect_y = %d %d %.8f %d %.16e\n", mu, nu, t_flow, y, O1_rect);
        double O1_imp = (m_c_plaq * corr_1x1[y] + m_c_rect * corr_1x2[y]);
        vout.general(m_vl, "  O1_clover_imp_y = %d %d %.8f %d %.16e\n", mu, nu, t_flow, y, O1_imp);
        result = O1_imp;
        vout.general(m_vl, "  O1_plaq_y = %d %d %.8f %d %.16e\n", mu, nu, t_flow, y, corr_plaq[y]);
      }

      /*
      double O1_imp = m_c_plaq * O1_1x1 + m_c_rect * O1_1x2;
      double O1_rect = l_c_rect * O1_1x2;
      vout.general(m_vl, "  O1_clover_plaq_sumy = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_1x1/Ly);
      vout.general(m_vl, "  O1_clover_rect_sumy = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_rect/Ly);
      vout.general(m_vl, "  O1_clover_imp_sumy = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_imp/Ly);
      vout.general(m_vl, "  O1_plaq_sumy = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_plaq/Ly);
      */
    }
  }

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  return result;
}


//====================================================================
double EnergyMomentumTensor::measure_EMT_at_y_FT(const double t_flow)
{
  assert(m_flag_field_set == 1);

  const int Ly = CommonParameters::Ly();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int Ndim = CommonParameters::Ndim();

  const double normalization = double(Ly) / Nvol / NPE;
  double       result        = 0.0; // dummy initialization

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  const int   Np = (2 * m_max_mom + 1);
  vector<int> source_position(4, 0);
  vector<int> momentum_sink(3);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        for (int mu = 0; mu < Ndim; ++mu) {
          for (int nu = 0; nu < Ndim; ++nu) {
            std::vector<double> corr_plaq(Ly, 0);
            std::vector<double> corr_1x1(Ly, 0);
            std::vector<double> corr_1x2(Ly, 0);
            std::vector<double> corr_scr(Ly);
            for (int rho = 0; rho < Ndim; rho++) {
              if ((mu != rho) && (nu != rho)) {
                int ind1 = index_munu2i(mu, rho);
                int ind2 = index_munu2i(nu, rho);
                int fac1 = factor(mu, rho);
                int fac2 = factor(nu, rho);
                m_field_strength.contract_at_y(corr_scr, m_Fmunu_plaq[ind1], m_Fmunu_plaq[ind2], momentum_sink, source_position);
                for (int t = 0; t < Ly; ++t) {
                  corr_plaq[t] += (fac1 * fac2 * corr_scr[t]);
                }
                m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x1[ind1], m_Fmunu_1x1[ind2], momentum_sink, source_position);
                for (int t = 0; t < Ly; ++t) {
                  corr_1x1[t] += (fac1 * fac2 * corr_scr[t]);
                }
                m_field_strength.contract_at_y(corr_scr, m_Fmunu_1x2[ind1], m_Fmunu_1x2[ind2], momentum_sink, source_position);
                for (int t = 0; t < Ly; ++t) {
                  corr_1x2[t] += (fac1 * fac2 * corr_scr[t]);
                }
              }
            }
            //double O1_plaq = 0.0;
            //double O1_1x1 = 0.0;
            //double O1_1x2 = 0.0;
            for (int t = 0; t < Ly; ++t) {
              corr_plaq[t] *= 2 * normalization;
              //O1_plaq += corr_plaq[t];
              corr_1x1[t] *= 2 * normalization;
              //O1_1x1 += corr_1x1[t];
              corr_1x2[t] *= 4 * normalization;
              //O1_1x2 += corr_1x2[t];
              vout.general(m_vl, "  O1_clover_plaq_y_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, corr_1x1[t]);
              double O1_rect = l_c_rect * corr_1x2[t];
              vout.general(m_vl, "  O1_clover_rect_y_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, O1_rect);
              double O1_imp = (m_c_plaq * corr_1x1[t] + m_c_rect * corr_1x2[t]);
              vout.general(m_vl, "  O1_clover_imp_y_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, O1_imp);
              result = O1_imp;
              vout.general(m_vl, "  O1_plaq_y_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, corr_plaq[t]);
            }

            /*
            double O1_imp = m_c_plaq * O1_1x1 + m_c_rect * O1_1x2;
            double O1_rect = l_c_rect * O1_1x2;
            vout.general(m_vl, "  O1_clover_plaq_sumy_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_1x1/Ly);
            //vout.general(m_vl, "  O1_clover_1x2_sumy_FT = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_1x2/Ly);
            vout.general(m_vl, "  O1_clover_rect_sumy_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_rect/Ly);
            vout.general(m_vl, "  O1_clover_imp_sumy_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_imp/Ly);
            vout.general(m_vl, "  O1_plaq_sumy_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_plaq/Ly);
            */
          }
        }
      }
    }
  }

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  return result;
}


//====================================================================
double EnergyMomentumTensor::measure_EMT_at_z(const double t_flow)
{
  assert(m_flag_field_set == 1);

  const int Lz = CommonParameters::Lz();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int Ndim = CommonParameters::Ndim();

  const double normalization = double(Lz) / Nvol / NPE;
  double       result        = 0.0; // dummy initialization

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = 0; nu < Ndim; ++nu) {
      std::vector<double> corr_plaq(Lz, 0);
      std::vector<double> corr_1x1(Lz, 0);
      std::vector<double> corr_1x2(Lz, 0);
      std::vector<double> corr_scr(Lz);
      for (int rho = 0; rho < Ndim; rho++) {
        if ((mu != rho) && (nu != rho)) {
          int ind1 = index_munu2i(mu, rho);
          int ind2 = index_munu2i(nu, rho);
          int fac1 = factor(mu, rho);
          int fac2 = factor(nu, rho);
          m_field_strength.contract_at_z(corr_scr, m_Fmunu_plaq[ind1], m_Fmunu_plaq[ind2]);
          for (int z = 0; z < Lz; ++z) {
            corr_plaq[z] += (fac1 * fac2 * corr_scr[z]);
          }
          m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x1[ind1], m_Fmunu_1x1[ind2]);
          for (int z = 0; z < Lz; ++z) {
            corr_1x1[z] += (fac1 * fac2 * corr_scr[z]);
          }
          m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x2[ind1], m_Fmunu_1x2[ind2]);
          for (int z = 0; z < Lz; ++z) {
            corr_1x2[z] += (fac1 * fac2 * corr_scr[z]);
          }
        }
      }
      //double O1_plaq = 0.0;
      //double O1_1x1 = 0.0;
      //double O1_1x2 = 0.0;
      for (int z = 0; z < Lz; ++z) {
        corr_plaq[z] *= 2 * normalization;
        //O1_plaq += corr_plaq[z];
        corr_1x1[z] *= 2 * normalization;
        //O1_1x1 += corr_1x1[z];
        corr_1x2[z] *= 4 * normalization;
        //O1_1x2 += corr_1x2[z];

        vout.general(m_vl, "  O1_clover_plaq_z = %d %d %.8f %d %.16e\n", mu, nu, t_flow, z, corr_1x1[z]);
        double O1_rect = l_c_rect * corr_1x2[z];
        vout.general(m_vl, "  O1_clover_rect_z = %d %d %.8f %d %.16e\n", mu, nu, t_flow, z, O1_rect);
        double O1_imp = (m_c_plaq * corr_1x1[z] + m_c_rect * corr_1x2[z]);
        vout.general(m_vl, "  O1_clover_imp_z = %d %d %.8f %d %.16e\n", mu, nu, t_flow, z, O1_imp);
        result = O1_imp;
        vout.general(m_vl, "  O1_plaq_z = %d %d %.8f %d %.16e\n", mu, nu, t_flow, z, corr_plaq[z]);
      }

      /*
      double O1_imp = m_c_plaq * O1_1x1 + m_c_rect * O1_1x2;
      double O1_rect = l_c_rect * O1_1x2;
      vout.general(m_vl, "  O1_clover_plaq_sumz = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_1x1/Lz);
      vout.general(m_vl, "  O1_clover_rect_sumz = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_rect/Lz);
      vout.general(m_vl, "  O1_clover_imp_sumz = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_imp/Lz);
      vout.general(m_vl, "  O1_plaq_sumz = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_plaq/Lz);
      */
    }
  }

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  return result;
}


//====================================================================
double EnergyMomentumTensor::measure_EMT_at_z_FT(const double t_flow)
{
  assert(m_flag_field_set == 1);

  const int Lz = CommonParameters::Lz();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int Ndim = CommonParameters::Ndim();

  const double normalization = double(Lz) / Nvol / NPE;
  double       result        = 0.0; // dummy initialization

  static const double l_c_rect = 1.0 / 8.0;
  m_c_rect = -1.0 / 12.0;
  m_c_plaq = 1.0 - 8 * m_c_rect;

  const int   Np = (2 * m_max_mom + 1);
  vector<int> source_position(4, 0);
  vector<int> momentum_sink(3);

  //- output
  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  for (int ipx = 0; ipx < m_max_mom + 1; ipx++) {
    for (int ipy = 0; ipy < Np; ipy++) {
      for (int ipz = 0; ipz < Np; ipz++) {
        momentum_sink[0] = ipx;
        momentum_sink[1] = ipy - m_max_mom;
        momentum_sink[2] = ipz - m_max_mom;
        for (int mu = 0; mu < Ndim; ++mu) {
          for (int nu = 0; nu < Ndim; ++nu) {
            std::vector<double> corr_plaq(Lz, 0);
            std::vector<double> corr_1x1(Lz, 0);
            std::vector<double> corr_1x2(Lz, 0);
            std::vector<double> corr_scr(Lz);
            for (int rho = 0; rho < Ndim; rho++) {
              if ((mu != rho) && (nu != rho)) {
                int ind1 = index_munu2i(mu, rho);
                int ind2 = index_munu2i(nu, rho);
                int fac1 = factor(mu, rho);
                int fac2 = factor(nu, rho);
                m_field_strength.contract_at_z(corr_scr, m_Fmunu_plaq[ind1], m_Fmunu_plaq[ind2], momentum_sink, source_position);
                for (int t = 0; t < Lz; ++t) {
                  corr_plaq[t] += (fac1 * fac2 * corr_scr[t]);
                }
                m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x1[ind1], m_Fmunu_1x1[ind2], momentum_sink, source_position);
                for (int t = 0; t < Lz; ++t) {
                  corr_1x1[t] += (fac1 * fac2 * corr_scr[t]);
                }
                m_field_strength.contract_at_z(corr_scr, m_Fmunu_1x2[ind1], m_Fmunu_1x2[ind2], momentum_sink, source_position);
                for (int t = 0; t < Lz; ++t) {
                  corr_1x2[t] += (fac1 * fac2 * corr_scr[t]);
                }
              }
            }
            //double O1_plaq = 0.0;
            //double O1_1x1 = 0.0;
            //double O1_1x2 = 0.0;
            for (int t = 0; t < Lz; ++t) {
              corr_plaq[t] *= 2 * normalization;
              //O1_plaq += corr_plaq[t];
              corr_1x1[t] *= 2 * normalization;
              //O1_1x1 += corr_1x1[t];
              corr_1x2[t] *= 4 * normalization;
              //O1_1x2 += corr_1x2[t];
              vout.general(m_vl, "  O1_clover_plaq_z_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, corr_1x1[t]);
              double O1_rect = l_c_rect * corr_1x2[t];
              vout.general(m_vl, "  O1_clover_rect_z_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, O1_rect);
              double O1_imp = (m_c_plaq * corr_1x1[t] + m_c_rect * corr_1x2[t]);
              vout.general(m_vl, "  O1_clover_imp_z_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, O1_imp);
              result = O1_imp;
              vout.general(m_vl, "  O1_plaq_z_FT = %d %d %d %d %d %.8f %d %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, t, corr_plaq[t]);
            }

            /*
            double O1_imp = m_c_plaq * O1_1x1 + m_c_rect * O1_1x2;
            double O1_rect = l_c_rect * O1_1x2;
            vout.general(m_vl, "  O1_clover_plaq_sumz_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_1x1/Lz);
            //vout.general(m_vl, "  O1_clover_1x2_sumz_FT = %d %d %.8f %.16e\n", mu, nu, t_flow, O1_1x2/Lz);
            vout.general(m_vl, "  O1_clover_rect_sumz_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_rect/Lz);
            vout.general(m_vl, "  O1_clover_imp_sumz_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_imp/Lz);
            vout.general(m_vl, "  O1_plaq_sumz_FT = %d %d %d %d %d %.8f %.16e\n", momentum_sink[0], momentum_sink[1], momentum_sink[2], mu, nu, t_flow, O1_plaq/Lz);
            */
          }
        }
      }
    }
  }

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }
  return result;
}


//====================================================================
void EnergyMomentumTensor::set_field_strength(const Field_G& U)
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
  m_flag_field_set = 1;
}


//====================================================================
int EnergyMomentumTensor::factor(const int mu, const int nu)
{
  if (mu < nu) {
    return 1;
  } else if (mu > nu) {
    return -1;
  } else {
    return 0;
  }
}


//====================================================================
int EnergyMomentumTensor::index_munu2i(int mu, int nu)
{
  const int Ndim = CommonParameters::Ndim();

  assert(mu < Ndim);
  assert(nu < Ndim);

  if (mu > nu) {
    int scr = mu;

    mu = nu;
    nu = scr;
  }
  if ((mu == 0) && (nu == 1)) return 0;
  else if ((mu == 0) && (nu == 2)) return 1;
  else if ((mu == 0) && (nu == 3)) return 2;
  else if ((mu == 1) && (nu == 2)) return 3;
  else if ((mu == 1) && (nu == 3)) return 4;
  else if ((mu == 2) && (nu == 3)) return 5;
  else return 0;
}


//====================================================================
//============================================================END=====
