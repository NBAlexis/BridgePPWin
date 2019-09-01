/*!
        @file    staple_lex.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "staple_lex.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Staple_lex::register_factory();
}
#endif

const std::string Staple_lex::class_name = "Staple_lex";

//====================================================================
void Staple_lex::set_parameters(const Parameters& params)
{
  m_filename_output = params.get_string("filename_output");
  if (m_filename_output.empty()) {
    m_filename_output = "stdout";
  }

  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
double Staple_lex::plaquette(const Field_G& U)
{
  const double result = (plaq_s(U) + plaq_t(U)) / 2;

  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  vout.general(m_vl, "plaq = %20.16e\n", result);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }

  return result;
}


//====================================================================
double Staple_lex::plaq_s(const Field_G& U)
{
  const int Nc       = CommonParameters::Nc();
  const int Ndim     = CommonParameters::Ndim();
  const int Ndim_spc = Ndim - 1;

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  double plaq = 0.0;

  for (int mu = 0; mu < Ndim_spc; ++mu) {
    int nu = (mu + 1) % Ndim_spc;

    Field_G Umu;
    Umu.setpart_ex(0, U, mu);

    upper(m_staple, U, mu, nu);

    //    #pragma omp parallel
    {
      double plaq1 = real(dotc(m_staple, Umu));

#pragma omp master
      plaq += plaq1;
    }
  }

  return plaq / Nvol / NPE / Nc / Ndim_spc;
}


//====================================================================
double Staple_lex::plaq_t(const Field_G& U)
{
  const int Nc       = CommonParameters::Nc();
  const int Ndim     = CommonParameters::Ndim();
  const int Ndim_spc = Ndim - 1;

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int mu = Ndim_spc;

  Field_G Umu;

  Umu.setpart_ex(0, U, mu);

  double plaq = 0.0;
  for (int nu = 0; nu < Ndim_spc; ++nu) {
    upper(m_staple, U, mu, nu);

    //    #pragma omp parallel
    {
      double plaq1 = real(dotc(m_staple, Umu));

#pragma omp master
      plaq += plaq1;
    }
  }

  return plaq / Nvol / NPE / Nc / Ndim_spc;
}


//====================================================================
void Staple_lex::staple(Field_G& W, const Field_G& U, const int mu)
{
  const int Ndim = CommonParameters::Ndim();

  W.set(0.0);

  for (int nu = 0; nu < Ndim; ++nu) {
    if (nu != mu) {
      //W += upper(U, mu, nu);
      //W += lower(U, mu, nu);
      Field_G Vud(W.nvol(), 1);

      upper(Vud, U, mu, nu);
      axpy(W, 1.0, Vud);

      lower(Vud, U, mu, nu);
      axpy(W, 1.0, Vud);
    }
  }
}


//====================================================================
void Staple_lex::upper(Field_G& c, const Field_G& U,
                       const int mu, const int nu)
{
  // (1)  mu (2)
  //    +-->--+
  // nu |     |
  //   i+     +

  //#pragma omp parallel
  {
    Field_G Umu;
    Umu.setpart_ex(0, U, mu);

    Field_G Unu;
    Unu.setpart_ex(0, U, nu);

    m_shift->backward(m_v, Unu, mu);
    m_shift->backward(c, Umu, nu);

    mult_Field_Gnd(m_w, 0, c, 0, m_v, 0);
    mult_Field_Gnn(c, 0, U, nu, m_w, 0);
  }
}


//====================================================================
void Staple_lex::lower(Field_G& c, const Field_G& U,
                       const int mu, const int nu)
{
  //    +     +
  // nu |     |
  //   i+-->--+
  //  (1)  mu (2)

  //#pragma omp parallel
  {
    Field_G Unu;
    Unu.setpart_ex(0, U, nu);
    m_shift->backward(m_w, Unu, mu);

    mult_Field_Gnn(m_v, 0, U, mu, m_w, 0);
    mult_Field_Gdn(m_w, 0, U, nu, m_v, 0);

    m_shift->forward(c, m_w, nu);
  }
}


//====================================================================
//============================================================END=====
