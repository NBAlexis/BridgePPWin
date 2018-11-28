#include "BridgeLib_Private.h"

/*!
        @file    $Id:: staple_eo.cpp #$

        @brief

        @author  UEDA, Satoru
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "staple_eo.h"


#ifdef USE_FACTORY
namespace {
  Staple *create_object()
  {
    return new Staple_eo();
  }


  bool init = Staple::Factory::Register("EvenOdd", create_object);
}
#endif


const std::string Staple_eo::class_name = "Staple_eo";

//====================================================================
void Staple_eo::set_parameters(const Parameters& params)
{
  m_filename_output = params.get_string("filename_output");
  if (m_filename_output.empty()) {
    m_filename_output = "stdout";
  }

  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
double Staple_eo::plaquette(const Field_G& U)
{
  Field_G  Ueo(U);
  Index_eo idx;

  idx.convertField(Ueo, U);

  return (plaq_s(Ueo) + plaq_t(Ueo)) / 2;
}


//====================================================================
double Staple_eo::plaq_s(const Field_G& Ueo)
{
  int Nc       = CommonParameters::Nc();
  int Lvol     = CommonParameters::Lvol();
  int Nvol     = CommonParameters::Nvol();
  int Ndim     = CommonParameters::Ndim();
  int Ndim_spc = Ndim - 1;

  double         plaq = 0.0;
  static Field_G staple;

  upper(staple, Ueo, 0, 1);
  for (int site = 0; site < Nvol; site++) {
    plaq += ReTr(Ueo.mat(site, 0) * staple.mat_dag(site));   // P_xy
  }

  upper(staple, Ueo, 1, 2);
  for (int site = 0; site < Nvol; site++) {
    plaq += ReTr(Ueo.mat(site, 1) * staple.mat_dag(site));   // P_yz
  }

  upper(staple, Ueo, 2, 0);
  for (int site = 0; site < Nvol; site++) {
    plaq += ReTr(Ueo.mat(site, 2) * staple.mat_dag(site));   // P_zx
  }

  plaq = Communicator::reduce_sum(plaq);

  return plaq / (Lvol * Nc * Ndim_spc);
}


//====================================================================
double Staple_eo::plaq_t(const Field_G& Ueo)
{
  int Nc       = CommonParameters::Nc();
  int Lvol     = CommonParameters::Lvol();
  int Nvol     = CommonParameters::Nvol();
  int Ndim     = CommonParameters::Ndim();
  int Ndim_spc = Ndim - 1;

  double         plaq = 0.0;
  static Field_G staple;

  for (int nu = 0; nu < Ndim - 1; nu++) {
    lower(staple, Ueo, 3, nu);
    for (int site = 0; site < Nvol; site++) {
      plaq += ReTr(Ueo.mat(site, 3) * staple.mat_dag(site));   // P_zx
    }
  }

  plaq = Communicator::reduce_sum(plaq);

  return plaq / (Lvol * Nc * Ndim_spc);
}


//====================================================================
void Staple_eo::staple(Field_G& W, const Field_G& Ueo, const int mu)
{
  int Ndim = CommonParameters::Ndim();

  W.set(0.0);
  Field_G u_tmp(W.nvol(), 1);
  for (int nu = 0; nu < Ndim; nu++) {
    if (nu != mu) {
      upper(u_tmp, Ueo, mu, nu);
      axpy(W, 1.0, u_tmp);
      lower(u_tmp, Ueo, mu, nu);
      axpy(W, 1.0, u_tmp);
    }
  }
}


//====================================================================
void Staple_eo::upper(Field_G& c, const Field_G& Ueo, const int mu, const int nu)
{
  // (1)  mu (2)
  //    +-->--+
  // nu |     |
  //   i+     +

  Field_G Umu, Unu;

  Umu.setpart_ex(0, Ueo, mu);
  Unu.setpart_ex(0, Ueo, nu);

  m_shift.backward(m_v, Unu, mu);
  m_shift.backward(c, Umu, nu);

  mult_Field_Gnd(m_w, 0, c, 0, m_v, 0);
  mult_Field_Gnn(c, 0, Unu, 0, m_w, 0);
}


//====================================================================
void Staple_eo::lower(Field_G& c, const Field_G& Ueo, const int mu, const int nu)
{
  //    +     +
  // nu |     |
  //   i+-->--+
  //  (1)  mu (2)

  Field_G Umu, Unu;

  Umu.setpart_ex(0, Ueo, mu);
  Unu.setpart_ex(0, Ueo, nu);

  m_shift.backward(m_w, Unu, mu);
  mult_Field_Gnn(m_v, 0, Umu, 0, m_w, 0);
  mult_Field_Gdn(m_w, 0, Unu, 0, m_v, 0);
  m_shift.forward(c, m_w, nu);
}


//====================================================================
//============================================================END=====
