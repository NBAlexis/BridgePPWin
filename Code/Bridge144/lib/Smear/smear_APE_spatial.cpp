/*!
        @file    $Id:: smear_APE_spatial.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "smear_APE_spatial.h"


#ifdef USE_FACTORY
namespace {
  Smear *create_object(Projection *proj)
  {
    return new Smear_APE_spatial(proj);
  }


  bool init = Smear::Factory::Register("APE_spatial", create_object);
}
#endif



const std::string Smear_APE_spatial::class_name = "Smear_APE_spatial";

//====================================================================
void Smear_APE_spatial::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double rho;

  int err = 0;
  err += params.fetch_double("rho", rho);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho);
}


//====================================================================
void Smear_APE_spatial::set_parameters(const double rho)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  rho = %10.6F\n", rho);

  //- range check
  // NB. rho == 0 is allowed.

  //- store values
  m_rho = rho;
}


//====================================================================
void Smear_APE_spatial::smear(Field_G& Usmear, const Field_G& U)
{
  int Nvol = CommonParameters::Nvol();

  assert(U.nvol() == Nvol);
  assert(U.nex() == m_Ndim);

  assert(Usmear.nvol() == Nvol);
  assert(Usmear.nex() == m_Ndim);

  int Ndim_spc = m_Ndim - 1;

  Field_G c_tmp(Nvol, 1), u_tmp(Nvol, 1), u_tmp2(Nvol, 1);

  Staple_lex staple;

  double plaq = staple.plaq_s(U);
  vout.general(m_vl, "  plaq_s(org  ) = %12.8f\n", plaq);
  plaq = staple.plaq_t(U);
  vout.general(m_vl, "  plaq_t(org  ) = %12.8f\n", plaq);

  Usmear.set(0.0);

  for (int mu = 0; mu < Ndim_spc; ++mu) {
    c_tmp.set(0.0);
    u_tmp.setpart_ex(0, U, mu);

    for (int nu = 0; nu < Ndim_spc; ++nu) {
      if (nu != mu) {
        staple.upper(u_tmp2, U, mu, nu);
        c_tmp.addpart_ex(0, u_tmp2, 0, m_rho);

        staple.lower(u_tmp2, U, mu, nu);
        c_tmp.addpart_ex(0, u_tmp2, 0, m_rho);
      }
    }

    m_proj->project(u_tmp2, m_rho, c_tmp, u_tmp);
    Usmear.setpart_ex(mu, u_tmp2, 0);
  }

  int mu = m_Ndim - 1; // temporal link: unsmeared.
  Usmear.setpart_ex(mu, U, mu);

  plaq = staple.plaq_s(Usmear);
  vout.general(m_vl, "  plaq_s(smear) = %12.8f\n", plaq);
  plaq = staple.plaq_t(Usmear);
  vout.general(m_vl, "  plaq_t(smear) = %12.8f\n", plaq);
}


//====================================================================
//============================================================END=====
