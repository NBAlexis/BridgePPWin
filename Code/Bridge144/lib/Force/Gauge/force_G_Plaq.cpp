#include "BridgeLib_Private.h"

/*!
        @file    $Id:: force_G_Plaq.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "force_G_Plaq.h"


#ifdef USE_FACTORY
namespace {
  Force_G *create_object()
  {
    return new Force_G_Plaq();
  }


  bool init = Force_G::Factory::Register("Force_G_Plaq", create_object);
}
#endif


const std::string Force_G_Plaq::class_name = "Force_G_Plaq";

//====================================================================
void Force_G_Plaq::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double beta;

  int err = 0;
  err += params.fetch_double("beta", beta);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(beta);
}


//====================================================================
void Force_G_Plaq::set_parameters(const double beta)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  beta = %8.4f\n", beta);

  //- range check
  // NB. kappa == 0 is allowed.

  //- store values
  m_beta = beta;
}


//====================================================================
void Force_G_Plaq::force_core(Field& force)
{
  int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  assert(m_U->nin() == Nc * Nc * 2);
  assert(m_U->nvol() == Nvol);
  assert(m_U->nex() == Ndim);

  assert(force.nin() == Nc * Nc * 2);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Ndim);

  Mat_SU_N ut(Nc); // Nc x Nc complex matrix (SA)
  Field_G  st(Nvol, 1);
  Field_G  force1(Nvol, 1);

  for (int mu = 0; mu < Ndim; ++mu) {
    m_staple.staple(st, *m_U, mu);
    // Calculate staple for m_U(all site,\mu) (SA)

    /* -->--
       |     |
       |     |
       U_mu
    */
    for (int site = 0; site < Nvol; ++site) {
      ut = m_U->mat(site, mu) * st.mat_dag(site); // U_\mu * (staple)^\dagger (SA)
      ut.at();                                    // anti-hermitian and traceless (SA)
      //ut *= -beta/Nc;
      // force = -beta*{U_\mu *staple^\dagger}_{traceless & anti-hermitian) (SA)
      force1.set_mat(site, 0, ut);
    }

    axpy(force, mu, -(m_beta / Nc), force1, 0);
  }

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);

  //- calculate average, max, deviation of force over (site, mu)
  vout.general(m_vl, "    Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//====================================================================
//============================================================END=====
