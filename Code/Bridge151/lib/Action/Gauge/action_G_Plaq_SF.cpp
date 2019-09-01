/*!
        @file    action_G_Plaq_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "action_G_Plaq_SF.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Action_G_Plaq_SF::register_factory();
}
#endif

const std::string Action_G_Plaq_SF::class_name = "Action_G_Plaq_SF";

//====================================================================
void Action_G_Plaq_SF::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double              beta, ct0, ct1, ct2;
  std::vector<double> phi, phipr;

  int err = 0;
  err += params.fetch_double("beta", beta);
  err += params.fetch_double("ct0", ct0);
  err += params.fetch_double("ct1", ct1);
  err += params.fetch_double("ct2", ct2);
  err += params.fetch_double_vector("phi", phi);
  err += params.fetch_double_vector("phipr", phipr);

  if (err) {
    vout.crucial(m_vl, "Error ar %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  const double gg = 6.0 / beta;
  const double ct = ct0 + ct1 * gg + ct2 * gg * gg;

  set_parameters(beta, &phi[0], &phipr[0], ct);

  //- post-process
  m_force_G->set_parameters(params);
}


//====================================================================

/*!
  Set parameters for the Wilson plaquette action with the SF boundary.
  <ul>
  <li>m_beta
  <li>m_phi, m_phipr: boundary spatial link
  <li>m_ct: improvement factor for the boundary temporal plaquette.
  </ul>
*/
void Action_G_Plaq_SF::set_parameters(const double beta,
                                      double *phi, double *phipr,
                                      const double ct)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  beta  = %12.6f\n", beta);
  vout.general(m_vl, "  phi1  = %12.6f\n", phi[0]);
  vout.general(m_vl, "  phi2  = %12.6f\n", phi[1]);
  vout.general(m_vl, "  phi3  = %12.6f\n", phi[2]);
  vout.general(m_vl, "  phipr1= %12.6f\n", phipr[0]);
  vout.general(m_vl, "  phipr2= %12.6f\n", phipr[1]);
  vout.general(m_vl, "  phipr3= %12.6f\n", phipr[2]);
  vout.general(m_vl, "  ct    = %12.6f\n", ct);

  //- range check
  // NB. beta,phi,phipr,ct = 0 is allowed.

  //- store values
  m_beta = beta;

  m_phi   = phi;
  m_phipr = phipr;
  m_ct    = ct;

  //- post-process
  m_staple.set_parameters(m_phi, m_phipr);
}


//====================================================================
double Action_G_Plaq_SF::langevin(RandomNumbers *rand)
{
  const double H_U = calcH();

  return H_U;
}


//====================================================================

/*!
  <ul>
  <li>The Wilson's plaquette gauge action:
\f[
 S_G=-\frac{\beta}{N_c}\sum_p{\rm Re}{\rm Tr} U_p
\f]
  <li>Calculate plaquette using Staple_SF.plaquette() function.
  <li>plaquette is NOT normalized to unity in Staple_SF since we have SF boundary and non-active spatial plaquette at the boundary.
  <li>The boundary improvement factor ct is implemented.
  <ul>
  <li>It is available by using plaquette_ct() for the plaquette sum.
  </ul>
  </ul>
 */
double Action_G_Plaq_SF::calcH()
{
  const int Nc = CommonParameters::Nc();

  const double plaq = m_staple.plaquette_ct(*m_U, m_ct);

  const double H_U = -(1.0 / Nc) * m_beta * plaq;

  vout.general(m_vl, "H_Gplaq    = %18.8f\n", H_U);

  return H_U;
}


//====================================================================

/*!
  The force for the Wilson plaquette action with the SF boundary.
  <ul>
  <li>The boundary improvement factor ct is implemented.
  <ul>
  <li>It is available by using m_staple.staple_ct for the staple.
  </ul>
  <li>Force for the boundary spatial link is set to zero.
  </ul>
*/
void Action_G_Plaq_SF::force(Field& force)
{
  //- check of argument element_type
  assert(force.nin() == m_U->nin());
  assert(force.nvol() == m_U->nvol());
  assert(force.nex() == m_U->nex());

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  force.set(0.0);

  m_force_G->force_core(force, m_U);
}


//====================================================================
//============================================================END=====
