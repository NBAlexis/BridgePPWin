/*!
        @file    force_G_Plaq_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "force_G_Plaq_SF.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Force_G_Plaq_SF::register_factory();
}
#endif

const std::string Force_G_Plaq_SF::class_name = "Force_G_Plaq_SF";

//====================================================================
void Force_G_Plaq_SF::set_parameters(const Parameters& params)
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
}


//====================================================================

/*!
  Set parameters for the Wilson plaquette force with the SF boundary.
  <ul>
  <li>m_beta
  <li>m_phi, m_phipr: boundary spatial link
  <li>m_ct: improvement factor for the boundary temporal plaquette.
  </ul>
*/
void Force_G_Plaq_SF::set_parameters(const double beta,
                                     double *phi, double *phipr, const double ct)
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

  m_ct = ct;

  m_phi   = phi;
  m_phipr = phipr;

  //- post-process
  m_staple.set_parameters(m_phi, m_phipr);
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
void Force_G_Plaq_SF::force_core(Field& force)
{
  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();
  const int Nc   = CommonParameters::Nc();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  for (int mu = 0; mu < Nex; ++mu) {
    Field_G st;
    m_staple.staple_ct(st, *m_U, mu, m_ct);

    Field_G force1;

    for (int site = 0; site < Nvol; ++site) {
      Mat_SU_N ut(Nc);
      ut = m_U->mat(site, mu) * st.mat_dag(site);
      ut.at();
      force1.set_mat(site, 0, ut);
    }
    axpy(force, mu, -(m_beta / Nc), force1, 0);
  }

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "    Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//====================================================================

/*!
  <ul>
  <li>The result is (x,y,z,t,mu,color,real part,imaginary part)
  </ul>
 */
void Force_G_Plaq_SF::print_force(const Field_G *U)
{
  const int Nc = CommonParameters::Nc();
  const int Nx = CommonParameters::Nx();
  const int Ny = CommonParameters::Ny();
  const int Nz = CommonParameters::Nz();
  const int Nt = CommonParameters::Nt();

  //  if(comm->ipe(3)==0){
  for (int t = 0; t < Nt; t++) {
    for (int x = 0; x < Nx; x++) {
      for (int y = 0; y < Ny; y++) {
        for (int z = 0; z < Nz; z++) {
          Index_lex index;
          int       ix = index.site(x, y, z, t);

          for (int mu = 0; mu < 4; mu++) {
            Mat_SU_N Tmp(Nc);
            Tmp = U->mat(ix, mu);

            for (int c = 0; c < Nc * Nc; ++c) {
              vout.general(m_vl, "%d %d %d %d %d %d %0.16e %0.16e\n",
                           x, y, z, t, mu, c, Tmp.r(c), Tmp.i(c));
            }
          }
        }
      }
    }
  }
}


//====================================================================
//============================================================END=====
