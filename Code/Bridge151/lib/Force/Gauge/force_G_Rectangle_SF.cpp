/*!
        @file    force_G_Rectangle_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "force_G_Rectangle_SF.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Force_G_Rectangle_SF::register_factory();
}
#endif

const std::string Force_G_Rectangle_SF::class_name = "Force_G_Rectangle_SF";

//====================================================================
void Force_G_Rectangle_SF::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double              beta, c_plaq, c_rect;
  double              ct0, ct1, ct2, ctr0, ctr1, ctr2;
  std::vector<double> phi, phipr;

  int err = 0;
  err += params.fetch_double("beta", beta);
  err += params.fetch_double("c_plaq", c_plaq);
  err += params.fetch_double("c_rect", c_rect);

  err += params.fetch_double("ct0", ct0);
  err += params.fetch_double("ct1", ct1);
  err += params.fetch_double("ct2", ct2);

  err += params.fetch_double("ctr0", ctr0);
  err += params.fetch_double("ctr1", ctr1);
  err += params.fetch_double("ctr2", ctr2);

  err += params.fetch_double_vector("phi", phi);
  err += params.fetch_double_vector("phipr", phipr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  const double gg = 6.0 / beta;

  const double ct  = ct0 + ct1 * gg + ct2 * gg * gg;
  const double ctr = ctr0 + ctr1 * gg + ctr2 * gg * gg;

  set_parameters(beta, c_plaq, c_rect, &phi[0], &phipr[0], ct, ctr);
}


//====================================================================

/*!
  Set parameters for the improved gauge action with the SF boundary.
  <ul>
    <li> m_beta
    <li> m_c_plaq, m_c_rect: plaquette and rectangle factor.
    <ul>
      <li> Iwasaki action: c_plaq =  3.648, c_rect = -0.331
    </ul>
    <li> m_phi, m_phipr: boundary spatial link
    <li> m_ct: improvement factor for the boundary temporal plaquette.
    <li> m_ctr: improvement factor for the boundary temporal rectangle with two links attached to the boundary.
  </ul>
*/
void Force_G_Rectangle_SF::set_parameters(const double beta, const double c_plaq, const double c_rect,
                                          double *phi, double *phipr, const double ct, const double ctr)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  beta   = %12.6f\n", beta);
  vout.general(m_vl, "  c_plaq = %12.6f\n", c_plaq);
  vout.general(m_vl, "  c_rect = %12.6f\n", c_rect);
  vout.general(m_vl, "  phi1   = %12.6f\n", phi[0]);
  vout.general(m_vl, "  phi2   = %12.6f\n", phi[1]);
  vout.general(m_vl, "  phi3   = %12.6f\n", phi[2]);
  vout.general(m_vl, "  phipr1 = %12.6f\n", phipr[0]);
  vout.general(m_vl, "  phipr2 = %12.6f\n", phipr[1]);
  vout.general(m_vl, "  phipr3 = %12.6f\n", phipr[2]);
  vout.general(m_vl, "  ct     = %12.6f\n", ct);
  vout.general(m_vl, "  ctr    = %12.6f\n", ctr);

  //- range check
  // NB. beta,c_plaq,c_rect,phi,phipr,ct,ctr = 0 is allowed.

  //- store values
  m_beta   = beta;
  m_c_plaq = c_plaq;
  m_c_rect = c_rect;

  m_ct  = ct;
  m_ctr = ctr;
  //  m_phi   = phi  ;
  //  m_phipr = phipr;

  //- post-process
  m_staple.set_parameters(phi, phipr);

  const int Lx = CommonParameters::Lx();

  double c0r = cos(phi[0] / Lx);
  double c0i = sin(phi[0] / Lx);
  double c1r = cos(phi[1] / Lx);
  double c1i = sin(phi[1] / Lx);
  double c2r = cos(phi[2] / Lx);
  double c2i = sin(phi[2] / Lx);
  m_wk.zero();
  m_wk.set(0, 0, c0r, c0i);
  m_wk.set(1, 1, c1r, c1i);
  m_wk.set(2, 2, c2r, c2i);

  c0r = cos(phipr[0] / Lx);
  c0i = sin(phipr[0] / Lx);
  c1r = cos(phipr[1] / Lx);
  c1i = sin(phipr[1] / Lx);
  c2r = cos(phipr[2] / Lx);
  c2i = sin(phipr[2] / Lx);
  m_wkpr.zero();
  m_wkpr.set(0, 0, c0r, c0i);
  m_wkpr.set(1, 1, c1r, c1i);
  m_wkpr.set(2, 2, c2r, c2i);
}


//====================================================================

/*!
  The force for the rectangle improved gauge action with the SF boundary.
  <ul>
  <li>We use Wk, Wk' for the boundary spatial link.
  <li>The boundary improvement factor ct, ctr is implemented.
  <li>ctr is multiplied to the following temporal rectangle staple
  <pre>
        +---+---+        +---+---+
        |  ctr              ctr  |
    t=0 +---+---x        x---+---+

        x   <---+        +---x   ^
        |  ctr  |        |  ctr  |
    t=0 +---+---+        +---+---+

        +---+---+        +---+---+
        |  ctr              ctr  |
 t=Nt-1 +---+---x        x---+---+

        +---+---+        +---+---+
        |  ctr  |        |  ctr  |
 t=Nt-1 x   <---+        +---x   v
  </pre>
  <li>Force for the boundary spatial link is set to zero.
  <pre>
        +---+---+             +---+---+
        |       |  --> 0      |       |  --> 0
    t=0 x   <---+         t=0 +---x   v

    t=0 x   <---+         t=0 +---x   ^
        |       |  --> 0      |       |  --> 0
        +---+---+             +---+---+
  </pre>
  <ul>
  <li>We notice that the upper and lower staple accompanied with the boundary spatial link is set to zero by Staple_SF::upper() and Staple_SF::lower().
  Corresponding contributions to the boundary spatial link are automaticaly zero.
  </ul>
  <li>Contribution from the non existing rectangle is automatically zero by Staple_SF::upper() and Staple_SF::lower().
  <pre>
        +---+        +---+
            |        |
    t=0 ^   +    t=0 +   ^  --> 0
        |   |        |   |
        +---+        +---+

        +---+        +---+
        |   |        |   |
   t=Nt +   +   t=Nt +   +  --> 0
            |        |
        <---+        +--->
  </pre>

  </ul>
*/
void Force_G_Rectangle_SF::force_core(Field& force)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int Nt   = CommonParameters::Nt();
  const int NPEt = CommonParameters::NPEt();

  assert(m_U->nin() == m_Nc * m_Nc * 2);
  assert(m_U->nvol() == Nvol);
  assert(m_U->nex() == Ndim);

  assert(force.nin() == m_Nc * m_Nc * 2);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Ndim);

  Mat_SU_N wzero(m_Nc);
  wzero.zero();

  for (int mu = 0; mu < Ndim; ++mu) {
    Field_G force1;
    force1.set(0.0);

    for (int nu = 0; nu < Ndim; ++nu) {
      if (nu == mu) continue;

      Field_G_SF Cup1;
      m_staple.upper(Cup1, *m_U, mu, nu);

      Field_G_SF Cup2;
      m_staple.upper(Cup2, *m_U, nu, mu);

      Field_G_SF Cdn1;
      m_staple.lower(Cdn1, *m_U, mu, nu);

      Field_G_SF Cdn2;
      m_staple.lower(Cdn2, *m_U, nu, mu);

      // plaquette term
      Field_G_SF Umu = Cup1;
      Field_G_SF Unu = Cdn1;
      if (Communicator::ipe(3) == 0) {
        if (mu == 3) {
          Umu.mult_ct_boundary(0, m_ct);
          Unu.mult_ct_boundary(0, m_ct);
        }
        if (nu == 3) {
          Unu.mult_ct_boundary(1, m_ct);
        }
      }
      if (Communicator::ipe(3) == NPEt - 1) {
        if (mu == 3) {
          Umu.mult_ct_boundary(Nt - 1, m_ct);
          Unu.mult_ct_boundary(Nt - 1, m_ct);
        }
        if (nu == 3) {
          Umu.mult_ct_boundary(Nt - 1, m_ct);
        }
      }
      force1.addpart_ex(0, Umu, 0, m_c_plaq);
      force1.addpart_ex(0, Unu, 0, m_c_plaq);

      // rectangular term
      Umu.setpart_ex(0, *m_U, mu);
      Unu.setpart_ex(0, *m_U, nu);
      // For the boundary spatial link, use Wk.
      if ((Communicator::ipe(3) == 0) && (nu == 3)) Umu.set_boundary_wk(m_wk);
      if ((Communicator::ipe(3) == 0) && (mu == 3)) Unu.set_boundary_wk(m_wk);

      //      +---+---+
      //      |       |   term
      //      x   <---+
      Field_G_SF v;
      m_shift.backward(v, Cup2, mu);

      Field_G_SF c;
      m_shift.backward(c, Umu, nu);
      // The force for the spatial link near the boundary. Multiplied with ctr.
      if ((Communicator::ipe(3) == NPEt - 1) && (nu == 3)) {
        c.set_boundary_wkpr(m_wkpr);
        c.mult_ct_boundary(Nt - 1, m_ctr);
      }

      Field_G_SF w;
      mult_Field_Gnd(w, 0, c, 0, v, 0);
      mult_Field_Gnn(c, 0, Unu, 0, w, 0);
      // The force for the boundary spatial link is set to zero.
      if ((Communicator::ipe(3) == 0) && (nu == 3)) c.set_boundary_zero();
      force1.addpart_ex(0, c, 0, m_c_rect);

      //      +---+
      //      |   |
      //      +   +   term
      //      |   |
      //      x   v

      m_shift.backward(v, Unu, mu);
      m_shift.backward(c, Cup1, nu);
      // The force for the temporal link near the boundary. Multiplied with ctr.
      if ((Communicator::ipe(3) == NPEt - 1) && (mu == 3)) {
        v.set_boundary_wkpr(m_wkpr);
        v.mult_ct_boundary(Nt - 1, m_ctr);
      }
      mult_Field_Gnd(w, 0, c, 0, v, 0);
      mult_Field_Gnn(c, 0, Unu, 0, w, 0);
      // The force for the boundary spatial link is set to zero.
      if ((Communicator::ipe(3) == 0) && (nu == 3)) c.set_boundary_zero();
      // The force for the boundary temporal link.
      if ((Communicator::ipe(3) == 0) && (mu == 3)) c.mult_ct_boundary(0, m_ctr);
      force1.addpart_ex(0, c, 0, m_c_rect);

      //      +---+---+
      //      |       |   term
      //      +---x   v

      m_shift.backward(v, Unu, mu);
      m_shift.backward(c, Umu, nu);
      if ((Communicator::ipe(3) == NPEt - 1) && (nu == 3)) {
        c.set_boundary_wkpr(m_wkpr);
        c.mult_ct_boundary(Nt - 1, m_ctr);
      }
      if ((Communicator::ipe(3) == NPEt - 1) && (mu == 3)) v.set_boundary_wkpr(m_wkpr);
      mult_Field_Gnd(w, 0, c, 0, v, 0);
      mult_Field_Gnn(c, 0, Cdn2, 0, w, 0);
      if ((Communicator::ipe(3) == 0) && (nu == 3)) c.set_boundary_zero();
      force1.addpart_ex(0, c, 0, m_c_rect);

      //      x   <---+
      //      |       |   term
      //      +---+---+

      m_shift.backward(v, Cup2, mu);
      mult_Field_Gnn(w, 0, Umu, 0, v, 0);
      mult_Field_Gdn(v, 0, Unu, 0, w, 0);
      if ((Communicator::ipe(3) == 0) && (nu == 3)) v.mult_ct_boundary(0, m_ctr);
      m_shift.forward(c, v, nu);
      if ((Communicator::ipe(3) == 0) && (nu == 3)) c.set_boundary_zero();
      force1.addpart_ex(0, c, 0, m_c_rect);

      //      x   ^
      //      |   |
      //      +   +   term
      //      |   |
      //      +---+

      m_shift.backward(v, Unu, mu);
      if ((Communicator::ipe(3) == NPEt - 1) && (mu == 3)) {
        v.set_boundary_wkpr(m_wkpr);
        v.mult_ct_boundary(Nt - 1, m_ctr);
      }
      mult_Field_Gnn(w, 0, Cdn1, 0, v, 0);
      mult_Field_Gdn(v, 0, Unu, 0, w, 0);
      if ((Communicator::ipe(3) == 0) && (mu == 3)) v.mult_ct_boundary(0, m_ctr);
      m_shift.forward(c, v, nu);
      if ((Communicator::ipe(3) == 0) && (nu == 3)) c.set_boundary_zero();
      force1.addpart_ex(0, c, 0, m_c_rect);

      //      +---x   ^
      //      |       |   term
      //      +---+---+

      m_shift.backward(v, Unu, mu);
      if ((Communicator::ipe(3) == NPEt - 1) && (mu == 3)) v.set_boundary_wkpr(m_wkpr);
      mult_Field_Gnn(w, 0, Umu, 0, v, 0);
      mult_Field_Gdn(v, 0, Cdn2, 0, w, 0);
      if ((Communicator::ipe(3) == 0) && (nu == 3)) v.mult_ct_boundary(0, m_ctr);
      m_shift.forward(c, v, nu);
      if ((Communicator::ipe(3) == 0) && (nu == 3)) c.set_boundary_zero();
      force1.addpart_ex(0, c, 0, m_c_rect);
    }

    Field_G force2;
    mult_Field_Gnd(force2, 0, *m_U, mu, force1, 0);
    at_Field_G(force2, 0);

    axpy(force, mu, -(m_beta / m_Nc), force2, 0);
  }

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "    Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//====================================================================
//============================================================END=====
