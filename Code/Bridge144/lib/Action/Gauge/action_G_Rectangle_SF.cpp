#include "BridgeLib_Private.h"

/*!
        @file    $Id:: action_G_Rectangle_SF.cpp #$

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "action_G_Rectangle_SF.h"


#ifdef USE_FACTORY
namespace {
  Action *create_object()
  {
    return new Action_G_Rectangle_SF();
  }


  bool init = Action::Factory::Register("Action_G_Rectangle_SF", create_object);
}
#endif



const std::string Action_G_Rectangle_SF::class_name = "Action_G_Rectangle_SF";

//====================================================================
void Action_G_Rectangle_SF::set_parameters(const Parameters& params)
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


  double gg = 6.0 / beta;

  double ct  = ct0 + ct1 * gg + ct2 * gg * gg;
  double ctr = ctr0 + ctr1 * gg + ctr2 * gg * gg;

  set_parameters(beta, c_plaq, c_rect, &phi[0], &phipr[0], ct, ctr);

  //- post-process
  m_force_G->set_parameters(params);
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
void Action_G_Rectangle_SF::set_parameters(double beta, double c_plaq, double c_rect,
                                           double *phi, double *phipr, double ct, double ctr)
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

  int    Lx = CommonParameters::Lx();
  double c0r, c0i, c1r, c1i, c2r, c2i;
  c0r = cos(phi[0] / Lx);
  c0i = sin(phi[0] / Lx);
  c1r = cos(phi[1] / Lx);
  c1i = sin(phi[1] / Lx);
  c2r = cos(phi[2] / Lx);
  c2i = sin(phi[2] / Lx);
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
double Action_G_Rectangle_SF::langevin(RandomNumbers *rand)
{
  double H_U = calcH();

  return H_U;
}


//====================================================================

/*!
  The improved gauge action with rectangle term:
\f[
 S_G=-\frac{\beta}{N_c}\left(c_0\sum_p{\rm Re}{\rm Tr} U_p
 +c_1\sum_r{\rm Re}{\rm Tr} U_r\right)
\f]
  <ul>
  <li>one plaquette term is added.
  <li>Two rectangular terms are added:
  <pre>
                               +---+
           +---+---+           |   |
           |       |           +   +
           x   <---+           |   |
                               x   v
  </pre>
  <li>We use Wk, Wk' for the boundary spatial link.
  <li>Contributions from the boundary spatial plaquettes and rectangles are set to zero.
  <li>The temporal rectangle that cross the boundary is set to zero.
<pre>
     +---+
     |   |
 t=0 +   +  --> 0
     |   |
     x   v
</pre>
  <ul>
  <li>These are automaticaly expected by a property of Staple_SF::upper()
  </ul>
  <li>Tree level improvement factor ct, ctr is implemented.
<pre>
      +---+       +---+---+
   ct |   |   ctr |       |
  t=0 x---+   t=0 x---+---+

      +---+       +---+---+
   ct |   |   ctr |       |
 Nt-1 x---+  Nt-1 x---+---+
</pre>
  </ul>
 */
double Action_G_Rectangle_SF::calcH()
{
  int Ndim = CommonParameters::Ndim();
  int Nvol = CommonParameters::Nvol();

  int Nt   = CommonParameters::Nt();
  int NPEt = CommonParameters::NPEt();

  // int Ndim2  = Ndim * (Ndim - 1) / 2;
  // int size_U = Lvol * Ndim2;

  //Field_G Cup1(Nvol,1), Cup2(Nvol,1);
  //Field_G Cdn1(Nvol,1), Cdn2(Nvol,1);
  //  Field_G Umu(Nvol,1), Unu(Nvol,1);
  //  Field_G v(Nvol,1), w(Nvol,1), c(Nvol,1);
  Field_G_SF Cup1, Cup2;
  Field_G_SF Cdn1, Cdn2;
  Field_G_SF Umu, Unu;
  Field_G_SF v, w, c;

  Mat_SU_N vmat(m_Nc), wmat(m_Nc), cmat(m_Nc);

  double plaqF = 0.0;
  double rectF = 0.0;

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = mu + 1; nu < Ndim; ++nu) {
      m_staple.upper(Cup1, *m_U, mu, nu);
      m_staple.upper(Cup2, *m_U, nu, mu);

      // plaquette term
      Unu = Cup2;
      // If the node is at the boundary the temporal plaquette is multiplied with ct.
      if ((nu == 3) && (Communicator::ipe(3) == 0)) {
        Unu.mult_ct_boundary(0, m_ct);
      }
      if ((nu == 3) && (Communicator::ipe(3) == NPEt - 1)) {
        Unu.mult_ct_boundary(Nt - 1, m_ct);
      }
      for (int site = 0; site < Nvol; ++site) {
        plaqF += ReTr(m_U->mat(site, nu) * Unu.mat_dag(site));
      }

      // rectangular terms

      //      +---+---+
      //      |       |   term
      //      x   <---+

      Umu.setpart_ex(0, *m_U, mu);
      Unu.setpart_ex(0, *m_U, nu);
      if ((Communicator::ipe(3) == 0) && (nu == 3)) {
        Umu.set_boundary_wk(m_wk);
      }

      m_shift.backward(v, Cup2, mu);
      m_shift.backward(c, Umu, nu);

      if ((Communicator::ipe(3) == 0) && (nu == 3)) {
        c.mult_ct_boundary(0, m_ctr);
      }
      if ((Communicator::ipe(3) == NPEt - 1) && (nu == 3)) {
        c.set_boundary_wkpr(m_wkpr);
        c.mult_ct_boundary(Nt - 1, m_ctr);
      }

      mult_Field_Gnd(w, 0, c, 0, v, 0);
      mult_Field_Gnn(c, 0, Unu, 0, w, 0);
      for (int site = 0; site < Nvol; ++site) {
        rectF += ReTr(Umu.mat(site) * c.mat_dag(site));
      }

      //      +---+
      //      |   |
      //      +   +   term
      //      |   |
      //      x   v

      m_shift.backward(v, Unu, mu);
      m_shift.backward(c, Cup1, nu);

      mult_Field_Gnd(w, 0, c, 0, v, 0);
      mult_Field_Gnn(c, 0, Unu, 0, w, 0);
      for (int site = 0; site < Nvol; ++site) {
        rectF += ReTr(Umu.mat(site) * c.mat_dag(site));
      }
    }
  }

  plaqF = Communicator::reduce_sum(plaqF);
  rectF = Communicator::reduce_sum(rectF);

  //  double plaq = plaqF/Nc;
  //  vout.general(m_vl,"    Plaquette    = %18.8f\n",plaq/size_U);

  //  double H_U =  m_c_plaq * (Ndim2*Lvol - plaqF/Nc)
  //              + m_c_rect * (Ndim2*Lvol*2 - rectF/Nc);
  double H_U = m_c_plaq * (-plaqF / m_Nc)
               + m_c_rect * (-rectF / m_Nc);

  H_U = m_beta * H_U;

  vout.general(m_vl, "    H_Grectangle = %18.8f\n", H_U);
  //  vout.general(m_vl,"    H_G/dof      = %18.8f\n",H_U/size_U);

  return H_U;
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
void Action_G_Rectangle_SF::force(Field& force)
{
  //- check of argument type
  assert(force.nin() == m_U->nin());
  assert(force.nvol() == m_U->nvol());
  assert(force.nex() == m_U->nex());

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  force.set(0.0);

  m_force_G->force_core(force, m_U);
}


//====================================================================
//============================================================END=====
