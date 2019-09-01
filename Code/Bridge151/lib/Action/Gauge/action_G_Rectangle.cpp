/*!
        @file    action_G_Rectangle.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "action_G_Rectangle.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Action_G_Rectangle::register_factory();
}
#endif

const std::string Action_G_Rectangle::class_name = "Action_G_Rectangle";

//====================================================================
void Action_G_Rectangle::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double beta, c_plaq, c_rect;

  int err = 0;
  err += params.fetch_double("beta", beta);
  err += params.fetch_double("c_plaq", c_plaq);
  err += params.fetch_double("c_rect", c_rect);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(beta, c_plaq, c_rect);

  //- post-process
  m_force_G->set_parameters(params);
}


//====================================================================
void Action_G_Rectangle::set_parameters(const double beta,
                                        const double c_plaq, const double c_rect)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  beta   = %12.6f\n", beta);
  vout.general(m_vl, "  c_plaq = %12.6f\n", c_plaq);
  vout.general(m_vl, "  c_rect = %12.6f\n", c_rect);

  //- range check
  // NB. beta,c_plaq,c_rect == 0 is allowed.

  //- store values
  m_beta   = beta;
  m_c_plaq = c_plaq;
  m_c_rect = c_rect;

  //- post-process
}


//====================================================================
double Action_G_Rectangle::langevin(RandomNumbers *rand)
{
  const double H_U = calcH();

  return H_U;
}


//====================================================================
double Action_G_Rectangle::calcH()
{
  const int Nc    = CommonParameters::Nc();
  const int Ndim  = CommonParameters::Ndim();
  const int Ndim2 = Ndim * (Ndim - 1) / 2;

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const double eps = CommonParameters::epsilon_criterion();


  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  double plaqF = 0.0;
  double rectF = 0.0;

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = mu + 1; nu < Ndim; ++nu) {
      Field_G Cup1;
      m_staple.upper(Cup1, *m_U, mu, nu);

      //- plaquette term
      for (int site = 0; site < Nvol; ++site) {
        plaqF += ReTr(m_U->mat(site, mu) * Cup1.mat_dag(site));
      }

      //- rectangular terms
      // NB. skip this part, if m_c_rect = 0.0
      if (fabs(m_c_rect) > eps) {
        Field_G Cup2;
        m_staple.upper(Cup2, *m_U, nu, mu);

        //      +---+---+
        //      |       |   term
        //      x   <---+

        Field_G Umu;
        copy(Umu, 0, *m_U, mu);

        Field_G Unu;
        copy(Unu, 0, *m_U, nu);

        Field_G v;
        m_shift.backward(v, Cup2, mu);

        Field_G c;
        m_shift.backward(c, Umu, nu);

        Field_G w;
        mult_Field_Gnd(w, 0, c, 0, v, 0);

        mult_Field_Gnn(c, 0, Unu, 0, w, 0);

        for (int site = 0; site < Nvol; ++site) {
          rectF += ReTr(m_U->mat(site, mu) * c.mat_dag(site));
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
          rectF += ReTr(m_U->mat(site, mu) * c.mat_dag(site));
        }
      }
    }
  }

  plaqF = Communicator::reduce_sum(plaqF);
  rectF = Communicator::reduce_sum(rectF);

  const double plaq = plaqF / Nc;
  vout.general(m_vl, "    Plaquette    = %18.8f\n", plaq / Nvol / NPE / Ndim2);

  double H_U = m_c_plaq * Ndim2 * Nvol * NPE - m_c_plaq * plaqF / Nc
               + m_c_rect * Ndim2 * Nvol * NPE * 2 - m_c_rect * rectF / Nc;

  H_U = m_beta * H_U;

  vout.general(m_vl, "    H_Grectangle = %18.8f\n", H_U);
  vout.general(m_vl, "    H_G/dof      = %18.8f\n", H_U / Nvol / NPE / Ndim2);

  return H_U;
}


//====================================================================
void Action_G_Rectangle::force(Field& force)
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
