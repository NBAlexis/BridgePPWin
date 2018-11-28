#include "BridgeLib_Private.h"

/*!
        @file    $Id:: action_G_Rectangle.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "action_G_Rectangle.h"


#ifdef USE_FACTORY
namespace {
  Action *create_object()
  {
    return new Action_G_Rectangle();
  }


  bool init = Action::Factory::Register("Action_G_Rectangle", create_object);
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
void Action_G_Rectangle::set_parameters(double beta,
                                        double c_plaq, double c_rect)
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
double Action_G_Rectangle::langevin(RandomNumbers *)
{
  double H_U = calcH();

  return H_U;
}


//====================================================================
double Action_G_Rectangle::calcH()
{
  int          Nc   = CommonParameters::Nc();
  int          Ndim = CommonParameters::Ndim();
  int          Nvol = CommonParameters::Nvol();
  int          Lvol = CommonParameters::Lvol();
  const double eps  = CommonParameters::epsilon_criterion();


  int Ndim2  = Ndim * (Ndim - 1) / 2;
  int size_U = Lvol * Ndim2;

  Field_G Cup1(Nvol, 1), Cup2(Nvol, 1);
  Field_G Cdn1(Nvol, 1), Cdn2(Nvol, 1);
  Field_G Umu(Nvol, 1), Unu(Nvol, 1);
  Field_G v(Nvol, 1), w(Nvol, 1), c(Nvol, 1);

  double plaqF = 0.0;
  double rectF = 0.0;

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());


  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = mu + 1; nu < Ndim; ++nu) {
      m_staple.upper(Cup1, *m_U, mu, nu);

      //- plaquette term
      for (int site = 0; site < Nvol; ++site) {
        plaqF += ReTr(m_U->mat(site, mu) * Cup1.mat_dag(site));
      }

      //- rectangular terms
      // NB. skip this part, if m_c_rect = 0.0
      if (fabs(m_c_rect) > eps) {
        m_staple.upper(Cup2, *m_U, nu, mu);

        //      +---+---+
        //      |       |   term
        //      x   <---+

        copy(Umu, 0, *m_U, mu);
        copy(Unu, 0, *m_U, nu);

        m_shift.backward(v, Cup2, mu);
        m_shift.backward(c, Umu, nu);

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

  double plaq = plaqF / Nc;
  vout.general(m_vl, "    Plaquette    = %18.8f\n", plaq / size_U);

  double H_U = m_c_plaq * (Ndim2 * Lvol - plaqF / Nc)
               + m_c_rect * (Ndim2 * Lvol * 2 - rectF / Nc);

  H_U = m_beta * H_U;

  vout.general(m_vl, "    H_Grectangle = %18.8f\n", H_U);
  vout.general(m_vl, "    H_G/dof      = %18.8f\n", H_U / size_U);

  return H_U;
}


//====================================================================
void Action_G_Rectangle::force(Field& force)
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
