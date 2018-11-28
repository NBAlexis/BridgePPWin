#include "BridgeLib_Private.h"

/*!
        @file    $Id:: force_G_Rectangle.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "force_G_Rectangle.h"


#ifdef USE_FACTORY
namespace {
  Force_G *create_object()
  {
    return new Force_G_Rectangle();
  }


  bool init = Force_G::Factory::Register("Force_G_Rectangle", create_object);
}
#endif


const std::string Force_G_Rectangle::class_name = "Force_G_Rectangle";

//====================================================================
void Force_G_Rectangle::set_parameters(const Parameters& params)
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
}


//====================================================================
void Force_G_Rectangle::set_parameters(double beta, double c_plaq, double c_rect)
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
void Force_G_Rectangle::force_core(Field& force)
{
  int          Nc   = CommonParameters::Nc();
  int          Nvol = CommonParameters::Nvol();
  int          Ndim = CommonParameters::Ndim();
  const double eps  = CommonParameters::epsilon_criterion();

  assert(m_U->nin() == Nc * Nc * 2);
  assert(m_U->nvol() == Nvol);
  assert(m_U->nex() == Ndim);

  assert(force.nin() == Nc * Nc * 2);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Ndim);

  Mat_SU_N ut(Nc);

  Field_G force1(Nvol, 1), force2(Nvol, 1);

  Field_G Cup1(Nvol, 1), Cup2(Nvol, 1);
  Field_G Cdn1(Nvol, 1), Cdn2(Nvol, 1);
  Field_G Umu(Nvol, 1), Unu(Nvol, 1);
  Field_G v(Nvol, 1), w(Nvol, 1), c(Nvol, 1);

  for (int mu = 0; mu < Ndim; ++mu) {
    force1.set(0.0);
    for (int nu = 0; nu < Ndim; ++nu) {
      if (nu == mu) continue;

      m_staple.upper(Cup1, *m_U, mu, nu);
      m_staple.upper(Cup2, *m_U, nu, mu);
      m_staple.lower(Cdn1, *m_U, mu, nu);
      m_staple.lower(Cdn2, *m_U, nu, mu);

      //- plaquette term
      force1.addpart_ex(0, Cup1, 0, m_c_plaq);
      force1.addpart_ex(0, Cdn1, 0, m_c_plaq);


      //- rectangular term
      // NB. skip this part, if m_c_rect = 0.0
      if (fabs(m_c_rect) > eps) {
        Umu.setpart_ex(0, *m_U, mu);
        Unu.setpart_ex(0, *m_U, nu);

        //      +---+---+
        //      |       |   term
        //      x   <---+

        m_shift.backward(v, Cup2, mu);
        m_shift.backward(c, Umu, nu);

        mult_Field_Gnd(w, 0, c, 0, v, 0);
        mult_Field_Gnn(c, 0, Unu, 0, w, 0);

        force1.addpart_ex(0, c, 0, m_c_rect);

        //      +---+
        //      |   |
        //      +   +   term
        //      |   |
        //      x   v

        m_shift.backward(v, Unu, mu);
        m_shift.backward(c, Cup1, nu);

        mult_Field_Gnd(w, 0, c, 0, v, 0);
        mult_Field_Gnn(c, 0, Unu, 0, w, 0);

        force1.addpart_ex(0, c, 0, m_c_rect);

        //      +---+---+
        //      |       |   term
        //      +---x   v

        m_shift.backward(v, Unu, mu);
        m_shift.backward(c, Umu, nu);

        mult_Field_Gnd(w, 0, c, 0, v, 0);
        mult_Field_Gnn(c, 0, Cdn2, 0, w, 0);

        force1.addpart_ex(0, c, 0, m_c_rect);

        //      x   <---+
        //      |       |   term
        //      +---+---+

        m_shift.backward(v, Cup2, mu);

        mult_Field_Gnn(w, 0, Umu, 0, v, 0);
        mult_Field_Gdn(v, 0, Unu, 0, w, 0);

        m_shift.forward(c, v, nu);

        force1.addpart_ex(0, c, 0, m_c_rect);

        //      x   ^
        //      |   |
        //      +   +   term
        //      |   |
        //      +---+

        m_shift.backward(v, Unu, mu);

        mult_Field_Gnn(w, 0, Cdn1, 0, v, 0);
        mult_Field_Gdn(v, 0, Unu, 0, w, 0);

        m_shift.forward(c, v, nu);

        force1.addpart_ex(0, c, 0, m_c_rect);

        //      +---x   ^
        //      |       |   term
        //      +---+---+

        m_shift.backward(v, Unu, mu);

        mult_Field_Gnn(w, 0, Umu, 0, v, 0);
        mult_Field_Gdn(v, 0, Cdn2, 0, w, 0);

        m_shift.forward(c, v, nu);

        force1.addpart_ex(0, c, 0, m_c_rect);
      }
    }

    mult_Field_Gnd(force2, 0, *m_U, mu, force1, 0);
    at_Field_G(force2, 0);

    axpy(force, mu, -(m_beta / Nc), force2, 0);
  }

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "    Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//====================================================================
//============================================================END=====
