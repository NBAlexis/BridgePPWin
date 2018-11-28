#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fopr_Wilson_Isochemical.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_Wilson_Isochemical.h"



#ifdef USE_FACTORY
namespace {
  Fopr *create_object()
  {
    return new Fopr_Wilson_Isochemical();
  }


  Fopr *create_object_with_repr(const std::string& repr)
  {
    return new Fopr_Wilson_Isochemical(repr);
  }


  bool init1 = Fopr::Factory_noarg::Register("Wilson_Isochemical", create_object);
  bool init2 = Fopr::Factory_string::Register("Wilson_Isochemical", create_object_with_repr);
}
#endif



const std::string Fopr_Wilson_Isochemical::class_name = "Fopr_Wilson_Isochemical";

//====================================================================
void Fopr_Wilson_Isochemical::set_parameters(const Parameters& params)
{
  const std::string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double           kappa, mu;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("isospin_chemical_potential", mu);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(kappa, mu, bc);
}


//====================================================================
void Fopr_Wilson_Isochemical::set_parameters(const double kappa,
                                             const double mu,
                                             const std::vector<int> bc)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  vout.general(m_vl, "  mu    = %12.8f\n", mu);
  for (int dir = 0; dir < m_Ndim; ++dir) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", dir, bc[dir]);
  }

  //- range check
  // NB. kappa,mu = 0 is allowed.
  assert(bc.size() == m_Ndim);

  //- store values
  m_kappa  = kappa;
  m_mu     = mu;
  m_exp_mu = exp(mu);

  m_boundary.resize(m_Ndim);
  for (int dir = 0; dir < m_Ndim; ++dir) {
    m_boundary[dir] = bc[dir];
  }

  //- propagate parameters
  m_fopr_w->set_parameters(kappa, bc);
}


//====================================================================
void Fopr_Wilson_Isochemical::init(std::string repr)
{
  m_Ndim = CommonParameters::Ndim();

  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int NinF = 2 * Nc * Nd;
  int Nvol = CommonParameters::Nvol();

  m_boundary.resize(m_Ndim);
  m_U = 0;

  m_fopr_w = new Fopr_Wilson(repr);

  m_vt.reset(NinF, Nvol, 1);
  m_w2.reset(NinF, Nvol, 1);
}


//====================================================================
void Fopr_Wilson_Isochemical::tidyup()
{
  delete m_fopr_w;
}


//====================================================================
void Fopr_Wilson_Isochemical::D(Field& w, const Field& f)
{
  int   Nvol = field_nvol();
  int   Nin  = field_nin();
  int   Nex  = field_nex();
  Field v(Nin, Nvol, Nex);

  Dspc(w, f);

  v.set(0.0);  // v = 0.0;
  m_fopr_w->mult_up(3, v, f);
  w.addpart_ex(0, v, 0, -m_kappa * m_exp_mu);

  v.set(0.0);  // v = 0.0;
  m_fopr_w->mult_dn(3, v, f);
  w.addpart_ex(0, v, 0, -m_kappa / m_exp_mu);

#pragma omp barrier
}


//====================================================================
void Fopr_Wilson_Isochemical::Dminmu(Field& w, const Field& f)
{
  int   Nvol = field_nvol();
  int   Nin  = field_nin();
  int   Nex  = field_nex();
  Field v(Nin, Nvol, Nex);

  Dspc(w, f);

  v.set(0.0);  // v = 0.0;
  m_fopr_w->mult_up(3, v, f);
  w.addpart_ex(0, v, 0, -m_kappa / m_exp_mu);

  v.set(0.0);  // v = 0.0;
  m_fopr_w->mult_dn(3, v, f);
  w.addpart_ex(0, v, 0, -m_kappa * m_exp_mu);

#pragma omp barrier
}


//====================================================================
void Fopr_Wilson_Isochemical::Dspc(Field& w, const Field& f)
{
  w.set(0.0);  // w = 0.0;

  m_fopr_w->mult_up(0, w, f);
  m_fopr_w->mult_dn(0, w, f);

  m_fopr_w->mult_up(1, w, f);
  m_fopr_w->mult_dn(1, w, f);

  m_fopr_w->mult_up(2, w, f);
  m_fopr_w->mult_dn(2, w, f);

  scal(w, -m_kappa); // w *= -m_kappa;
  axpy(w, 1.0, f);   // w += f;

#pragma omp barrier
}


//====================================================================
void Fopr_Wilson_Isochemical::mult_gm5(Field& v, const Field& f)
{
  m_fopr_w->mult_gm5(v, f);
}


//====================================================================
void Fopr_Wilson_Isochemical::mult_gm5p(int mu, Field_F& v, const Field_F& w)
{
  m_vt.set(0.0);  // vt = 0.0;

  m_fopr_w->mult_up(mu, m_vt, (Field)w);
  m_fopr_w->mult_gm5(v, m_vt);
}


//====================================================================
void Fopr_Wilson_Isochemical::DdagD(Field& w, const Field& f)
{
  D(m_w2, f);
  mult_gm5(w, m_w2);
  Dminmu(m_w2, w);
  mult_gm5(w, m_w2);
}


//====================================================================
void Fopr_Wilson_Isochemical::Ddag(Field& w, const Field& f)
{
  mult_gm5(w, f);
  Dminmu(m_w2, w);
  mult_gm5(w, m_w2);
}


//====================================================================
void Fopr_Wilson_Isochemical::H(Field& w, const Field& f)
{
  D(m_w2, f);
  mult_gm5(w, m_w2);
}


//====================================================================
void Fopr_Wilson_Isochemical::Hdag(Field& w, const Field& f)
{
  Dminmu(m_w2, f);
  mult_gm5(w, m_w2);
}


//====================================================================
double Fopr_Wilson_Isochemical::flop_count()
{
  //- Counting of floating point operations.
  //  not implemented, yet.

  double flop = 0.0;

  return flop;
}


//====================================================================
//============================================================END=====
