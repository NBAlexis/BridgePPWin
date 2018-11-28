#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fopr_Chebyshev.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_Chebyshev.h"


#ifdef USE_FACTORY
namespace {
  Fopr *create_object(Fopr *fopr)
  {
    return new Fopr_Chebyshev(fopr);
  }


  bool init = Fopr::Factory_fopr::Register("Chevyshev", create_object);
}
#endif



const std::string Fopr_Chebyshev::class_name = "Fopr_Chebyshev";

//====================================================================
void Fopr_Chebyshev::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int    Np;
  double v_threshold, v_max;

  int err = 0;
  err += params.fetch_int("degree_of_polynomial", Np);
  err += params.fetch_double("threshold_value", v_threshold);
  err += params.fetch_double("upper_bound", v_max);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Np, v_threshold, v_max);
}


//====================================================================
void Fopr_Chebyshev::set_parameters(int Np, double v_threshold, double v_max)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Np          = %d\n", Np);
  vout.general(m_vl, "  v_threshold = %16.8e\n", v_threshold);
  vout.general(m_vl, "  v_max       = %16.8e\n", v_max);


  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Np);
  // NB. v_threshold,v_max == 0 is allowed.

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Npcb = Np;

  double b_max = v_max / v_threshold;
  double r     = 2.0 / (b_max * b_max - 1.0);
  double s     = v_threshold / sqrt(0.5 * r);

  m_Fcb1 = 2.0 / (s * s);
  m_Fcb2 = -(1.0 + r);

  vout.general(m_vl, "  Fcb1        = %16.8e\n", m_Fcb1);
  vout.general(m_vl, "  Fcb2        = %16.8e\n", m_Fcb2);
}


//====================================================================
void Fopr_Chebyshev::mult(Field& v, const Field& w)
{
  std::vector<Field> dj(3);
  int                Nin  = w.nin();
  int                Nvol = w.nvol();
  int                Nex  = w.nex();

  assert(v.nin() == Nin);
  assert(v.nvol() == Nvol);
  assert(v.nex() == Nex);

  for (int k = 0; k < 3; ++k) {
    dj[k].reset(Nin, Nvol, Nex);
  }

  dj[0] = w;
  scal(dj[0], -1.0);
  dj[1].set(0.0);

  int jn  = 2;
  int jp1 = 1;
  int jp2 = 0;

  for (int j = m_Npcb; j >= 2; --j) {
    m_fopr->mult(dj[jn], dj[jp1]);
    scal(dj[jn], m_Fcb1);
    axpy(dj[jn], m_Fcb2, dj[jp1]);

    scal(dj[jn], 2.0);
    axpy(dj[jn], -1.0, dj[jp2]);

    jn  = (jn + 1) % 3;
    jp1 = (jp1 + 1) % 3;
    jp2 = (jp2 + 1) % 3;
  }

  m_fopr->mult(v, dj[jp1]);
  scal(v, m_Fcb1);
  axpy(v, m_Fcb2, dj[jp1]);
  axpy(v, -1.0, dj[jp2]);
}


//====================================================================
double Fopr_Chebyshev::mult(double x)
{
  std::vector<double> dj(3);

  dj[0] = -1.0;
  dj[1] = 0.0;

  int jn  = 2;
  int jp1 = 1;
  int jp2 = 0;

  for (int j = m_Npcb; j >= 2; --j) {
    dj[jn]  = x * dj[jp1];
    dj[jn] *= m_Fcb1;
    dj[jn] += m_Fcb2 * dj[jp1];

    dj[jn] *= 2.0;
    dj[jn] -= 1.0 * dj[jp2];

    jn  = (jn + 1) % 3;
    jp1 = (jp1 + 1) % 3;
    jp2 = (jp2 + 1) % 3;
  }

  double v = x * dj[jp1];
  v *= m_Fcb1;
  v += m_Fcb2 * dj[jp1];
  v -= dj[jp2];

  return v;
}


//====================================================================
//============================================================END=====
