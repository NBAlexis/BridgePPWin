/*!
        @file    source_Wall_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "source_Wall_SF.h"

const std::string Source_Wall_SF::class_name = "Source_Wall_SF";

//====================================================================
void Source_Wall_SF::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double ct_tilde;

  int err = 0;
  err += params.fetch_double("ct_tilde", ct_tilde);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(ct_tilde);
}


//====================================================================
void Source_Wall_SF::set_parameters(const double ct_tilde)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  ct_tilde = %10.6f\n", ct_tilde);

  //- range check
  // NB. ct_tilde=0 is allowed.

  //- store values
  m_ct_tilde = ct_tilde;
}


//====================================================================
void Source_Wall_SF::set_parameters(Field_G *U, const double ct_tilde)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  ct_tilde = %10.6f\n", ct_tilde);

  //- range check
  // NB. ct_tilde=0 is allowed.

  //- store values
  m_U        = U;
  m_ct_tilde = ct_tilde;
}


//====================================================================

/*!
  Use the smeared link for the SF Wall source.
*/
void Source_Wall_SF::set_parameters(Field_G *U, Director_Smear *dr_smear,
                                    const double ct_tilde)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  ct_tilde = %10.6f\n", ct_tilde);

  //- range check
  // NB. ct_tilde=0 is allowed.

  //- store values
  dr_smear->set_config(U);
  const int Nsmear = dr_smear->get_Nsmear();
  m_U = (Field_G *)dr_smear->getptr_smearedConfig(Nsmear);

  m_ct_tilde = ct_tilde;
}


//====================================================================

/*!
Set the wall source for the SF boundary propagator.
\f[
b_1\left({y},\alpha,a\right)_{\beta,b}=\widetilde{c}_t
U_0^\dagger\left(\vec{y},0\right)_{ab}
\delta_{y_0,1}\left(P_+\right)_{\alpha\beta}.
\f]
<ul>
<li>u0dag is used as \f$U_0^\dagger\left(\vec{y},0\right)_{ab}\f$.
<li>a=ac, b=ic
<li>\f$\alpha=\beta=\f$id
<li>[14 Apr 2012 Y.Taniguchi]
</ul>
 */
void Source_Wall_SF::set_t0(Field_F& src, const int ic, const int id)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  std::vector<int> Nsize(Ndim);

  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();
  Nsize[3] = CommonParameters::Nt();

  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();

  assert(ic < Nc);
  assert(id < Nd / 2);
  assert(src.nvol() == Nvol);
  assert(src.nex() == 1);

  Mat_SU_N u0dag(Nc);
  src.set(0.0);

  if (Communicator::ipe(3) == 0) {
    const int t = 1;

    for (int z = 0; z < Nsize[2]; ++z) {
      for (int y = 0; y < Nsize[1]; ++y) {
        for (int x = 0; x < Nsize[0]; ++x) {
          int site  = m_index.site(x, y, z, t);
          int site0 = m_index.site(x, y, z, t - 1);

          u0dag = m_U->mat_dag(site0, 3);

          for (int ac = 0; ac < Nc; ++ac) {
            src.set_ri(ac, id, site, 0, u0dag.r(ac, ic), u0dag.i(ac, ic));
          }
        }
      }
    }
    scal(src, m_ct_tilde);
  }
}


//====================================================================

/*!
Set the wall source for the SF boundary propagator.
\f[
b_{T-1}\left({y},\alpha,a\right)_{\beta,b}=\widetilde{c}_tU_0(\vec{y},T-1)_{ab}
\left(P_-\right)_{\alpha\beta}\delta_{y_0,T-1}
\f]
<ul>
<li>u0 is used as \f$U_0\left(\vec{y},T-a\right)_{ab}\f$.
<li>a=ac, b=ic
<li>\f$\alpha=\beta=\f$id
<li>[14 Apr 2012 Y.Taniguchi]
</ul>
 */
void Source_Wall_SF::set_tT(Field_F& src, const int ic, const int id)
{
  const int NPEt = CommonParameters::NPEt();
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  std::vector<int> Nsize(Ndim);

  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();
  Nsize[3] = CommonParameters::Nt();

  const int Nc = CommonParameters::Nc();
  const int Nd = CommonParameters::Nd();

  assert(ic < Nc);
  assert(id > Nd / 2 - 1);
  assert(id < Nd);
  assert(src.nvol() == Nvol);
  assert(src.nex() == 1);

  Mat_SU_N u0(Nc);
  src.set(0.0);

  if (Communicator::ipe(3) == (NPEt - 1)) {
    const int t = Nsize[3] - 1;

    for (int z = 0; z < Nsize[2]; ++z) {
      for (int y = 0; y < Nsize[1]; ++y) {
        for (int x = 0; x < Nsize[0]; ++x) {
          int site = m_index.site(x, y, z, t);
          u0 = m_U->mat(site, 3);

          for (int ac = 0; ac < Nc; ++ac) {
            src.set_ri(ac, id, site, 0, u0.r(ac, ic), u0.i(ac, ic));
          }
        }
      }
    }
    scal(src, m_ct_tilde);
  }
}


//====================================================================
//============================================================END=====
