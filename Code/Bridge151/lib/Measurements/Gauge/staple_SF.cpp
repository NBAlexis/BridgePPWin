/*!
        @file    staple_SF.cpp

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "staple_SF.h"

/*!
  Set the SF BC for wk and wkpr.
  \f[
  U_k(x)|_{t=0}=W_k(\vec{x})=\exp\left(a C_k\right),\quad
  C_k=\frac{i}{L}\pmatrix{\phi_1\cr&\phi_2\cr&&\phi_3\cr},
  \f]
  \f[
  U_k(x)|_{t=T}=W_k'(\vec{x})=\exp\left(a C_k'\right),\quad
  C_k'=\frac{i}{L}\pmatrix{\phi'_1\cr&\phi'_2\cr&&\phi'_3\cr},
  \f]
  omega0 is set to the default value
  \f[
  \verb|iomega0|=i\Omega_0=i\pmatrix{1\cr&-\frac{1}{2}\cr&&-\frac{1}{2}\cr}
  \f]
*/

const std::string Staple_SF::class_name = "Staple_SF";

//====================================================================
void Staple_SF::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  std::vector<double> phi, phipr, p_omega;

  int err = 0;
  err += params.fetch_double_vector("phi", phi);
  err += params.fetch_double_vector("phipr", phipr);
  err += params.fetch_double_vector("p_omega", p_omega);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(phi, phipr, p_omega);  // call std::vector version
}


//====================================================================
void Staple_SF::set_parameters(const std::vector<double>& phi, const std::vector<double>& phipr,
                               const std::vector<double>& p_omega)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  phi1     = %12.6f\n", phi[0]);
  vout.general(m_vl, "  phi2     = %12.6f\n", phi[1]);
  vout.general(m_vl, "  phi3     = %12.6f\n", phi[2]);
  vout.general(m_vl, "  phipr1   = %12.6f\n", phipr[0]);
  vout.general(m_vl, "  phipr2   = %12.6f\n", phipr[1]);
  vout.general(m_vl, "  phipr3   = %12.6f\n", phipr[2]);
  vout.general(m_vl, "  p_omega1 = %12.6f\n", p_omega[0]);
  vout.general(m_vl, "  p_omega2 = %12.6f\n", p_omega[1]);
  vout.general(m_vl, "  p_omega3 = %12.6f\n", p_omega[2]);

  //- range check
  // NB. phi,phipr,p_omega == 0 is allowed.
  assert(phi.size() == 3);
  assert(phipr.size() == 3);
  assert(p_omega.size() == 3);

  //- store values
  set_parameters(&phi[0], &phipr[0], &p_omega[0]);  // call double[] version
}


//====================================================================
void Staple_SF::set_parameters(const double *phi, const double *phipr)
{
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

  m_i_omega0.zero();
  m_i_omega0.set(0, 0, 0.0, 1.0);
  m_i_omega0.set(1, 1, 0.0, -0.5);
  m_i_omega0.set(2, 2, 0.0, -0.5);

  m_initialized = 1;
}


//====================================================================

/*!
  Set the SF BC for wk, wkpr and omega0.
  \f[
  U_k(x)|_{t=0}=W_k(\vec{x})=\exp\left(a C_k\right),\quad
  C_k=\frac{i}{L}\pmatrix{\phi_1\cr&\phi_2\cr&&\phi_3\cr},
  \f]
  \f[
  U_k(x)|_{t=T}=W_k'(\vec{x})=\exp\left(a C_k'\right),\quad
  C_k'=\frac{i}{L}\pmatrix{\phi'_1\cr&\phi'_2\cr&&\phi'_3\cr},
  \f]
  omega0 is set to the input value
  \f[
  \verb|omega0|=i\Omega_0
  \f]
*/
void Staple_SF::set_parameters(const double *phi, const double *phipr, const double *p_omega)
{
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

  m_i_omega0.zero();
  m_i_omega0.set(0, 0, 0.0, p_omega[0]);
  m_i_omega0.set(1, 1, 0.0, p_omega[1]);
  m_i_omega0.set(2, 2, 0.0, p_omega[2]);

  m_initialized = 1;
}


//====================================================================

/*!
  Evaluate boudary plaqutte VEV for the SF coupling.
  The SF running coupling is defined as
  \f[
  \frac{k}{\overline{g}^2}=
 \left\langle\frac{\partial S_g}{\partial\eta}\right\rangle
+\left\langle\frac{\partial S_f}{\partial\eta}\right\rangle
  \f]
  \f[
  k=12\left(\frac{L}{a}\right)^2
  \left(\sin\theta+\sin\left(2\theta\right)\right),
  \quad
  \theta=\frac{1}{3}\pi\left(\frac{a^2}{TL}\right)
  \f]
  \f[
  \frac{\partial S_g}{\partial\eta}=
  -\frac{2}{g_0^2}\sum_{\vec{x}}
  c_t\sum_{k=1}^{3}\frac{a}{L}{\rm Re\ }{\rm tr}\Bigl(
  i\Omega_0W_k(\vec{x})U_0(\vec{x}+\hat{k},0)U_k^\dagger(\vec{x},1)
  U_0^\dagger(\vec{x},0)
  -i\Omega_0W_k'U_0^\dagger(\vec{x}+\hat{k},T-1)U_k^\dagger(\vec{x},T-1)
  U_0(\vec{x},T-1)
  \Bigr)
  \f]
  \f[
  \frac{\partial S_g}{\partial\eta}=
  -\frac{2}{g_0^2}\sum_{\vec{x}}
  c_t\sum_{k=1}^{3}\frac{a}{L}{\rm Re\ }{\rm tr}\Bigl(
  i\Omega_0\verb|upper|(3,k)(\vec{x},0)U_0^\dagger(\vec{x},0)
  -i\Omega_0\verb|upper|(3,k)(\vec{x},T-1)^\dagger U_0(\vec{x},T-1)
  \Bigr)
  \f]
  <ul>
  <li>Evaluete the following quantity and print the result
  \f[
  -c_t\sum_{\vec{x}}\sum_{k=1}^{3}\frac{a}{L}{\rm Re\ }{\rm tr}\Bigl(
  i\Omega_0\verb|upper|(3,k)(\vec{x},0)U_0^\dagger(\vec{x},0)
  -i\Omega_0\verb|upper|(3,k)(\vec{x},T-1)^\dagger U_0(\vec{x},T-1)
  \Bigr)
  \f]
  <li>Clover term contribution \f$\partial S_f/\partial\eta\f$ is not evaluated here.
  </ul>
*/
double Staple_SF::sf_coupling_plaq(const Field_G& U, const double ct)
{
  if (!m_initialized) {
    vout.crucial(m_vl, "Error at %s: Parameter is not initialized.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  const int Nx   = CommonParameters::Nx();
  const int Ny   = CommonParameters::Ny();
  const int Nz   = CommonParameters::Nz();
  const int Nt   = CommonParameters::Nt();
  const int NPEt = CommonParameters::NPEt();

  const int Lx = CommonParameters::Lx();

  double plaq   = 0.0;
  double plaqt0 = 0.0;
  double plaqtT = 0.0;

  for (int nu = 0; nu < Ndim - 1; nu++) {
    Field_G_SF staple;
    upper(staple, U, 3, nu);

    for (int z = 0; z < Nz; z++) {
      for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
          // boundary
          if (Communicator::ipe(3) == 0) {
            int t    = 0;
            int site = m_index.site(x, y, z, t);

            Mat_SU_N up(Nc);
            up = staple.mat(site) * U.mat_dag(site, 3);
            double scr = ReTr(m_i_omega0 * up);

            /*
            up.unit();
            up *= staple.mat(site);
            up *= U->mat_dag(site,3);
            up *= m_i_omega0;
            scr = ReTr( up );
            */
            plaq   -= scr;
            plaqt0 += scr;
          }
          // boundary
          if (Communicator::ipe(3) == NPEt - 1) {
            int t    = Nt - 1;
            int site = m_index.site(x, y, z, t);

            Mat_SU_N up(Nc);
            up = staple.mat_dag(site) * U.mat(site, 3);
            double scr = ReTr(m_i_omega0 * up);

            /*
            up.unit();
            up *= staple.mat_dag(site);
            up *= U->mat(site,3);
            up *= m_i_omega0;
            scr = ReTr( up );
            */
            plaq   += scr;
            plaqtT += scr;
          }
        }
      }
    }
  }
  plaq   = Communicator::reduce_sum(plaq);
  plaqt0 = Communicator::reduce_sum(plaqt0);
  plaqtT = Communicator::reduce_sum(plaqtT);

  plaq   *= ct / Lx;
  plaqt0 *= ct / Lx;
  plaqtT *= ct / Lx;

  vout.general(m_vl, "SF_delSg_plaq, from 0, from T = %.8f %.8f %.8f\n", plaq, plaqt0, plaqtT);

  return plaq;
}


//====================================================================

/*!
  Evaluate the temporal rectangle VEV at the boundary for the SF running coupling.

  The SF running coupling is given by
\f[
  \frac{k}{\overline{g}^2}=
 \left\langle\frac{\partial S_g^{\rm plaq}}{\partial\eta}\right\rangle
+\left\langle\frac{\partial S_g^{\rm rect}}{\partial\eta}\right\rangle
+\left\langle\frac{\partial S_f}{\partial\eta}\right\rangle
\f]
\f[
\frac{\partial S_g}{\partial\eta}=
-\frac{2}{g_0^2}c_1\sum_{\vec{x}}\sum_{k=1}^{3}
\frac{ia}{L}{\rm Re}{\rm Tr}\Bigl(
c_t^R\left(
\Omega S_{kk0}(\vec{x},0)
+\Omega S_{0kk}(\vec{x},0)\right)
+\Omega S_{00k}(\vec{x},0)
\Bigr)
\f]
\f[
+\frac{2}{g_0^2}c_1\sum_{\vec{x}}\sum_{k=1}^{3}
\frac{ia}{L}{\rm Re}{\rm Tr}\Bigl(
c_t^R\left(
\Omega'S_{kk0}(\vec{x},T)
+\Omega'S_{0kk}(\vec{x},T)\right)
+\Omega'S_{00k}(\vec{x},T)
\Bigr)
\f]
\f[
S_{kk0}(\vec{x},0)=W_kW_kU_0(\vec{x}+2\hat{k},0)
U_k^\dagger(\vec{x}+\hat{k},1)U_k^\dagger(\vec{x},1)U_0^\dagger(\vec{x},0)
\f]
\f[
S_{0kk}(\vec{x},0)=W_kU_0(\vec{x}+\hat{k},0)U_k^\dagger(\vec{x},1)
U_k^\dagger(\vec{x}-\hat{k},1)U_0^\dagger(\vec{x}-\hat{k},0)W_k
\f]
\f[
S_{00k}(\vec{x},0)=W_kU_0(\vec{x}+\hat{k},0)U_0(\vec{x}+\hat{k},1)
U_k^\dagger(\vec{x},2)U_0^\dagger(\vec{x},1)U_0^\dagger(\vec{x},0)
\f]
\f[
S_{kk0}(\vec{x},T)=W_k'W_k'U_0^\dagger(\vec{x}+2\hat{k},T)
U_k^\dagger(\vec{x}+\hat{k},T-1)U_k^\dagger(\vec{x},T-1)U_0(\vec{x},T-1)
\f]
\f[
S_{0kk}(\vec{x},T)=W_k'U_0^\dagger(\vec{x}+\hat{k},T-1)
U_k^\dagger(\vec{x},T-1)U_k^\dagger(\vec{x}-\hat{k},T-1)
U_0(\vec{x}-\hat{k},T-1)W_k'
\f]
\f[
S_{00k}(\vec{x},T)=W_k'U_0^\dagger(\vec{x}+\hat{k},T-1)
U_0^\dagger(\vec{x}+\hat{k},T-2)U_k^\dagger(\vec{x},T-2)U_0(\vec{x},T-2)
U_0(\vec{x},T-1)
\f]
In this function we evaluate
\f[
-\sum_{\vec{x}}\sum_{k=1}^{3}
\frac{ia}{L}{\rm Re}{\rm Tr}\Bigl(
c_t^R\left(
\Omega S_{kk0}(\vec{x},0)
+\Omega S_{0kk}(\vec{x},0)\right)
+\Omega S_{00k}(\vec{x},0)
\Bigr)
\f]
\f[
+\sum_{\vec{x}}\sum_{k=1}^{3}
\frac{ia}{L}{\rm Re}{\rm Tr}\Bigl(
c_t^R\left(
\Omega'S_{kk0}(\vec{x},T)
+\Omega'S_{0kk}(\vec{x},T)\right)
+\Omega'S_{00k}(\vec{x},T)
\Bigr)
\f]
The following quantities are also printed out.
<pre>
  rect01         rect02
      <---<---+      <---<---+
      |       |      |       |
  t=0 x--->---+  t=0 +---x---+
  omega0               omega0

  rectt1         rectt2
  omega0              omega0
 t=Nt x--->---+ t=Nt +---x---+
      |       |      |       |
      +---<---+      <---<---+

  rect03      rectt3
              omega0
      +---+  t=Nt x--->
      |   |       |   |
      v   ^       ^   v
      |   |       |   |
  t=0 x--->       +---+
   omega0
</pre>
*/
double Staple_SF::sf_coupling_rect(const Field_G& U, const double ctr)
{
  if (!m_initialized) {
    vout.crucial(m_vl, "Error at %s: Parameter is not initialized.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();

  const int Nx   = CommonParameters::Nx();
  const int Ny   = CommonParameters::Ny();
  const int Nz   = CommonParameters::Nz();
  const int Nt   = CommonParameters::Nt();
  const int NPEt = CommonParameters::NPEt();

  const int Lx = CommonParameters::Lx();

  const int nu = 3;

  double rect01 = 0.0;
  double rect02 = 0.0;
  double rect03 = 0.0;
  double rectt1 = 0.0;
  double rectt2 = 0.0;
  double rectt3 = 0.0;

  for (int mu = 0; mu < Ndim - 1; mu++) {
    // rect01
    //      <---<---+
    //      |       |
    //  t=0 x--->---+
    //  omega0

    // rect02
    //      <---<---+
    //      |       |
    //  t=0 +---x---+
    //       omega0

    // rectt1
    // omega0
    // t=Nt x--->---+
    //      |       |
    //      +---<---+

    // rectt2
    //       omega0
    // t=Nt +---x---+
    //      |       |
    //      <---<---+
    Field_G_SF Cup2;
    upper(Cup2, U, nu, mu);

    Field_G Umu;
    Umu.setpart_ex(0, U, mu);

    Field_G Unu;
    Unu.setpart_ex(0, U, nu);

    Field_G v;
    m_shift.backward(v, Cup2, mu);

    Field_G c;
    m_shift.backward(c, Umu, nu);

    for (int z = 0; z < Nz; z++) {
      for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
          int t = 0;
          if (Communicator::ipe(3) == 0) {
            int site = m_index.site(x, y, z, t);

            Mat_SU_N wmat(Nc);
            wmat = m_wk * v.mat(site);

            Mat_SU_N cmat(Nc);
            cmat = wmat * c.mat_dag(site);

            wmat    = cmat * Unu.mat_dag(site);
            rect01 += ReTr(m_i_omega0 * wmat);

            wmat    = m_i_omega0 * v.mat(site);
            cmat    = wmat * c.mat_dag(site);
            wmat    = cmat * Unu.mat_dag(site);
            rect02 += ReTr(m_wk * wmat);
          }

          t = Nt - 1;
          if (Communicator::ipe(3) == NPEt - 1) {
            int site = m_index.site(x, y, z, t);

            Mat_SU_N cmat(Nc);
            cmat = m_i_omega0 * m_wkpr;

            Mat_SU_N wmat(Nc);
            wmat = cmat * v.mat_dag(site);

            cmat    = Unu.mat(site) * wmat;
            rectt1 += ReTr(cmat * U.mat_dag(site, mu));

            cmat    = m_i_omega0 * v.mat_dag(site);
            wmat    = m_wkpr * cmat;
            cmat    = Unu.mat(site) * wmat;
            rectt2 += ReTr(cmat * U.mat_dag(site, mu));
          }
        }
      }
    }

    // rect03
    //      +---+
    //      |   |
    //      v   ^
    //      |   |
    //  t=0 x--->
    //   omega0
    Field_G_SF Cup1;
    upper(Cup1, U, mu, nu);

    m_shift.backward(v, Unu, mu);
    m_shift.backward(c, Cup1, nu);

    for (int z = 0; z < Nz; z++) {
      for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
          int t = 0;
          if (Communicator::ipe(3) == 0) {
            int site = m_index.site(x, y, z, t);

            Mat_SU_N wmat(Nc);
            wmat = c.mat(site) * v.mat_dag(site);

            Mat_SU_N cmat(Nc);
            cmat = Unu.mat(site) * wmat;

            wmat    = m_wk * cmat.dag();
            rect03 += ReTr(m_i_omega0 * wmat);
          }
        }
      }
    }

    // rectt3
    //   omega0
    // t=Nt x--->
    //      |   |
    //      ^   v
    //      |   |
    //      +---+
    Field_G_SF Cdn1;
    lower(Cdn1, U, mu, nu);

    m_shift.backward(v, Unu, mu);

    for (int z = 0; z < Nz; z++) {
      for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
          int t = Nt - 1;
          if (Communicator::ipe(3) == NPEt - 1) {
            int site = m_index.site(x, y, z, t);

            Mat_SU_N wmat(Nc);
            wmat = m_i_omega0 * m_wkpr;

            Mat_SU_N cmat(Nc);
            cmat = wmat * v.mat_dag(site);

            wmat    = cmat * Cdn1.mat_dag(site);
            rectt3 += ReTr(Unu.mat(site) * wmat);
          }
        }
      }
    }
  }
  rect01 = Communicator::reduce_sum(rect01);
  rect02 = Communicator::reduce_sum(rect02);
  rect03 = Communicator::reduce_sum(rect03);
  rectt1 = Communicator::reduce_sum(rectt1);
  rectt2 = Communicator::reduce_sum(rectt2);
  rectt3 = Communicator::reduce_sum(rectt3);

  rect01 *= ctr / Lx;
  rect02 *= ctr / Lx;
  rect03 /= Lx;
  rectt1 *= ctr / Lx;
  rectt2 *= ctr / Lx;
  rectt3 /= Lx;

  double rect = -rect01 - rect02 - rect03 + rectt1 + rectt2 + rectt3;

  vout.general(m_vl, "SF_delSg_rect, at 01, 02, 03, at T1, T2, T3 = %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",
               rect, rect01, rect02, rect03, rectt1, rectt2, rectt3);

  return rect;
}


//====================================================================

/*!
  Evaluate summed ReTr plaquette with SF BC.
  <ul>
  <li>Plaquette is NOT normalized at all.
  <li>Contribution from the spatial plaquettes at t=0 and t=T are not included since they are always unity.
  </ul>
*/
double Staple_SF::plaquette(const Field_G& U)
{
  if (!m_initialized) {
    vout.crucial(m_vl, "Error at %s: Parameter is not initialized.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
  //  return (plaq_s(U) +plaq_t(U))/2;
  return(plaq_s(U) + plaq_t(U));
}


//====================================================================

/*!
  Evaluate summed ReTr plaquette with SF BC.
  <ul>
  <li>The temporal plaquette attached to the boundary is multiplied with the improvement factor ct.
<pre>
     +---+
  ct |   |
 t=0 x---+
</pre>
  <li>Plaquette is NOT normalized at all.
  <li>Contribution from the spatial plaquettes at t=0 and t=T are not included since they are always unity.
  </ul>
*/
double Staple_SF::plaquette_ct(const Field_G& U, const double ct)
{
  if (!m_initialized) {
    vout.crucial(m_vl, "Error at %s: Parameter is not initialized.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
  return(plaq_s(U) + plaq_t_ct(U, ct));
}


//====================================================================

/*!
  Evaluate summed spatial ReTr plaquette with SF BC.
  <ul>
  <li>Plaquette is NOT normalized at all.
  <li>Contribution from the spatial plaquettes at t=0 and t=T are not included since they are always unity.
  <li>The function upper(U,mu=k,nu) is modified to give zero.
  </ul>
*/
double Staple_SF::plaq_s(const Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();

  Field_G_SF staple;
  double     plaq = 0.0;

  upper(staple, U, 0, 1);
  for (int site = 0; site < Nvol; site++) {
    plaq += ReTr(U.mat(site, 0) * staple.mat_dag(site));   // P_xy
  }

  upper(staple, U, 1, 2);
  for (int site = 0; site < Nvol; site++) {
    plaq += ReTr(U.mat(site, 1) * staple.mat_dag(site));   // P_yz
  }

  upper(staple, U, 2, 0);
  for (int site = 0; site < Nvol; site++) {
    plaq += ReTr(U.mat(site, 2) * staple.mat_dag(site));   // P_zx
  }

  plaq = Communicator::reduce_sum(plaq);

  return plaq;
}


//====================================================================

/*!
  Evaluate summed temporal ReTr plaquette with SF BC.
  <ul>
  <li>Plaquette is NOT normalized at all.
  <li>The temporal plaquettes attached to t=0 and t=T boundary are evaluated with Wk and Wk'.
  <li>The function lower(U,3,nu) is modified to use Wk and Wk'.
  </ul>
*/
double Staple_SF::plaq_t(const Field_G& U)
{
  if (!m_initialized) {
    vout.crucial(m_vl, "Error at %s: Parameter is not initialized.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  const int Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();

  Field_G_SF staple;
  double     plaq = 0.0;

  for (int nu = 0; nu < Ndim - 1; nu++) {
    lower(staple, U, 3, nu);
    //    staple = upper(U,3,nu);
    for (int site = 0; site < Nvol; site++) {
      plaq += ReTr(U.mat(site, 3) * staple.mat_dag(site));   // P_tk
    }
  }

  plaq = Communicator::reduce_sum(plaq);

  return plaq;
}


//====================================================================

/*!
  Evaluate summed temporal ReTr plaquette with SF BC.
  <ul>
  <li>The temporal plaquette attached to the boundary is multiplied with ct.
<pre>
     +---+
  ct |   |
 t=0 x---+
</pre>
  <li>The temporal plaquettes attached to t=0 and t=T boundary are evaluated with Wk and Wk'.
  <li>The function lower(U,3,nu) is modified to use Wk and Wk'.
  <li>Plaquette is NOT normalized at all.
  </ul>
*/
double Staple_SF::plaq_t_ct(const Field_G& U, const double ct)
{
  if (!m_initialized) {
    vout.crucial(m_vl, "Error at %s: Parameter is not initialized.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  const int Ndim = CommonParameters::Ndim();
  const int Nt   = CommonParameters::Nt();
  const int Nvol = CommonParameters::Nvol();
  const int NPEt = CommonParameters::NPEt();

  Field_G_SF staple;
  double     plaq = 0.0;

  for (int nu = 0; nu < Ndim - 1; nu++) {
    lower(staple, U, 3, nu);
    // If the node is at the boundary the temporal plaquette is multiplied with ct.
    if (Communicator::ipe(3) == 0) {
      staple.mult_ct_boundary(0, ct);
    }
    if (Communicator::ipe(3) == NPEt - 1) {
      staple.mult_ct_boundary(Nt - 1, ct);
    }
    for (int site = 0; site < Nvol; site++) {
      plaq += ReTr(U.mat(site, 3) * staple.mat_dag(site));   // P_tk
    }
  }

  plaq = Communicator::reduce_sum(plaq);

  return plaq;
}


//====================================================================

/*!
  Evaluate staple for all the links in mu direction with SF BC.
<pre>
   (1)  mu (2)
      +-->--+
   nu |     |
     i+     +

      +     +
   nu |     |
     i+-->--+
    (1)  mu (2)
</pre>
  <ul>
  <li>We should notice that all the staple attached to the boundary spatial link is set to zero.
  <li>This staple at the boundary shall not be used since the corresponding force for the boundary spatial link is always set to zero.
  </ul>
*/
void Staple_SF::staple(Field_G& W, const Field_G& U, const int mu)
{
  if (!m_initialized) {
    vout.crucial(m_vl, "Error at %s: Parameter is not initialized.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  const int Ndim = CommonParameters::Ndim();

  W.set(0.0);

  for (int nu = 0; nu < Ndim; nu++) {
    if (nu != mu) {
      Field_G_SF c_tmp;

      upper(c_tmp, U, mu, nu);
      axpy(W, 1.0, c_tmp);

      lower(c_tmp, U, mu, nu);
      axpy(W, 1.0, c_tmp);
    }
  }
}


//====================================================================

/*!
  Evaluate staple for all the links in mu direction with SF BC and boundary improvement factor ct.
<pre>
   (1)  mu (2)
      +-->--+
   nu |     |
     i+     +

      +     +
   nu |     |
     i+-->--+
    (1)  mu (2)
</pre>
  <ul>
  <li>staple attached to spatial boundary link is set to zero.
  <li>contributions from the temporal plaquette attached to the boundary is multiplied with ct.
  </ul>
*/
void Staple_SF::staple_ct(Field_G& W, const Field_G& U, const int mu, const double ct)
{
  const int Ndim = CommonParameters::Ndim();
  const int Nt   = CommonParameters::Nt();
  const int NPEt = CommonParameters::NPEt();

  W.set(0.0);

  for (int nu = 0; nu < Ndim; nu++) {
    if (nu != mu) {
      Field_G_SF staple_upper;
      Field_G_SF staple_lower;

      upper(staple_upper, U, mu, nu);
      lower(staple_lower, U, mu, nu);

      if (Communicator::ipe(3) == 0) {
        if (mu == 3) {
          staple_upper.mult_ct_boundary(0, ct);
          staple_lower.mult_ct_boundary(0, ct);
        }
        if (nu == 3) {
          staple_lower.mult_ct_boundary(1, ct);
        }
      }

      if (Communicator::ipe(3) == NPEt - 1) {
        if (mu == 3) {
          staple_upper.mult_ct_boundary(Nt - 1, ct);
          staple_lower.mult_ct_boundary(Nt - 1, ct);
        }
        if (nu == 3) {
          staple_upper.mult_ct_boundary(Nt - 1, ct);
        }
      }

      axpy(W, 1.0, staple_upper);
      axpy(W, 1.0, staple_lower);
    }
  }
}


//====================================================================
void Staple_SF::upper(Field_G_SF& c, const Field_G& U, const int mu, const int nu)
{
  if (!m_initialized) {
    vout.crucial(m_vl, "Error at %s: Parameter is not initialized.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  const int Nvol = CommonParameters::Nvol();

  // (1)  mu (2)
  //    +-->--+
  // nu |     |
  //   i+     +
  Field_G_SF Umu;
  Umu.setpart_ex(0, U, mu);

  Field_G_SF Unu;
  Unu.setpart_ex(0, U, nu);
  if (mu != 3) Umu.set_boundary_wk(m_wk);
  if (nu != 3) Unu.set_boundary_wk(m_wk);

  Field_G_SF v;
  m_shift.backward(v, Unu, mu);
  m_shift.backward(c, Umu, nu);
  if (mu == 3) v.set_boundary_wkpr(m_wkpr);
  if (nu == 3) c.set_boundary_wkpr(m_wkpr);

  Field_G_SF w;
  mult_Field_Gnd(w, 0, c, 0, v, 0);
  mult_Field_Gnn(c, 0, Unu, 0, w, 0);
  if (mu != 3) c.set_boundary_zero();
}


//====================================================================
void Staple_SF::lower(Field_G_SF& c, const Field_G& U, const int mu, const int nu)
{
  if (!m_initialized) {
    vout.crucial(m_vl, "Error at %s: Parameter is not initialized.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  const int Nvol = CommonParameters::Nvol();

  //    +     +
  // nu |     |
  //   i+-->--+
  //  (1)  mu (2)
  Field_G_SF Umu;
  Umu.setpart_ex(0, U, mu);

  Field_G_SF Unu;
  Unu.setpart_ex(0, U, nu);
  if (mu != 3) Umu.set_boundary_wk(m_wk);
  if (nu != 3) Unu.set_boundary_wk(m_wk);

  Field_G_SF w;
  m_shift.backward(w, Unu, mu);
  if (mu == 3) w.set_boundary_wkpr(m_wkpr);

  Field_G_SF v;
  mult_Field_Gnn(v, 0, Umu, 0, w, 0);
  mult_Field_Gdn(w, 0, Unu, 0, v, 0);

  m_shift.forward(c, w, nu);
  if (mu != 3) c.set_boundary_zero();
}


//====================================================================

/*!
  Print out the plaquette for two definitions.
\f[
\frac{1}{3}\frac{1}{6L_xL_yL_zL_t-3L_xL_yL_z}\sum_{p\notin{\rm boundary}}U_p
\f]
  <ul>
  <li>Contributions from the boundary spatial plaquettes is not included.
  <li>Normalized with the number of plaquettes Lx*Ly*Lz*(6Lt-3) and color factor 3.
  </ul>
\f[
\frac{1}{3}\frac{1}{6L_xL_yL_zL_t}\sum_{p}U_p
\f]
  <ul>
  <li>Contributions from the boundary spatial plaquettes is included.
  <li>Normalized with the number of plaquettes 6*Lx*Ly*Lz*Lt and color factor 3.
  </ul>
 */
void Staple_SF::print_plaquette(const Field_G& U)
{
  const int Lx = CommonParameters::Lx();
  const int Ly = CommonParameters::Ly();
  const int Lz = CommonParameters::Lz();
  const int Lt = CommonParameters::Lt();

  const double plaq  = plaquette(U);
  const double plaq2 = plaq + 3 * 3 * Lx * Ly * Lz;

  vout.general(m_vl, "plaq_SF without boundary spatial plaq = %.8f\n",
               plaq / (3 * Lx * Ly * Lz * (6 * Lt - 3)));
  vout.general(m_vl, "plaq_SF with boundary spatial plaq = %.8f\n",
               plaq2 / (3 * 6 * Lx * Ly * Lz * Lt));
}


//============================================================END=====
