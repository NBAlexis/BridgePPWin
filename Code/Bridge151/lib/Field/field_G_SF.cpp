/*!
        @file    field_G_SF.cpp

        @brief

        @author  "Taniguchi Yusuke"  (tanigchi)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "field_G_SF.h"

//====================================================================

/*!
  Set the boundary spatial link \f$U_k(t=0)\f$ to its proper Dirichlet bounday \f$W_k\f$.
  <ul>
  <li>Supposed to be used for u_mu given by u_mu.setpart_ex(0, *g,mu) when mu is a spatial.
  <li>Lorentz index is set to mn=0.
  <li>Wk is given as an argument U.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
void Field_G_SF::set_boundary_wk(const Mat_SU_N& U)
{
  const int mn = 0;

  if (Communicator::ipe(3) == 0) {
    for (int site = 0; site < m_Svol; ++site) {
      for (int cc = 0; cc < m_Nc * m_Nc; ++cc) {
        field[myindex(2 * cc, site, mn)]     = U.r(cc);
        field[myindex(2 * cc + 1, site, mn)] = U.i(cc);
      }
    }
  }
}


//====================================================================

/*!
  Set the boundary spatial link \f$U_k(t=N_t-1)\f$ to its proper Dirichlet bounday \f$W_k'\f$.
  <ul>
  <li>Supposed to be used for v given by shift.backward(v,u_nu,mu) when mu is a temporal.
  <li>Lorentz index is set to mn=0.
  <li>Wk' is given as an argument U.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
void Field_G_SF::set_boundary_wkpr(const Mat_SU_N& U)
{
  const int mn = 0;

  if (Communicator::ipe(3) == m_NPEt - 1) {
    for (int site = m_Nvol - m_Svol; site < m_Nvol; ++site) {
      for (int cc = 0; cc < m_Nc * m_Nc; ++cc) {
        field[myindex(2 * cc, site, mn)]     = U.r(cc);
        field[myindex(2 * cc + 1, site, mn)] = U.i(cc);
      }
    }
  }
}


//====================================================================

/*!
  Set the boundary matrix to zero.
  <ul>
  <li>Supposed to be used for a staple for the boundary spatial plaquette in order to keep the corresponding force to zero.
  <li>Lorentz index is set to mn=0.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
void Field_G_SF::set_boundary_zero()
{
  const int mn = 0;

  if (Communicator::ipe(3) == 0) {
    for (int site = 0; site < m_Svol; ++site) {
      for (int cc = 0; cc < 2 * m_Nc * m_Nc; ++cc) {
        field[myindex(cc, site, mn)] = 0.0;
      }
    }
  }
}


//====================================================================

/*!
  Set the boundary spatial link to zero.
  <ul>
  <li>Supposed to be used for a force for the boundary spatial link in order to keep it zero.
  <li>Lorentz index mn runs for spatial value.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
void Field_G_SF::set_boundary_spatial_link_zero()
{
  if (Communicator::ipe(3) == 0) {
    for (int mn = 0; mn < 3; ++mn) {
      for (int site = 0; site < m_Svol; ++site) {
        for (int cc = 0; cc < 2 * m_Nc * m_Nc; ++cc) {
          field[myindex(cc, site, mn)] = 0.0;
        }
      }
    }
  }
}


//====================================================================

/*!
  Multiply the boundary improvement factor ct or ctr to an SU(N) matrix object attached to the boundary t=0 or t=Nt.
  <ul>
  <li>Supposed to be used for a matrix given by upper(g,mu,nu) or lower(g,mu,nu) when it is attached to the t=0 boundary.
  <li>The SU(N) matrix object may belong to a site at t=0 or t=1 or t=Nt.
  <li>This function performs no check if the object is a proper one at the boundary.
  <li>Lorentz index is set to mn=0.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
void Field_G_SF::mult_ct_boundary(const int t, const double ct)
{
  const int mn  = 0;
  const int ini = m_Nx * m_Ny * m_Nz * t;
  const int fin = m_Nx * m_Ny * m_Nz * (t + 1);

  if (((Communicator::ipe(3) == 0) && ((t == 0) || (t == 1))) || ((Communicator::ipe(3) == m_NPEt - 1) && (t == m_Nt - 1))) {
    for (int site = ini; site < fin; ++site) {
      for (int cc = 0; cc < 2 * m_Nc * m_Nc; ++cc) {
        field[myindex(cc, site, mn)] *= ct;
      }
    }
  }
}


//====================================================================

/*!
  Set the boundary spatial link \f$U_k(t=0)\f$ to its proper Dirichlet bounday \f$W_k\f$.
  <ul>
  <li>Supposed to be used for u_mu given by u_mu.setpart_ex(0, *g,mu) when mu is a spatial.
  <li>Lorentz index is set to mn=0.
  <li>Field is given as an argument.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
void Field_G_SF::set_boundary_wk(Field_G& f)
{
  const int mn = 0;

  assert(f.nex() == 1);
  if (Communicator::ipe(3) == 0) {
    for (int site = 0; site < m_Svol; ++site) {
      f.set_mat(site, mn, m_wk);
    }
  }
}


//====================================================================

/*!
  Set the boundary spatial link \f$U_k(t=N_t-1)\f$ to its proper Dirichlet bounday \f$W_k'\f$.
  <ul>
  <li>Supposed to be used for v given by shift.backward(v,u_nu,mu) when mu is a temporal.
  <li>Lorentz index is set to mn=0.
  <li>Field is given as an argument.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
void Field_G_SF::set_boundary_wkpr(Field_G& f)
{
  const int mn = 0;

  assert(f.nex() == 1);
  if (Communicator::ipe(3) == m_NPEt - 1) {
    for (int site = m_Nvol - m_Svol; site < m_Nvol; ++site) {
      f.set_mat(site, mn, m_wkpr);
    }
  }
}


//====================================================================

/*!
  Set the boundary matrix to zero.
  <ul>
  <li>Supposed to be used for a staple for the boundary spatial plaquette in order to keep the corresponding force to zero.
  <li>Lorentz index is set to mn=0.
  <li>Field is given as an argument.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
void Field_G_SF::set_boundary_zero(Field_G& f)
{
  const int mn = 0;

  assert(f.nex() == 1);
  if (Communicator::ipe(3) == 0) {
    for (int site = 0; site < m_Svol; ++site) {
      for (int cc = 0; cc < m_Nc * m_Nc; ++cc) {
        f.set_ri(cc, site, mn, 0.0, 0.0);
      }
    }
  }
}


//====================================================================

/*!
  Set the boundary spatial link to zero.
  <ul>
  <li>Supposed to be used for a force for the boundary spatial link in order to keep it zero.
  <li>Lorentz index mn runs for spatial value.
  <li>Field is given as an argument.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
void Field_G_SF::set_boundary_spatial_link_zero(Field_G& f)
{
  assert(f.nex() >= 3);
  if (Communicator::ipe(3) == 0) {
    for (int mn = 0; mn < 3; ++mn) {
      for (int site = 0; site < m_Svol; ++site) {
        for (int cc = 0; cc < m_Nc * m_Nc; ++cc) {
          f.set_ri(cc, site, mn, 0.0, 0.0);
        }
      }
    }
  }
}


//====================================================================
void Field_G_SF::set_parameters(const std::vector<double>& phi,
                                const std::vector<double>& phipr)
{
  double aphi[3];
  double aphipr[3];

  for (int i = 0; i < 3; ++i) {
    aphi[i]   = phi[i];
    aphipr[i] = phipr[i];
  }

  set_parameters(aphi, aphipr);
}


//====================================================================
void Field_G_SF::set_parameters(double *phi, double *phipr)
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
}


//====================================================================
