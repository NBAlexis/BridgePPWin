/*!
        @file    field_G_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/


#ifndef FIELD_G_SF_INCLUDED
#define FIELD_G_SF_INCLUDED

#include "Parameters/parameters.h"
#include "field_G.h"

//! SU(N) gauge field class BAPI in which a few functions are added for the SF.

/*!
  This class BAPI defines SU(N) gauge field, which is used such as gauge configuration.
  <ul>
  <li>A derived class BAPI from Field_G in order to add a few functions to manipulate the boundary link variables.
  <li>An inheritance was adopted since we need to manipulate the Field contents and upcast into Filed_G object.
  <li>[23 Mar 2012 Y.Taniguchi]
  </ul>
 */


class BAPI Field_G_SF : public Field_G
{
 private:
  //! number of color elements
  int m_Nc;
  //! number of components as real values
  int m_Ndf;
  //! number of sites
  int m_Nvol;
  //! extra degree of freedom, such as mu, nu.
  int m_Nex;

  int m_Nx, m_Ny, m_Nz, m_Nt, m_Svol;
  int m_NPEt;

  //! SF boundary condition at t=0
  Mat_SU_N m_wk;
  //! SF boundary condition at t=Nt
  Mat_SU_N m_wkpr;

 public:
  //  Field_G_SF(){}
  Field_G_SF(const int Nvol = CommonParameters::Nvol(), const int Nex = 1) :
    m_Nc(CommonParameters::Nc()), m_Nvol(Nvol), m_Nex(Nex),
    m_wk(m_Nc), m_wkpr(m_Nc)
  {
    m_Nx   = CommonParameters::Nx();
    m_Ny   = CommonParameters::Ny();
    m_Nz   = CommonParameters::Nz();
    m_Nt   = CommonParameters::Nt();
    m_NPEt = CommonParameters::NPEt();
    m_Svol = m_Nx * m_Ny * m_Nz;
    m_Ndf  = 2 * m_Nc * m_Nc;
    Field::reset(m_Ndf, m_Nvol, m_Nex);
  }

  Field_G_SF(const Field& x) :
    Field_G(x),
    m_Nc(CommonParameters::Nc()), m_Ndf(x.nin()), m_Nvol(x.nvol()), m_Nex(x.nex()),
    m_wk(m_Nc), m_wkpr(m_Nc)
  {
    assert(m_Ndf == 2 * m_Nc * m_Nc);
    m_Nx   = CommonParameters::Nx();
    m_Ny   = CommonParameters::Ny();
    m_Nz   = CommonParameters::Nz();
    m_Nt   = CommonParameters::Nt();
    m_NPEt = CommonParameters::NPEt();
    m_Svol = m_Nx * m_Ny * m_Nz;
  }

  Field_G_SF(double *phi, double *phipr) :
    m_Nc(CommonParameters::Nc()), m_wk(m_Nc), m_wkpr(m_Nc)
  {
    const int Lx = CommonParameters::Lx();

    m_Nex  = 1;
    m_Nvol = CommonParameters::Nvol();
    m_Nx   = CommonParameters::Nx();
    m_Ny   = CommonParameters::Ny();
    m_Nz   = CommonParameters::Nz();
    m_Nt   = CommonParameters::Nt();
    m_NPEt = CommonParameters::NPEt();
    m_Svol = m_Nx * m_Ny * m_Nz;
    m_Ndf  = 2 * m_Nc * m_Nc;

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

  Field_G_SF& operator=(const Field_G_SF& v) { copy(*this, v); return *this; }

  //! Set the parameter by giving vector objects.
  void set_parameters(const std::vector<double>& phi,
                      const std::vector<double>& phipr);

  //! Set the parameter for the boundary links Wk, Wk'.
  void set_parameters(double *phi, double *phipr);

  //! Set the boundary spatial link at t=0 for SF bc.
  void set_boundary_wk(const Mat_SU_N& U);

  //! Set the boundary spatial link at t=Nt-1 for SF bc.
  void set_boundary_wkpr(const Mat_SU_N& U);

  //! Set the boundary matrix to 0 for SF bc.
  void set_boundary_zero();

  //! Set the boundary spatial link to 0 for SF bc.
  void set_boundary_spatial_link_zero();

  //! Multiply the boundary improvement factor ct or ctr to an SU(N) matrix object which belongs to a site at t.
  void mult_ct_boundary(const int t, const double ct);


  //! Set the boundary spatial link at t=0 for SF bc.
  void set_boundary_wk(Field_G& f);

  //! Set the boundary spatial link at t=Nt-1 for SF bc.
  void set_boundary_wkpr(Field_G& f);

  //! Set the boundary matrix to 0 for SF bc.
  void set_boundary_zero(Field_G& f);

  //! Set the boundary spatial link to 0 for SF bc.
  void set_boundary_spatial_link_zero(Field_G& f);
};
#endif
