#include "BridgeLib_Private.h"

/*!
        @file    $Id:: fopr_Clover_SF.cpp #$

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fopr_Clover_SF.h"



#ifdef USE_FACTORY
namespace {
  Fopr *create_object()
  {
    return new Fopr_Clover_SF();
  }


  bool init = Fopr::Factory_noarg::Register("Clover_SF", create_object);
}
#endif



const std::string Fopr_Clover_SF::class_name = "Fopr_Clover_SF";

//====================================================================
void Fopr_Clover_SF::set_parameters(const Parameters& params)
{
  const std::string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  std::string         str_gmset_type;
  double              kappa, cSW;
  std::vector<int>    bc;
  std::vector<double> phi, phipr;

  int err          = 0;
  int err_optional = 0;
  err_optional += params.fetch_string("gamma_matrix_type", str_gmset_type);
  err          += params.fetch_double("hopping_parameter", kappa);
  err          += params.fetch_double("clover_coefficient", cSW);
  err          += params.fetch_int_vector("boundary_condition", bc);
  err          += params.fetch_double_vector("phi", phi);
  err          += params.fetch_double_vector("phipr", phipr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  m_repr = str_gmset_type;

  set_parameters(kappa, cSW, bc, &phi[0], &phipr[0]);
}


//====================================================================
void Fopr_Clover_SF::set_parameters(double kappa, double cSW, std::vector<int> bc,
                                    double *phi, double *phipr)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }
  vout.general(m_vl, "  phi1  = %12.8f\n", phi[0]);
  vout.general(m_vl, "  phi2  = %12.8f\n", phi[1]);
  vout.general(m_vl, "  phi3  = %12.8f\n", phi[2]);
  vout.general(m_vl, "  phipr1= %12.8f\n", phipr[0]);
  vout.general(m_vl, "  phipr2= %12.8f\n", phipr[1]);
  vout.general(m_vl, "  phipr3= %12.8f\n", phipr[2]);

  //- range check
  // NB. kappa,cSW == 0 is allowed.
  // NB. phi,phipr == 0 is allowed.
  assert(bc.size() == m_Ndim);

  //- store values
  m_kappa = kappa;
  m_cSW   = cSW;

  for (int i = 0; i < 3; ++i) {
    m_phi[i]   = phi[i];
    m_phipr[i] = phipr[i];
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
  }

  //- propagate parameters
  m_fopr_w->set_parameters(m_kappa, m_boundary);
}


//====================================================================
namespace {
  inline double mult_uv_r(const double *u, const double *v)
  {
    return u[0] * v[0] - u[1] * v[1]
           + u[2] * v[2] - u[3] * v[3]
           + u[4] * v[4] - u[5] * v[5];
  }


  inline double mult_uv_i(const double *u, const double *v)
  {
    return u[0] * v[1] + u[1] * v[0]
           + u[2] * v[3] + u[3] * v[2]
           + u[4] * v[5] + u[5] * v[4];
  }
}

//====================================================================
void Fopr_Clover_SF::init(std::string repr)
{
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nc   = CommonParameters::Nc();
  m_Nd   = CommonParameters::Nd();
  m_NinF = 2 * m_Nc * m_Nd;

  m_U = 0;

  m_repr = repr;

  m_boundary.resize(m_Ndim);
  m_GM.resize(m_Ndim + 1);
  m_SG.resize(m_Ndim * m_Ndim);

  GammaMatrixSet *gmset = GammaMatrixSet::New(m_repr);

  m_GM[0] = gmset->get_GM(gmset->GAMMA1);
  m_GM[1] = gmset->get_GM(gmset->GAMMA2);
  m_GM[2] = gmset->get_GM(gmset->GAMMA3);
  m_GM[3] = gmset->get_GM(gmset->GAMMA4);
  m_GM[4] = gmset->get_GM(gmset->GAMMA5);

  m_SG[sg_index(0, 1)] = gmset->get_GM(gmset->SIGMA12);
  m_SG[sg_index(1, 2)] = gmset->get_GM(gmset->SIGMA23);
  m_SG[sg_index(2, 0)] = gmset->get_GM(gmset->SIGMA31);
  m_SG[sg_index(3, 0)] = gmset->get_GM(gmset->SIGMA41);
  m_SG[sg_index(3, 1)] = gmset->get_GM(gmset->SIGMA42);
  m_SG[sg_index(3, 2)] = gmset->get_GM(gmset->SIGMA43);

  m_SG[sg_index(1, 0)] = m_SG[sg_index(0, 1)].mult(-1);
  m_SG[sg_index(2, 1)] = m_SG[sg_index(1, 2)].mult(-1);
  m_SG[sg_index(0, 2)] = m_SG[sg_index(2, 0)].mult(-1);
  m_SG[sg_index(0, 3)] = m_SG[sg_index(3, 0)].mult(-1);
  m_SG[sg_index(1, 3)] = m_SG[sg_index(3, 1)].mult(-1);
  m_SG[sg_index(2, 3)] = m_SG[sg_index(3, 2)].mult(-1);

  m_SG[sg_index(0, 0)] = gmset->get_GM(gmset->UNITY);
  m_SG[sg_index(1, 1)] = gmset->get_GM(gmset->UNITY);
  m_SG[sg_index(2, 2)] = gmset->get_GM(gmset->UNITY);
  m_SG[sg_index(3, 3)] = gmset->get_GM(gmset->UNITY);
  // these 4 gamma matrices are actually not used.


  //  m_fopr_w = new Fopr_Wilson_SF(repr);
  // Dirac reps. only!
  m_fopr_w = new Fopr_Wilson_SF();

  if (m_repr == "Dirac") {
    m_csw = &Fopr_Clover_SF::mult_csw_dirac;
    // }else if(m_repr == "Chiral"){
    // m_csw = &Fopr_Clover_SF::mult_csw_chiral;
  }

  delete gmset;
}


//====================================================================
void Fopr_Clover_SF::tidyup()
{
  delete m_fopr_w;
}


//====================================================================
void Fopr_Clover_SF::DdagD(Field& w, const Field& f)
{
  assert(f.nex() == 1);
  Field w2(f.nin(), f.nvol(), 1);

  D(w2, f);
  mult_gm5(w, w2);
  D(w2, w);
  mult_gm5(w, w2);
}


//====================================================================
void Fopr_Clover_SF::Ddag(Field& w, const Field& f)
{
  assert(f.nex() == 1);
  Field w2(f.nin(), f.nvol(), 1);

  mult_gm5(w, f);
  D(w2, w);
  mult_gm5(w, w2);
}


//====================================================================
void Fopr_Clover_SF::H(Field& w, const Field& f)
{
  assert(f.nex() == 1);
  Field w2(f.nin(), f.nvol(), 1);

  D(w2, f);
  mult_gm5(w, w2);
}


//====================================================================
void Fopr_Clover_SF::D(Field& w, const Field& f)
{
  assert(f.nex() == 1);
  int Nvol = f.nvol();

  Field_F w2(Nvol, 1);

  //  Field_F w3(f);
  //  setzero.set_boundary_zero(w3);
  //  m_fopr_w->D(w,w3);
  //  mult_csw(w2,w3);

  m_fopr_w->D(w, f);
  mult_csw(w2, (Field_F)f);
  axpy(w, -1.0, w2); // w -= (Field)w2;
  setzero.set_boundary_zero(w);

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_SF::mult_isigma(Field_F& v, const Field_F& w,
                                 const int mu, const int nu)
{
  assert(mu != nu);

  mult_iGM(v, m_SG[sg_index(mu, nu)], w);
}


//====================================================================
void Fopr_Clover_SF::mult_csw(Field_F& v, const Field_F& w)
{
  (this->*m_csw)(v, w);
}


//====================================================================
void Fopr_Clover_SF::mult_csw_chiral(Field_F& v, const Field_F& w)
{
  assert(w.nex() == 1);

  int Nc   = CommonParameters::Nc();
  int Nvc  = 2 * Nc;
  int Ndf  = 2 * Nc * Nc;
  int Nd   = CommonParameters::Nd();
  int Nvol = w.nvol();

  int id1 = 0;
  int id2 = Nvc;
  int id3 = Nvc * 2;
  int id4 = Nvc * 3;

  const double *w2 = w.ptr(0);
  double       *v2 = v.ptr(0);

  double *Bx = m_Bx.ptr(0);
  double *By = m_By.ptr(0);
  double *Bz = m_Bz.ptr(0);
  double *Ex = m_Ex.ptr(0);
  double *Ey = m_Ey.ptr(0);
  double *Ez = m_Ez.ptr(0);

  v.set(0.0);

  for (int site = 0; site < Nvol; ++site) {
    int iv = Nvc * Nd * site;
    int ig = Ndf * site;

    for (int ic = 0; ic < Nc; ++ic) {
      int icr = 2 * ic;
      int ici = icr + 1;
      int icg = ic * Nvc + ig;

      // isigma_23 * Bx
      v2[icr + id1 + iv] -= mult_uv_i(&Bx[icg], &w2[id2 + iv]);
      v2[ici + id1 + iv] += mult_uv_r(&Bx[icg], &w2[id2 + iv]);
      v2[icr + id2 + iv] -= mult_uv_i(&Bx[icg], &w2[id1 + iv]);
      v2[ici + id2 + iv] += mult_uv_r(&Bx[icg], &w2[id1 + iv]);

      v2[icr + id3 + iv] -= mult_uv_i(&Bx[icg], &w2[id4 + iv]);
      v2[ici + id3 + iv] += mult_uv_r(&Bx[icg], &w2[id4 + iv]);
      v2[icr + id4 + iv] -= mult_uv_i(&Bx[icg], &w2[id3 + iv]);
      v2[ici + id4 + iv] += mult_uv_r(&Bx[icg], &w2[id3 + iv]);

      // isigma_31 * By
      v2[icr + id1 + iv] += mult_uv_r(&By[icg], &w2[id2 + iv]);
      v2[ici + id1 + iv] += mult_uv_i(&By[icg], &w2[id2 + iv]);
      v2[icr + id2 + iv] -= mult_uv_r(&By[icg], &w2[id1 + iv]);
      v2[ici + id2 + iv] -= mult_uv_i(&By[icg], &w2[id1 + iv]);

      v2[icr + id3 + iv] += mult_uv_r(&By[icg], &w2[id4 + iv]);
      v2[ici + id3 + iv] += mult_uv_i(&By[icg], &w2[id4 + iv]);
      v2[icr + id4 + iv] -= mult_uv_r(&By[icg], &w2[id3 + iv]);
      v2[ici + id4 + iv] -= mult_uv_i(&By[icg], &w2[id3 + iv]);

      // isigma_12 * Bz
      v2[icr + id1 + iv] -= mult_uv_i(&Bz[icg], &w2[id1 + iv]);
      v2[ici + id1 + iv] += mult_uv_r(&Bz[icg], &w2[id1 + iv]);
      v2[icr + id2 + iv] += mult_uv_i(&Bz[icg], &w2[id2 + iv]);
      v2[ici + id2 + iv] -= mult_uv_r(&Bz[icg], &w2[id2 + iv]);

      v2[icr + id3 + iv] -= mult_uv_i(&Bz[icg], &w2[id3 + iv]);
      v2[ici + id3 + iv] += mult_uv_r(&Bz[icg], &w2[id3 + iv]);
      v2[icr + id4 + iv] += mult_uv_i(&Bz[icg], &w2[id4 + iv]);
      v2[ici + id4 + iv] -= mult_uv_r(&Bz[icg], &w2[id4 + iv]);

      // isigma_41 * Ex
      v2[icr + id1 + iv] += mult_uv_i(&Ex[icg], &w2[id2 + iv]);
      v2[ici + id1 + iv] -= mult_uv_r(&Ex[icg], &w2[id2 + iv]);
      v2[icr + id2 + iv] += mult_uv_i(&Ex[icg], &w2[id1 + iv]);
      v2[ici + id2 + iv] -= mult_uv_r(&Ex[icg], &w2[id1 + iv]);

      v2[icr + id3 + iv] -= mult_uv_i(&Ex[icg], &w2[id4 + iv]);
      v2[ici + id3 + iv] += mult_uv_r(&Ex[icg], &w2[id4 + iv]);
      v2[icr + id4 + iv] -= mult_uv_i(&Ex[icg], &w2[id3 + iv]);
      v2[ici + id4 + iv] += mult_uv_r(&Ex[icg], &w2[id3 + iv]);

      // isigma_42 * Ey
      v2[icr + id1 + iv] -= mult_uv_r(&Ey[icg], &w2[id2 + iv]);
      v2[ici + id1 + iv] -= mult_uv_i(&Ey[icg], &w2[id2 + iv]);
      v2[icr + id2 + iv] += mult_uv_r(&Ey[icg], &w2[id1 + iv]);
      v2[ici + id2 + iv] += mult_uv_i(&Ey[icg], &w2[id1 + iv]);

      v2[icr + id3 + iv] += mult_uv_r(&Ey[icg], &w2[id4 + iv]);
      v2[ici + id3 + iv] += mult_uv_i(&Ey[icg], &w2[id4 + iv]);
      v2[icr + id4 + iv] -= mult_uv_r(&Ey[icg], &w2[id3 + iv]);
      v2[ici + id4 + iv] -= mult_uv_i(&Ey[icg], &w2[id3 + iv]);

      // isigma_43 * Ez
      v2[icr + id1 + iv] += mult_uv_i(&Ez[icg], &w2[id1 + iv]);
      v2[ici + id1 + iv] -= mult_uv_r(&Ez[icg], &w2[id1 + iv]);
      v2[icr + id2 + iv] -= mult_uv_i(&Ez[icg], &w2[id2 + iv]);
      v2[ici + id2 + iv] += mult_uv_r(&Ez[icg], &w2[id2 + iv]);

      v2[icr + id3 + iv] -= mult_uv_i(&Ez[icg], &w2[id3 + iv]);
      v2[ici + id3 + iv] += mult_uv_r(&Ez[icg], &w2[id3 + iv]);
      v2[icr + id4 + iv] += mult_uv_i(&Ez[icg], &w2[id4 + iv]);
      v2[ici + id4 + iv] -= mult_uv_r(&Ez[icg], &w2[id4 + iv]);
    }
  }

  scal(v, m_kappa * m_cSW); // v *= m_kappa * m_cSW;
}


//====================================================================
void Fopr_Clover_SF::mult_csw_dirac(Field_F& v, const Field_F& w)
{
  assert(w.nex() == 1);

  int Nc   = CommonParameters::Nc();
  int Nvc  = 2 * Nc;
  int Ndf  = 2 * Nc * Nc;
  int Nd   = CommonParameters::Nd();
  int Nvol = w.nvol();

  int id1 = 0;
  int id2 = Nvc;
  int id3 = Nvc * 2;
  int id4 = Nvc * 3;

  const double *w2 = w.ptr(0);
  double       *v2 = v.ptr(0);

  double *Bx = m_Bx.ptr(0);
  double *By = m_By.ptr(0);
  double *Bz = m_Bz.ptr(0);
  double *Ex = m_Ex.ptr(0);
  double *Ey = m_Ey.ptr(0);
  double *Ez = m_Ez.ptr(0);

  v.set(0.0);

  for (int site = 0; site < Nvol; ++site) {
    int iv = Nvc * Nd * site;
    int ig = Ndf * site;

    for (int ic = 0; ic < Nc; ++ic) {
      int icr = 2 * ic;
      int ici = icr + 1;
      int icg = ic * Nvc + ig;

      // isigma_23 * Bx
      v2[icr + id1 + iv] -= mult_uv_i(&Bx[icg], &w2[id2 + iv]);
      v2[ici + id1 + iv] += mult_uv_r(&Bx[icg], &w2[id2 + iv]);
      v2[icr + id2 + iv] -= mult_uv_i(&Bx[icg], &w2[id1 + iv]);
      v2[ici + id2 + iv] += mult_uv_r(&Bx[icg], &w2[id1 + iv]);

      v2[icr + id3 + iv] -= mult_uv_i(&Bx[icg], &w2[id4 + iv]);
      v2[ici + id3 + iv] += mult_uv_r(&Bx[icg], &w2[id4 + iv]);
      v2[icr + id4 + iv] -= mult_uv_i(&Bx[icg], &w2[id3 + iv]);
      v2[ici + id4 + iv] += mult_uv_r(&Bx[icg], &w2[id3 + iv]);

      // isigma_31 * By
      v2[icr + id1 + iv] += mult_uv_r(&By[icg], &w2[id2 + iv]);
      v2[ici + id1 + iv] += mult_uv_i(&By[icg], &w2[id2 + iv]);
      v2[icr + id2 + iv] -= mult_uv_r(&By[icg], &w2[id1 + iv]);
      v2[ici + id2 + iv] -= mult_uv_i(&By[icg], &w2[id1 + iv]);

      v2[icr + id3 + iv] += mult_uv_r(&By[icg], &w2[id4 + iv]);
      v2[ici + id3 + iv] += mult_uv_i(&By[icg], &w2[id4 + iv]);
      v2[icr + id4 + iv] -= mult_uv_r(&By[icg], &w2[id3 + iv]);
      v2[ici + id4 + iv] -= mult_uv_i(&By[icg], &w2[id3 + iv]);

      // isigma_12 * Bz
      v2[icr + id1 + iv] -= mult_uv_i(&Bz[icg], &w2[id1 + iv]);
      v2[ici + id1 + iv] += mult_uv_r(&Bz[icg], &w2[id1 + iv]);
      v2[icr + id2 + iv] += mult_uv_i(&Bz[icg], &w2[id2 + iv]);
      v2[ici + id2 + iv] -= mult_uv_r(&Bz[icg], &w2[id2 + iv]);

      v2[icr + id3 + iv] -= mult_uv_i(&Bz[icg], &w2[id3 + iv]);
      v2[ici + id3 + iv] += mult_uv_r(&Bz[icg], &w2[id3 + iv]);
      v2[icr + id4 + iv] += mult_uv_i(&Bz[icg], &w2[id4 + iv]);
      v2[ici + id4 + iv] -= mult_uv_r(&Bz[icg], &w2[id4 + iv]);

      // isigma_41 * Ex
      v2[icr + id1 + iv] += mult_uv_i(&Ex[icg], &w2[id4 + iv]);
      v2[ici + id1 + iv] -= mult_uv_r(&Ex[icg], &w2[id4 + iv]);
      v2[icr + id2 + iv] += mult_uv_i(&Ex[icg], &w2[id3 + iv]);
      v2[ici + id2 + iv] -= mult_uv_r(&Ex[icg], &w2[id3 + iv]);

      v2[icr + id3 + iv] += mult_uv_i(&Ex[icg], &w2[id2 + iv]);
      v2[ici + id3 + iv] -= mult_uv_r(&Ex[icg], &w2[id2 + iv]);
      v2[icr + id4 + iv] += mult_uv_i(&Ex[icg], &w2[id1 + iv]);
      v2[ici + id4 + iv] -= mult_uv_r(&Ex[icg], &w2[id1 + iv]);

      // isigma_42 * Ey
      v2[icr + id1 + iv] -= mult_uv_r(&Ey[icg], &w2[id4 + iv]);
      v2[ici + id1 + iv] -= mult_uv_i(&Ey[icg], &w2[id4 + iv]);
      v2[icr + id2 + iv] += mult_uv_r(&Ey[icg], &w2[id3 + iv]);
      v2[ici + id2 + iv] += mult_uv_i(&Ey[icg], &w2[id3 + iv]);

      v2[icr + id3 + iv] -= mult_uv_r(&Ey[icg], &w2[id2 + iv]);
      v2[ici + id3 + iv] -= mult_uv_i(&Ey[icg], &w2[id2 + iv]);
      v2[icr + id4 + iv] += mult_uv_r(&Ey[icg], &w2[id1 + iv]);
      v2[ici + id4 + iv] += mult_uv_i(&Ey[icg], &w2[id1 + iv]);

      // isigma_43 * Ez
      v2[icr + id1 + iv] += mult_uv_i(&Ez[icg], &w2[id3 + iv]);
      v2[ici + id1 + iv] -= mult_uv_r(&Ez[icg], &w2[id3 + iv]);
      v2[icr + id2 + iv] -= mult_uv_i(&Ez[icg], &w2[id4 + iv]);
      v2[ici + id2 + iv] += mult_uv_r(&Ez[icg], &w2[id4 + iv]);

      v2[icr + id3 + iv] += mult_uv_i(&Ez[icg], &w2[id1 + iv]);
      v2[ici + id3 + iv] -= mult_uv_r(&Ez[icg], &w2[id1 + iv]);
      v2[icr + id4 + iv] -= mult_uv_i(&Ez[icg], &w2[id2 + iv]);
      v2[ici + id4 + iv] += mult_uv_r(&Ez[icg], &w2[id2 + iv]);
    }
  }

  scal(v, m_kappa * m_cSW); // v *= m_kappa * m_cSW;
}


//====================================================================
void Fopr_Clover_SF::set_csw()
{
  set_fieldstrength(m_Bx, 1, 2);
  set_fieldstrength(m_By, 2, 0);
  set_fieldstrength(m_Bz, 0, 1);
  set_fieldstrength(m_Ex, 3, 0);
  set_fieldstrength(m_Ey, 3, 1);
  set_fieldstrength(m_Ez, 3, 2);
}


/*!
  The field strength defined by clover with the SF BC.
  <ul>
  <li>The field strength is set to zero at t=0 boundary.
  <li>This is performed automatically by a use of Staple_SF.
  </ul>
*/
//====================================================================
void Fopr_Clover_SF::set_fieldstrength(Field_G& Fst,
                                       const int mu, const int nu)
{
  int Nvol = CommonParameters::Nvol();

  Staple_SF staple;

  staple.set_parameters(m_phi, m_phipr);

  Field_G_SF Cup(Nvol, 1), Cdn(Nvol, 1);
  Field_G_SF Umu(Nvol, 1);
  Field_G    v(Nvol, 1), v2(Nvol, 1);

  staple.upper(Cup, *m_U, mu, nu);
  staple.lower(Cdn, *m_U, mu, nu);
  Umu.setpart_ex(0, *m_U, mu);

  mult_Field_Gnd(Fst, 0, Umu, 0, Cup, 0);
  multadd_Field_Gnd(Fst, 0, Umu, 0, Cdn, 0, -1.0);

  mult_Field_Gdn(v, 0, Cup, 0, Umu, 0);
  multadd_Field_Gdn(v, 0, Cdn, 0, Umu, 0, -1.0);

  m_shift.forward(v2, v, mu);

  axpy(Fst, 1.0, v2); // Fst += v2;

  ah_Field_G(Fst, 0);
  scal(Fst, 0.25); // Fst *= 0.25;
}


//====================================================================
double Fopr_Clover_SF::flop_count()
{
  //- Counting of floating point operations.
  //  not implemented, yet.

  double flop = 0.0;

  return flop;
}


//====================================================================
//============================================================END=====
