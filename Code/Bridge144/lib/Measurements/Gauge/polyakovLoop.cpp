#include "BridgeLib_Private.h"

/*!
        @file    $Id:: polyakovLoop.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2017-03-03 15:12:38 #$

        @version $LastChangedRevision: 1582 $
*/

#include "polyakovLoop.h"

const std::string PolyakovLoop::class_name = "PolyakovLoop";

//====================================================================
void PolyakovLoop::set_parameters(const Parameters& params)
{
  m_filename_output = params.get_string("filename_output");
  if (m_filename_output.empty()) {
    m_filename_output = "stdout";
  }

  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

#if 0
  //- fetch and check input parameters
  int Nspc_size, Ntype;

  int err = 0;
  err += params.fetch_int("spatial_correlator_size", Nspc_size);
  err += params.fetch_int("number_of_correlator_type", Ntype);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Nspc_size, Ntype);
#endif
}


//====================================================================
void PolyakovLoop::set_parameters(int Nspc_size, int Ntype)
{
#if 0
  //- print input parameters
  vout.general(m_vl, "Polyakov loop measurement:\n");
  vout.general(m_vl, "  Nspc_size = %d\n", Nspc_size);
  vout.general(m_vl, "  Ntype     = %d\n", Ntype);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Nspc_size);
  err += ParameterCheck::non_negative(Ntype);

  //! The following setting explicitly depends on the definition
  //! of unit vectors.
  if (Ntype > 6) ++err;

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Nspc_size = Nspc_size;
  m_Ntype     = Ntype;


  //- post-process
  m_Nx_ext   = CommonParameters::Nx() + m_Nspc_size + 1;
  m_Ny_ext   = CommonParameters::Ny() + m_Nspc_size + 1;
  m_Nz_ext   = CommonParameters::Nz() + m_Nspc_size + 1;
  m_Nt_ext   = 1;
  m_Nvol_ext = m_Nx_ext * m_Ny_ext * m_Nz_ext * m_Nt_ext;

  //! The following setting explicitly depends on the definition
  //! of unit vectors.
  m_Nmax[0] = Nspc_size;
  m_Nmax[1] = Nspc_size;
  m_Nmax[2] = Nspc_size / 2;
  m_Nmax[3] = Nspc_size;
  m_Nmax[4] = Nspc_size / 2;
  m_Nmax[5] = Nspc_size / 2;
#endif
}


//====================================================================
void PolyakovLoop::init()
{
  int Ndim = CommonParameters::Ndim();

  assert(Ndim == 4);

  m_filename_output = "stdout";

#if 0
  m_Ntype_max = 6;
  int Ndim_spc = Ndim - 1;

  m_Nunit.resize(m_Ntype_max);
  m_Nmax.resize(m_Ntype_max);

  for (int i = 0; i < m_Ntype_max; ++i) {
    m_Nunit[i].resize(Ndim_spc);
  }

  // The following setting explicitly depends on the definition
  // of unit vectors.
  assert(m_Ntype_max >= 6);

  m_Nunit[0][0] = 1;
  m_Nunit[0][1] = 0;
  m_Nunit[0][2] = 0;

  m_Nunit[1][0] = 1;
  m_Nunit[1][1] = 1;
  m_Nunit[1][2] = 0;

  m_Nunit[2][0] = 2;
  m_Nunit[2][1] = 1;
  m_Nunit[2][2] = 0;

  m_Nunit[3][0] = 1;
  m_Nunit[3][1] = 1;
  m_Nunit[3][2] = 1;

  m_Nunit[4][0] = 2;
  m_Nunit[4][1] = 1;
  m_Nunit[4][2] = 1;

  m_Nunit[5][0] = 2;
  m_Nunit[5][1] = 2;
  m_Nunit[5][2] = 1;
#endif
}


//====================================================================
dcomplex PolyakovLoop::measure_ploop(Field_G& U)
{
  int Ndim     = CommonParameters::Ndim();
  int Nx       = CommonParameters::Nx();
  int Ny       = CommonParameters::Ny();
  int Nz       = CommonParameters::Nz();
  int Nc       = CommonParameters::Nc();
  int Npet     = Communicator::npe(Ndim - 1);
  int Lvol_spc = CommonParameters::Lvol() / CommonParameters::Lt();
  int Nvol_spc = Nx * Ny * Nz;

  Field_G P(Nvol_spc, 1);

  calc_ploop(P, U);

  Mat_SU_N utmp(Nc);
  double   re_tr = 0.0;
  double   im_tr = 0.0;
  for (int site = 0; site < Nvol_spc; ++site) {
    P.mat(utmp, site, 0);

    re_tr += ReTr(utmp);
    im_tr += ImTr(utmp);
  }
  re_tr /= (double(Nc));
  im_tr /= (double(Nc));

  double re_ploop = Communicator::reduce_sum(re_tr);
  double im_ploop = Communicator::reduce_sum(im_tr);

  re_ploop /= double(Lvol_spc * Npet);
  im_ploop /= double(Lvol_spc * Npet);


  std::ostream& log_file_previous = vout.getStream();
  std::ofstream log_file;

  if (m_filename_output != "stdout") {
    log_file.open(m_filename_output.c_str(), std::ios::app);
    vout.init(log_file);
  }

  vout.general(m_vl, "PolyakovLoop = %20.16e %20.16e\n", re_ploop, im_ploop);

  if (m_filename_output != "stdout") {
    log_file.close();
    vout.init(log_file_previous);
  }


  dcomplex ploop = cmplx(re_ploop, im_ploop);

  return ploop;
}


//====================================================================
void PolyakovLoop::calc_ploop(Field_G& P, Field_G& U)
{
  int Ndim = CommonParameters::Ndim();
  // int Ndim_spc = Ndim - 1;

  int Nx = CommonParameters::Nx();
  int Ny = CommonParameters::Ny();
  int Nz = CommonParameters::Nz();
  int Nt = CommonParameters::Nt();
  int Nc = CommonParameters::Nc();

  int Nvol = U.nvol();

  assert(Nvol == Nx * Ny * Nz * Nt);

  int Nvol_spc = P.nvol();
  assert(Nvol_spc == Nx * Ny * Nz);

  Index_lex index;
  Index_lex index_spc(Nx, Ny, Nz, 1);

  Field_G  Ut(Nvol, 1);
  Mat_SU_N utmp1(Nc), utmp2(Nc), utmp3(Nc);

  Ut.setpart_ex(0, U, Ndim - 1);

  //- Definition of local Polyakov loops
  int t = 0;
  for (int z = 0; z < Nz; ++z) {
    for (int y = 0; y < Ny; ++y) {
      for (int x = 0; x < Nx; ++x) {
        int site  = index.site(x, y, z, t);
        int site2 = index_spc.site(x, y, z, 0);
        Ut.mat(utmp1, site, 0);
        P.set_mat(site2, 0, utmp1);
      }
    }
  }

  for (int t = 1; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          int site  = index.site(x, y, z, t);
          int site2 = index_spc.site(x, y, z, 0);
          Ut.mat(utmp1, site, 0);
          P.mat(utmp2, site2, 0);
          utmp3.mult_nn(utmp2, utmp1);
          P.set_mat(site2, 0, utmp3);
        }
      }
    }
  }

  //- global Polyakov loops
  int Npet = Communicator::npe(Ndim - 1);

  if (Npet > 1) {
    Field_G Pcp1(Nvol_spc, 1), Pcp2(Nvol_spc, 1);
    int     size_cp = P.nin() * Nvol_spc;

    for (int ipet = 1; ipet < Npet; ++ipet) {
      if (ipet == 1) {
        Pcp1.setpart_ex(0, P, 0);
      } else {
        Pcp1.setpart_ex(0, Pcp2, 0);
      }

      Communicator::exchange(size_cp, Pcp2.ptr(0), Pcp1.ptr(0), Ndim - 1, 1, 0);

      mult_Field_Gnn(Pcp1, 0, P, 0, Pcp2, 0);
      P.setpart_ex(0, Pcp1, 0);
    }
  }
}


//====================================================================
//============================================================END=====
