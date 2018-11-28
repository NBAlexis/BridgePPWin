#include "BridgeLib_Private.h"
#if USE_ORG

/*!
        @file    $Id:: field_G.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "Field/field_G.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//====================================================================
void Field_G::check()
{
  //  int Nthread = ThreadManager_OpenMP::get_num_threads();
  //  int i_thread = ThreadManager_OpenMP::get_thread_id();
  //  assert(Nthread == 1);
  //  assert(i_thread == 0);
  // Org-version only support single thread environment.
}


//====================================================================
void Field_G::set_unit()
{
  Mat_SU_N ut(m_Nc);

  ut.unit();

  for (int mu = 0; mu < m_Nex; ++mu) {
    for (int site = 0; site < m_Nvol; ++site) {
      set_mat(site, mu, ut);
    }
  }
}


//====================================================================
void Field_G::set_random(RandomNumbers *rand)
{
  rand->gauss_lex_global(*this);

  Mat_SU_N ut(m_Nc);

  for (int mu = 0; mu < m_Nex; ++mu) {
    for (int site = 0; site < m_Nvol; ++site) {
      this->mat(ut, site, mu);
      ut.reunit();
      this->set_mat(site, mu, ut);
    }
  }
}


//====================================================================
void Field_G::set_random(unique_ptr<RandomNumbers>& rand)
{
  set_random(rand.get());
}


//====================================================================
void Field_G::reunit()
{
  Mat_SU_N ut(m_Nc);

  for (int mu = 0; mu < m_Nex; ++mu) {
    for (int site = 0; site < m_Nvol; ++site) {
      mat(ut, site, mu);
      ut.reunit();
      set_mat(site, mu, ut);
    }
  }
}


//====================================================================
void mult_Field_Gnn(Field_G& W, const int ex,
                    const Field_G& U1, const int ex1,
                    const Field_G& U2, const int ex2)
{
  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  const int Nvol = W.nvol();

  for (int site = 0; site < Nvol; ++site) {
    W.set_mat(site, ex, U1.mat(site, ex1) * U2.mat(site, ex2));
  }
}


//====================================================================
void mult_Field_Gdn(Field_G& W, const int ex,
                    const Field_G& U1, const int ex1,
                    const Field_G& U2, const int ex2)
{
  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  const int Nvol = W.nvol();

  for (int site = 0; site < Nvol; ++site) {
    W.set_mat(site, ex, U1.mat_dag(site, ex1) * U2.mat(site, ex2));
  }
}


//====================================================================
void mult_Field_Gnd(Field_G& W, const int ex,
                    const Field_G& U1, const int ex1,
                    const Field_G& U2, const int ex2)
{
  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  const int Nvol = W.nvol();

  for (int site = 0; site < Nvol; ++site) {
    W.set_mat(site, ex, U1.mat(site, ex1) * U2.mat_dag(site, ex2));
  }
}


//====================================================================
void mult_Field_Gdd(Field_G& W, const int ex,
                    const Field_G& U1, const int ex1,
                    const Field_G& U2, const int ex2)
{
  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  const int Nvol = W.nvol();

  for (int site = 0; site < Nvol; ++site) {
    W.set_mat(site, ex, U1.mat_dag(site, ex1) * U2.mat_dag(site, ex2));
  }
}


//====================================================================
void multadd_Field_Gnn(Field_G& W, const int ex,
                       const Field_G& U1, const int ex1,
                       const Field_G& U2, const int ex2,
                       const double ff)
{
  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  const int Nvol = W.nvol();

  for (int site = 0; site < Nvol; ++site) {
    W.add_mat(site, ex, U1.mat(site, ex1) * U2.mat(site, ex2) * ff);
  }
}


//====================================================================
void multadd_Field_Gdn(Field_G& W, const int ex,
                       const Field_G& U1, const int ex1,
                       const Field_G& U2, const int ex2,
                       const double ff)
{
  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  const int Nvol = W.nvol();

  for (int site = 0; site < Nvol; ++site) {
    W.add_mat(site, ex, U1.mat_dag(site, ex1) * U2.mat(site, ex2) * ff);
  }
}


//====================================================================
void multadd_Field_Gnd(Field_G& W, const int ex,
                       const Field_G& U1, const int ex1,
                       const Field_G& U2, const int ex2,
                       const double ff)
{
  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  const int Nvol = W.nvol();

  for (int site = 0; site < Nvol; ++site) {
    W.add_mat(site, ex, U1.mat(site, ex1) * U2.mat_dag(site, ex2) * ff);
  }
}


//====================================================================
void multadd_Field_Gdd(Field_G& W, const int ex,
                       const Field_G& U1, const int ex1,
                       const Field_G& U2, const int ex2,
                       const double ff)
{
  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  const int Nvol = W.nvol();

  for (int site = 0; site < Nvol; ++site) {
    W.add_mat(site, ex, U1.mat_dag(site, ex1) * U2.mat_dag(site, ex2) * ff);
  }
}


//====================================================================
void at_Field_G(Field_G& W, const int ex)
{
  assert(ex < W.nex());

  const int Nvol = W.nvol();

  for (int site = 0; site < Nvol; ++site) {
    W.set_mat(site, ex, W.mat(site, ex).at());
  }
}


//====================================================================
void ah_Field_G(Field_G& W, const int ex)
{
  assert(ex < W.nex());

  const int Nvol = W.nvol();

  for (int site = 0; site < Nvol; ++site) {
    W.set_mat(site, ex, W.mat(site, ex).ah());
  }
}


//====================================================================
void mult_exp_Field_G(Field_G& W, const double alpha, const Field_G& iP, const Field_G& U, const int Nprec)
{
  // W = exp(alpha * iP) * U
  //   = (U + alpha * iP * (U + alpha/2 * iP * ( ... (U+ alpha/n * iP * U) ...
  int Nvol = U.nvol();
  int Nex  = U.nex();
  int Nc   = CommonParameters::Nc();

  Mat_SU_N u0(Nc), v0(Nc), w0(Nc);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = 0; site < Nvol; ++site) {
      u0 = U.mat(site, ex);
      v0 = iP.mat(site, ex);

      w0 = SU_N::mat_exp(alpha, v0, u0, Nprec);

      W.set_mat(site, ex, w0);
    }
  }
}


//====================================================================
//============================================================END=====
#endif
