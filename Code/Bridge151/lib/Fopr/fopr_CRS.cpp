/*!
        @file    fopr_CRS.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "fopr_CRS.h"

#include <fstream>

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_CRS::register_factory();
}
#endif

const std::string Fopr_CRS::class_name = "Fopr_CRS";

//====================================================================
void Fopr_CRS::set_parameters(const Parameters& params)
{
  const std::string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);
}


//====================================================================
void Fopr_CRS::DdagD(Field& v, const Field& w)
{
  Field w2(w.nin(), w.nvol(), w.nex());

  D(w2, w);
  Ddag(v, w2);
}


//====================================================================
void Fopr_CRS::D(Field& v, const Field& w)
{
  const int Nsize = w.size() / 2;

  assert(Nsize == m_Nsize);

  int Ncolumn = 0;

  // v = 0.0;
  v.set(0.0);

  for (int j = 0; j < Nsize; ++j) {
    if (j == Nsize - 1) {
      Ncolumn = m_Nnz - m_rowidx_nz[j];
    } else {
      Ncolumn = m_rowidx_nz[j + 1] - m_rowidx_nz[j];
    }

    double vr = 0.0;
    double vi = 0.0;
    for (int i = 0; i < Ncolumn; ++i) {
      int    i2 = i + m_rowidx_nz[j];
      int    k  = m_column_nz[i2];
      double wr = w.cmp(2 * k);
      double wi = w.cmp(2 * k + 1);

      vr += m_elem_nz[2 * i2] * wr - m_elem_nz[2 * i2 + 1] * wi;
      vi += m_elem_nz[2 * i2] * wi + m_elem_nz[2 * i2 + 1] * wr;

      v.set(2 * j, vr);
      v.set(2 * j + 1, vi);
    }
  }
}


//====================================================================
void Fopr_CRS::Ddag(Field& v, const Field& w)
{
  const int Nsize = w.size() / 2;

  assert(Nsize == m_Nsize);

  int Ncolumn = 0;

  // v = 0.0;
  v.set(0.0);

  for (int j = 0; j < Nsize; ++j) {
    if (j == Nsize - 1) {
      Ncolumn = m_Nnz - m_rowidx_nz[j];
    } else {
      Ncolumn = m_rowidx_nz[j + 1] - m_rowidx_nz[j];
    }

    double wr = w.cmp(2 * j);
    double wi = w.cmp(2 * j + 1);

    for (int i = 0; i < Ncolumn; ++i) {
      int    i2 = i + m_rowidx_nz[j];
      int    k  = m_column_nz[i2];
      double vr = m_elem_nz[2 * i2] * wr + m_elem_nz[2 * i2 + 1] * wi;
      double vi = m_elem_nz[2 * i2] * wi - m_elem_nz[2 * i2 + 1] * wr;

      v.add(2 * k, vr);
      v.add(2 * k + 1, vi);
    }
  }
}


//====================================================================
void Fopr_CRS::set_matrix()
{
  // the implementation of this function assumes that the field
  // is complex valued.

  assert(CommonParameters::NPE() == 1);

  m_Nin  = m_fopr->field_nin();
  m_Nvol = m_fopr->field_nvol();
  m_Nex  = m_fopr->field_nex();

  m_Nsize = (m_Nin / 2) * m_Nvol * m_Nex;

  Field w(m_Nin, m_Nvol, m_Nex), v(m_Nin, m_Nvol, m_Nex);

  int                 Nnz1;
  std::vector<int>    index_nz1(m_Nsize);
  std::vector<double> elem_nz1(2 * m_Nsize);

  vout.general(m_vl, "Setting SRC matrix format.\n");
  vout.general(m_vl, "  Nex = %4d\n", m_Nex);

  //- estimate of matrix size
  int Nnz_tot = 0;

  {
    int in = 0;
    {
      int site = 0;
      for (int ex = 0; ex < m_Nex; ++ex) {
        w.set(0.0);
        w.set(in, site, ex, 1.0);
        m_fopr->mult_dag(v, w);
        set_matrix_1row(Nnz1, index_nz1, elem_nz1, v);
        Nnz_tot += Nnz1;
        vout.general(m_vl, "  ex = %4d  Nnz1 = %d\n", ex, Nnz1);
      }
    }
  }
  Nnz_tot *= (m_Nin / 2) * m_Nvol;
  vout.general(m_vl, "  Nnz_tot     = %d\n", Nnz_tot);

  //- resize data array
  m_Nnz = Nnz_tot;
  m_rowidx_nz.resize(m_Nsize);
  m_column_nz.resize(m_Nnz);
  m_elem_nz.resize(2 * m_Nnz);

  //- setting CRS matrix
  Nnz_tot = 0;
  for (int ex = 0; ex < m_Nex; ++ex) {
    vout.general(m_vl, "  ex = %4d starts.\n", ex);
    for (int site = 0; site < m_Nvol; ++site) {
      for (int in2 = 0; in2 < (m_Nin / 2); ++in2) {
        int j = in2 + (m_Nin / 2) * (site + m_Nvol * ex);

        w.set(0.0);
        w.set(2 * in2, site, ex, 1.0);
        m_fopr->mult_dag(v, w);
        set_matrix_1row(Nnz1, index_nz1, elem_nz1, v);
        if (Nnz_tot + Nnz1 > m_Nnz) {
          vout.crucial(m_vl, "Error at %s: unexpected data size ini set_matrix\n", class_name.c_str());
          exit(EXIT_FAILURE);
        }

        m_rowidx_nz[j] = Nnz_tot;
        for (int i = 0; i < Nnz1; ++i) {
          m_column_nz[Nnz_tot + i]         = index_nz1[i];
          m_elem_nz[2 * (Nnz_tot + i)]     = elem_nz1[2 * i];
          m_elem_nz[2 * (Nnz_tot + i) + 1] = elem_nz1[2 * i + 1];
        }
        Nnz_tot += Nnz1;
      }
    }
  }

  m_Nnz = Nnz_tot;

  vout.general(m_vl, "  Nnz(recalc) = %d\n", m_Nnz);
}


//====================================================================
void Fopr_CRS::set_matrix_1row(int& Nnz, std::vector<int>& index_nz,
                               std::vector<double>& elem_nz, const Field& v)
{
  // the implementation of this function assumes that the field
  // is complex valued.

  const int Nsize = v.size() / 2;

  int j = 0;

  for (int i = 0; i < Nsize; ++i) {
    if ((v.cmp(2 * i) != 0.0) || (v.cmp(2 * i + 1) != 0.0)) {
      index_nz[j]        = i;
      elem_nz[2 * j]     = v.cmp(2 * i);
      elem_nz[2 * j + 1] = -v.cmp(2 * i + 1);

      ++j;
    }
  }

  Nnz = j;
}


//====================================================================
void Fopr_CRS::set_matrix(const std::string fname)
{
  std::fstream config;

  config.open(fname.c_str(), std::ios::in);
  if (!config.is_open()) {
    vout.crucial(m_vl, "Error at %s: Failed to open the next file %s in %s(%d)",
                 class_name.c_str(), fname.c_str(), __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "Reading CRS matrix from %s\n", fname.c_str());

  config >> m_Nin;
  config >> m_Nvol;
  config >> m_Nex;
  config >> m_Nsize;
  config >> m_Nnz;

  m_rowidx_nz.resize(m_Nsize);
  m_column_nz.resize(m_Nnz);
  m_elem_nz.resize(2 * m_Nnz);

  for (int j = 0; j < m_Nsize; ++j) {
    config >> m_rowidx_nz[j];
    m_rowidx_nz[j] -= 1;
  }

  for (int j = 0; j < m_Nnz; ++j) {
    config >> m_column_nz[j];
    m_column_nz[j] -= 1;
  }

  for (int j = 0; j < m_Nnz; ++j) {
    config >> m_elem_nz[2 * j];
    config >> m_elem_nz[2 * j + 1];
  }

  config.close();

  vout.general(m_vl, "Reading CRS matrix finished.\n");
}


//====================================================================
void Fopr_CRS::write_matrix(const std::string fname)
{
  std::fstream config;

  config.open(fname.c_str(), std::ios::out);
  if (!config.is_open()) {
    vout.crucial(m_vl, "Error at %s: Failed to open the text file %s(%d)\n",
                 class_name.c_str(), __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "Writing CRS matrix to %s\n", fname.c_str());

  config << m_Nin << std::endl;
  config << m_Nvol << std::endl;
  config << m_Nex << std::endl;
  config << m_Nsize << std::endl;
  config << m_Nnz << std::endl;

  for (int j = 0; j < m_Nsize; ++j) {
    config << m_rowidx_nz[j] + 1 << std::endl;
  }

  for (int j = 0; j < m_Nnz; ++j) {
    config << m_column_nz[j] + 1 << std::endl;
  }

  config.setf(std::ios_base::scientific, std::ios_base::floatfield);
  config.precision(14);

  for (int j = 0; j < m_Nnz; ++j) {
    config << m_elem_nz[2 * j] << std::endl;
    config << m_elem_nz[2 * j + 1] << std::endl;
  }

  config.close();

  vout.general(m_vl, "Writing CRS matrix finished.\n");
}


//====================================================================
//============================================================END=====
