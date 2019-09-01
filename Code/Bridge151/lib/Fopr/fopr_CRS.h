/*!
        @file    fopr_CRS.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/


#ifndef FOPR_CRS_INCLUDED
#define FOPR_CRS_INCLUDED

#include <string>

#include "fopr.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Fermion operator with CRS matrix format.

/*!
   This fermion operator class BAPI is defined by CRS matrix format
   which is widely used in studies of linear algorithms.
   The matrix is defined by giving fermion operator, or reading
   from a file (filename is provided at construction).
   In present implementation, this class BAPI works only on single
   node.  Parallel version should be implemented.
                                      [07 Dec 2011  H.Matsufuru]
*/

class BAPI Fopr_CRS : public Fopr
{
 public:
  static const std::string class_name;

 private:
  int m_Nin, m_Nvol, m_Nex;
  int m_Nsize, m_Nnz;
  std::vector<int> m_rowidx_nz;
  std::vector<int> m_column_nz;
  std::vector<double> m_elem_nz;
  std::string m_mode;
  Fopr *m_fopr;

 public:

  Fopr_CRS(Fopr *fopr)
    : Fopr(), m_fopr(fopr)
  {
    set_matrix();
  }

  Fopr_CRS(unique_ptr<Fopr>& fopr)
    : Fopr(), m_fopr(fopr.get())
  {
    set_matrix();
  }

  Fopr_CRS(const std::string fname)
    : Fopr(), m_fopr(0)
  {
    set_matrix(fname);
  }

  void set_parameters(const Parameters&);

  void write_matrix(const std::string);

  void set_config(Field *U)
  {
    if (m_fopr == 0) {
      vout.crucial(m_vl, "Error at %s: fopr is not set.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    } else {
      m_fopr->set_config(U);
    }
  }

  void set_config(unique_ptr<Field_G>& U)
  {
    if (m_fopr == 0) {
      vout.crucial(m_vl, "Error at %s: fopr is not set.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    } else {
      m_fopr->set_config(U.get());
    }
  }

  void set_mode(const std::string mode)
  {
    m_mode = mode;
  }

  std::string get_mode() const
  {
    return m_mode;
  }

  void mult(Field& v, const Field& f)
  {
    if (m_mode == "D") {
      D(v, f);
    } else if (m_mode == "DdagD") {
      DdagD(v, f);
    } else if (m_mode == "Ddag") {
      Ddag(v, f);
    } else {
      vout.crucial(m_vl, "Error at %s: mode unknown: '%s'.\n", class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }

  void mult_dag(Field& v, const Field& f)
  {
    if (m_mode == "D") {
      Ddag(v, f);
    } else if (m_mode == "DdagD") {
      DdagD(v, f);
    } else if (m_mode == "Ddag") {
      D(v, f);
    } else {
      vout.crucial(m_vl, "Error at %s: mode unknown: '%s'.\n", class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }

  void DdagD(Field&, const Field&);
  void D(Field&, const Field&);
  void Ddag(Field&, const Field&);
  void H(Field&, const Field&);

  int field_nvol() { return m_Nvol; }
  int field_nin() { return m_Nin; }
  int field_nex() { return m_Nex; }

 private:
  void set_matrix();
  void set_matrix(const std::string);

  void set_matrix_1row(int&, std::vector<int>&,
                       std::vector<double>&, const Field&);

#ifdef USE_FACTORY
 private:
  static Fopr *create_object_with_fopr(Fopr *fopr)
  {
    return new Fopr_CRS(fopr);
  }

  static Fopr *create_object_with_filename(const std::string& fname)
  {
    return new Fopr_CRS(fname);
  }

 public:
  static bool register_factory()
  {
    bool init1 = Fopr::Factory_fopr::Register("CRS", create_object_with_fopr);
    bool init2 = Fopr::Factory_string::Register("CRS", create_object_with_filename);

    return init1 && init2;
  }
#endif
};
#endif
