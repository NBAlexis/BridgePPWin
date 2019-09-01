/*!
        @file    index_eo.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef INDEX_EO_INCLUDED
#define INDEX_EO_INCLUDED

#include "field.h"
#include "index_lex.h"

#include "Communicator/communicator.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Even-odd site index.

/*!
    This class BAPI defines even-odd site index.
    Only the site degree of freedom is concerned.
    Nx (x-extent inside a node) must be even in the present
    implementation.
    Coverting from and reverting to the lexical site index are
    implemented as member functions of this class BAPI.
    In present implementation, there is no superclass structure,
    and thus polymorphism is not available.
    Some of method names might be confusing; restructuring may
    be helpful.
                                      [25 Dec 2011 H.Matsufuru]
*/
class BAPI Index_eo {
 private:
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nvol;
  int m_Nx2, m_Nvol2;
  std::valarray<int> m_yzt_eo;
  std::valarray<int> m_Site_up;
  std::valarray<int> m_Site_dn;
  Index_lex m_index_lex;
  int m_node_eo;                 // {0, 1} -- local origin is on even/odd side

  Bridge::VerboseLevel m_vl;
 public:
  Index_eo() :
    m_Nx(CommonParameters::Nx()),
    m_Ny(CommonParameters::Ny()),
    m_Nz(CommonParameters::Nz()),
    m_Nt(CommonParameters::Nt()),
    m_Nvol(CommonParameters::Nvol()),
    m_Nx2(CommonParameters::Nx() / 2),
    m_Nvol2(CommonParameters::Nvol() / 2),
    m_vl(CommonParameters::Vlevel())
  {
    if ((m_Nx % 2) == 1) {
      vout.crucial(m_vl, "Error at Index_eo: Nx = %d, which must be even.\n", m_Nx);
      exit(EXIT_FAILURE);
    }

    //- grid coordinate
    const int ipe_x = Communicator::ipe(0);
    const int ipe_y = Communicator::ipe(1);
    const int ipe_z = Communicator::ipe(2);
    const int ipe_t = Communicator::ipe(3);

    //- node even/odd
    m_node_eo = (ipe_x * m_Nx + ipe_y * m_Ny + ipe_z * m_Nz + ipe_t * m_Nt) % 2;

    m_yzt_eo.resize(m_Ny * m_Nz * m_Nt);
    m_Site_up.resize(m_Nx2 * m_Ny * m_Nz * m_Nt * 2);
    m_Site_dn.resize(m_Nx2 * m_Ny * m_Nz * m_Nt * 2);

    for (int t = 0; t < m_Nt; ++t) {
      for (int z = 0; z < m_Nz; ++z) {
        for (int y = 0; y < m_Ny; ++y) {
          int t_global = t + ipe_t * m_Nt;
          int z_global = z + ipe_z * m_Nz;
          int y_global = y + ipe_y * m_Ny;

          m_yzt_eo[y + m_Ny * (z + m_Nz * t)] = (y_global + z_global + t_global) % 2;
        }
      }
    }

    for (int t = 0; t < m_Nt; ++t) {
      for (int z = 0; z < m_Nz; ++z) {
        for (int y = 0; y < m_Ny; ++y) {
          int yzt = y + m_Ny * (z + m_Nz * t);

          for (int x2 = 0; x2 < m_Nx2; ++x2) {
            int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t));
            m_Site_up[s]           = ((x2 + m_yzt_eo[yzt]) % m_Nx2) + m_Nx2 * (y + m_Ny * (z + m_Nz * t));
            m_Site_up[s + m_Nvol2] = ((x2 + 1 - m_yzt_eo[yzt]) % m_Nx2) + m_Nx2 * (y + m_Ny * (z + m_Nz * t));
            m_Site_dn[s]           = ((x2 - 1 + m_yzt_eo[yzt] + m_Nx2) % m_Nx2) + m_Nx2 * (y + m_Ny * (z + m_Nz * t));
            m_Site_dn[s + m_Nvol2] = ((x2 - m_yzt_eo[yzt] + m_Nx2) % m_Nx2) + m_Nx2 * (y + m_Ny * (z + m_Nz * t));
          }
        }
      }
    }
  }

  int leo(const int y, const int z, const int t) const
  {
    return m_yzt_eo[y + m_Ny * (z + m_Nz * t)];
  }

  int site(const int x2, const int y, const int z, const int t,
           const int ieo) const
  {
    return x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;
  }

  int site(const int is, const int ieo) const
  {
    return is + m_Nvol2 * ieo;
  }

  int site_up(const int x2, const int y, const int z, const int t,
              const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;

    return m_Site_up[s] + m_Nvol2 * (1 - ieo);
  }

  int site_xup(const int x2, const int y, const int z, const int t,
               const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;

    return m_Site_up[s] + m_Nvol2 * (1 - ieo);
  }

  int site_dn(const int x2, const int y, const int z, const int t,
              const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;

    return m_Site_dn[s] + m_Nvol2 * (1 - ieo);
  }

  int site_xdn(const int x2, const int y, const int z, const int t,
               const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;

    return m_Site_dn[s] + m_Nvol2 * (1 - ieo);
  }

  int siteh(const int x2, const int y, const int z, const int t)
  const
  {
    return x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t));
  }

  int siteh_up(const int x2, const int y, const int z, const int t,
               const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;

    return m_Site_up[s];
  }

  int siteh_xup(const int x2, const int y, const int z, const int t,
                const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;

    return m_Site_up[s];
  }

  int siteh_dn(const int x2, const int y, const int z, const int t,
               const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;

    return m_Site_dn[s];
  }

  int siteh_xdn(const int x2, const int y, const int z, const int t,
                const int ieo) const
  {
    int s = x2 + m_Nx2 * (y + m_Ny * (z + m_Nz * t)) + m_Nvol2 * ieo;

    return m_Site_dn[s];
  }

  void convertField(Field& eo, const Field& lex);
  void convertField(Field& eo, const Field& lex, const int ieo);

  void reverseField(Field& lex, const Field& eo);
  void reverseField(Field& lex, const Field& eo, const int ieo);

  void splitField(Field& e, Field& o, const Field& eo);

  void mergeField(Field& eo, const Field& e, const Field& o);
};
#endif
