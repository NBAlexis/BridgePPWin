/*!
        @file    $Id:: index_eo.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2017-05-08 18:39:54 #$

        @version $LastChangedRevision: 1623 $
*/

#ifndef INDEX_EO_INCLUDED
#define INDEX_EO_INCLUDED

#include "field_G.h"
#include "Communicator/communicator.h"
#include "index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Even-odd site index.

/*!
    This class BAPI defines even-odd site index.
    Only the site degree of freedom is concerned.
    Nx (x-extent inside a node) must be even in the present
    implementation.
    Coverting from and reverting to the lexical site index are
    implemented as member functions of this class.
    In present implementation, there is no superclass structure,
    and thus polymorphism is not available.
    Some of method names might be confusing; restructuring may
    be helpful.
                                      [25 Dec 2011 H.Matsufuru]
*/
class BAPI Index_eo {
 private:
  int                Nx, Ny, Nz, Nt, Nvol;
  int                Nx2, Nvol2;
  std::valarray<int> Leo;
  std::valarray<int> Site_up;
  std::valarray<int> Site_dn;
  Index_lex          m_index_lex;
  int                m_node_eo;  // {0, 1} -- local origin is on even/odd side

  Bridge::VerboseLevel m_vl;
 public:
  Index_eo() :
    Nx(CommonParameters::Nx()),
    Ny(CommonParameters::Ny()),
    Nz(CommonParameters::Nz()),
    Nt(CommonParameters::Nt()),
    Nvol(CommonParameters::Nvol()),
    Nx2(CommonParameters::Nx() / 2),
    Nvol2(CommonParameters::Nvol() / 2),
    m_vl(CommonParameters::Vlevel())
  {
    if ((Nx % 2) == 1) {
      vout.crucial(m_vl, "Error at Index_eo: Nx is not even.\n");
      exit(EXIT_FAILURE);
    }

    // grid coordinate
    int gx = Communicator::ipe(0);
    int gy = Communicator::ipe(1);
    int gz = Communicator::ipe(2);
    int gt = Communicator::ipe(3);

    // node even/odd
    m_node_eo = (gx * Nx + gy * Ny + gz * Nz + gt * Nt) % 2;

    Leo.resize(Ny * Nz * Nt);
    Site_up.resize(Nx2 * Ny * Nz * Nt * 2);
    Site_dn.resize(Nx2 * Ny * Nz * Nt * 2);

    for (int t = 0; t < Nt; ++t) {
      for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
          int t2 = t + gt * Nt;
          int z2 = z + gz * Nz;
          int y2 = y + gy * Ny;
          // int t2 = Communicator::ipe(3) * Nt + t;
          // int z2 = Communicator::ipe(2) * Nz + z;
          // int y2 = Communicator::ipe(1) * Ny + y;
          Leo[y + Ny * (z + Nz * t)] = (y2 + z2 + t2) % 2;
        }
      }
    }

    for (int t = 0; t < Nt; ++t) {
      for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
          int yzt = y + Ny * (z + Nz * t);
          // int t2  = t;
          // int z2  = z;
          // int y2  = y;
          for (int x2 = 0; x2 < Nx2; ++x2) {
            int s = x2 + Nx2 * (y + Ny * (z + Nz * t));
            Site_up[s]         = ((x2 + Leo[yzt]) % Nx2) + Nx2 * (y + Ny * (z + Nz * t));
            Site_up[s + Nvol2] = ((x2 + 1 - Leo[yzt]) % Nx2) + Nx2 * (y + Ny * (z + Nz * t));
            Site_dn[s]         = ((x2 - 1 + Leo[yzt] + Nx2) % Nx2) + Nx2 * (y + Ny * (z + Nz * t));
            Site_dn[s + Nvol2] = ((x2 - Leo[yzt] + Nx2) % Nx2) + Nx2 * (y + Ny * (z + Nz * t));
          }
        }
      }
    }
  }

  int leo(const int y, const int z, const int t) const
  {
    return Leo[y + Ny * (z + Nz * t)];
  }

  int site(const int x2, const int y, const int z, const int t,
           const int ieo) const
  {
    return x2 + Nx2 * (y + Ny * (z + Nz * t)) + Nvol2 * ieo;
  }

  int site(const int is, const int ieo) const
  {
    return is + Nvol2 * ieo;
  }

  int site_up(const int x2, const int y, const int z, const int t,
              const int ieo) const
  {
    int s = x2 + Nx2 * (y + Ny * (z + Nz * t)) + Nvol2 * ieo;

    return Site_up[s] + Nvol2 * (1 - ieo);
  }

  int site_xup(const int x2, const int y, const int z, const int t,
               const int ieo) const
  {
    int s = x2 + Nx2 * (y + Ny * (z + Nz * t)) + Nvol2 * ieo;

    return Site_up[s] + Nvol2 * (1 - ieo);
  }

  int site_dn(const int x2, const int y, const int z, const int t,
              const int ieo) const
  {
    int s = x2 + Nx2 * (y + Ny * (z + Nz * t)) + Nvol2 * ieo;

    return Site_dn[s] + Nvol2 * (1 - ieo);
  }

  int site_xdn(const int x2, const int y, const int z, const int t,
               const int ieo) const
  {
    int s = x2 + Nx2 * (y + Ny * (z + Nz * t)) + Nvol2 * ieo;

    return Site_dn[s] + Nvol2 * (1 - ieo);
  }

  int siteh(const int x2, const int y, const int z, const int t)
  const
  {
    return x2 + Nx2 * (y + Ny * (z + Nz * t));
  }

  int siteh_up(const int x2, const int y, const int z, const int t,
               const int ieo) const
  {
    int s = x2 + Nx2 * (y + Ny * (z + Nz * t)) + Nvol2 * ieo;

    return Site_up[s];
  }

  int siteh_xup(const int x2, const int y, const int z, const int t,
                const int ieo) const
  {
    int s = x2 + Nx2 * (y + Ny * (z + Nz * t)) + Nvol2 * ieo;

    return Site_up[s];
  }

  int siteh_dn(const int x2, const int y, const int z, const int t,
               const int ieo) const
  {
    int s = x2 + Nx2 * (y + Ny * (z + Nz * t)) + Nvol2 * ieo;

    return Site_dn[s];
  }

  int siteh_xdn(const int x2, const int y, const int z, const int t,
                const int ieo) const
  {
    int s = x2 + Nx2 * (y + Ny * (z + Nz * t)) + Nvol2 * ieo;

    return Site_dn[s];
  }

  void convertField(Field& eo, const Field& lex);
  void convertField(Field& eo, const Field& lex, const int ieo);

  void reverseField(Field& lex, const Field& eo);
  void reverseField(Field& lex, const Field& eo, const int ieo);

  void splitField(Field& e, Field& o, const Field& eo);

  void mergeField(Field& eo, const Field& e, const Field& o);
};
#endif
