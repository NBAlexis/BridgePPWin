#include "BridgeLib_Private.h"

/*!
        @file    $Id:: index_eo.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2017-05-08 18:39:54 #$

        @version $LastChangedRevision: 1623 $
*/


#include "index_eo.h"

#include <assert.h>

//====================================================================
void Index_eo::convertField(Field& field_eo, const Field& field_lex)
{
  assert(field_eo.nin() == field_lex.nin());
  assert(field_eo.nex() == field_lex.nex());
  assert(Nvol == field_lex.nvol());
  assert(Nvol == field_eo.nvol());

  int Nin = field_lex.nin();
  int Nex = field_lex.nex();

  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          int x2  = x / 2;
          int ieo = (x + y + z + t + m_node_eo) % 2;
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              field_eo.set(in, site(x2, y, z, t, ieo), ex,
                           field_lex.cmp(in, m_index_lex.site(x, y, z, t), ex));
            }
          }
        }
      }
    }
  }
}


//====================================================================
void Index_eo::convertField(Field& field_eo, const Field& field_lex,
                            const int ieo)
{
  assert(field_eo.nin() == field_lex.nin());
  assert(field_eo.nex() == field_lex.nex());
  assert(Nvol == field_lex.nvol());
  assert(Nvol2 == field_eo.nvol());

  int Nin = field_lex.nin();
  int Nex = field_lex.nex();

  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        for (int x2 = 0; x2 < Nx2; ++x2) {
          int x = 2 * x2 + ((y + z + t + m_node_eo + ieo) % 2);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              field_eo.set(in, siteh(x2, y, z, t), ex,
                           field_lex.cmp(in, m_index_lex.site(x, y, z, t), ex));
            }
          }
        }
      }
    }
  }
}


//====================================================================
void Index_eo::reverseField(Field& field_lex, const Field& field_eo,
                            const int ieo)
{
  assert(field_eo.nin() == field_lex.nin());
  assert(field_eo.nex() == field_lex.nex());
  assert(Nvol == field_lex.nvol());
  assert(Nvol2 == field_eo.nvol());

  int Nin = field_lex.nin();
  int Nex = field_lex.nex();

  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        for (int x2 = 0; x2 < Nx2; ++x2) {
          int x = 2 * x2 + ((y + z + t + m_node_eo + ieo) % 2);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              field_lex.set(in, m_index_lex.site(x, y, z, t), ex,
                            field_eo.cmp(in, siteh(x2, y, z, t), ex));
            }
          }
        }
      }
    }
  }
}


//====================================================================
void Index_eo::reverseField(Field& field_lex, const Field& field_eo)
{
  assert(field_eo.nin() == field_lex.nin());
  assert(field_eo.nex() == field_lex.nex());
  assert(Nvol == field_lex.nvol());
  assert(Nvol == field_eo.nvol());

  int Nin = field_lex.nin();
  int Nex = field_lex.nex();

  for (int t = 0; t < Nt; ++t) {
    for (int z = 0; z < Nz; ++z) {
      for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
          int x2  = x / 2;
          int ieo = (x + y + z + t + m_node_eo) % 2;
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              field_lex.set(in, m_index_lex.site(x, y, z, t), ex,
                            field_eo.cmp(in, site(x2, y, z, t, ieo), ex));
            }
          }
        }
      }
    }
  }
}


//====================================================================
void Index_eo::splitField(Field& field_e, Field& field_o,
                          const Field& field_eo)
{
  int Nin = field_eo.nin();
  int Nex = field_eo.nex();

  assert(field_e.nin() == Nin);
  assert(field_e.nex() == Nex);

  assert(field_o.nin() == Nin);
  assert(field_o.nex() == Nex);

  int Nvol2 = field_eo.nvol() / 2;

  for (int iex = 0; iex < Nex; ++iex) {
    for (int ivol = 0; ivol < Nvol2; ++ivol) {
      for (int iin = 0; iin < Nin; ++iin) {
        field_e.set(iin, ivol, iex,
                    field_eo.cmp(iin, ivol, iex));
        field_o.set(iin, ivol, iex,
                    field_eo.cmp(iin, ivol + Nvol2, iex));
      }
    }
  }
}


//====================================================================
void Index_eo::mergeField(Field& field_eo,
                          const Field& field_e, const Field& field_o)
{
  int Nin = field_eo.nin();
  int Nex = field_eo.nex();

  assert(field_e.nin() == Nin);
  assert(field_e.nex() == Nex);

  assert(field_o.nin() == Nin);
  assert(field_o.nex() == Nex);

  int Nvol2 = field_eo.nvol() / 2;

  for (int iex = 0; iex < Nex; ++iex) {
    for (int ivol = 0; ivol < Nvol2; ++ivol) {
      for (int iin = 0; iin < Nin; ++iin) {
        field_eo.set(iin, ivol, iex,
                     field_e.cmp(iin, ivol, iex));
        field_eo.set(iin, ivol + Nvol2, iex,
                     field_o.cmp(iin, ivol, iex));
      }
    }
  }
}


//====================================================================
//============================================================END=====
