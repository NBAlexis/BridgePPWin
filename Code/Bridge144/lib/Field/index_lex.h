/*!
        @file    $Id:: index_lex.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/


#ifndef INDEX_LEX_INCLUDED
#define INDEX_LEX_INCLUDED

#include "Parameters/commonParameters.h"

//! Lexical site index.

/*!
  This class defines lexicographycal site index.
  Only the site degree of freedom is concerned.
  In present implementation, there is no superclass structure,
  and thus polymorphism is not available.
  Is it better to be renamed Index_lex and derived from
  generic Index class ?
                                       [25 Dec 2011 H.Matsufuru]

  Nx,Ny,Nz,Nt are enabled to be given at the construction.
                                       [26 May 2012 H.Matsufuru]
*/
class Index_lex {
 protected:
  int m_Nx, m_Ny, m_Nz, m_Nt;

 public:
  Index_lex() :
    m_Nx(CommonParameters::Nx()),
    m_Ny(CommonParameters::Ny()),
    m_Nz(CommonParameters::Nz()),
    m_Nt(CommonParameters::Nt()) { }

  Index_lex(int Nx, int Ny, int Nz, int Nt)
  {
    m_Nx = Nx;
    m_Ny = Ny;
    m_Nz = Nz;
    m_Nt = Nt;
  }

  int site(const int& x, const int& y, const int& z, const int& t)
  const
  {
    return m_Nx * (m_Ny * (m_Nz * t + z) + y) + x;
  }
};
#endif
