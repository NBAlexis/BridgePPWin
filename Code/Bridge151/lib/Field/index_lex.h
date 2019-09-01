/*!
        @file    index_lex.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/


#ifndef INDEX_LEX_INCLUDED
#define INDEX_LEX_INCLUDED

#include "Parameters/commonParameters.h"

//! Lexical site index.

/*!
  This class BAPI defines lexicographycal site index.
  Only the site degree of freedom is concerned.
  In present implementation, there is no superclass structure,
  and thus polymorphism is not available.
  Is it better to be renamed Index_lex and derived from
  generic Index class BAPI ?
                                       [25 Dec 2011 H.Matsufuru]

  Nx,Ny,Nz,Nt are enabled to be given at the construction.
                                       [26 May 2012 H.Matsufuru]
*/
class BAPI Index_lex {
 protected:
  int m_Nx, m_Ny, m_Nz, m_Nt;

 public:
  Index_lex() :
    m_Nx(CommonParameters::Nx()),
    m_Ny(CommonParameters::Ny()),
    m_Nz(CommonParameters::Nz()),
    m_Nt(CommonParameters::Nt()) { }

  Index_lex(const int Nx, const int Ny, const int Nz, const int Nt)
  {
    m_Nx = Nx;
    m_Ny = Ny;
    m_Nz = Nz;
    m_Nt = Nt;
  }

  int site(const int& x, const int& y, const int& z, const int& t)
  const
  {
    return x + m_Nx * (y + m_Ny * (z + m_Nz * t));
  }
};
#endif
