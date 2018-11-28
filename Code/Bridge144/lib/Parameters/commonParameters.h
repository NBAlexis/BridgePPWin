/*!
        @file    $Id:: commonParameters.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef COMMONPARAMETERS_INCLUDED
#define COMMONPARAMETERS_INCLUDED

#include <string>
using std::string;

#include <cmath>

//#include <valarray>
#include <vector>

#include "IO/bridgeIO.h"

//! Common parameter class: provides parameters as singleton.

/*!
    At present stage, several sets of parameters are explicitly
    given in the implementation file corresponding to some
    strings, for convenience of test.
                                      [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    epsilon_criterion is implemented. [14 Nov 2012 Y.Namekawa]
    General value of Nc become available.
                                     [07 May 2014 H.Matsufuru]
 */

class CommonParameters {
 public:
  static const std::string class_name;

 private:
  //- color, spinor and space-time dimension
  static int       m_Nc;
  static const int m_Nd;
  static const int m_Ndim;

  static const double m_epsilon_criterion;

  //- global lattice size
  static int m_Lx, m_Ly, m_Lz, m_Lt;
  static int m_Lvol;

  //- Number of processors assigined in each direction
  static int  m_NPEx, m_NPEy, m_NPEz, m_NPEt, m_NPE;
  static char m_map_grid[16];

  //- local lattice size
  static int m_Nx, m_Ny, m_Nz, m_Nt, m_Nvol;

  static bool m_initialized;

  CommonParameters() {}
  CommonParameters(const CommonParameters&) {}
  CommonParameters& operator=(const CommonParameters&);

  // static void setup(const string&);
  static bool check_parameters();

  static Bridge::VerboseLevel m_vlevel;

 public:
  //! initialization (Nc=3 is assumed).
  static void init(const std::vector<int>& lattice_size,
                   const std::vector<int>& grid_size);

  //! initialization with a given value of Nc.
  static void init(const std::vector<int>& lattice_size,
                   const std::vector<int>& grid_size,
                   const int Nc);

  //! initialization for default verbose level.
  static void init_Vlevel(Bridge::VerboseLevel vlevel);

  static void print_parameters();

  static int Lx() { return m_Lx; }
  static int Ly() { return m_Ly; }
  static int Lz() { return m_Lz; }
  static int Lt() { return m_Lt; }
  static int Lvol() { return m_Lvol; }

  static int NPEx() { return m_NPEx; }
  static int NPEy() { return m_NPEy; }
  static int NPEz() { return m_NPEz; }
  static int NPEt() { return m_NPEt; }
  static int NPE() { return m_NPE; }

  static char *Grid_map() { return m_map_grid; }

  static int Nx() { return m_Nx; }
  static int Ny() { return m_Ny; }
  static int Nz() { return m_Nz; }
  static int Nt() { return m_Nt; }
  static int Nvol() { return m_Nvol; }

  static int Lsize(const int dir);
  static int Nsize(const int dir);
  static int NPEsize(const int dir);

  static int Nc() { return m_Nc; }
  static int Nd() { return m_Nd; }
  static int Ndim() { return m_Ndim; }

  static double epsilon_criterion() { return m_epsilon_criterion; }
  static double epsilon_criterion2() { return pow(m_epsilon_criterion, 2); }

  static Bridge::VerboseLevel Vlevel() { return m_vlevel; }
};
#endif
