/*!
        @file    shiftField_lex_imp.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"

#if USE_IMP

#include "Field/shiftField_lex.h"
#include <vector>

const std::string ShiftField_lex::class_name = "ShiftField_lex";

//====================================================================
void ShiftField_lex::backward(Field& v,
                              const Field& w, const int mu)
{
  const int boundary_condition = 1;

  if (mu == 0) {  // x-direction
    up_x(&v, &w, boundary_condition);
  } else if (mu == 1) {
    up_y(&v, &w, boundary_condition);
  } else if (mu == 2) {
    up_z(&v, &w, boundary_condition);
  } else if (mu == 3) {
    up_t(&v, &w, boundary_condition);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_lex::forward(Field& v,
                             const Field& w, const int mu)
{
  const int boundary_condition = 1;

  if (mu == 0) {
    dn_x(&v, &w, boundary_condition);
  } else if (mu == 1) {
    dn_y(&v, &w, boundary_condition);
  } else if (mu == 2) {
    dn_z(&v, &w, boundary_condition);
  } else if (mu == 3) {
    dn_t(&v, &w, boundary_condition);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_lex::backward(Field& v,
                              const Field& w, const int boundary_condition, const int mu)
{
  if (mu == 0) {  // x-direction
    up_x(&v, &w, boundary_condition);
  } else if (mu == 1) {
    up_y(&v, &w, boundary_condition);
  } else if (mu == 2) {
    up_z(&v, &w, boundary_condition);
  } else if (mu == 3) {
    up_t(&v, &w, boundary_condition);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_lex::forward(Field& v,
                             const Field& w, const int boundary_condition, const int mu)
{
  if (mu == 0) {
    dn_x(&v, &w, boundary_condition);
  } else if (mu == 1) {
    dn_y(&v, &w, boundary_condition);
  } else if (mu == 2) {
    dn_z(&v, &w, boundary_condition);
  } else if (mu == 3) {
    dn_t(&v, &w, boundary_condition);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_lex::up_x(Field *v,
                          const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(0) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin   = w->nin();
  const int Nex   = w->nex();
  const int Nvol  = w->nvol();
  const int Nvol2 = Nvol / m_Nx;

  double       *vp = v->ptr(0);
  const double *wp = w->ptr(0);

  std::vector<double> wt(Nin * Nvol2 * Nex);
  std::vector<double> vt(Nin * Nvol2 * Nex);

  for (int s2 = 0; s2 < m_Ny * m_Nz * m_Nt; s2++) {
    // bulk
    for (int x = 0; x < m_Nx - 1; x++) {
      int ix = x + m_Nx * s2;
      int px = ix + 1;
      for (int ex = 0; ex < Nex; ex++) {
        for (int in = 0; in < Nin; in++) {
          vp[in + Nin * (ix + Nvol * ex)] = wp[in + Nin * (px + Nvol * ex)];
        }
      }
    }

    // boundary (x=m_Nx-1)
    int px = m_Nx * s2;
    // int ix = m_Nx - 1 + px;
    for (int ex = 0; ex < Nex; ex++) {
      for (int in = 0; in < Nin; in++) {
        wt[in + Nin * (s2 + Nvol2 * ex)] = bc2 * wp[in + Nin * (px + Nvol * ex)];
      }
    }
  }

  const int size = Nin * Nvol2 * Nex;
  Communicator::exchange(size, &vt[0], &wt[0], 0, 1, 0);

  for (int s2 = 0; s2 < m_Ny * m_Nz * m_Nt; s2++) {
    // boundary (x=m_Nx-1)
    int px = m_Nx * s2;
    int ix = m_Nx - 1 + px;
    for (int in = 0; in < Nin; in++) {
      for (int ex = 0; ex < Nex; ex++) {
        vp[in + Nin * (ix + Nvol * ex)] = vt[in + Nin * (s2 + Nvol2 * ex)];
      }
    }
  }
}


//====================================================================
void ShiftField_lex::dn_x(Field *v,
                          const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(0) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin   = w->nin();
  const int Nex   = w->nex();
  const int Nvol  = w->nvol();
  const int Nvol2 = Nvol / m_Nx;

  double       *vp = v->ptr(0);
  const double *wp = w->ptr(0);

  std::vector<double> wt(Nin * Nvol2 * Nex);
  std::vector<double> vt(Nin * Nvol2 * Nex);

  for (int s2 = 0; s2 < m_Ny * m_Nz * m_Nt; s2++) {
    // bulk
    for (int x = 1; x < m_Nx; x++) {
      int ix = x + m_Nx * s2;
      int px = ix - 1;
      for (int ex = 0; ex < Nex; ex++) {
        for (int in = 0; in < Nin; in++) {
          vp[in + Nin * (ix + Nvol * ex)] = wp[in + Nin * (px + Nvol * ex)];
        }
      }
    }

    // boundary (x=0)
    int px = m_Nx - 1 + m_Nx * s2;
    for (int ex = 0; ex < Nex; ex++) {
      for (int in = 0; in < Nin; in++) {
        wt[in + Nin * (s2 + Nvol2 * ex)] = wp[in + Nin * (px + Nvol * ex)];
      }
    }
  }

  const int size = Nin * Nvol2 * Nex;
  Communicator::exchange(size, &vt[0], &wt[0], 0, -1, 4);

  for (int s2 = 0; s2 < m_Ny * m_Nz * m_Nt; s2++) {
    // boundary (x=0)
    int ix = m_Nx * s2;
    for (int ex = 0; ex < Nex; ex++) {
      for (int in = 0; in < Nin; in++) {
        vp[in + Nin * (ix + Nvol * ex)] = bc2 * vt[in + Nin * (s2 + Nvol2 * ex)];
      }
    }
  }
}


//====================================================================
void ShiftField_lex::up_y(Field *v,
                          const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(1) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin   = w->nin();
  const int Nex   = w->nex();
  const int Nvol  = w->nvol();
  const int Nvol2 = Nvol / m_Ny;

  double       *vp = v->ptr(0);
  const double *wp = w->ptr(0);

  std::vector<double> wt(Nin * Nvol2 * Nex);
  std::vector<double> vt(Nin * Nvol2 * Nex);

  for (int zt = 0; zt < m_Nz * m_Nt; zt++) {
    for (int x = 0; x < m_Nx; x++) {
      int s1 = x + m_Nx * m_Ny * zt;
      int s2 = x + m_Nx * zt;

      // bulk
      for (int y = 0; y < m_Ny - 1; y++) {
        int ix = s1 + m_Nx * y;
        int px = ix + m_Nx;
        for (int ex = 0; ex < Nex; ex++) {
          for (int in = 0; in < Nin; in++) {
            vp[in + Nin * (ix + Nvol * ex)] = wp[in + Nin * (px + Nvol * ex)];
          }
        }
      }

      // boundary (y=m_Ny-1)
      int px = s1;
      for (int ex = 0; ex < Nex; ex++) {
        for (int in = 0; in < Nin; in++) {
          wt[in + Nin * (s2 + Nvol2 * ex)] = bc2 * wp[in + Nin * (px + Nvol * ex)];
        }
      }
    }
  }

  const int size = Nin * Nvol2 * Nex;
  Communicator::exchange(size, &vt[0], &wt[0], 1, 1, 1);

  for (int zt = 0; zt < m_Nz * m_Nt; zt++) {
    for (int x = 0; x < m_Nx; x++) {
      int s1 = x + m_Nx * m_Ny * zt;
      int s2 = x + m_Nx * zt;

      // boundary (y=m_Ny-1)
      int ix = s1 + m_Nx * (m_Ny - 1);
      for (int ex = 0; ex < Nex; ex++) {
        for (int in = 0; in < Nin; in++) {
          vp[in + Nin * (ix + Nvol * ex)] = vt[in + Nin * (s2 + Nvol2 * ex)];
        }
      }
    }
  }
}


//====================================================================
void ShiftField_lex::dn_y(Field *v,
                          const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(1) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin   = w->nin();
  const int Nex   = w->nex();
  const int Nvol  = w->nvol();
  const int Nvol2 = Nvol / m_Ny;

  double       *vp = v->ptr(0);
  const double *wp = w->ptr(0);

  std::vector<double> wt(Nin * Nvol2 * Nex);
  std::vector<double> vt(Nin * Nvol2 * Nex);

  for (int zt = 0; zt < m_Nz * m_Nt; zt++) {
    for (int x = 0; x < m_Nx; x++) {
      int s1 = x + m_Nx * m_Ny * zt;
      int s2 = x + m_Nx * zt;

      // bulk
      for (int y = 1; y < m_Ny; y++) {
        int ix = s1 + m_Nx * y;
        int px = ix - m_Nx;
        for (int ex = 0; ex < Nex; ex++) {
          for (int in = 0; in < Nin; in++) {
            vp[in + Nin * (ix + Nvol * ex)] = wp[in + Nin * (px + Nvol * ex)];
          }
        }
      }

      // boundary (y=0)
      int px = s1 + m_Nx * (m_Ny - 1);
      for (int ex = 0; ex < Nex; ex++) {
        for (int in = 0; in < Nin; in++) {
          wt[in + Nin * (s2 + Nvol2 * ex)] = wp[in + Nin * (px + Nvol * ex)];
        }
      }
    }
  }

  const int size = Nin * Nvol2 * Nex;
  Communicator::exchange(size, &vt[0], &wt[0], 1, -1, 5);

  for (int zt = 0; zt < m_Nz * m_Nt; zt++) {
    for (int x = 0; x < m_Nx; x++) {
      int s1 = x + m_Nx * m_Ny * zt;
      int s2 = x + m_Nx * zt;

      // boundary (y=0)
      int ix = s1;
      for (int ex = 0; ex < Nex; ex++) {
        for (int in = 0; in < Nin; in++) {
          vp[in + Nin * (ix + Nvol * ex)] = bc2 * vt[in + Nin * (s2 + Nvol2 * ex)];
        }
      }
    }
  }
}


//====================================================================
void ShiftField_lex::up_z(Field *v,
                          const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(2) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin   = w->nin();
  const int Nex   = w->nex();
  const int Nvol  = w->nvol();
  const int Nvol2 = Nvol / m_Nz;

  double       *vp = v->ptr(0);
  const double *wp = w->ptr(0);

  std::vector<double> wt(Nin * Nvol2 * Nex);
  std::vector<double> vt(Nin * Nvol2 * Nex);

  for (int t = 0; t < m_Nt; t++) {
    for (int xy = 0; xy < m_Nx * m_Ny; xy++) {
      int s1 = xy + m_Nx * m_Ny * m_Nz * t;
      int s2 = xy + m_Nx * m_Ny * t;

      // bulk
      for (int z = 0; z < m_Nz - 1; z++) {
        int ix = s1 + m_Nx * m_Ny * z;
        int px = ix + m_Nx * m_Ny;
        for (int ex = 0; ex < Nex; ex++) {
          for (int in = 0; in < Nin; in++) {
            vp[in + Nin * (ix + Nvol * ex)] = wp[in + Nin * (px + Nvol * ex)];
          }
        }
      }

      // boundary (z=m_Nz-1)
      int px = s1;
      for (int ex = 0; ex < Nex; ex++) {
        for (int in = 0; in < Nin; in++) {
          wt[in + Nin * (s2 + Nvol2 * ex)] = bc2 * wp[in + Nin * (px + Nvol * ex)];
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Nz) * Nex;
  Communicator::exchange(size, &vt[0], &wt[0], 2, 1, 2);

  for (int t = 0; t < m_Nt; t++) {
    for (int xy = 0; xy < m_Nx * m_Ny; xy++) {
      int s1 = xy + m_Nx * m_Ny * m_Nz * t;
      int s2 = xy + m_Nx * m_Ny * t;

      // boundary (z=m_Nz-1)
      int ix = s1 + m_Nx * m_Ny * (m_Nz - 1);
      for (int ex = 0; ex < Nex; ex++) {
        for (int in = 0; in < Nin; in++) {
          vp[in + Nin * (ix + Nvol * ex)] = vt[in + Nin * (s2 + Nvol2 * ex)];
        }
      }
    }
  }
}


//====================================================================
void ShiftField_lex::dn_z(Field *v,
                          const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(2) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin   = w->nin();
  const int Nex   = w->nex();
  const int Nvol  = w->nvol();
  const int Nvol2 = Nvol / m_Nz;

  double       *vp = v->ptr(0);
  const double *wp = w->ptr(0);

  std::vector<double> wt(Nin * Nvol2 * Nex);
  std::vector<double> vt(Nin * Nvol2 * Nex);

  for (int t = 0; t < m_Nt; t++) {
    for (int xy = 0; xy < m_Nx * m_Ny; xy++) {
      int s1 = xy + m_Nx * m_Ny * m_Nz * t;
      int s2 = xy + m_Nx * m_Ny * t;

      // bulk
      for (int z = 1; z < m_Nz; z++) {
        int ix = s1 + m_Nx * m_Ny * z;
        int px = ix - m_Nx * m_Ny;
        for (int ex = 0; ex < Nex; ex++) {
          for (int in = 0; in < Nin; in++) {
            vp[in + Nin * (ix + Nvol * ex)] = wp[in + Nin * (px + Nvol * ex)];
          }
        }
      }
      // boundary (z=0)
      int px = s1 + m_Nx * m_Ny * (m_Nz - 1);
      for (int in = 0; in < Nin; in++) {
        for (int ex = 0; ex < Nex; ex++) {
          wt[in + Nin * (s2 + Nvol2 * ex)] = wp[in + Nin * (px + Nvol * ex)];
        }
      }
    }
  }

  const int size = Nin * Nvol2 * Nex;
  Communicator::exchange(size, &vt[0], &wt[0], 2, -1, 6);

  for (int t = 0; t < m_Nt; t++) {
    for (int xy = 0; xy < m_Nx * m_Ny; xy++) {
      int s1 = xy + m_Nx * m_Ny * m_Nz * t;
      int s2 = xy + m_Nx * m_Ny * t;

      // boundary (z=0)
      int ix = s1;
      for (int in = 0; in < Nin; in++) {
        for (int ex = 0; ex < Nex; ex++) {
          vp[in + Nin * (ix + Nvol * ex)] = bc2 * vt[in + Nin * (s2 + Nvol2 * ex)];
        }
      }
    }
  }
}


//====================================================================
void ShiftField_lex::up_t(Field *v,
                          const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(3) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin   = w->nin();
  const int Nex   = w->nex();
  const int Nvol  = w->nvol();
  const int Nvol2 = Nvol / m_Nt;

  double       *vp = v->ptr(0);
  const double *wp = w->ptr(0);

  std::vector<double> wt(Nin * Nvol2 * Nex);
  std::vector<double> vt(Nin * Nvol2 * Nex);

  for (int s2 = 0; s2 < m_Nx * m_Ny * m_Nz; s2++) {
    // bulk
    for (int t = 0; t < m_Nt - 1; t++) {
      int ix = s2 + m_Nx * m_Ny * m_Nz * t;
      int px = ix + m_Nx * m_Ny * m_Nz;
      for (int ex = 0; ex < Nex; ex++) {
        for (int in = 0; in < Nin; in++) {
          vp[in + Nin * (ix + Nvol * ex)] = wp[in + Nin * (px + Nvol * ex)];
        }
      }
    }

    // boundary (t=m_Nt-1)
    int px = s2;
    for (int ex = 0; ex < Nex; ex++) {
      for (int in = 0; in < Nin; in++) {
        wt[in + Nin * (s2 + Nvol2 * ex)] = bc2 * wp[in + Nin * (px + Nvol * ex)];
      }
    }
  }

  const int size = Nin * Nvol2 * Nex;
  Communicator::exchange(size, &vt[0], &wt[0], 3, 1, 3);

  for (int s2 = 0; s2 < m_Nx * m_Ny * m_Nz; s2++) {
    // boundary (t=m_Nt-1)
    int ix = s2 + m_Nx * m_Ny * m_Nz * (m_Nt - 1);
    for (int ex = 0; ex < Nex; ex++) {
      for (int in = 0; in < Nin; in++) {
        vp[in + Nin * (ix + Nvol * ex)] = vt[in + Nin * (s2 + Nvol2 * ex)];
      }
    }
  }
}


//====================================================================
void ShiftField_lex::dn_t(Field *v,
                          const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(3) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin   = w->nin();
  const int Nex   = w->nex();
  const int Nvol  = w->nvol();
  const int Nvol2 = Nvol / m_Nt;

  double       *vp = v->ptr(0);
  const double *wp = w->ptr(0);

  std::vector<double> wt(Nin * Nvol2 * Nex);
  std::vector<double> vt(Nin * Nvol2 * Nex);

  for (int s2 = 0; s2 < m_Nx * m_Ny * m_Nz; s2++) {
    // bulk
    for (int t = 1; t < m_Nt; t++) {
      int ix = s2 + m_Nx * m_Ny * m_Nz * t;
      int px = ix - m_Nx * m_Ny * m_Nz;
      for (int ex = 0; ex < Nex; ex++) {
        for (int in = 0; in < Nin; in++) {
          vp[in + Nin * (ix + Nvol * ex)] = wp[in + Nin * (px + Nvol * ex)];
        }
      }
    }
    // boundary (t=0)
    int px = s2 + m_Nx * m_Ny * m_Nz * (m_Nt - 1);
    for (int ex = 0; ex < Nex; ex++) {
      for (int in = 0; in < Nin; in++) {
        wt[in + Nin * (s2 + Nvol2 * ex)] = wp[in + Nin * (px + Nvol * ex)];
      }
    }
  }

  const int size = Nin * Nvol2 * Nex;
  Communicator::exchange(size, &vt[0], &wt[0], 3, -1, 7);

  for (int s2 = 0; s2 < m_Nx * m_Ny * m_Nz; s2++) {
    // boundary (t=0)
    int ix = s2;
    for (int ex = 0; ex < Nex; ex++) {
      for (int in = 0; in < Nin; in++) {
        vp[in + Nin * (ix + Nvol * ex)] = bc2 * vt[in + Nin * (s2 + Nvol2 * ex)];
      }
    }
  }
}
#endif

//====================================================================
//============================================================END=====
