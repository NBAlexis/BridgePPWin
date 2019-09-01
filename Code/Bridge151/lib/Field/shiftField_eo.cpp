/*!
        @file    shiftField_eo.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "shiftField_eo.h"

const std::string ShiftField_eo::class_name = "ShiftField_eo";

//====================================================================
void ShiftField_eo::backward_h(Field& v, const Field& w,
                               const int mu, const int ieo)
{
  const int boundary_condition = 1;

  if (mu == 0) {        // x-direction
    up_xh(v, w, boundary_condition, ieo);
  } else if (mu == 1) { // y-direction
    up_yh(v, w, boundary_condition, ieo);
  } else if (mu == 2) { // z-direction
    up_zh(v, w, boundary_condition, ieo);
  } else if (mu == 3) { // t-direction
    up_th(v, w, boundary_condition, ieo);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_eo::forward_h(Field& v, const Field& w,
                              const int mu, const int ieo)
{
  const int boundary_condition = 1;

  if (mu == 0) {        // x-direction
    dn_xh(v, w, boundary_condition, ieo);
  } else if (mu == 1) { // y-direction
    dn_yh(v, w, boundary_condition, ieo);
  } else if (mu == 2) { // z-direction
    dn_zh(v, w, boundary_condition, ieo);
  } else if (mu == 3) { // t-direction
    dn_th(v, w, boundary_condition, ieo);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_eo::backward_h(Field& v, const Field& w,
                               const int boundary_condition, const int mu, const int ieo)
{
  if (mu == 0) {        // x-direction
    up_xh(v, w, boundary_condition, ieo);
  } else if (mu == 1) { // y-direction
    up_yh(v, w, boundary_condition, ieo);
  } else if (mu == 2) { // z-direction
    up_zh(v, w, boundary_condition, ieo);
  } else if (mu == 3) { // t-direction
    up_th(v, w, boundary_condition, ieo);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_eo::forward_h(Field& v, const Field& w,
                              const int boundary_condition, const int mu, const int ieo)
{
  if (mu == 0) {        // x-direction
    dn_xh(v, w, boundary_condition, ieo);
  } else if (mu == 1) { // y-direction
    dn_yh(v, w, boundary_condition, ieo);
  } else if (mu == 2) { // z-direction
    dn_zh(v, w, boundary_condition, ieo);
  } else if (mu == 3) { // t-direction
    dn_th(v, w, boundary_condition, ieo);
  } else {
    vout.crucial(m_vl, "Error at %s: wrong mu = %d\n", class_name.c_str(), mu);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void ShiftField_eo::backward(Field& v, const Field& w, const int mu)
{
  const int Nin   = w.nin();
  const int Nvol2 = w.nvol() / 2;
  const int Nex   = w.nex();

  Field field_e(Nin, Nvol2, Nex);
  Field field_o(Nin, Nvol2, Nex);

  Field field_se(Nin, Nvol2, Nex);
  Field field_so(Nin, Nvol2, Nex);

  m_index_eo.splitField(field_e, field_o, w);

  backward_h(field_se, field_o, mu, 0);
  backward_h(field_so, field_e, mu, 1);

  m_index_eo.mergeField(v, field_se, field_so);
}


//====================================================================
void ShiftField_eo::forward(Field& v, const Field& w, const int mu)
{
  const int Nin   = w.nin();
  const int Nvol2 = w.nvol() / 2;
  const int Nex   = w.nex();

  Field field_e(Nin, Nvol2, Nex);
  Field field_o(Nin, Nvol2, Nex);

  Field field_se(Nin, Nvol2, Nex);
  Field field_so(Nin, Nvol2, Nex);

  m_index_eo.splitField(field_e, field_o, w);

  forward_h(field_se, field_o, mu, 0);
  forward_h(field_so, field_e, mu, 1);

  m_index_eo.mergeField(v, field_se, field_so);
}


//====================================================================
void ShiftField_eo::backward(Field& v, const Field& w,
                             const int boundary_condition, const int mu)
{
  const int Nin   = w.nin();
  const int Nvol2 = w.nvol() / 2;
  const int Nex   = w.nex();

  Field field_e(Nin, Nvol2, Nex);
  Field field_o(Nin, Nvol2, Nex);

  Field field_se(Nin, Nvol2, Nex);
  Field field_so(Nin, Nvol2, Nex);

  m_index_eo.splitField(field_e, field_o, w);

  backward_h(field_se, field_o, boundary_condition, mu, 0);
  backward_h(field_so, field_e, boundary_condition, mu, 1);

  m_index_eo.mergeField(v, field_se, field_so);
}


//====================================================================
void ShiftField_eo::forward(Field& v, const Field& w,
                            const int boundary_condition, const int mu)
{
  const int Nin   = w.nin();
  const int Nvol2 = w.nvol() / 2;
  const int Nex   = w.nex();

  Field field_e(Nin, Nvol2, Nex);
  Field field_o(Nin, Nvol2, Nex);

  Field field_se(Nin, Nvol2, Nex);
  Field field_so(Nin, Nvol2, Nex);

  m_index_eo.splitField(field_e, field_o, w);

  forward_h(field_se, field_o, boundary_condition, mu, 0);
  forward_h(field_so, field_e, boundary_condition, mu, 1);

  m_index_eo.mergeField(v, field_se, field_so);
}


//====================================================================
void ShiftField_eo::up_xh(Field& v, const Field& w, const int boundary_condition,
                          const int ieo)
{
  double bc2;

  if (Communicator::ipe(0) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin   = w.nin();
  const int Nex   = w.nex();
  const int Nvol  = w.nvol();
  const int Nvol2 = (1 + Nvol / m_Nx2) / 2;

  Field wt(Nin, Nvol2, Nex);
  Field vt(Nin, Nvol2, Nex);

  int s2 = 0;
  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int y = 0; y < m_Ny; y++) {
        if (m_index_eo.leo(y, z, t) == ieo) {
          // bulk
          for (int x = 0; x < m_Nx2; x++) {
            int ix = m_index_eo.siteh(x, y, z, t);
            for (int ex = 0; ex < Nex; ex++) {
              for (int in = 0; in < Nin; in++) {
                v.set(in, ix, ex, w.cmp(in, ix, ex));
              }
            }
          }
        } else {
          // bulk
          for (int x = 0; x < m_Nx2 - 1; x++) {
            int ix = m_index_eo.siteh(x, y, z, t);
            int px = m_index_eo.siteh_xup(x, y, z, t, ieo);
            for (int ex = 0; ex < Nex; ex++) {
              for (int in = 0; in < Nin; in++) {
                v.set(in, ix, ex, w.cmp(in, px, ex));
              }
            }
          }
          // boundary (x=m_Nx2-1)
          int px = m_index_eo.siteh(0, y, z, t);
          for (int in = 0; in < Nin; in++) {
            for (int ex = 0; ex < Nex; ex++) {
              wt.set(in, s2, ex, bc2 * w.cmp(in, px, ex));
            }
          }

          s2++;
        }
      }
    }
  }

  if (s2 > Nvol2) {
    vout.crucial(m_vl, "Error at %s: invalid size\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  const int size = Nin * Nvol2 * Nex;
  exchange(size, &vt, &wt, 0, 1, 0);

  s2 = 0;
  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int y = 0; y < m_Ny; y++) {
        // boundary (x=m_Nx-1)
        if (m_index_eo.leo(y, z, t) != ieo) {
          int ix = m_index_eo.siteh(m_Nx2 - 1, y, z, t);

          for (int in = 0; in < Nin; in++) {
            for (int ex = 0; ex < Nex; ex++) {
              v.set(in, ix, ex, vt.cmp(in, s2, ex));
            }
          }

          s2++;
        }
      }
    }
  }
}


//====================================================================
void ShiftField_eo::dn_xh(Field& v, const Field& w, const int boundary_condition,
                          const int ieo)
{
  double bc2;

  if (Communicator::ipe(0) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin   = w.nin();
  const int Nex   = w.nex();
  const int Nvol  = w.nvol();
  const int Nvol2 = (1 + Nvol / m_Nx2) / 2;

  Field wt(Nin, Nvol2, Nex);
  Field vt(Nin, Nvol2, Nex);

  int s2 = 0;
  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int y = 0; y < m_Ny; y++) {
        if (m_index_eo.leo(y, z, t) == (1 - ieo)) {
          // bulk
          for (int x = 0; x < m_Nx2; x++) {
            int ix = m_index_eo.siteh(x, y, z, t);
            for (int ex = 0; ex < Nex; ex++) {
              for (int in = 0; in < Nin; in++) {
                v.set(in, ix, ex, w.cmp(in, ix, ex));
              }
            }
          }
        } else {
          // bulk
          for (int x = 1; x < m_Nx2; x++) {
            int ix = m_index_eo.siteh(x, y, z, t);
            int mx = m_index_eo.siteh_xdn(x, y, z, t, ieo);
            for (int ex = 0; ex < Nex; ex++) {
              for (int in = 0; in < Nin; in++) {
                v.set(in, ix, ex, w.cmp(in, mx, ex));
              }
            }
          }
          // boundary (x=0)
          int mx = m_index_eo.siteh(m_Nx2 - 1, y, z, t);
          for (int in = 0; in < Nin; in++) {
            for (int ex = 0; ex < Nex; ex++) {
              wt.set(in, s2, ex, w.cmp(in, mx, ex));
            }
          }

          s2++;
        }
      }
    }
  }

  const int size = Nin * Nvol2 * Nex;
  exchange(size, &vt, &wt, 0, -1, 4);

  s2 = 0;
  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int y = 0; y < m_Ny; y++) {
        // boundary (x=0)
        if (m_index_eo.leo(y, z, t) != (1 - ieo)) {
          int ix = m_index_eo.siteh(0, y, z, t);

          for (int in = 0; in < Nin; in++) {
            for (int ex = 0; ex < Nex; ex++) {
              v.set(in, ix, ex, bc2 * vt.cmp(in, s2, ex));
            }
          }

          s2++;
        }
      }
    }
  }
}


//====================================================================
void ShiftField_eo::up_yh(Field& v, const Field& w, const int boundary_condition,
                          const int ieo)
{
  double bc2;

  if (Communicator::ipe(1) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w.nin();
  const int Nex  = w.nex();
  const int Nvol = w.nvol();

  Field wt(Nin, Nvol / m_Ny, Nex);
  Field vt(Nin, Nvol / m_Ny, Nex);

  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int x = 0; x < m_Nx2; x++) {
        // bulk
        for (int y = 0; y < m_Ny - 1; y++) {
          int ix = m_index_eo.siteh(x, y, z, t);
          int px = m_index_eo.siteh(x, y + 1, z, t);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v.set(in, ix, ex, w.cmp(in, px, ex));
            }
          }
        }

        // boundary (y=m_Ny-1)
        int s2 = x + m_Nx2 * (z + m_Nz * t);
        int px = m_index_eo.siteh(x, 0, z, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            wt.set(in, s2, ex, bc2 * w.cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Ny) * Nex;
  exchange(size, &vt, &wt, 1, 1, 1);

  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int x = 0; x < m_Nx2; x++) {
        // boundary (y=m_Ny-1)
        int s2 = x + m_Nx2 * (z + m_Nz * t);
        int ix = m_index_eo.siteh(x, m_Ny - 1, z, t);

        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            v.set(in, ix, ex, vt.cmp(in, s2, ex));
          }
        }
      }
    }
  }
}


//====================================================================
void ShiftField_eo::dn_yh(Field& v, const Field& w, const int boundary_condition,
                          const int ieo)
{
  double bc2;

  if (Communicator::ipe(1) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w.nin();
  const int Nex  = w.nex();
  const int Nvol = w.nvol();

  Field wt(Nin, Nvol / m_Ny, Nex);
  Field vt(Nin, Nvol / m_Ny, Nex);

  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int x = 0; x < m_Nx2; x++) {
        // bulk
        for (int y = 1; y < m_Ny; y++) {
          int ix = m_index_eo.siteh(x, y, z, t);
          int px = m_index_eo.siteh(x, y - 1, z, t);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v.set(in, ix, ex, w.cmp(in, px, ex));
            }
          }
        }
        // boundary (y=0)
        int s2 = x + m_Nx2 * (z + m_Nz * t);
        int px = m_index_eo.siteh(x, m_Ny - 1, z, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            wt.set(in, s2, ex, w.cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Ny) * Nex;
  exchange(size, &vt, &wt, 1, -1, 5);

  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int x = 0; x < m_Nx2; x++) {
        int s2 = x + m_Nx2 * (z + m_Nz * t);

        // boundary (y=0)
        int ix = m_index_eo.siteh(x, 0, z, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            v.set(in, ix, ex, bc2 * vt.cmp(in, s2, ex));
          }
        }
      }
    }
  }
}


//====================================================================
void ShiftField_eo::up_zh(Field& v, const Field& w, const int boundary_condition,
                          const int ieo)
{
  double bc2;

  if (Communicator::ipe(2) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w.nin();
  const int Nex  = w.nex();
  const int Nvol = w.nvol();

  Field wt(Nin, Nvol / m_Nz, Nex);
  Field vt(Nin, Nvol / m_Nz, Nex);

  for (int t = 0; t < m_Nt; t++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx2; x++) {
        // bulk
        for (int z = 0; z < m_Nz - 1; z++) {
          int ix = m_index_eo.siteh(x, y, z, t);
          int px = m_index_eo.siteh(x, y, z + 1, t);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v.set(in, ix, ex, w.cmp(in, px, ex));
            }
          }
        }

        // boundary (z=m_Nz-1)
        int s2 = x + m_Nx2 * (y + m_Ny * t);
        int px = m_index_eo.siteh(x, y, 0, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            wt.set(in, s2, ex, bc2 * w.cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Nz) * Nex;
  exchange(size, &vt, &wt, 2, 1, 2);

  for (int t = 0; t < m_Nt; t++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx2; x++) {
        int s2 = x + m_Nx2 * (y + m_Ny * t);

        // boundary (z=m_Nz-1)
        int ix = m_index_eo.siteh(x, y, m_Nz - 1, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            v.set(in, ix, ex, vt.cmp(in, s2, ex));
          }
        }
      }
    }
  }
}


//====================================================================
void ShiftField_eo::dn_zh(Field& v, const Field& w, const int boundary_condition,
                          const int ieo)
{
  double bc2;

  if (Communicator::ipe(2) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w.nin();
  const int Nex  = w.nex();
  const int Nvol = w.nvol();

  Field wt(Nin, Nvol / m_Nz, Nex);
  Field vt(Nin, Nvol / m_Nz, Nex);

  for (int t = 0; t < m_Nt; t++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx2; x++) {
        int s2 = x + m_Nx2 * (y + m_Ny * t);

        // bulk
        for (int z = 1; z < m_Nz; z++) {
          int ix = m_index_eo.siteh(x, y, z, t);
          int px = m_index_eo.siteh(x, y, z - 1, t);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v.set(in, ix, ex, w.cmp(in, px, ex));
            }
          }
        }
        // boundary (z=0)
        int px = m_index_eo.siteh(x, y, m_Nz - 1, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            wt.set(in, s2, ex, w.cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Nz) * Nex;
  exchange(size, &vt, &wt, 2, -1, 6);

  for (int t = 0; t < m_Nt; t++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx2; x++) {
        int s2 = x + m_Nx2 * (y + m_Ny * t);

        // boundary (z=0)
        int ix = m_index_eo.siteh(x, y, 0, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            v.set(in, ix, ex, bc2 * vt.cmp(in, s2, ex));
          }
        }
      }
    }
  }
}


//====================================================================
void ShiftField_eo::up_th(Field& v, const Field& w, const int boundary_condition,
                          const int ieo)
{
  double bc2;

  if (Communicator::ipe(3) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w.nin();
  const int Nex  = w.nex();
  const int Nvol = w.nvol();

  Field wt(Nin, Nvol / m_Nt, Nex);
  Field vt(Nin, Nvol / m_Nt, Nex);

  for (int z = 0; z < m_Nz; z++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx2; x++) {
        int s2 = x + m_Nx2 * (y + m_Ny * z);

        // bulk
        for (int t = 0; t < m_Nt - 1; t++) {
          int ix = m_index_eo.siteh(x, y, z, t);
          int px = m_index_eo.siteh(x, y, z, t + 1);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v.set(in, ix, ex, w.cmp(in, px, ex));
            }
          }
        }

        // boundary (t=m_Nt-1)
        int px = m_index_eo.siteh(x, y, z, 0);
        for (int ex = 0; ex < Nex; ex++) {
          for (int in = 0; in < Nin; in++) {
            wt.set(in, s2, ex, bc2 * w.cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Nt) * Nex;
  exchange(size, &vt, &wt, 3, 1, 3);

  for (int z = 0; z < m_Nz; z++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx2; x++) {
        int s2 = x + m_Nx2 * (y + m_Ny * z);

        // boundary (t=m_Nt-1)
        int ix = m_index_eo.siteh(x, y, z, m_Nt - 1);
        for (int ex = 0; ex < Nex; ex++) {
          for (int in = 0; in < Nin; in++) {
            v.set(in, ix, ex, vt.cmp(in, s2, ex));
          }
        }
      }
    }
  }
}


//====================================================================
void ShiftField_eo::dn_th(Field& v, const Field& w, const int boundary_condition,
                          const int ieo)
{
  double bc2;

  if (Communicator::ipe(3) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w.nin();
  const int Nex  = w.nex();
  const int Nvol = w.nvol();

  Field wt(Nin, Nvol / m_Nt, Nex);
  Field vt(Nin, Nvol / m_Nt, Nex);

  for (int z = 0; z < m_Nz; z++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx2; x++) {
        int s2 = x + m_Nx2 * (y + m_Ny * z);

        //- bulk
        for (int t = 1; t < m_Nt; t++) {
          int ix = m_index_eo.siteh(x, y, z, t);
          int px = m_index_eo.siteh(x, y, z, t - 1);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v.set(in, ix, ex, w.cmp(in, px, ex));
            }
          }
        }
        //- boundary (t=0)
        // int ix = m_index_eo.siteh(x, y, z, 0);
        int px = m_index_eo.siteh(x, y, z, m_Nt - 1);
        for (int ex = 0; ex < Nex; ex++) {
          for (int in = 0; in < Nin; in++) {
            wt.set(in, s2, ex, w.cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Nt) * Nex;
  exchange(size, &vt, &wt, 3, -1, 7);

  for (int z = 0; z < m_Nz; z++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx2; x++) {
        int s2 = x + m_Nx2 * (y + m_Ny * z);

        //- boundary (t=0)
        int ix = m_index_eo.siteh(x, y, z, 0);
        // int px = m_index_eo.siteh(x, y, z, m_Nt - 1);
        for (int ex = 0; ex < Nex; ex++) {
          for (int in = 0; in < Nin; in++) {
            v.set(in, ix, ex, bc2 * vt.cmp(in, s2, ex));
          }
        }
      }
    }
  }
}


//====================================================================
//============================================================END=====
