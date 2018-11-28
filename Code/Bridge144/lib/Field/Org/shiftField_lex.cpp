#include "BridgeLib_Private.h"
#if USE_ORG

/*!
        @file    $Id:: shiftField_lex.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2017-04-05 18:26:17 #$

        @version $LastChangedRevision: 1608 $
*/

#include "Field/shiftField_lex.h"

const std::string ShiftField_lex::class_name = "ShiftField_lex";

//====================================================================
void ShiftField_lex::backward(Field& v, const Field& w,
                              const int mu)
{
  int boundary_condition = 1;

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
void ShiftField_lex::forward(Field& v, const Field& w,
                             const int mu)
{
  int boundary_condition = 1;

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
void ShiftField_lex::backward(Field& v, const Field& w,
                              const int boundary_condition, const int mu)
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
void ShiftField_lex::forward(Field& v, const Field& w,
                             const int boundary_condition, const int mu)
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
void ShiftField_lex::up_x(Field *v, const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(0) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w->nin();
  const int Nex  = w->nex();
  const int Nvol = w->nvol();

  unique_ptr<Field> wt(new Field(Nin, Nvol / m_Nx, Nex));
  unique_ptr<Field> vt(new Field(Nin, Nvol / m_Nx, Nex));

  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int y = 0; y < m_Ny; y++) {
        int s2 = y + m_Ny * (z + m_Nz * t);
        // bulk
        for (int x = 0; x < m_Nx - 1; x++) {
          int ix = m_index_lex.site(x, y, z, t);
          int px = m_index_lex.site(x + 1, y, z, t);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v->set(in, ix, ex, w->cmp(in, px, ex));
            }
          }
        }
        // boundary (x=m_Nx-1)
        int ix = m_index_lex.site(m_Nx - 1, y, z, t);
        int px = m_index_lex.site(0, y, z, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            wt->set(in, s2, ex, bc2 * w->cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Nx) * Nex;
  Communicator::exchange(size, vt->ptr(0), wt->ptr(0), 0, 1, 0);

  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int y = 0; y < m_Ny; y++) {
        int s2 = y + m_Ny * (z + m_Nz * t);
        // boundary (x=m_Nx-1)
        int ix = m_index_lex.site(m_Nx - 1, y, z, t);
        int px = m_index_lex.site(0, y, z, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            v->set(in, ix, ex, vt->cmp(in, s2, ex));
          }
        }
      }
    }
  }

}


//====================================================================
void ShiftField_lex::up_y(Field *v, const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(1) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w->nin();
  const int Nex  = w->nex();
  const int Nvol = w->nvol();

  unique_ptr<Field> wt(new Field(Nin, Nvol / m_Ny, Nex));
  unique_ptr<Field> vt(new Field(Nin, Nvol / m_Ny, Nex));

  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int x = 0; x < m_Nx; x++) {
        int s2 = x + m_Nx * (z + m_Nz * t);

        // bulk
        for (int y = 0; y < m_Ny - 1; y++) {
          int ix = m_index_lex.site(x, y, z, t);
          int px = m_index_lex.site(x, y + 1, z, t);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v->set(in, ix, ex, w->cmp(in, px, ex));
            }
          }
        }

        // boundary (y=m_Ny-1)
        int ix = m_index_lex.site(x, m_Ny - 1, z, t);
        int px = m_index_lex.site(x, 0, z, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            wt->set(in, s2, ex, bc2 * w->cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Ny) * Nex;
  Communicator::exchange(size, vt->ptr(0), wt->ptr(0), 1, 1, 1);

  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int x = 0; x < m_Nx; x++) {
        int s2 = x + m_Nx * (z + m_Nz * t);

        // boundary (y=m_Ny-1)
        int ix = m_index_lex.site(x, m_Ny - 1, z, t);
        int px = m_index_lex.site(x, 0, z, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            v->set(in, ix, ex, vt->cmp(in, s2, ex));
          }
        }
      }
    }
  }

}


//====================================================================
void ShiftField_lex::up_z(Field *v, const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(2) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w->nin();
  const int Nex  = w->nex();
  const int Nvol = w->nvol();

  unique_ptr<Field> wt(new Field(Nin, Nvol / m_Nz, Nex));
  unique_ptr<Field> vt(new Field(Nin, Nvol / m_Nz, Nex));

  for (int t = 0; t < m_Nt; t++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx; x++) {
        int s2 = x + m_Nx * (y + m_Ny * t);

        // bulk
        for (int z = 0; z < m_Nz - 1; z++) {
          int ix = m_index_lex.site(x, y, z, t);
          int px = m_index_lex.site(x, y, z + 1, t);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v->set(in, ix, ex, w->cmp(in, px, ex));
            }
          }
        }

        // boundary (z=m_Nz-1)
        int ix = m_index_lex.site(x, y, m_Nz - 1, t);
        int px = m_index_lex.site(x, y, 0, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            wt->set(in, s2, ex, bc2 * w->cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Nz) * Nex;
  Communicator::exchange(size, vt->ptr(0), wt->ptr(0), 2, 1, 2);

  for (int t = 0; t < m_Nt; t++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx; x++) {
        int s2 = x + m_Nx * (y + m_Ny * t);

        // boundary (z=m_Nz-1)
        int ix = m_index_lex.site(x, y, m_Nz - 1, t);
        int px = m_index_lex.site(x, y, 0, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            v->set(in, ix, ex, vt->cmp(in, s2, ex));
          }
        }
      }
    }
  }

}


//====================================================================
void ShiftField_lex::up_t(Field *v, const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(3) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w->nin();
  const int Nex  = w->nex();
  const int Nvol = w->nvol();

  unique_ptr<Field> wt(new Field(Nin, Nvol / m_Nt, Nex));
  unique_ptr<Field> vt(new Field(Nin, Nvol / m_Nt, Nex));

  for (int z = 0; z < m_Nz; z++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx; x++) {
        int s2 = x + m_Nx * (y + m_Ny * z);

        // bulk
        for (int t = 0; t < m_Nt - 1; t++) {
          int ix = m_index_lex.site(x, y, z, t);
          int px = m_index_lex.site(x, y, z, t + 1);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v->set(in, ix, ex, w->cmp(in, px, ex));
            }
          }
        }

        // boundary (t=m_Nt-1)
        int ix = m_index_lex.site(x, y, z, m_Nt - 1);
        int px = m_index_lex.site(x, y, z, 0);
        for (int ex = 0; ex < Nex; ex++) {
          for (int in = 0; in < Nin; in++) {
            wt->set(in, s2, ex, bc2 * w->cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Nt) * Nex;
  Communicator::exchange(size, vt->ptr(0), wt->ptr(0), 3, 1, 3);

  for (int z = 0; z < m_Nz; z++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx; x++) {
        int s2 = x + m_Nx * (y + m_Ny * z);

        // boundary (t=m_Nt-1)
        int ix = m_index_lex.site(x, y, z, m_Nt - 1);
        int px = m_index_lex.site(x, y, z, 0);
        for (int ex = 0; ex < Nex; ex++) {
          for (int in = 0; in < Nin; in++) {
            v->set(in, ix, ex, vt->cmp(in, s2, ex));
          }
        }
      }
    }
  }

}


//====================================================================
void ShiftField_lex::dn_x(Field *v, const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(0) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w->nin();
  const int Nex  = w->nex();
  const int Nvol = w->nvol();

  unique_ptr<Field> wt(new Field(Nin, Nvol / m_Nx, Nex));
  unique_ptr<Field> vt(new Field(Nin, Nvol / m_Nx, Nex));

  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int y = 0; y < m_Ny; y++) {
        int s2 = y + m_Ny * (z + m_Nz * t);

        // bulk
        for (int x = 1; x < m_Nx; x++) {
          int ix = m_index_lex.site(x, y, z, t);
          int px = m_index_lex.site(x - 1, y, z, t);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v->set(in, ix, ex, w->cmp(in, px, ex));
            }
          }
        }
        // boundary (x=0)
        //int ix = m_index_lex.site(0   ,y,z,t);
        int px = m_index_lex.site(m_Nx - 1, y, z, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            wt->set(in, s2, ex, w->cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Nx) * Nex;
  Communicator::exchange(size, vt->ptr(0), wt->ptr(0), 0, -1, 4);

  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int y = 0; y < m_Ny; y++) {
        int s2 = y + m_Ny * (z + m_Nz * t);
        // boundary (x=0)
        int ix = m_index_lex.site(0, y, z, t);
        //int px = m_index_lex.site(m_Nx-1,y,z,t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            v->set(in, ix, ex, bc2 * vt->cmp(in, s2, ex));
          }
        }
      }
    }
  }

}


//====================================================================
void ShiftField_lex::dn_y(Field *v, const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(1) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w->nin();
  const int Nex  = w->nex();
  const int Nvol = w->nvol();

  unique_ptr<Field> wt(new Field(Nin, Nvol / m_Ny, Nex));
  unique_ptr<Field> vt(new Field(Nin, Nvol / m_Ny, Nex));

  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int x = 0; x < m_Nx; x++) {
        int s2 = x + m_Nx * (z + m_Nz * t);

        // bulk
        for (int y = 1; y < m_Ny; y++) {
          int ix = m_index_lex.site(x, y, z, t);
          int px = m_index_lex.site(x, y - 1, z, t);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v->set(in, ix, ex, w->cmp(in, px, ex));
            }
          }
        }
        // boundary (y=0)
        int ix = m_index_lex.site(x, 0, z, t);
        int px = m_index_lex.site(x, m_Ny - 1, z, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            wt->set(in, s2, ex, w->cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Ny) * Nex;
  Communicator::exchange(size, vt->ptr(0), wt->ptr(0), 1, -1, 5);

  for (int t = 0; t < m_Nt; t++) {
    for (int z = 0; z < m_Nz; z++) {
      for (int x = 0; x < m_Nx; x++) {
        int s2 = x + m_Nx * (z + m_Nz * t);

        // boundary (y=0)
        int ix = m_index_lex.site(x, 0, z, t);
        int px = m_index_lex.site(x, m_Ny - 1, z, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            v->set(in, ix, ex, bc2 * vt->cmp(in, s2, ex));
          }
        }
      }
    }
  }

}


//====================================================================
void ShiftField_lex::dn_z(Field *v, const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(2) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w->nin();
  const int Nex  = w->nex();
  const int Nvol = w->nvol();

  unique_ptr<Field> wt(new Field(Nin, Nvol / m_Nz, Nex));
  unique_ptr<Field> vt(new Field(Nin, Nvol / m_Nz, Nex));

  for (int t = 0; t < m_Nt; t++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx; x++) {
        int s2 = x + m_Nx * (y + m_Ny * t);

        // bulk
        for (int z = 1; z < m_Nz; z++) {
          int ix = m_index_lex.site(x, y, z, t);
          int px = m_index_lex.site(x, y, z - 1, t);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v->set(in, ix, ex, w->cmp(in, px, ex));
            }
          }
        }
        // boundary (z=0)
        int ix = m_index_lex.site(x, y, 0, t);
        int px = m_index_lex.site(x, y, m_Nz - 1, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            wt->set(in, s2, ex, w->cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Nz) * Nex;
  Communicator::exchange(size, vt->ptr(0), wt->ptr(0), 2, -1, 6);

  for (int t = 0; t < m_Nt; t++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx; x++) {
        int s2 = x + m_Nx * (y + m_Ny * t);

        // boundary (z=0)
        int ix = m_index_lex.site(x, y, 0, t);
        int px = m_index_lex.site(x, y, m_Nz - 1, t);
        for (int in = 0; in < Nin; in++) {
          for (int ex = 0; ex < Nex; ex++) {
            v->set(in, ix, ex, bc2 * vt->cmp(in, s2, ex));
          }
        }
      }
    }
  }

}


//====================================================================
void ShiftField_lex::dn_t(Field *v, const Field *w, const int boundary_condition)
{
  double bc2;

  if (Communicator::ipe(3) == 0) {
    bc2 = boundary_condition;
  } else {
    bc2 = 1.0;
  }

  const int Nin  = w->nin();
  const int Nex  = w->nex();
  const int Nvol = w->nvol();

  unique_ptr<Field> wt(new Field(Nin, Nvol / m_Nt, Nex));
  unique_ptr<Field> vt(new Field(Nin, Nvol / m_Nt, Nex));

  for (int z = 0; z < m_Nz; z++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx; x++) {
        int s2 = x + m_Nx * (y + m_Ny * z);

        // bulk
        for (int t = 1; t < m_Nt; t++) {
          int ix = m_index_lex.site(x, y, z, t);
          int px = m_index_lex.site(x, y, z, t - 1);
          for (int ex = 0; ex < Nex; ex++) {
            for (int in = 0; in < Nin; in++) {
              v->set(in, ix, ex, w->cmp(in, px, ex));
            }
          }
        }
        // boundary (t=0)
        int ix = m_index_lex.site(x, y, z, 0);
        int px = m_index_lex.site(x, y, z, m_Nt - 1);
        for (int ex = 0; ex < Nex; ex++) {
          for (int in = 0; in < Nin; in++) {
            wt->set(in, s2, ex, w->cmp(in, px, ex));
          }
        }
      }
    }
  }

  const int size = Nin * (Nvol / m_Nt) * Nex;
  Communicator::exchange(size, vt->ptr(0), wt->ptr(0), 3, -1, 7);

  for (int z = 0; z < m_Nz; z++) {
    for (int y = 0; y < m_Ny; y++) {
      for (int x = 0; x < m_Nx; x++) {
        int s2 = x + m_Nx * (y + m_Ny * z);

        // boundary (t=0)
        int ix = m_index_lex.site(x, y, z, 0);
        int px = m_index_lex.site(x, y, z, m_Nt - 1);
        for (int ex = 0; ex < Nex; ex++) {
          for (int in = 0; in < Nin; in++) {
            v->set(in, ix, ex, bc2 * vt->cmp(in, s2, ex));
          }
        }
      }
    }
  }

}


//====================================================================
//============================================================END=====
#endif
