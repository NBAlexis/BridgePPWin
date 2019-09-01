/*!
        @file    field.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "field.h"

#include <cstring>
#include "ResourceManager/threadManager_OpenMP.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

using std::string;

//====================================================================
namespace {
  inline
  void set_threadtask(int& i_thread, int& Nthread, int& is, int& ns,
                      const int size)
  {
    Nthread  = ThreadManager_OpenMP::get_num_threads();
    i_thread = ThreadManager_OpenMP::get_thread_id();

    is = static_cast<long_t>(size) * i_thread / Nthread;
    ns = static_cast<long_t>(size) * (i_thread + 1) / Nthread;
  }
}

//====================================================================
void Field::check()
{
  //    vout.general("Field was constructed.\n");
}


//====================================================================
double dot(const Field& y, const Field& x)
{
  const double *yp = y.ptr(0);
  const double *xp = x.ptr(0);

  //  int size = x.ntot();
  int i_thread, Nthread, is, ns;

  set_threadtask(i_thread, Nthread, is, ns, x.ntot());

  double a = 0.0;

  for (int k = is; k < ns; ++k) {
    a += yp[k] * xp[k];
  }
  ThreadManager_OpenMP::reduce_sum_global(a, i_thread, Nthread);

  return a;
}


//====================================================================
double dot(const Field& y, const int exy, const Field& x, const int exx)
{
  assert(x.nin() == y.nin());
  assert(x.nvol() == y.nvol());

  const double *yp = y.ptr(0, 0, exy);
  const double *xp = x.ptr(0, 0, exx);

  //  int size = x.nin() * x.nvol();
  int i_thread, Nthread, is, ns;
  set_threadtask(i_thread, Nthread, is, ns, x.ntot());

  double a = 0.0;

  for (int k = is; k < ns; ++k) {
    a += yp[k] * xp[k];
  }
  ThreadManager_OpenMP::reduce_sum_global(a, i_thread, Nthread);

  return a;
}


//====================================================================
dcomplex dotc(const Field& y, const Field& x)
{
  const double *yp = y.ptr(0);
  const double *xp = x.ptr(0);

  assert(x.ntot() == y.ntot());

  //  int ntot = x.ntot();
  int i_thread, Nthread, is, ns;
  set_threadtask(i_thread, Nthread, is, ns, x.ntot());

  if ((y.field_element_type() == Element_type::COMPLEX) &&
      (x.field_element_type() == Element_type::COMPLEX)) {
    double prd_r = 0.0;
    double prd_i = 0.0;

    //    for (int k = 0; k < ntot; k += 2) {
    for (int k = is; k < ns; k += 2) {
      prd_r += yp[k] * xp[k] + yp[k + 1] * xp[k + 1];
      prd_i += yp[k] * xp[k + 1] - yp[k + 1] * xp[k];
    }
    ThreadManager_OpenMP::reduce_sum_global(prd_r, i_thread, Nthread);
    ThreadManager_OpenMP::reduce_sum_global(prd_i, i_thread, Nthread);

    return cmplx(prd_r, prd_i);
  } else if ((y.field_element_type() == Element_type::REAL) &&
             (y.field_element_type() == Element_type::REAL)) {
    return cmplx(dot(y, x), 0.0);
  } else {
    vout.crucial("Error at %s: unsupported arg types.\n", __func__);
    exit(EXIT_FAILURE);

    //return cmplx(0.0, 0.0);  // never reached.
  }
}


//====================================================================
dcomplex dotc(const Field& y, const int exy, const Field& x, const int exx)
{
  const double *yp = y.ptr(0, 0, exy);
  const double *xp = x.ptr(0, 0, exx);

  assert(x.nin() == y.nin());
  assert(x.nvol() == y.nvol());

  //  int size = x.nin() * x.nvol();
  int i_thread, Nthread, is, ns;
  set_threadtask(i_thread, Nthread, is, ns, x.nin() * x.nvol());

  if ((y.field_element_type() == Element_type::COMPLEX) &&
      (x.field_element_type() == Element_type::COMPLEX)) {
    double prd_r = 0.0;
    double prd_i = 0.0;

    for (int k = is; k < ns; k += 2) {
      prd_r += yp[k] * xp[k] + yp[k + 1] * xp[k + 1];
      prd_i += yp[k] * xp[k + 1] - yp[k + 1] * xp[k];
    }
    ThreadManager_OpenMP::reduce_sum_global(prd_r, i_thread, Nthread);
    ThreadManager_OpenMP::reduce_sum_global(prd_i, i_thread, Nthread);

    return cmplx(prd_r, prd_i);
  } else if ((y.field_element_type() == Element_type::REAL) &&
             (y.field_element_type() == Element_type::REAL)) {
    return cmplx(dot(y, exy, x, exx), 0.0);
  } else {
    vout.crucial("Error at %s: unsupported arg types.\n", __func__);
    exit(EXIT_FAILURE);

    //return cmplx(0.0, 0.0);  // never reached.
  }
}


//====================================================================
void axpy(Field& y, const double a, const Field& x)
{
  double       *yp = y.ptr(0);
  const double *xp = x.ptr(0);

  //  int size = x.ntot();
  int i_thread, Nthread, is, ns;

  set_threadtask(i_thread, Nthread, is, ns, x.ntot());

  for (int k = is; k < ns; ++k) {
    yp[k] += a * xp[k];
  }
}


//====================================================================
void axpy(Field& y, const int exy, const double a, const Field& x, const int exx)
{
  assert(x.nin() == y.nin());
  assert(x.nvol() == y.nvol());

  double       *yp = y.ptr(0, 0, exy);
  const double *xp = x.ptr(0, 0, exx);

  //  int size = x.nin() * x.nvol();
  int i_thread, Nthread, is, ns;
  set_threadtask(i_thread, Nthread, is, ns, x.nin() * x.nvol());

  for (int k = is; k < ns; ++k) {
    yp[k] += a * xp[k];
  }
}


//====================================================================
void axpy(Field& y, const dcomplex a, const Field& x)
{
  if (imag(a) == 0.0) {
    return axpy(y, real(a), x);
  } else if ((y.field_element_type() == Element_type::REAL) &&
             (x.field_element_type() == Element_type::REAL)) {
    vout.crucial("Error at %s: real vector and complex parameter.\n", __func__);
    exit(EXIT_FAILURE);

    //return axpy(y, real(a), x);  // ignore imaginary part of a
  } else if ((y.field_element_type() == Element_type::COMPLEX) &&
             (x.field_element_type() == Element_type::COMPLEX)) {
    double       *yp = y.ptr(0);
    const double *xp = x.ptr(0);

    assert(x.ntot() == y.ntot());

    //    int ntot = x.ntot();
    int i_thread, Nthread, is, ns;
    set_threadtask(i_thread, Nthread, is, ns, x.ntot());

    double ar = real(a);
    double ai = imag(a);

    //    for (int k = 0; k < ntot; k += 2) {
    for (int k = is; k < ns; k += 2) {
      yp[k]     += ar * xp[k] - ai * xp[k + 1];
      yp[k + 1] += ar * xp[k + 1] + ai * xp[k];
    }

    return;
  } else {
    vout.crucial("Error at %s: unsupported types.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void axpy(Field& y, const int exy, const dcomplex a, const Field& x, const int exx)
{
  if (imag(a) == 0.0) {
    return axpy(y, exy, real(a), x, exx);
  } else if ((y.field_element_type() == Element_type::REAL) &&
             (x.field_element_type() == Element_type::REAL)) {
    vout.crucial("Error at %s: real vector and complex parameter.\n", __func__);
    exit(EXIT_FAILURE);

    //return axpy(y, real(a), x);  // ignore imaginary part of a
  } else if ((y.field_element_type() == Element_type::COMPLEX) &&
             (x.field_element_type() == Element_type::COMPLEX)) {
    double       *yp = y.ptr(0, 0, exy);
    const double *xp = x.ptr(0, 0, exx);

    assert(x.nin() == y.nin());
    assert(x.nvol() == y.nvol());

    //    int size = x.nin() * x.nvol();
    int i_thread, Nthread, is, ns;
    set_threadtask(i_thread, Nthread, is, ns, x.nin() * x.nvol());

    double ar = real(a);
    double ai = imag(a);

    for (int k = is; k < ns; k += 2) {
      yp[k]     += ar * xp[k] - ai * xp[k + 1];
      yp[k + 1] += ar * xp[k + 1] + ai * xp[k];
    }

    return;
  } else {
    vout.crucial("Error at %s: unsupported types.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void scal(Field& x, const double a)
{
  //  x.field *= a;

  double *xp = x.ptr(0, 0, 0);
  int    i_thread, Nthread, is, ns;

  set_threadtask(i_thread, Nthread, is, ns, x.ntot());

  for (int k = is; k < ns; ++k) {
    xp[k] *= a;
  }
}


//====================================================================
void scal(Field& x, const int exx, const double a)
{
  double *xp = x.ptr(0, 0, exx);
  //  int size = x.nin() * x.nvol();
  int i_thread, Nthread, is, ns;

  set_threadtask(i_thread, Nthread, is, ns, x.nin() * x.nvol());

  //  for (int k = 0; k < size; ++k) {
  for (int k = is; k < ns; ++k) {
    //    x.field[k] *= a;
    xp[k] *= a;
  }
}


//====================================================================
void scal(Field& x, const dcomplex a)
{
  if (x.field_element_type() == Element_type::REAL) {
    vout.crucial("Error at %s: real vector and complex parameter.\n", __func__);
    exit(EXIT_FAILURE);

    //    x.field *= real(a);  // ignore imaginary part of a
    //scal(x, a);
  } else if (x.field_element_type() == Element_type::COMPLEX) {
    double *xp = x.ptr(0);
    //int ntot = x.ntot();
    int i_thread, Nthread, is, ns;
    set_threadtask(i_thread, Nthread, is, ns, x.ntot());

    double ar = real(a);
    double ai = imag(a);

    //    for (int k = 0; k < ntot; k += 2) {
    for (int k = is; k < ns; k += 2) {
      double xr = xp[k];
      double xi = xp[k + 1];

      xp[k]     = ar * xr - ai * xi;
      xp[k + 1] = ar * xi + ai * xr;
    }
  } else {
    vout.crucial("Error at %s: unsupported field type.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void scal(Field& x, const int exx, const dcomplex a)
{
  if (x.field_element_type() == Element_type::REAL) {
    vout.crucial("Error at %s: real vector and complex parameter.\n", __func__);
    exit(EXIT_FAILURE);

    //    x.field *= real(a);  // ignore imaginary part of a
    //scal(x, exx, a);
  } else if (x.field_element_type() == Element_type::COMPLEX) {
    double *xp = x.ptr(0, 0, exx);
    //    int size = x.nin() * x.nvol();
    int i_thread, Nthread, is, ns;
    set_threadtask(i_thread, Nthread, is, ns, x.nin() * x.nvol());

    double ar = real(a);
    double ai = imag(a);

    //    for (int k = 0; k < size; k += 2) {
    for (int k = is; k < ns; k += 2) {
      double xr = xp[k];
      double xi = xp[k + 1];

      xp[k]     = ar * xr - ai * xi;
      xp[k + 1] = ar * xi + ai * xr;
    }
  } else {
    vout.crucial("Error ar %s: unsupported field type.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void copy(Field& y, const Field& x)
{
  //  y.field = x.field;

  assert(x.nin() == y.nin());
  assert(x.nvol() == y.nvol());
  assert(x.nex() == y.nex());

  double       *yp = y.ptr(0, 0, 0);
  const double *xp = x.ptr(0, 0, 0);

  int i_thread, Nthread, is, ns;
  set_threadtask(i_thread, Nthread, is, ns, x.ntot());

  for (int k = is; k < ns; ++k) {
    yp[k] = xp[k];
  }
}


//====================================================================
void copy(Field& y, const int exy, const Field& x, const int exx)
{
  assert(x.nin() == y.nin());
  assert(x.nvol() == y.nvol());

  double       *yp = y.ptr(0, 0, exy);
  const double *xp = x.ptr(0, 0, exx);

  //  int size = x.nin() * x.nvol();
  int i_thread, Nthread, is, ns;
  set_threadtask(i_thread, Nthread, is, ns, x.nin() * x.nvol());

  for (int k = is; k < ns; ++k) {
    yp[k] = xp[k];
  }

  ////  for (int k = 0; k < size; ++k) {
  ////    yp[k] = xp[k];
  ////  }
  //  memcpy(yp, xp, sizeof(double) * size);
}


//====================================================================
void Field::set(double a)
{
  int i_thread, Nthread, is, ns;

  set_threadtask(i_thread, Nthread, is, ns, ntot());

  double *yp = this->ptr(0);

  for (int k = is; k < ns; ++k) {
    yp[k] = a;
  }
}


//====================================================================
double Field::norm2() const
{
  int i_thread, Nthread, is, ns;

  set_threadtask(i_thread, Nthread, is, ns, ntot());

  const double *yp = this->ptr(0);

  double a = 0.0;

  for (int k = is; k < ns; ++k) {
    a += yp[k] * yp[k];
  }
  ThreadManager_OpenMP::reduce_sum_global(a, i_thread, Nthread);

  return a;
}


//====================================================================
void aypx(const double a, Field& y, const Field& x)
{
  int i_thread, Nthread, is, ns;

  set_threadtask(i_thread, Nthread, is, ns, x.ntot());

  double       *yp = y.ptr(0);
  const double *xp = x.ptr(0);

  for (int k = is; k < ns; ++k) {
    yp[k] = a * yp[k] + xp[k];
  }
}


//====================================================================
void aypx(const dcomplex a, Field& y, const Field& x)
{
  if (imag(a) == 0.0) {
    return aypx(real(a), y, x);
  } else if ((y.field_element_type() == Element_type::REAL) &&
             (x.field_element_type() == Element_type::REAL)) {
    vout.crucial("Error at %s: real vector and complex parameter.\n", __func__);
    exit(EXIT_FAILURE);

    //return aypx(real(a), y, x);  // ignore imaginary part of a
  } else if ((y.field_element_type() == Element_type::COMPLEX) &&
             (x.field_element_type() == Element_type::COMPLEX)) {
    double       *yp = y.ptr(0);
    const double *xp = x.ptr(0);

    assert(x.ntot() == y.ntot());

    int i_thread, Nthread, is, ns;
    set_threadtask(i_thread, Nthread, is, ns, x.ntot());

    double ar = real(a);
    double ai = imag(a);

    for (int k = is; k < ns; k += 2) {
      double ypr = yp[k];
      double ypi = yp[k + 1];
      yp[k]     = ar * ypr - ai * ypi + xp[k];
      yp[k + 1] = ar * ypi + ai * ypr + xp[k + 1];
    }

    return;
  } else {
    vout.crucial("Error at %s: unsupported types.\n", __func__);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Field::stat(double& Fave, double& Fmax, double& Fdev) const
{
  assert(ThreadManager_OpenMP::get_num_threads() == 1);

  double sum  = 0.0;
  double sum2 = 0.0;

  Fmax = 0.0;
  for (int ex = 0, nnex = nex(); ex < nnex; ++ex) {
    for (int site = 0, nnvol = nvol(); site < nnvol; ++site) {
      double fst = 0.0;
      for (int in = 0, nnin = nin(); in < nnin; ++in) {
        double fv = field[myindex(in, site, ex)];
        fst += fv * fv;
      }
      sum2 += fst;
      fst   = sqrt(fst);
      sum  += fst;
      if (fst > Fmax) Fmax = fst;
    }
  }

  sum  = Communicator::reduce_sum(sum);
  sum2 = Communicator::reduce_sum(sum2);
  Fmax = Communicator::reduce_max(Fmax);

  double perdeg = 1.0 / ((double)CommonParameters::Lvol() * nex());
  Fave = sum * perdeg;

  double fval = sum2 * perdeg;
  fval -= Fave * Fave;
  Fdev  = sqrt(fval);
}


//====================================================================
void report_field_stat(const Bridge::VerboseLevel vl,
                       const std::string& msg, const Field& f)
{
  assert(ThreadManager_OpenMP::get_num_threads() == 1);

  double favg, fmax, fdev;
  f.stat(favg, fmax, fdev);

  vout.general(vl,
               "    %s: avg = %12.6f  max = %12.6f  dev = %12.6f\n",
               msg.c_str(),
               favg, fmax, fdev);
}


//====================================================================
//============================================================END=====
