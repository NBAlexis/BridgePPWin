#include "BridgeLib_Private.h"

#ifdef COMMUNICATOR_USE_SINGLE

/*!
        @file    $Id: communicator_single.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-03-21 15:28:34 #$

        @version $LastChangedRevision: 1679 $
 */

#include "Communicator/communicator.h"

#include <cstdarg>
#include <cstring>
#include <cassert>

#include <sys/time.h>
#include <time.h>

#include "Parameters/commonParameters.h"
//#include "Field/field.h"

//namespace Communicator {
static int m_Ndim = CommonParameters::Ndim();
//}

//====================================================================
int Communicator::init(int *pargc, char ***pargv)
{
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::finalize()
{
  return EXIT_SUCCESS;
}


//====================================================================
void Communicator::abort()
{
  ::abort();  // call system function.
}


//====================================================================
int Communicator::setup(int ninstance)
{
  m_Ndim = CommonParameters::Ndim();
  return EXIT_SUCCESS;
}


//====================================================================
bool Communicator::is_primary()
{
  return true;
}


//====================================================================
bool Communicator::is_primary_master()
{
  return true;
}


//====================================================================
int Communicator::self()
{
  return 0;
}


//====================================================================
int Communicator::size()
{
  return 1;
}


//====================================================================
int Communicator::self_global()
{
  return 0;
}


//====================================================================
int Communicator::ipe(const int dir)
{
  return 0;
}


//====================================================================
int Communicator::npe(const int dir)
{
  return 1;
}


//====================================================================
int Communicator::grid_rank(int *rank, const int *gcoord)
{
  *rank = 0;
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::grid_coord(int *gcoord, const int rank)
{
  for (int i = 0; i < m_Ndim; ++i) {
    gcoord[i] = 0;
  }
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::grid_dims(int *gdims)
{
  for (int i = 0; i < m_Ndim; ++i) {
    gdims[i] = 1;
  }
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::sync()
{
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::sync_global()
{
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::Base::broadcast(size_t size, void *data, int sender)
{
  // stay intact
  return EXIT_SUCCESS;
}


int Communicator::broadcast(int count, double *data, int sender)
{
  // stay intact
  return EXIT_SUCCESS;
}


int Communicator::broadcast(int count, int *data, int sender)
{
  // stay intact
  return EXIT_SUCCESS;
}


int Communicator::broadcast(int count, string& data, int sender)
{
  // stay intact
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::Base::exchange(size_t size, void *recv_buf, void *send_buf, int idir, int ipm, int tag)
{
  memcpy(recv_buf, send_buf, size);
  return EXIT_SUCCESS;
}


int Communicator::exchange(int count, double *recv_buf, double *send_buf, int idir, int ipm, int tag)
{
  memcpy(recv_buf, send_buf, sizeof(double) * count);
  return EXIT_SUCCESS;
}


int Communicator::exchange(int count, int *recv_buf, int *send_buf, int idir, int ipm, int tag)
{
  memcpy(recv_buf, send_buf, sizeof(int) * count);
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::Base::send_1to1(size_t size, void *recv_buf, void *send_buf, int send_to, int recv_from, int tag)
{
  memcpy(recv_buf, send_buf, size);
  return EXIT_SUCCESS;
}


int Communicator::send_1to1(int count, double *recv_buf, double *send_buf, int send_to, int recv_from, int tag)
{
  memcpy(recv_buf, send_buf, sizeof(double) * count);
  return EXIT_SUCCESS;
}


int Communicator::send_1to1(int count, int *recv_buf, int *send_buf, int send_to, int recv_from, int tag)
{
  memcpy(recv_buf, send_buf, sizeof(int) * count);
  return EXIT_SUCCESS;
}


//====================================================================
int Communicator::reduce_sum(int count, double *recv_buf, double *send_buf, int pattern)
{
  memcpy(recv_buf, send_buf, sizeof(double) * count);
  return EXIT_SUCCESS;
}


int Communicator::reduce_sum(int count, int *recv_buf, int *send_buf, int pattern)
{
  memcpy(recv_buf, send_buf, sizeof(int) * count);
  return EXIT_SUCCESS;
}


double Communicator::reduce_sum(double v)
{
  return v;
}


//====================================================================
double Communicator::reduce_max(double v)
{
  return v;
}


//====================================================================
double Communicator::reduce_min(double v)
{
  return v;
}


//====================================================================
int Communicator::status()
{
#ifdef DEBUG
  printf("Communicator Single\n");
#endif
  return EXIT_SUCCESS;
}


//====================================================================
double Communicator::get_time()
{
  struct timeval now;

  if (gettimeofday(&now, (struct timezone *)0) != 0) {
    return double();
  }

  double sec = (double)now.tv_sec + ((double)now.tv_usec) * 1.0e-6;

  return sec;
}


//====================================================================
//============================================================END=====
#endif //COMMUNICATOR_USE_SINGLE