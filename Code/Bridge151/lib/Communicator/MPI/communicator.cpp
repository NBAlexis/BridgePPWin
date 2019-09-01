/*!
        @file    communicator.cpp

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
 */
#include "BridgeLib_Private.h"
//#include "layout.h"

//====================================================================

/**
   MPI implementation of Communicator class BAPI
   which delegates to Communicator_impl class BAPI.
 */
int Communicator::init(int *pargc, char ***pargv)
{
  return Communicator_impl::init(pargc, pargv);
}


//====================================================================
int Communicator::finalize()
{
  return Communicator_impl::finalize();
}


//====================================================================
void Communicator::abort()
{
  return Communicator_impl::abort();
}


//====================================================================
int Communicator::setup(int ninstance)
{
  return Communicator_impl::setup(ninstance);
}


//====================================================================
int Communicator::setup(
  const std::vector<int>& lattice_size,
  std::vector<int>& grid_size,
  int ninstance)
{
  return Communicator_impl::setup(lattice_size, grid_size, ninstance);
}


//====================================================================
bool Communicator::is_primary()
{
  return Communicator_impl::is_primary();
}


//====================================================================
bool Communicator::is_primary_master()
{
  return Communicator_impl::is_primary_master();
}


//====================================================================
int Communicator::self()
{
  return Communicator_impl::self();
}


//====================================================================
int Communicator::size()
{
  return Communicator_impl::size();
}


//====================================================================
#ifdef ENABLE_MULTI_INSTANCE
int Communicator::self_global()
{
  return Communicator_impl::self_global();
}


//====================================================================
int Communicator::world_id()
{
  return Communicator_impl::world_id();
}


#endif

//====================================================================
int Communicator::ipe(const int idir)
{
  return Communicator_impl::Layout::ipe(idir);
}


//====================================================================
int Communicator::npe(const int idir)
{
  return Communicator_impl::Layout::npe(idir);
}


//====================================================================
int Communicator::grid_rank(int *rank, const int *grid_coord)
{
  return Communicator_impl::Layout::grid_rank(rank, grid_coord);
}


//====================================================================
int Communicator::grid_coord(int *grid_coord, const int rank)
{
  return Communicator_impl::Layout::grid_coord(grid_coord, rank);
}


//====================================================================
int Communicator::grid_dims(int *grid_dims)
{
  return Communicator_impl::Layout::grid_dims(grid_dims);
}


//====================================================================
int Communicator::sync()
{
  return Communicator_impl::sync();
}


//====================================================================
#ifdef ENABLE_MULTI_INSTANCE
int Communicator::sync_global()
{
  return Communicator_impl::sync_global();
}


#endif

//====================================================================
int Communicator::Base::broadcast(size_t size, void *data, int sender)
{
  return Communicator_impl::Base::broadcast(size, data, sender);
}


int Communicator::broadcast(int count, double *data, int sender)
{
  return Communicator_impl::Base::broadcast(sizeof(double) * count, (void *)data, sender);
}


int Communicator::broadcast(int count, float *data, int sender)
{
  return Communicator_impl::Base::broadcast(sizeof(float) * count, (void *)data, sender);
}


int Communicator::broadcast(int count, int *data, int sender)
{
  return Communicator_impl::Base::broadcast(sizeof(int) * count, (void *)data, sender);
}


int Communicator::broadcast(int count, string& data, int sender)
{
  return Communicator_impl::broadcast_string(count, data, sender);
}


//====================================================================
int Communicator::Base::exchange(size_t size, void *recv_buf, void *send_buf, int idir, int ipm, int itag)
{
  return Communicator_impl::Base::exchange(size, recv_buf, send_buf, idir, ipm, itag);
}


int Communicator::exchange(int count, double *recv_buf, double *send_buf, int idir, int ipm, int itag)
{
  return Communicator_impl::Base::exchange(sizeof(double) * count, (void *)recv_buf, (void *)send_buf, idir, ipm, itag);
}


int Communicator::exchange(int count, float *recv_buf, float *send_buf, int idir, int ipm, int itag)
{
  return Communicator_impl::Base::exchange(sizeof(float) * count, (void *)recv_buf, (void *)send_buf, idir, ipm, itag);
}


int Communicator::exchange(int count, int *recv_buf, int *send_buf, int idir, int ipm, int itag)
{
  return Communicator_impl::Base::exchange(sizeof(int) * count, (void *)recv_buf, (void *)send_buf, idir, ipm, itag);
}


//====================================================================
int Communicator::Base::send_1to1(size_t size, void *recv_buf, void *send_buf, int send_to, int recv_from, int tag)
{
  return Communicator_impl::Base::send_1to1(size, recv_buf, send_buf, send_to, recv_from, tag);
}


int Communicator::send_1to1(int count, double *recv_buf, double *send_buf, int send_to, int recv_from, int tag)
{
  return Communicator_impl::Base::send_1to1(sizeof(double) * count, (void *)recv_buf, (void *)send_buf, send_to, recv_from, tag);
}


int Communicator::send_1to1(int count, float *recv_buf, float *send_buf, int send_to, int recv_from, int tag)
{
  return Communicator_impl::Base::send_1to1(sizeof(float) * count, (void *)recv_buf, (void *)send_buf, send_to, recv_from, tag);
}


int Communicator::send_1to1(int count, int *recv_buf, int *send_buf, int send_to, int recv_from, int tag)
{
  return Communicator_impl::Base::send_1to1(sizeof(int) * count, (void *)recv_buf, (void *)send_buf, send_to, recv_from, tag);
}


//====================================================================
int Communicator::reduce_sum(int count, double *recv_buf, double *send_buf, int pattern)
{
  return Communicator_impl::Base::reduce(count, (void *)recv_buf, (void *)send_buf, MPI_DOUBLE, MPI_SUM, pattern);
}


int Communicator::reduce_sum(int count, float *recv_buf, float *send_buf, int pattern)
{
  return Communicator_impl::Base::reduce(count, (void *)recv_buf, (void *)send_buf, MPI_FLOAT, MPI_SUM, pattern);
}


int Communicator::reduce_sum(int count, int *recv_buf, int *send_buf, int pattern)
{
  return Communicator_impl::Base::reduce(count, (void *)recv_buf, (void *)send_buf, MPI_INT, MPI_SUM, pattern);
}


//====================================================================
int Communicator::reduce_max(int count, double *recv_buf, double *send_buf, int pattern)
{
  return Communicator_impl::Base::reduce(count, (void *)recv_buf, (void *)send_buf, MPI_DOUBLE, MPI_MAX, pattern);
}


int Communicator::reduce_max(int count, float *recv_buf, float *send_buf, int pattern)
{
  return Communicator_impl::Base::reduce(count, (void *)recv_buf, (void *)send_buf, MPI_FLOAT, MPI_MAX, pattern);
}


int Communicator::reduce_max(int count, int *recv_buf, int *send_buf, int pattern)
{
  return Communicator_impl::Base::reduce(count, (void *)recv_buf, (void *)send_buf, MPI_INT, MPI_MAX, pattern);
}


//====================================================================
int Communicator::reduce_min(int count, double *recv_buf, double *send_buf, int pattern)
{
  return Communicator_impl::Base::reduce(count, (void *)recv_buf, (void *)send_buf, MPI_DOUBLE, MPI_MIN, pattern);
}


int Communicator::reduce_min(int count, float *recv_buf, float *send_buf, int pattern)
{
  return Communicator_impl::Base::reduce(count, (void *)recv_buf, (void *)send_buf, MPI_FLOAT, MPI_MIN, pattern);
}


int Communicator::reduce_min(int count, int *recv_buf, int *send_buf, int pattern)
{
  return Communicator_impl::Base::reduce(count, (void *)recv_buf, (void *)send_buf, MPI_INT, MPI_MIN, pattern);
}


//====================================================================
double Communicator::reduce_sum(double a)
{
  double ar = double();

  reduce_sum(1, &ar, &a, 0);
  return ar;
}


//====================================================================
double Communicator::reduce_max(double a)
{
  double ar = double();

  reduce_max(1, &ar, &a, 0);
  return ar;
}


//====================================================================
double Communicator::reduce_min(double a)
{
  double ar = double();

  reduce_min(1, &ar, &a, 0);
  return ar;
}


//====================================================================
float Communicator::reduce_sum(float a)
{
  float ar = float();

  reduce_sum(1, &ar, &a, 0);
  return ar;
}


//====================================================================
float Communicator::reduce_max(float a)
{
  float ar = float();

  reduce_max(1, &ar, &a, 0);
  return ar;
}


//====================================================================
float Communicator::reduce_min(float a)
{
  float ar = float();

  reduce_min(1, &ar, &a, 0);
  return ar;
}


//====================================================================
double Communicator::get_time()
{
  return Communicator_impl::get_time();
}


//====================================================================
int Communicator::status()
{
  return Communicator_impl::status();
}


//====================================================================
Channel *Communicator::send_init(int count, int idir, int ipm)
{
  return Communicator_impl::send_init(count, idir, ipm);
}


//====================================================================
Channel *Communicator::recv_init(int count, int idir, int ipm)
{
  return Communicator_impl::recv_init(count, idir, ipm);
}


//====================================================================
//============================================================END=====
