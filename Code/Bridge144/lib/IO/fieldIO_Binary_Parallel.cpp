/*!
        @file    $Id: fieldIO_Binary_Parallel.cpp #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-01-22 13:51:53 #$

        @version $LastChangedRevision: 1646 $
*/

// this code only makes sense in MPI environment.
#ifdef USE_MPI

#include "fieldIO_Binary_Parallel.h"

const std::string FieldIO_Binary_Parallel::class_name = "FieldIO_Binary_Parallel";

//====================================================================
FieldIO_Binary_Parallel::FieldIO_Binary_Parallel(const IO_Format::Format *format)
  : FieldIO(format),
    m_is_initialized(false),
    m_nvol(0), m_nin_file(0), m_nex_file(0)
{}

//====================================================================
FieldIO_Binary_Parallel::~FieldIO_Binary_Parallel()
{
  finalize();
}

//====================================================================
void FieldIO_Binary_Parallel::read_file(Field *u, string filename)
{
  static const char _function_name[] = "FieldIO_Binary_Parallel::read_file";

  if (!u) {
    vout.crucial(m_vl, "Warning at %s: field is null.\n", _function_name);
    return;
  }

  initialize(u);

  MPI_File fh;
  int      ret;

  // fetch data from file into buffer: buffer is in file order.
  double *buf = new double [m_nin_file * m_nvol * m_nex_file];
  if (!buf) {
    vout.crucial(m_vl, "Error at %s: allocate buffer failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

#ifdef USE_BGNET
  ret = MPI_File_open(MPI_COMM_WORLD, const_cast<char *>(filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
#else
  ret = MPI_File_open(Communicator_impl::world(), const_cast<char *>(filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
#endif
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_open failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_File_set_view(fh, 0, m_type_vector, m_type_tiled, const_cast<char *>("native"), MPI_INFO_NULL);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_set_view failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_File_read_all(fh, (void *)buf, m_nvol * m_nex_file, m_type_vector, MPI_STATUS_IGNORE);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_read_all failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_File_close(&fh);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_close failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  if (!is_bigendian()) {
    convert_endian(buf, sizeof(double), m_nvol * m_nin_file * m_nex_file);
  }

  // unpack buffer
  double *p = buf;

  for (int j = 0; j < m_nex_file; ++j) {
    for (int isite = 0; isite < m_nvol; ++isite) {
      for (int i = 0; i < m_nin_file; ++i) {
        int s, t;
        m_format->file_to_field(s, t, i, j);

        u->set(s, isite, t, *p++);
      }
    }
  }

  delete [] buf;

  finalize();
}


//====================================================================
void FieldIO_Binary_Parallel::write_file(Field *u, string filename)
{
  static const char _function_name[] = "FieldIO_Binary_Parallel::write_file";

  if (!u) {
    vout.crucial(m_vl, "Warning at %s: field is null.\n", _function_name);
  }

  initialize(u);

  double *buf = new double [m_nin_file * m_nvol * m_nex_file];
  if (!buf) {
    vout.crucial(m_vl, "Error at %s: allocate buffer failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  // pack buffer
  double *p = buf;

  for (int j = 0; j < m_nex_file; ++j) {
    for (int isite = 0; isite < m_nvol; ++isite) {
      for (int i = 0; i < m_nin_file; ++i) {
        int s, t;
        m_format->file_to_field(s, t, i, j);

        *p++ = u->cmp(s, isite, t);
      }
    }
  }

  if (!is_bigendian()) {
    convert_endian(buf, sizeof(double), m_nin_file * m_nvol * m_nex_file);
  }

  MPI_File fh;
  int      ret;

#ifdef USE_BGNET
  ret = MPI_File_open(MPI_COMM_WORLD, const_cast<char *>(filename.c_str()), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
#else
  ret = MPI_File_open(Communicator_impl::world(), const_cast<char *>(filename.c_str()), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
#endif
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_open failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_File_set_view(fh, 0, m_type_vector, m_type_tiled, const_cast<char *>("native"), MPI_INFO_NULL);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_set_view failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_File_write_all(fh, (void *)buf, m_nvol * m_nex_file, m_type_vector, MPI_STATUS_IGNORE);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_write_all failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_File_close(&fh);
  if (ret) {
    vout.crucial(m_vl, "Error at %s: MPI_File_close failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  delete [] buf;

  finalize();
}


//====================================================================
int FieldIO_Binary_Parallel::initialize(const Field* v)
{
  static const char _function_name[] = "FieldIO_Binary_Parallel::initialize";

  // N.B. trivial layout returns 0 for nin/nex_file.
  int nin_file = m_format->nin() ? m_format->nin() : v->nin();
  int nex_file = m_format->nex() ? m_format->nex() : v->nex();
  int nvol = v->nvol();
  
  // check if layout is already generated and recyclable.
  if (m_is_initialized
      && nin_file == m_nin_file
      && nex_file == m_nex_file
      && nvol == m_nvol)
  {
    vout.detailed(m_vl, "%s: layout recycled.\n", _function_name);
    return EXIT_SUCCESS;
  }
  
  // first, cleanup pre-existing layout, if any.
  clear_layout();

  // local parameters
  m_nin_file = nin_file;
  m_nex_file = nex_file;

  int ndim = CommonParameters::Ndim();

  int *global_dims = new int[ndim];
  global_dims[0] = CommonParameters::Lx();
  global_dims[1] = CommonParameters::Ly();
  global_dims[2] = CommonParameters::Lz();
  global_dims[3] = CommonParameters::Lt();

  int *local_dims = new int[ndim];
  local_dims[0] = CommonParameters::Nx();
  local_dims[1] = CommonParameters::Ny();
  local_dims[2] = CommonParameters::Nz();
  local_dims[3] = CommonParameters::Nt();

  m_nvol = 1;
  for (int i = 0; i < ndim; ++i) {
    m_nvol *= local_dims[i];
  }

  int *grid_pos = new int[ndim];
  for (int i = 0; i < ndim; ++i) {
    grid_pos[i] = Communicator::ipe(i);
  }

  int *starts = new int[ndim];
  for (int i = 0; i < ndim; ++i) {
    starts[i] = local_dims[i] * grid_pos[i];
  }

  int ret = 0;

  // MPI_Datatype m_type_vector;
  ret = MPI_Type_contiguous(sizeof(double) * m_nin_file, MPI_BYTE, &m_type_vector);
  if (ret) {
    vout.crucial(m_vl, "%s: MPI_Type_Contiguous failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_Type_commit(&m_type_vector);
  if (ret) {
    vout.crucial(m_vl, "%s: MPI_Type_commit failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  //  MPI_Datatype m_type_tiled;
  ret = MPI_Type_create_subarray(ndim, global_dims, local_dims, starts, MPI_ORDER_FORTRAN, m_type_vector, &m_type_tiled);
  if (ret) {
    vout.crucial(m_vl, "%s: MPI_Type_create_subarray failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  ret = MPI_Type_commit(&m_type_tiled);
  if (ret) {
    vout.crucial(m_vl, "%s: MPI_Type_commit failed.\n", _function_name);
    exit(EXIT_FAILURE);
  }

  delete [] starts;
  delete [] grid_pos;
  delete [] local_dims;
  delete [] global_dims;

  // initialization done.
  m_is_initialized = true;

  vout.detailed(m_vl, "%s: layout initialized.\n", _function_name);

  return EXIT_SUCCESS;
}

//====================================================================
int FieldIO_Binary_Parallel::clear_layout()
{
  const char _function_name[] = "FieldIO_Binary_Parallel::clear_layout";

  if (m_is_initialized) {
    m_nin_file = 0;
    m_nex_file = 0;
    m_nvol = 0;

    int ret = 0;
    
    ret = MPI_Type_free(&m_type_vector);
    if (ret) {
      vout.crucial(m_vl, "%s: MPI_Type_free for type_vector failed.\n", _function_name);
      exit(EXIT_FAILURE);
    }

    ret = MPI_Type_free(&m_type_tiled);
    if (ret) {
      vout.crucial(m_vl, "%s: MPI_Type_free for type_tiled failed.\n", _function_name);
      exit(EXIT_FAILURE);
    }

    m_is_initialized = false;
  }

  return EXIT_SUCCESS;
}


//====================================================================
int FieldIO_Binary_Parallel::finalize()
{
  static const char _function_name[] = "FieldIO_Binary_Parallel::finalize";

  clear_layout();

  vout.detailed(m_vl, "%s via MPI I/O finalize done.\n", class_name.c_str());

  return EXIT_SUCCESS;
}


//====================================================================
#endif

//====================================================================
//============================================================END=====
