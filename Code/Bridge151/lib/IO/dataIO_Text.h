/*!
        @file    dataIO_Text.h

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-07-22 15:50:06 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef DATAIO_TEXT_INCLUDED
#define DATAIO_TEXT_INCLUDED

#include "dataIO.h"
#include <string>
using std::string;

#include <fstream>
using std::fstream;

enum
{
  default_format_precision = 14,
};

//! DataIO_Text class BAPI for general file I/O in plain Text format.

/**
   DataIO_Text class BAPI provides file I/O in plain Text format.
   output format of numerical precision can be specified through
   set_parameter method. default is given as enum.

   Data are read from/written to file through the primary node
   (rank 0 usually).

   The read/write methods in various data formats are all derived from
   template functions whose implementation given in dataIO_Text_impl.h
   included at the bottom.
 */

class BAPI DataIO_Text : public DataIO
{
 public:
  DataIO_Text() : m_format_precision(default_format_precision) {}
  ~DataIO_Text() {}

  virtual void read_file(double *v, const size_t n, const string& f)
  { return read_file_base(v, n, f); }
  virtual void write_file(const double *v, const size_t n, const string& f, const bool append = true)
  { return write_file_base(v, n, f, append); }

  virtual void read_file(dcomplex *v, const size_t n, const string& f)
  { return read_file_base(v, n, f); }
  virtual void write_file(const dcomplex *v, const size_t n, const string& f, const bool append = true)
  { return write_file_base(v, n, f, append); }

  virtual void read_file(std::vector<double>& v, const string& f)
  { return read_file_base(&v[0], v.size(), f); }
  virtual void write_file(const std::vector<double>& v, const string& f, const bool append = true)
//    { return write_file_base(&v[0], v.size(), f, append); }
  {    // workaround for standard-compliant definition of operator[]
    std::vector<double>& p = const_cast<std::vector<double>&>(v);

    return write_file_base(&p[0], v.size(), f, append);
  }

  virtual void read_file(std::vector<dcomplex>& v, const string& f)
  { return read_file_base(&v[0], v.size(), f); }
  virtual void write_file(const std::vector<dcomplex>& v, const string& f, const bool append = true)
//    { return write_file_base(&v[0], v.size(), f, append); }
  {
    std::vector<dcomplex>& p = const_cast<std::vector<dcomplex>&>(v);

    return write_file_base(&p[0], v.size(), f, append);
  }

  void set_parameter(const int precision) { m_format_precision = precision; }

 private:
  int m_format_precision;

  template<typename T>
  void read_file_base(T *v, const size_t n, const string&);

  template<typename T>
  void write_file_base(const T *v, const size_t n, const string&, const bool append = true);

  template<typename T>
  inline void fetch_data(T& v, fstream& fs);

  template<typename T>
  inline void store_data(const T& v, fstream& fs);
};

#include "dataIO_Text_impl.h"
#endif /* DATAIO_TEXT_INCLUDED */
