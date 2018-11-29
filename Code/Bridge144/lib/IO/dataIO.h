/*!
        @file    $Id: dataIO.h #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2013-07-22 15:50:06 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef DATAIO_INCLUDED
#define DATAIO_INCLUDED

#include "bridge_complex.h"

#include "Parameters/commonParameters.h"
#include "bridgeIO.h"
using Bridge::vout;

#include <string>
using std::string;


//! DataIO class for file I/O of general collection of data.

/**
   DataIO class provides abstract base class for file I/O of
   general collection of data that do not have space-time index.

   interfaces defined for reading and writing array of double and complex
   of size n, and std::vector of double and complex.

 */

class BAPI DataIO
{
 public:
  DataIO() : m_vl(CommonParameters::Vlevel()) {}
  virtual ~DataIO() {}

 private:
  // non-copyable
  DataIO(const DataIO&);
  DataIO& operator=(const DataIO&);

 public:

  virtual void read_file(double *v, const size_t n, const string&) = 0;
  virtual void write_file(const double *v, const size_t n, const string&, const bool append = true) = 0;

  virtual void read_file(dcomplex *v, const size_t n, const string&) = 0;
  virtual void write_file(const dcomplex *v, const size_t n, const string&, const bool append = true) = 0;

  virtual void read_file(std::vector<double>&, const string&) = 0;
  virtual void write_file(const std::vector<double>&, const string&, const bool append = true) = 0;

  virtual void read_file(std::vector<dcomplex>&, const string&) = 0;
  virtual void write_file(const std::vector<dcomplex>&, const string&, const bool append = true) = 0;

 protected:
  Bridge::VerboseLevel m_vl;
};
#endif /* DATAIO_INCLUDED */
