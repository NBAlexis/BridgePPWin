/*!
        @file    gaugeConfig.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/


#ifndef GAUGECONFIG_INCLUDED
#define GAUGECONFIG_INCLUDED

#include <string>
using std::string;

#include "Field/field_G.h"

#include "fieldIO_Text.h"
#include "fieldIO_Text_4x4x4x8.h"
#include "fieldIO_Binary.h"
#include "fieldIO_Binary_Parallel.h"
#include "fieldIO_Binary_Distributed.h"
#include "fieldIO_Fortran.h"
#include "fieldIO_LIME.h"
#include "fieldIO_LIME_Parallel.h"
#include "fieldIO_None.h"

#include "io_format_gauge.h"

#include "bridgeIO.h"
using Bridge::vout;

//! GaugeConfig class BAPI for file I/O of gauge configuration.

/*!
   This class BAPI family is used to setup and output the gauge
   configuration.
   This is the base class BAPI, which implements common functions,
   and read/write methods are implemented in subclasses.
   At present, cutting off the gauge field for each node
   and deliver it to the node is implemented in this class BAPI.
   It may be better to separate that to other class BAPI for
   general usage for other field objects.
                                [28 Dec 2011 H.Matsufuru]

   GaugeConfig class BAPI provides file I/O of gauge configuration.
   It provides an interface to underlying FieldIO class BAPI family;

   The file format is specified by a string argument element_type to
   the constructor (in a somewhat similar manner as factory).
   Data layout is ILDG layout for most of the cases except
   for Fortran_JLQCD in which JLQCD layout is applied.

   "NO_OUTPUT" is added to GaugeConfig::write_file()
                                [22 Feb 2015 Y.Namekawa]
   unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
   Config types "Unit" and "Random" are introduced to generate
   unit and random gauge configurations for cold and hot start,
   respectively.
   The methods are renamed to read/write. read_file/write_file
   are left for compatibility. these methods now take args of
   Field_G* element_type instead of Field*.
                                [24 June 2016 T.Aoyama]
   "Null" config element_type is introduced that suppresses output.
   Specifying "NO_OUTPUT" for filename is also valid.
                                [8 July 2016 T.Aoyama]
   "Null" renamed by "None" in compliance with
   yaml specifications. (kept for backward compatilibity)
                                [21 November 2018 T.Aoyama]
 */

class BAPI GaugeConfig
{
 public:
  static const std::string class_name;

 public:
  GaugeConfig(const string& type);
  virtual ~GaugeConfig();

 private:
  // non-copyable
  GaugeConfig(const GaugeConfig&);
  GaugeConfig& operator=(const GaugeConfig&);

 public:
  void read(Field_G *U, const string& filename = string());

  void read(unique_ptr<Field_G>& U, const string& filename = string())
  { return read(U.get(), filename); }

  void write(Field_G *U, const string& filename = string());

  void write(unique_ptr<Field_G>& U, const string& filename = string())
  { return write(U.get(), filename); }


  void read_file(Field_G *U, const string& filename)
  { return read(U, filename); }
  void read_file(unique_ptr<Field_G>& U, const string& filename)
  { return read_file(U.get(), filename); }

  void write_file(Field_G *U, const string& filename)
  { return write(U, filename); }
  void write_file(unique_ptr<Field_G>& U, const string& filename)
  { return write_file(U.get(), filename); }

 protected:
  Bridge::VerboseLevel m_vl;
  FieldIO *m_fieldio;

 private:
  class BAPI DataSource;
  class BAPI DataSource_Unit;
  class BAPI DataSource_Random;

  DataSource *m_datasource;
};
#endif
