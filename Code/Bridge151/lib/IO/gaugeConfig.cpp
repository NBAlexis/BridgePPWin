/*!
        @file    gaugeConfig.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-02-13 20:10:32 #$

        @version $LastChangedRevision: 1942 $
*/
#include "BridgeLib_Private.h"
#include "gaugeConfig.h"
#include "Tools/randomNumberManager.h"

const std::string GaugeConfig::class_name = "GaugeConfig";

//====================================================================
class BAPI GaugeConfig::DataSource
{
 public:
  virtual ~DataSource() {}
  virtual void set(Field_G *U) = 0;
};

//--------------------------------------------------------------------
class BAPI GaugeConfig::DataSource_Unit : public DataSource
{
 public:
  virtual void set(Field_G *U)
  {
    if (U) U->set_unit();
  }
};

//--------------------------------------------------------------------
class BAPI GaugeConfig::DataSource_Random : public DataSource
{
 public:
  DataSource_Random()
  {
    m_rand = RandomNumberManager::getInstance();
  }

  ~DataSource_Random()
  {
  }

  virtual void set(Field_G *U)
  {
    if (U) U->set_random(m_rand);
  }

 private:
  RandomNumbers *m_rand;
};

//====================================================================
GaugeConfig::GaugeConfig(const string& type)
  : m_vl(CommonParameters::Vlevel()),
  m_fieldio(),
  m_datasource()
{
  if (type == "Text") {
    m_fieldio = new FieldIO_Text(IO_Format::Gauge::ILDG);
  } else if (type == "Text_4x4x4x8") {
    m_fieldio = new FieldIO_Text_4x4x4x8(IO_Format::Gauge::ILDG);
  } else if (type == "Binary") {
    m_fieldio = new FieldIO_Binary(IO_Format::Gauge::ILDG);
  } else if (type == "Fortran_JLQCD") {
    m_fieldio = new FieldIO_Fortran(IO_Format::Gauge::JLQCD);
  } else if (type == "Fortran_ILDG") {
    m_fieldio = new FieldIO_Fortran(IO_Format::Gauge::ILDG);
  } else if (type == "ILDG") {
    m_fieldio = new FieldIO_LIME(IO_Format::Gauge::ILDG);
  } else if (type == "Binary_Parallel") {
    m_fieldio = new FieldIO_Binary_Parallel(IO_Format::Gauge::ILDG);
  } else if (type == "Binary_Distributed") {
    m_fieldio = new FieldIO_Binary_Distributed(IO_Format::Gauge::ILDG);
  } else if (type == "ILDG_Parallel") {
    m_fieldio = new FieldIO_LIME_Parallel(IO_Format::Gauge::ILDG);
  } else if ((type == "None") || (type == "Null")) {  // for backward compatibility
    m_fieldio = new FieldIO_None(IO_Format::Trivial);
  } else if (type == "Unit") {
    m_datasource = new DataSource_Unit;
  } else if (type == "Random") {
    m_datasource = new DataSource_Random;
  } else {
    vout.crucial("Error at %s: unsupported element_type \"%s\".\n", class_name.c_str(), type.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
GaugeConfig::~GaugeConfig()
{
  if (m_fieldio) delete m_fieldio;
  if (m_datasource) delete m_datasource;
}


//====================================================================
void GaugeConfig::read(Field_G *U, const string& filename)
{
  if (m_fieldio) {
    return m_fieldio->read_file(U, filename);
  }

  if (m_datasource) {
    return m_datasource->set(U);
  }

  vout.crucial("Error at %s::read_file(): FieldIO or DataSource not set.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void GaugeConfig::write(Field_G *U, const string& filename)
{
  if (m_fieldio) {
    if (filename != "NO_OUTPUT") {
      return m_fieldio->write_file(U, filename);
    } else {
      return;  // do nothig
    }
  }

  vout.crucial("Error at %s::write_file(): FieldIO not set.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
//============================================================END=====
