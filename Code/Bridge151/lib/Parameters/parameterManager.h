/*!
        @file    parameterManager.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
 */

#ifndef PARAMETERMANAGER_INCLUDED
#define PARAMETERMANAGER_INCLUDED

//#include "configure.h"
//#include "bridge_defs.h"
//#include "parameters.h"
//#include "commonParameters.h"

//! Base class BAPI of parameter manager.

/*!
                          [17 Jun 2012 H.Matsufuru]
 */

class BAPI ParameterManager
{
 public:
  static const std::string class_name;

 protected:

  Bridge::VerboseLevel m_vl;

 public:

  ParameterManager() : m_vl(CommonParameters::Vlevel()) {}

  virtual ~ParameterManager() {}

 private:
  // non-copyable
  ParameterManager(const ParameterManager&);
  ParameterManager& operator=(const ParameterManager&);

 public:

  virtual void
  read_params(const std::string& params_file, Parameters& params) = 0;

  static
  void read(const std::string& params_file, Parameters& params);

  static
  Parameters read(const std::string& params_file);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }
};
#endif
