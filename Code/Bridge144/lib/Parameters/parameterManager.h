/*!
        @file    $Id:: parameterManager.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
 */

#ifndef PARAMETERMANAGER_INCLUDED
#define PARAMETERMANAGER_INCLUDED

#include "configure.h"
#include "defs.h"
#include "parameters.h"
#include "commonParameters.h"

//! Base class of parameter manager.

/*!
                          [17 Jun 2012 H.Matsufuru]
 */

class ParameterManager
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
