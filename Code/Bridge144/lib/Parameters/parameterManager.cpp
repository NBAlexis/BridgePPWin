#include "BridgeLib_Private.h"

/*!
        @file    $Id:: parameterManager.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "parameterManager_YAML.h"
#include "parameterManager_XML.h"

const std::string ParameterManager::class_name = "ParameterManager";

//====================================================================
Parameters ParameterManager::read(const std::string& params_file)
{
  Parameters params;

  read(params_file, params);
  return params;
}


//====================================================================
void ParameterManager::read(const std::string& params_file, Parameters& params)
{
  if (params_file.size() == 0) return;

  std::string ext = params_file.substr(params_file.find_last_of('.'));

  vout.paranoiac("ext = %s\n", ext.c_str());

  if (ext == ".yaml") {
    return ParameterManager_YAML().read_params(params_file, params);
  } else if (ext == ".xml") {
    return ParameterManager_XML().read_params(params_file, params);
  } else {
    vout.crucial("Error at %s: unrecognized file type: %s\n", class_name.c_str(), params_file.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
//============================================================END=====
