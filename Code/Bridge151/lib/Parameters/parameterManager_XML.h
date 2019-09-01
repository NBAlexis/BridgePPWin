/*!
        @file    parameterManager_XML.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
 */

#ifndef PARAMETERMANAGER_XML_INCLUDED
#define PARAMETERMANAGER_XML_INCLUDED

//#include "parameterManager.h"
//#include <string>

//! Parameter manager with YAML parser.

/*!
   This is a simple parser to read parameters from a file
   prepared with YAML format.
   Only simple cases were checked.
                                      [17 Jul 2012 H.Matsufuru]

   read and set parameters from XML file, using tinyxml-2 parser.
   [16 Mar 2015 T.Aoyama]

 */
class BAPI ParameterManager_XML : public ParameterManager
{
 public:
  static const std::string class_name;

  ParameterManager_XML() {}

  //! read parameters from file.
  void read_params(const std::string& params_file, Parameters& params);
};
#endif
