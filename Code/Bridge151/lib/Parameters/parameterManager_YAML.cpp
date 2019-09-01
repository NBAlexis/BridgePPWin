/*!
        @file    parameterManager_YAML.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
 */
#include "BridgeLib_Private.h"

//#include "parameterManager_YAML.h"
//
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <cstring>
//#include <map>
//#include <vector>
//#include <stack>
//
//#ifdef USE_YAMLCPPLIB
//#include <yaml-cpp/yaml.h>
//#endif
//
//#include "IO/bridgeIO.h"

using Bridge::vout;

using std::string;
using std::vector;
using std::map;
using std::pair;

const std::string ParameterManager_YAML::class_name = "ParameterManager_YAML";

//====================================================================
#ifdef USE_YAMLCPPLIB

//====================================================================
class BAPI Parser_yamlcpp
{
 public:
  static const string class_name;

 public:
  void parse(std::istream& fin, Parameters& params);

  Parser_yamlcpp()
    : m_vl(CommonParameters::Vlevel()) {}

 private:

  int store_map_data(const YAML::Node&, Parameters&);

  int store_vector_data(const YAML::Node&, vector<string>&);

  string type_string(const YAML::Node&) const;
  int traverse(const YAML::Node&, const string& indent = "") const;

  Bridge::VerboseLevel m_vl;
};

//--------------------------------------------------------------------
const string Parser_yamlcpp::class_name = "ParameterManager_YAML::Parser_yamlcpp";

//--------------------------------------------------------------------
void Parser_yamlcpp::parse(std::istream& fin, Parameters& params)
{
  YAML::Node node = YAML::Load(fin);

  int result = store_map_data(node, params);

  if (result != EXIT_SUCCESS) {
    vout.crucial(m_vl, "%s: parse failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//--------------------------------------------------------------------
int Parser_yamlcpp::store_map_data(const YAML::Node& node, Parameters& params)
{
  if (!node.IsMap()) {
    vout.general(m_vl, "Error: map expected. found %s\n", type_string(node).c_str());
    return EXIT_FAILURE;
  }

  int retv = EXIT_SUCCESS;

  for (YAML::const_iterator p = node.begin(); p != node.end(); ++p) {
    string key = p->first.as<string>();

    if (p->second.IsScalar()) {
      string value = p->second.as<string>();

      params.set_string(key, value);
    } else if (p->second.IsSequence()) {
      vector<string> v;

      if (store_vector_data(p->second, v) != EXIT_SUCCESS) {
        retv = EXIT_FAILURE;
      }

      params.set_string_vector(key, v);
    } else if (p->second.IsMap()) {
      Parameters pp;
      store_map_data(p->second, pp);

      params.set_Parameters(key, pp);
    } else {
      vout.general(m_vl, "Error: unexpected type %s\n", type_string(p->second).c_str());
      retv = EXIT_FAILURE;
    }
  }

  return retv;
}


//--------------------------------------------------------------------
int Parser_yamlcpp::store_vector_data(const YAML::Node& node, vector<string>& v)
{
  int retv = EXIT_SUCCESS;

  if (!node.IsSequence()) {
    vout.general(m_vl, "Error: vector expected. found %s\n", type_string(node).c_str());
    return EXIT_FAILURE;
  }

  for (size_t i = 0; i < node.size(); ++i) {
    if (!node[i].IsScalar()) {
      vout.general(m_vl, "Error: item[%zu] scalar expected. found %s\n", i, type_string(node).c_str());
      retv = EXIT_FAILURE;
    } else {
      v.push_back(node[i].as<string>());
    }
  }

  return retv;
}


//--------------------------------------------------------------------
string Parser_yamlcpp::type_string(const YAML::Node& node) const
{
  switch (node.Type())
  {
  case YAML::NodeType::Undefined:
    return "UndefinedType";

  case YAML::NodeType::Null:
    return "NullType";

  case YAML::NodeType::Scalar:
    return "Scalar";

  case YAML::NodeType::Sequence:
    return "Sequence";

  case YAML::NodeType::Map:
    return "Map";

  default:
    return "(Unknown)";
  }
  return "(ERROR)";
}


//--------------------------------------------------------------------
int Parser_yamlcpp::traverse(const YAML::Node& node, const string& indent) const
{
  switch (node.Type())
  {
  case YAML::NodeType::Null:
    vout.general("%sNull type\n", indent.c_str());
    break;

  case YAML::NodeType::Scalar:
    vout.general("%sScalar type \"%s\"\n", indent.c_str(), node.as<string>().c_str());
    break;

  case YAML::NodeType::Sequence:
    vout.general("%sSequence type\n", indent.c_str());
    for (size_t i = 0; i < node.size(); ++i) {
      vout.general("%s%zu:\n", indent.c_str(), i);
      traverse(node[i], indent + string("  "));
    }
    break;

  case YAML::NodeType::Map:
    vout.general("%sMap type\n", indent.c_str());
    for (YAML::const_iterator p = node.begin(); p != node.end(); ++p) {
      vout.general("%s%s:\n", indent.c_str(), p->first.as<string>().c_str());
      traverse(p->second, indent + string("  "));
    }
    break;

  default:
    vout.general("%sUnknown type\n", indent.c_str());
  }

  return 0;
}


//====================================================================
#else
//====================================================================

class BAPI Parser_bridge
{
 public:
  static const string class_name;

 public:
  Parser_bridge()
    : m_vl(CommonParameters::Vlevel()) {}

  void parse(std::istream& iss, Parameters& params);

 private:
  int parse_stream(std::istream& iss, Parameters& params);

  int parse_line(char *buf, string& key, string& value);
  int parse_vector(char *buf, vector<string>& vec);

  Bridge::VerboseLevel m_vl;
};

//--------------------------------------------------------------------
const string Parser_bridge::class_name = "ParameterManger_YAML::Parser_bridge";

//--------------------------------------------------------------------
void Parser_bridge::parse(std::istream& iss, Parameters& params)
{
  int result = parse_stream(iss, params);

  if (result != EXIT_SUCCESS) {
    vout.crucial(m_vl, "%s: parse failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//--------------------------------------------------------------------
int Parser_bridge::parse_vector(char *buf, vector<string>& vec)
{
  // N.B. buf modified on exit.
  const char sep = ',';

  int count = 0;

  if ((buf[0] != '[') || (buf[strlen(buf) - 1] != ']')) {
    return -1;
  }

  buf[strlen(buf) - 1] = '\0';

  char *p = buf + 1;

  while (p && *p)
  {
    // skip heading spaces
    while (*p == ' ')
    {
      ++p;
    }

    // find separator
    char *q = strchr(p, sep);

    if (q) {
      // eliminate separator
      *q = '\0';

      // eliminate trailing spaces of the item
      char *r = q - 1;
      while (*r == ' ')
      {
        *r-- = '\0';
      }

      vec.push_back(string(p));
      ++count;

      // go to next item
      p = q + 1;
    } else {
      // separator is not found; maybe last item in the sequence.

      // eliminate trailing spaces;
      char *r = p + strlen(p) - 1;
      while (r >= p && *r == ' ')
      {
        *r-- = '\0';
      }

      if (strlen(p) > 0) {
        vec.push_back(string(p));
        ++count;
      } else {
        // discard
      }

      p = q;
    }
  }

  return count;
}


//--------------------------------------------------------------------
int Parser_bridge::parse_line(char *buf, string& key, string& value)
{
  // N.B. buf modified on exit.

  const char delim = ':';

  // remove comments
  if (char *q = strchr(buf, '#')) { *q = '\0'; }

  // remove trailing spaces
  char *s = buf + strlen(buf) - 1;
  while (*s == ' ')
  {
    *s-- = '\0';
  }

  // find indent
  int indent = 0;

  char *p = buf;
  while (*p == ' ')
  {
    ++p;
    ++indent;
  }

  // find key-value separator
  char *q = strchr(buf, delim);

  if (!q) {
    key   = string();
    value = string();

    return -1;
  }

  // find key
  char *r = q;

  *r = '\0';
  --r;
  while (r >= p && *r == ' ')
  {
    *r-- = '\0';
  }

  key = string(p);

  // find value
  ++q;
  while (*q == ' ')
  {
    ++q;
  }

  value = string(q);

  // return indent
  return indent;
}


//--------------------------------------------------------------------
int Parser_bridge::parse_stream(std::istream& iss, Parameters& params)
{
  int retv = EXIT_SUCCESS;

  const size_t buf_size = 1024;
  char         buf[buf_size];

  typedef pair<string, Parameters *>   env_t;
  typedef pair<int, env_t>             level_t;
  typedef std::stack<level_t>          stack_t;

  stack_t levels;

  Parameters *current_params = &params;
  int        current_indent  = 0;

  bool expect_map = false;

  while (iss.getline(buf, buf_size))
  {
    string key, value;

    int indent = parse_line(buf, key, value);

    if (indent < 0) {
      vout.paranoiac(m_vl, "%s: empty line. skip.\n", class_name.c_str());
      continue;
    }

    // level up/down
    if (indent > current_indent) {
      if (!expect_map) {
        vout.general(m_vl, "%s: Error: unexpected nest level.\n", class_name.c_str());
        retv = EXIT_FAILURE;
        continue;
      }

      // start new level
      current_params = new Parameters;
      current_indent = indent;

      expect_map = false;
    } else {
      if (expect_map) {
        // open key in the previous line actually correspond to empty value

        level_t lv = levels.top();
        levels.pop();

        string     key_l            = lv.second.first;
        Parameters *stored_params = lv.second.second;

        stored_params->set_string(key_l, string()); // null string.

        // current_params is yet not newly created. no need to delete.

        // restore upper level -- maybe unnecessary.
        current_params = stored_params;
        current_indent = lv.first;

        expect_map = false;
      }

      if (indent < current_indent) {
        while (indent < current_indent)
        {
          level_t lv = levels.top();
          levels.pop();

          string     key_l            = lv.second.first;
          Parameters *stored_params = lv.second.second;

          stored_params->set_Parameters(key_l, *current_params);
          delete current_params;

          // restore upper level
          current_params = stored_params;
          current_indent = lv.first;
        }
      } else {  // indent == current_indent
        //
      }
    }

    // store key-value
    if (value.length() > 0) {
      if (value[0] == '[') {
        memset(buf, '\0', buf_size);
        value.copy(buf, buf_size);  // reuse buffer

        vector<string> v;

        int nvalues = parse_vector(buf, v);

        if (nvalues < 0) {
          vout.general(m_vl, "%s: ERROR: parse_vector failed.\n", class_name.c_str());
          continue;
        }

        // store key - vector value.
        current_params->set_string_vector(key, v);
      } else {
        // store key - scalar value.
        current_params->set_string(key, value);
      }
    } else {
      // key for a map in subsequent lines
      expect_map = true;

      // push current environment to stack
      levels.push(level_t(indent, env_t(key, current_params)));
    }
  }

  while (current_indent > 0)
  {
    level_t lv = levels.top();
    levels.pop();

    string     key            = lv.second.first;
    Parameters *stored_params = lv.second.second;

    stored_params->set_Parameters(key, *current_params);
    delete current_params;

    // restore upper level
    current_params = stored_params;
    current_indent = lv.first;
  }

  return retv;
}


#endif

//====================================================================
void ParameterManager_YAML::read_params(std::istream& fin, Parameters& params)
{
#ifdef USE_YAMLCPPLIB
  Parser_yamlcpp().parse(fin, params);
#else
  Parser_bridge().parse(fin, params);
#endif
}


//====================================================================
void ParameterManager_YAML::read_params(const std::string& params_file, Parameters& params)
{
  const int io_node = 0;  // node id for file i/o.

  int  filesize = 0;
  char *buf     = 0;

  if (Communicator::nodeid() == io_node) {
    // load and distribute
    std::ifstream fin(params_file.c_str());
    if (!fin) {
      vout.crucial(m_vl, "Error at %s: unable to read parameter file: %s.\n", class_name.c_str(), params_file.c_str());

      Communicator::abort();
    }

    fin.seekg(0, std::ios::end);
    filesize = static_cast<int>(fin.tellg());
    fin.seekg(0, std::ios::beg);

    int padding = 8 - (filesize % 8);

    vout.paranoiac(m_vl, "%s::read_params: filesize = %d, padding = %d\n", class_name.c_str(), filesize, padding);

    filesize += padding;

    Communicator::broadcast(1, &filesize, io_node);

    buf = new char [filesize];
    memset(buf, 0, filesize);

    fin.read(buf, filesize - padding);

    Communicator::Base::broadcast(filesize, buf, io_node);
  } else {
    // receive from io_node

    Communicator::broadcast(1, &filesize, io_node);

    buf = new char [filesize];
    memset(buf, 0, filesize);

    Communicator::Base::broadcast(filesize, buf, io_node);
  }

  std::istringstream iss(buf);
  read_params(iss, params);

  delete [] buf;
}


//====================================================================
//============================================================END=====
