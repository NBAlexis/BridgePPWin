#include "BridgeLib_Private.h"

/*!
        @file    $Id:: parameterManager_XML.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
 */

#include "parameterManager_XML.h"

#include "Tools/evalexpr.h"

#include <fstream>
#include <sstream>
#include <map>
#include <vector>

#include "tinyxml2.h"
using namespace tinyxml2;

#include "IO/bridgeIO.h"
using Bridge::vout;

using std::string;
using std::vector;
using std::map;

const std::string ParameterManager_XML::class_name = "ParameterManager_XML";

//====================================================================
class Parser_tinyxml
{
 public:
  static const string class_name;

 public:
  int parse(const char *buf, Parameters& params);

  Parser_tinyxml()
    : m_vl(CommonParameters::Vlevel()) {}

 private:

  int store_map_data(const XMLElement *, Parameters&);
  int store_vector_data(const XMLElement *, vector<string>&);

  void traverse(const XMLElement *, const string& indent = "") const;

 protected:
  Bridge::VerboseLevel m_vl;
};

//--------------------------------------------------------------------
const string Parser_tinyxml::class_name = "ParameterManager_XML::Parser_tinyxml";

//--------------------------------------------------------------------
int Parser_tinyxml::parse(const char *buf, Parameters& params)
{
  if (!buf) return EXIT_FAILURE;

  XMLDocument doc;

  XMLError err = doc.Parse(buf);

  if (err != XML_NO_ERROR) {
    vout.crucial(m_vl, "Error at %s: parse failed: %s\n", class_name.c_str(), doc.ErrorName());

    const char *err_str1 = doc.GetErrorStr1();
    if (err_str1) {
      vout.crucial(m_vl, "  %s\n", err_str1);
    }
    const char *err_str2 = doc.GetErrorStr2();
    if (err_str2) {
      vout.crucial(m_vl, "  %s\n", err_str2);
    }

    exit(EXIT_FAILURE);
  }

  XMLElement *root = doc.FirstChildElement("Parameters");
  if (!root) {
    vout.crucial(m_vl, "Error at %s: tag \"Parameters\" not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  // traverse(root);  // for debug

  for (const XMLElement *e = root->FirstChildElement(); e; e = e->NextSiblingElement()) {
    int retv = store_map_data(e, params);

    if (retv != EXIT_SUCCESS) {
      vout.crucial(m_vl, "Error at %s: parse failed.\n", class_name.c_str());
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}


//--------------------------------------------------------------------
int Parser_tinyxml::store_vector_data(const XMLElement *elem, vector<string>& v)
{
  if (!elem) return EXIT_FAILURE;

  const XMLElement *ee = elem->FirstChildElement();
  if (!ee) return EXIT_FAILURE;

  map<int, string> mm;
  int              max_index = 0;

  if (strcmp(ee->Name(), "value") == 0) {  // <value id="n">item</value>
    int count = 0;

    for (const XMLElement *e = elem->FirstChildElement(); e; e = e->NextSiblingElement()) {
      int id = 0;

      if (e->QueryIntAttribute("id", &id) == XML_SUCCESS) {
        --id;  // id = 1, 2, 3, ...
      } else {
        id = count;
      }

      if (id < 0) {
        fprintf(stderr, "Error. inappropriate id\n");
        return EXIT_FAILURE;
      }

      if (id > max_index) max_index = id;

      mm[id] = e->GetText();
      ++count;
    }
  } else {  // <x>item</x> <y>item</y> ...
    for (const XMLElement *e = elem->FirstChildElement(); e; e = e->NextSiblingElement()) {
      int id = 0;

      string tag = e->Name();

      if (tag == "x") {
        id = 0;
      } else if (tag == "y") {
        id = 1;
      } else if (tag == "z") {
        id = 2;
      } else if (tag == "t") {
        id = 3;
      } else if (tag == "w") {
        id = 4;
      } else {
        fprintf(stderr, "unknown tag.\n");
        return EXIT_FAILURE;
      }

      if (id > max_index) max_index = id;

      mm[id] = e->GetText();
    }
  }

  vector<string> vv(max_index + 1);

  for (map<int, string>::const_iterator p = mm.begin(); p != mm.end(); ++p) {
    vv[p->first] = p->second;
  }

  v = vv;

  return EXIT_SUCCESS;
}


//--------------------------------------------------------------------
int Parser_tinyxml::store_map_data(const XMLElement *elem, Parameters& params)
{
  if (!elem) return EXIT_FAILURE;

  const char *elem_name  = elem->Name();
  const char *elem_value = elem->GetText();

  if (elem_value != NULL) {  /* assume leaf node */
    // store to scalar map
    params.set_string(elem_name, elem_value);

    return EXIT_SUCCESS;
  }

  const char *elem_attr = elem->Attribute("type");

  if (elem_attr && (strcmp(elem_attr, "sequence") == 0)) {  /* sequence */
    vector<string> v;
    int            retv = store_vector_data(elem, v);

    if (retv == EXIT_SUCCESS) {
      // store to sequence map
      params.set_string_vector(elem_name, v);
    }

    return retv;
  }

  Parameters pp;

  for (const XMLElement *e = elem->FirstChildElement(); e; e = e->NextSiblingElement()) {
    int retv = store_map_data(e, pp);

    if (retv != EXIT_SUCCESS) {
      return EXIT_FAILURE;
    }
  }

  // store to map map
  params.set_Parameters(elem_name, pp);

  return EXIT_SUCCESS;
}


//--------------------------------------------------------------------
void Parser_tinyxml::traverse(const XMLElement *elem, const string& indent) const
{
  if (!elem) return;

  const char *elem_name  = elem->Name();
  const char *elem_value = elem->GetText();

  if (elem_value != NULL) {  /* assume leaf node */
    vout.general(m_vl, "%sElementName = %s, value = %s\n", indent.c_str(), elem_name, elem_value);
    return;
  } else {
    vout.general(m_vl, "%sElementName = %s\n", indent.c_str(), elem_name);

    const char *elem_attr = elem->Attribute("type");
    if (elem_attr && (strcmp(elem_attr, "sequence") == 0)) {
      vout.general(m_vl, "%s- is a sequence.\n", indent.c_str());
    }

    for (const XMLElement *e = elem->FirstChildElement(); e; e = e->NextSiblingElement()) {
      traverse(e, indent + "  ");
    }
  }
}


//--------------------------------------------------------------------

//====================================================================
void ParameterManager_XML::read_params(const std::string& params_file, Parameters& params)
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
    filesize = fin.tellg();
    fin.seekg(0, std::ios::beg);

    int padding = 8 - (filesize % 8);

    vout.paranoiac(m_vl, "%s::%s: filesize = %d, padding = %d\n", class_name.c_str(), __func__, filesize, padding);

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

  Parser_tinyxml().parse(buf, params);

  // params.dump();  // for debug

  delete [] buf;
}


//====================================================================
//============================================================END=====
