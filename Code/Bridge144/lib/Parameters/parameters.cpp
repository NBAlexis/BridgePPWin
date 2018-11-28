#include "BridgeLib_Private.h"

/*!
        @file    $Id:: parameters.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
 */

#include "parameters.h"
#include "Tools/evalexpr.h"

const std::string Parameters::class_name = "Parameters";

//====================================================================
Parameters::Parameters()
  : m_vlevel(CommonParameters::Vlevel())
{
  // default verbosity from common_parameters.
}


//====================================================================
void Parameters::set_double(const string& key, const double value)
{ m_map_double[key] = value; }

void Parameters::set_int(const string& key, const int value)
{ m_map_int[key] = value; }

void Parameters::set_string(const string& key, const string& value)
{ m_map_string[key] = value; }

void Parameters::set_double_vector(const string& key, const vector<double>& value)
{ m_map_double_vector[key] = value; }

void Parameters::set_int_vector(const string& key, const vector<int>& value)
{ m_map_int_vector[key] = value; }

void Parameters::set_string_vector(const string& key, const vector<string>& value)
{ m_map_string_vector[key] = value; }

void Parameters::set_Parameters(const string& key, const Parameters& value)
{ m_map_parameters[key] = value; }

void Parameters::set_VerboseLevel(const Bridge::VerboseLevel value)
{ m_vlevel = value; }

//====================================================================
double Parameters::convert_to_double(const string& value) const
{
  return EvalExpr(value).parse();
}


vector<double> Parameters::convert_to_double(const vector<string>& value) const
{
  vector<double> v;
  for (size_t i = 0; i < value.size(); ++i) {
    v.push_back(convert_to_double(value[i]));
  }
  return v;
}


//====================================================================
double Parameters::get_double(const string& key) const
{
  map<string, double>::const_iterator p = m_map_double.find(key);
  if (p != m_map_double.end()) {
    return p->second;
  }

  map<string, string>::const_iterator q = m_map_string.find(key);
  if (q != m_map_string.end()) {
    return convert_to_double(q->second);
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return double();
}


int Parameters::get_int(const string& key) const
{
  map<string, int>::const_iterator p = m_map_int.find(key);
  if (p != m_map_int.end()) {
    return p->second;
  }

  map<string, string>::const_iterator q = m_map_string.find(key);
  if (q != m_map_string.end()) {
    return atoi(q->second.c_str());
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return int();
}


unsigned long Parameters::get_unsigned_long(const string& key) const
{
  map<string, string>::const_iterator q = m_map_string.find(key);
  if (q != m_map_string.end()) {
    return strtoul(q->second.c_str(), NULL, 0);
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return 0;
}


string Parameters::get_string(const string& key) const
{
  map<string, string>::const_iterator p = m_map_string.find(key);
  if (p != m_map_string.end()) {
    return p->second;
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return string();
}


vector<double> Parameters::get_double_vector(const string& key) const
{
  map<string, vector<double> >::const_iterator p = m_map_double_vector.find(key);
  if (p != m_map_double_vector.end()) {
    return p->second;
  }

  map<string, vector<string> >::const_iterator q = m_map_string_vector.find(key);
  if (q != m_map_string_vector.end()) {
    return convert_to_double(q->second);
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return vector<double>();
}


vector<int> Parameters::get_int_vector(const string& key) const
{
  map<string, vector<int> >::const_iterator p = m_map_int_vector.find(key);
  if (p != m_map_int_vector.end()) {
    return p->second;
  }

  map<string, vector<string> >::const_iterator q = m_map_string_vector.find(key);
  if (q != m_map_string_vector.end()) {
    vector<int> v;

    for (size_t i = 0; i < q->second.size(); ++i) {
      v.push_back(atoi(q->second[i].c_str()));
    }

    return v;
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return vector<int>();
}


vector<string> Parameters::get_string_vector(const string& key) const
{
  map<string, vector<string> >::const_iterator p = m_map_string_vector.find(key);
  if (p != m_map_string_vector.end()) {
    return p->second;
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return vector<string>();
}


Parameters Parameters::get_Parameters(const string& key) const
{
  map<string, Parameters>::const_iterator p = m_map_parameters.find(key);
  if (p != m_map_parameters.end()) {
    return p->second;
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return Parameters();
}


Parameters& Parameters::get_Parameters(const string& key)
{
  map<string, Parameters>::iterator p = m_map_parameters.find(key);
  if (p != m_map_parameters.end()) {
    return p->second;
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return *this;
}


Bridge::VerboseLevel Parameters::get_VerboseLevel() const
{
  return m_vlevel;
}


//====================================================================
int Parameters::fetch_double(const string& key, double& value) const
{
  map<string, double>::const_iterator p = m_map_double.find(key);
  if (p != m_map_double.end()) {
    value = p->second;
    return EXIT_SUCCESS;
  }

  map<string, string>::const_iterator q = m_map_string.find(key);
  if (q != m_map_string.end()) {
    value = convert_to_double(q->second);
    return EXIT_SUCCESS;
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return EXIT_FAILURE;
}


int Parameters::fetch_int(const string& key, int& value) const
{
  map<string, int>::const_iterator p = m_map_int.find(key);
  if (p != m_map_int.end()) {
    value = p->second;
    return EXIT_SUCCESS;
  }

  map<string, string>::const_iterator q = m_map_string.find(key);
  if (q != m_map_string.end()) {
    value = atoi(q->second.c_str());
    return EXIT_SUCCESS;
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return EXIT_FAILURE;
}


int Parameters::fetch_unsigned_long(const string& key, unsigned long& value) const
{
  map<string, string>::const_iterator q = m_map_string.find(key);
  if (q != m_map_string.end()) {
    value = strtoul(q->second.c_str(), NULL, 0);
    return EXIT_SUCCESS;
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return EXIT_FAILURE;
}


int Parameters::fetch_string(const string& key, string& value) const
{
  map<string, string>::const_iterator p = m_map_string.find(key);
  if (p != m_map_string.end()) {
    value = p->second;
    return EXIT_SUCCESS;
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return EXIT_FAILURE;
}


int Parameters::fetch_double_vector(const string& key, vector<double>& value) const
{
  map<string, vector<double> >::const_iterator p = m_map_double_vector.find(key);
  if (p != m_map_double_vector.end()) {
    value = p->second;
    return EXIT_SUCCESS;
  }

  map<string, vector<string> >::const_iterator q = m_map_string_vector.find(key);
  if (q != m_map_string_vector.end()) {
    value = convert_to_double(q->second);
    return EXIT_SUCCESS;
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return EXIT_FAILURE;
}


int Parameters::fetch_int_vector(const string& key, vector<int>& value) const
{
  map<string, vector<int> >::const_iterator p = m_map_int_vector.find(key);
  if (p != m_map_int_vector.end()) {
    value = p->second;
    return EXIT_SUCCESS;
  }

  map<string, vector<string> >::const_iterator q = m_map_string_vector.find(key);
  if (q != m_map_string_vector.end()) {
    vector<int> v;

    for (size_t i = 0; i < q->second.size(); ++i) {
      v.push_back(atoi(q->second[i].c_str()));
    }

    value = v;
    return EXIT_SUCCESS;
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return EXIT_FAILURE;
}


int Parameters::fetch_string_vector(const string& key, vector<string>& value) const
{
  map<string, vector<string> >::const_iterator p = m_map_string_vector.find(key);
  if (p != m_map_string_vector.end()) {
    value = p->second;
    return EXIT_SUCCESS;
  }

  vout.crucial("%s: %s: key '%s' not found.\n", class_name.c_str(), __func__, key.c_str());
  return EXIT_FAILURE;
}


int Parameters::fetch_VerboseLevel(Bridge::VerboseLevel& value) const
{
  value = m_vlevel;
  return EXIT_SUCCESS;
}


//====================================================================
void Parameters::Register_double(const string& key, const double defvalue)
{ vout.crucial("%s: %s: unsupported.\n", class_name.c_str(), __func__); }

void Parameters::Register_int(const string& key, const int defvalue)
{ vout.crucial("%s: %s: unsupported.\n", class_name.c_str(), __func__); }

void Parameters::Register_string(const string& key, const string& defvalue)
{ vout.crucial("%s: %s: unsupported.\n", class_name.c_str(), __func__); }

void Parameters::Register_double_vector(const string& key, const vector<double>& defvalue)
{ vout.crucial("%s: %s: unsupported.\n", class_name.c_str(), __func__); }

void Parameters::Register_int_vector(const string& key, const vector<int>& defvalue)
{ vout.crucial("%s: %s: unsupported.\n", class_name.c_str(), __func__); }

void Parameters::Register_string_vector(const string& key, const vector<string>& defvalue)
{ vout.crucial("%s: %s: unsupported.\n", class_name.c_str(), __func__); }

void Parameters::Register_Parameters(const string& key, const Parameters& defvalue)
{ vout.crucial("%s: %s: unsupported.\n", class_name.c_str(), __func__); }

void Parameters::Register_Parameters(const string& key, Parameters *const defvalue)
{ vout.crucial("%s: %s: unsupported.\n", class_name.c_str(), __func__); }

void Parameters::Register_Parameters(const string& key, const unique_ptr<Parameters>& defvalue)
{ vout.crucial("%s: %s: unsupported.\n", class_name.c_str(), __func__); }

//====================================================================
bool Parameters::find_double(const string& key) const
{ return m_map_double.find(key) != m_map_double.end() || m_map_string.find(key) != m_map_string.end(); }

bool Parameters::find_int(const string& key) const
{ return m_map_int.find(key) != m_map_int.end() || m_map_string.find(key) != m_map_string.end(); }

bool Parameters::find_unsigned_long(const string& key) const
{
  // unsigned long value is stored as a string.
  return m_map_string.find(key) != m_map_string.end();
}


bool Parameters::find_string(const string& key) const
{ return m_map_string.find(key) != m_map_string.end(); }

bool Parameters::find_double_vector(const string& key) const
{ return m_map_double_vector.find(key) != m_map_double_vector.end() || m_map_string_vector.find(key) != m_map_string_vector.end(); }

bool Parameters::find_int_vector(const string& key) const
{ return m_map_int_vector.find(key) != m_map_int_vector.end() || m_map_string_vector.find(key) != m_map_string_vector.end(); }

bool Parameters::find_string_vector(const string& key) const
{ return m_map_string_vector.find(key) != m_map_string_vector.end(); }

bool Parameters::find_Parameters(const string& key) const
{ return m_map_parameters.find(key) != m_map_parameters.end(); }

bool Parameters::is_set(const string& key) const
{
  return find_double(key) ||
         find_int(key) ||
         find_string(key) ||
         find_double_vector(key) ||
         find_int_vector(key) ||
         find_string_vector(key) ||
         find_Parameters(key)
  ;
}


//====================================================================
void Parameters::dump(const string& indent) const
{
  printf("%sScalar<double>:\n", indent.c_str());
  if (m_map_double.size() == 0) {
    printf("%s  (none)\n", indent.c_str());
  } else {
    for (map<string, double>::const_iterator p = m_map_double.begin(); p != m_map_double.end(); ++p) {
      printf("%s  key = %s\tvalue = %e\n", indent.c_str(), p->first.c_str(), p->second);
    }
  }

  printf("%sScalar<int>:\n", indent.c_str());
  if (m_map_int.size() == 0) {
    printf("%s  (none)\n", indent.c_str());
  } else {
    for (map<string, int>::const_iterator p = m_map_int.begin(); p != m_map_int.end(); ++p) {
      printf("%s  key = %s\tvalue = %d\n", indent.c_str(), p->first.c_str(), p->second);
    }
  }

  printf("%sScalar<string>:\n", indent.c_str());
  if (m_map_string.size() == 0) {
    printf("%s  (none)\n", indent.c_str());
  } else {
    for (map<string, string>::const_iterator p = m_map_string.begin(); p != m_map_string.end(); ++p) {
      printf("%s  key = %s\tvalue = %s\n", indent.c_str(), p->first.c_str(), p->second.c_str());
    }
  }

  printf("%sVector<double>:\n", indent.c_str());
  if (m_map_double_vector.size() == 0) {
    printf("%s  (none)\n", indent.c_str());
  } else {
    for (map<string, vector<double> >::const_iterator p = m_map_double_vector.begin(); p != m_map_double_vector.end(); ++p) {
      printf("%s  key = %s\tvalue = [ \n", indent.c_str(), p->first.c_str());
      for (size_t i = 0; i < p->second.size(); ++i) {
        printf("%e, ", p->second[i]);
      }
      printf("]\n");
    }
  }

  printf("%sVector<int>:\n", indent.c_str());
  if (m_map_int_vector.size() == 0) {
    printf("%s  (none)\n", indent.c_str());
  } else {
    for (map<string, vector<int> >::const_iterator p = m_map_int_vector.begin(); p != m_map_int_vector.end(); ++p) {
      printf("%s  key = %s\tvalue = [ \n", indent.c_str(), p->first.c_str());
      for (size_t i = 0; i < p->second.size(); ++i) {
        printf("%d, ", p->second[i]);
      }
      printf("]\n");
    }
  }

  printf("%sVector<string>:\n", indent.c_str());
  if (m_map_string_vector.size() == 0) {
    printf("%s  (none)\n", indent.c_str());
  } else {
    for (map<string, vector<string> >::const_iterator p = m_map_string_vector.begin(); p != m_map_string_vector.end(); ++p) {
      printf("%s  key = %s\tvalue = [ ", indent.c_str(), p->first.c_str());
      for (size_t i = 0; i < p->second.size(); ++i) {
        printf("%s, ", p->second[i].c_str());
      }
      printf("]\n");
    }
  }

  printf("%sParameters:\n", indent.c_str());
  if (m_map_parameters.size() == 0) {
    printf("%s  (none)\n", indent.c_str());
  } else {
    for (map<string, Parameters>::const_iterator p = m_map_parameters.begin(); p != m_map_parameters.end(); ++p) {
      printf("%s  key = %s, value:\n", indent.c_str(), p->first.c_str());
      p->second.dump(indent + "    ");
    }
  }
}


//====================================================================
//============================================================END=====
