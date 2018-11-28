/*!
        @file    $Id:: parameters.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
 */

#ifndef PARAMETERS_INCLUDED
#define PARAMETERS_INCLUDED

#include <string>
#include <map>
#include <vector>
#include <sstream>

#include "Communicator/communicator.h"
#include "checker.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

using std::string;
using std::map;
using std::vector;

//! Class for parameters

/*!
    Base class of Parameters.
    This class gives general basis of structured parameter sets.
                                       [17 Jul 2012 H.Matsufuru]
    fetch is modified to return int, instead of bool in the original
    Aoyama-san's code.                 [16 Jun 2013 Y.Namekawa]

    Renewed to be a parameter container.
                                       [25 June 2016 T.Aoyama]

 */

class Parameters
{
 public:
  static const std::string class_name;

 public:
  Parameters();
  virtual ~Parameters() {}

  void set_double(const string& key, const double value);
  void set_int(const string& key, const int value);
  void set_string(const string& key, const string& value);
  void set_double_vector(const string& key, const vector<double>& value);
  void set_int_vector(const string& key, const vector<int>& value);
  void set_string_vector(const string& key, const vector<string>& value);
  void set_Parameters(const string& key, const Parameters& value);
  void set_VerboseLevel(const Bridge::VerboseLevel value);

  double get_double(const string& key) const;
  int get_int(const string& key) const;
  unsigned long get_unsigned_long(const string& key) const;
  string get_string(const string& key) const;

  vector<double> get_double_vector(const string& key) const;

  vector<int> get_int_vector(const string& key) const;

  vector<string> get_string_vector(const string& key) const;
  Parameters get_Parameters(const string& key) const;
  Parameters& get_Parameters(const string& key);
  Bridge::VerboseLevel get_VerboseLevel() const;

  Parameters lookup(const string& key) const { return get_Parameters(key); }
  Parameters& lookup(const string& key) { return get_Parameters(key); }

  int fetch_double(const string& key, double& value) const;
  int fetch_int(const string& key, int& value) const;
  int fetch_unsigned_long(const string& key, unsigned long& value) const;
  int fetch_string(const string& key, string& value) const;
  int fetch_double_vector(const string& key, vector<double>& value) const;
  int fetch_int_vector(const string& key, vector<int>& value) const;
  int fetch_string_vector(const string& key, vector<string>& value) const;
  int fetch_VerboseLevel(Bridge::VerboseLevel& value) const;

  bool find_double(const string& key) const;
  bool find_int(const string& key) const;
  bool find_unsigned_long(const string& key) const;
  bool find_string(const string& key) const;
  bool find_double_vector(const string& key) const;
  bool find_int_vector(const string& key) const;
  bool find_string_vector(const string& key) const;
  bool find_Parameters(const string& key) const;

  bool is_set(const string& key) const;

  // obsolete. just for compatibility.
  void Register_double(const string& key, const double defvalue);
  void Register_int(const string& key, const int defvalue);
  void Register_string(const string& key, const string& defvalue);
  void Register_double_vector(const string& key, const vector<double>& defvalue);
  void Register_int_vector(const string& key, const vector<int>& defvalue);
  void Register_string_vector(const string& key, const vector<string>& defvalue);
  void Register_Parameters(const string& key, const Parameters& defvalue);
  void Register_Parameters(const string& key, Parameters *const defvalue);
  void Register_Parameters(const string& key, const unique_ptr<Parameters>& defvalue);

 private:
  // scalar
  map<string, double> m_map_double;
  map<string, int>    m_map_int;
  map<string, string> m_map_string;
  // array
  map<string, vector<double> > m_map_double_vector;
  map<string, vector<int> >    m_map_int_vector;
  map<string, vector<string> > m_map_string_vector;
  // map
  map<string, Parameters> m_map_parameters;
  // verbose level
  Bridge::VerboseLevel m_vlevel;

  // utility
  double convert_to_double(const string&) const;

  vector<double> convert_to_double(const vector<string>&) const;

 public:
  void dump(const string& indent = "") const;

  // utility
  template<typename T>
  static
  string to_string(const vector<T>& v)
  {
    std::stringstream ss;

    ss << "[";
    for (size_t i = 0, n = v.size(); i < n; ++i) {
      if (i > 0) ss << ", ";
      ss << v[i];
    }
    ss << "]";

    return ss.str();
  }
};
#endif
