/*!
        @file    $Id:: filename.h #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FILENAME_INCLUDED
#define FILENAME_INCLUDED

//! Filename utility.

/*!
    Filename class generates a string by keyword substitution.

    usage:
        Filename filename(format_string);

        string fn = filename.format(param1, ...);
        string fn = filename.format(key_value_map);

    format_string contains keywords of a format ${keyword:n}
    where keyword is an alpha-numeric string ('_' and '.' are also
    permitted) and an optional number with preceding colon specifies
    the number of digits (zero padded).

    arguments param1, ... to format() are in the order of appearance
    in the format string. format() also takes key-value map.

                                   [25 March 2015 T.Aoyama]
 */

#include <string>
#include <map>
#include <vector>

class Filename
{
 public:
  typedef std::string                          keyword_type;
  typedef int                                  value_type;
  typedef std::map<keyword_type, value_type>   keyvalue_type;
  typedef std::vector<keyword_type>            keywordlist_type;

 public:
  Filename(const std::string& base);

 private:
  // non-copyable
  Filename(const Filename&);
  Filename& operator=(const Filename&);

 public:

  std::string format()
  {
    return format_(0);
  }

  std::string format(const int v1)
  {
    return format_(1, v1);
  }

  std::string format(const int v1, const int v2)
  {
    return format_(2, v1, v2);
  }

  std::string format(const int v1, const int v2, const int v3)
  {
    return format_(3, v1, v2, v3);
  }

  std::string format(const int v1, const int v2, const int v3, const int v4)
  {
    return format_(4, v1, v2, v3, v4);
  }

  std::string format(const int v1, const int v2, const int v3, const int v4, const int v5)
  {
    return format_(5, v1, v2, v3, v4, v5);
  }

  std::string format(const int v1, const int v2, const int v3, const int v4, const int v5, const int v6)
  {
    return format_(6, v1, v2, v3, v4, v5, v6);
  }

  std::string format(const keyvalue_type& kv);

 private:

  static const int format_str_size_ = 1024;
  char             format_str_[format_str_size_];

  keywordlist_type keywords_;

  std::string format_(int narg, ...);

  keywordlist_type generate_pattern(char *fmt, const size_t fmt_size, const char *msg);
};
#endif
