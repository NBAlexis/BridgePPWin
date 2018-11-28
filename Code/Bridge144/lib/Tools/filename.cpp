#include "BridgeLib_Private.h"

/*!
        @file    $Id:: filename.cpp #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include <cstdio>
#if BRIDGE_WIN
#include <regex>
#else
#include <regex.h>
#endif
#include <string.h>
#include <cassert>
#include "filename.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//====================================================================
Filename::Filename(const std::string& base)
{
  keywords_ = generate_pattern(format_str_, format_str_size_, base.c_str());
}


//====================================================================
std::string Filename::format(const Filename::keyvalue_type& kv)
{
  size_t nkeyword = keywords_.size();

  std::vector<int> v(nkeyword);

  for (size_t i = 0; i < nkeyword; ++i) {
    Filename::keyvalue_type::const_iterator p = kv.find(keywords_[i]);

    if (p != kv.end()) {
      v[i] = p->second;
    } else {
      vout.general("error: key \"%s\" not found.\n", keywords_[i].c_str());
      return std::string();
    }
  }

  switch (nkeyword)
  {
  case 0:
    return format_(0);

  case 1:
    return format_(1, v[0]);

  case 2:
    return format_(2, v[0], v[1]);

  case 3:
    return format_(3, v[0], v[1], v[2]);

  case 4:
    return format_(4, v[0], v[1], v[2], v[3]);

  case 5:
    return format_(5, v[0], v[1], v[2], v[3], v[4]);

  case 6:
    return format_(6, v[0], v[1], v[2], v[3], v[4], v[5]);

  default:
    vout.general("error: too many parameters.\n");
  }

  return std::string();
}


//====================================================================
std::string Filename::format_(int narg, ...)
{
  assert(narg == int(keywords_.size()));

  static const int buf_size = 4096;
  static char      buf[buf_size];

  va_list arg;

  va_start(arg, narg);
  vsnprintf(buf, buf_size, format_str_, arg);
  va_end(arg);

  return std::string(buf);
}

#if !BRIDGE_WIN
//====================================================================
Filename::keywordlist_type Filename::generate_pattern(char *fmt, const size_t fmt_size, const char *msg)
{
  keywordlist_type keyword_list;

  const char pattern[]    = "\\{([._[:alnum:]]+)(:([[:digit:]]+))?\\}";
  const int  num_regmatch = 4; // number of parentheses + 1

  regex_t reg;

  if (int retv = regcomp(&reg, pattern, REG_EXTENDED)) {
    vout.crucial("regcomp failed. %s(%d)\n", strerror(retv), retv);
    exit(1);
  }

  regmatch_t regmatch[num_regmatch];

  memset(fmt, '\0', fmt_size);
  char *bufp = fmt;

  int      nmatch = 0;
  regoff_t idx    = 0;

  while (true)
  {
    int v = regexec(&reg, msg + idx, num_regmatch, regmatch, 0);

    if (v != 0) break;

    ++nmatch;

    strncpy(bufp, msg + idx, regmatch[0].rm_so);
    bufp += regmatch[0].rm_so;

    char keyword[128];
    memset(keyword, '\0', 128);
    strncpy(keyword, msg + idx + regmatch[1].rm_so, regmatch[1].rm_eo - regmatch[1].rm_so);

    keyword_list.push_back(keyword_type(keyword));

    int width = 0;
    if (regmatch[3].rm_so != (regoff_t)-1) {
      width = atoi(msg + idx + regmatch[3].rm_so);
    }

    if (width > 0) {
      char tmp[128];
      memset(tmp, '\0', 128);

      snprintf(tmp, 128, "%%0%dd", width);

      strncpy(bufp, tmp, strlen(tmp));
      bufp += strlen(tmp);
    } else {
      strncpy(bufp, "%d", strlen("%d"));
      bufp += strlen("%d");
    }

    idx += regmatch[0].rm_eo;
  }

  strncpy(bufp, msg + idx, strlen(msg) - idx);

  regfree(&reg);

  return keyword_list;
}
#endif

#if BRIDGE_WIN

//What does this function do?????
Filename::keywordlist_type Filename::generate_pattern(char *fmt, const size_t fmt_size, const char *msg)
{
    keywordlist_type keyword_list;

    const char pattern[] = "\\{([._[:alnum:]]+)(:([[:digit:]]+))?\\}";
    const int  num_regmatch = 4; // number of parentheses + 1

    std::cmatch narrowMatch;
    std::regex rx(pattern);
    bool bFound = std::regex_search(msg, msg + strlen(msg), narrowMatch, rx);
    memset(fmt, '\0', fmt_size);
    if (bFound)
    {
        for (int i = 0; i < narrowMatch.size(); ++i)
        {
            keyword_list.push_back(keyword_type(narrowMatch[i].str()));
        }
    }

    return keyword_list;
}
#endif

//====================================================================
//============================================================END=====
