/*!
        @file    $Id: evalexpr_symbol.h #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef CALC_SYMBOL_H
#define CALC_SYMBOL_H

#include <map>
#include <string>
#include "IO/bridgeIO.h"

class SymbolTable
{
 private:

  enum ValueType { VARIABLE, FUNCTION, };

  struct SymbolRecord
  {
    ValueType type;
    union
    {
      double     val;
      function_t fptr;
    }
              value;
  };

  typedef std::map<std::string, SymbolRecord>   SymbolMap_t;

  SymbolMap_t table;

 public:

  bool find_symbol(const std::string& name);

  bool put_symbol(const std::string& name, const double value);
  bool put_symbol(const std::string& name, const function_t tptr);

  double get_symbol_value(const std::string& name) const;
  function_t get_symbol_function(const std::string& name) const;

  void dump() const;
};
#endif
