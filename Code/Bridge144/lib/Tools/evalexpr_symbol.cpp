#include "BridgeLib_Private.h"
#if USE_EVALEXPR

/*!
        @file    $Id: evalexpr_symbol.cpp #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "evalexpr_symbol.h"

using Bridge::vout;

//====================================================================
bool SymbolTable::find_symbol(const std::string& name)
{
  return table.find(name) != table.end();
}


//====================================================================
bool SymbolTable::put_symbol(const std::string& name, const double value)
{
  if (find_symbol(name)) {
    vout.detailed("key \"%s\" already exists. overwrite.\n", name.c_str());
  }

  SymbolRecord rec;
  rec.type      = VARIABLE;
  rec.value.val = value;
  table[name]   = rec;

  return true;
}


//====================================================================
bool SymbolTable::put_symbol(const std::string& name, const function_t tptr)
{
  if (find_symbol(name)) {
    vout.detailed("key \"%s\" already exists. overwrite.\n", name.c_str());
  }

  SymbolRecord rec;
  rec.type       = FUNCTION;
  rec.value.fptr = tptr;
  table[name]    = rec;

  return true;
}


//====================================================================
double SymbolTable::get_symbol_value(const std::string& name) const
{
  SymbolMap_t::const_iterator p = table.find(name);

  if (p != table.end()) {
    return p->second.value.val;
  } else {
    vout.detailed("key \"%s\" not found.\n", name.c_str());
    return double();
  }
}


//====================================================================
function_t SymbolTable::get_symbol_function(const std::string& name) const
{
  SymbolMap_t::const_iterator p = table.find(name);

  if (p != table.end()) {
    return p->second.value.fptr;
  } else {
    vout.detailed("key \"%s\" not found.\n", name.c_str());
    return (function_t)0;
  }
}


//====================================================================
void SymbolTable::dump() const
{
  for (SymbolMap_t::const_iterator p = table.begin(); p != table.end(); ++p) {
    vout.paranoiac("key = %s, ", (p->first).c_str());

    ValueType t = (p->second).type;
    if (t == VARIABLE) {
      vout.paranoiac("type = VARIABLE, value = %f", (p->second).value.val);
    } else if (t == FUNCTION) {
      vout.paranoiac("type = FUNCTION, value = %p", (p->second).value.fptr);
    } else {
      vout.paranoiac("type = UNKNOWN,  ");
    }

    vout.paranoiac("\n");
  }
}


//==========================================================
//==================================================END=====

#endif
