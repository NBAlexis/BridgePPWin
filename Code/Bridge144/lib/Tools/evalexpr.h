/*!
        @file    $Id: evalexpr.h #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef EVALEXPR_INCLUDED
#define EVALEXPR_INCLUDED

#include <map>
#include <string>
#include <cstring>

#include "evalexpr_parser.h"
#include "evalexpr_symbol.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! EvalExpr class for algebraic expression in parameter strings

/**
   EvalExpr class provides support to evaluating algebraic expressions
   in parameter strings. The expression is evaluated into a single double
   precision value.

   Supported expression includes arithmetic operations in in-fix notation,
   several transcendental functions including trigonometric and logarithmic
   functions, and numerical constant such as pi.
 */

class EvalExpr
{
 public:

  typedef yy::parser::semantic_type   semantic_type;
  typedef yy::parser::token           token;

  EvalExpr(const std::string& line, bool debug = false)
    : m_src(line), m_pos(0), m_result(0),
      m_trace(debug) {}

  // lexer
  int next_token(semantic_type& yylval);

  // parse
//  bool parse(double& result);
  double parse();

  // symbol_table manipulation
  double get_symbol_value(char const *name);
  function_t get_symbol_function(char const *name);

  // parser interface
  void set_result(double result);

  void error(const std::string& msg);

 private:
  const std::string m_src;
  unsigned int      m_pos;
  double            m_result;

  bool m_trace;

  static SymbolTable global_symbol_table;
};


// lexer interface for parser
int yylex(EvalExpr::semantic_type *yylval, EvalExpr& driver);
#endif
