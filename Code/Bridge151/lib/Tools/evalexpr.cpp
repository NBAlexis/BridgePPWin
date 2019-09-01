#include "BridgeLib_Private.h"
#if USE_EVALEXPR

/*!
        @file    $Id: evalexpr.cpp #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#include "evalexpr.h"

#ifdef __PGI
#define USE_STRDUP
#endif

#ifdef NECSX
#define USE_STRDUP
#endif

//====================================================================
// global symbol entries
//namespace {
//  struct _init_function
//  {
//    char const *name;
//    function_t func;
//  };
//
//  struct _init_variable
//  {
//    char const *name;
//    double     val;
//  };
//}

// import global symbol definitions
#include "evalexpr_global.h"

// register global symbols
namespace {
  SymbolTable init_global_symbol()
  {
    SymbolTable table;

    for (int i = 0; arithmetic_functions[i].name != 0; ++i) {
      table.put_symbol(arithmetic_functions[i].name, arithmetic_functions[i].func);
    }

    for (int i = 0; predefined_constants[i].name != 0; ++i) {
      table.put_symbol(predefined_constants[i].name, predefined_constants[i].val);
    }

    return table;
  }
}

SymbolTable EvalExpr::global_symbol_table = init_global_symbol();

//====================================================================
// lexer interface for parser
int yylex(EvalExpr::semantic_type *yylval, EvalExpr& driver)
{
  return driver.next_token(*yylval);
}


//====================================================================
int EvalExpr::next_token(semantic_type& yylval)
{
  if (m_trace) vout.paranoiac("%s: pos = %d, length = %zu\n", __func__, m_pos, m_src.length());

  char *p = const_cast<char *>(m_src.c_str()) + m_pos;

  int c = *p;

  // skip spaces
  while ((c == ' ') || (c == '\t'))
  {
    c = *++p;
    ++m_pos;
  }

  if (c == '\0') return 0;

  if (c == '\n') return token::EOL;

  // find number [0-9](.[0-9]*)(e[+-]?[0-9]+)
  if (isdigit(c)) {
    enum
    {
      IN_INT, IN_FRAC, IN_EXP, IN_EXP_DIGIT
    }
    state = IN_INT;

    int i;
    for (i = 0; ; c = p[++i]) {
      if (state == IN_INT) {
        if (isdigit(c)) {
        } else if (c == '.') {
          state = IN_FRAC;
        } else if (c == 'e') {
          state = IN_EXP;
        } else {
          break;  // accept
        }
      } else if (state == IN_FRAC) {
        if (isdigit(c)) {
        } else if (c == 'e') {
          state = IN_EXP;
        } else {
          break;  // accept
        }
      } else if ((state == IN_EXP) || (state == IN_EXP_DIGIT)) {
        if (isdigit(c)) {
        } else if ((state == IN_EXP) && ((c == '+') || (c == '-'))) {
          state = IN_EXP_DIGIT;
        } else {
          break;  // accept
        }
      } else {
        exit(EXIT_FAILURE);
      }
    }

    double x = atof(p);

#ifdef USE_STRDUP
    char *symbol = _strdup(p);
#else
    char *symbol = strndup(p, i + 1);
#endif
    symbol[i] = '\0';

    if (m_trace) vout.paranoiac("%s: accept: number: \"%s\"(%d), %f\n", __func__, symbol, i, x);

    free(symbol);

    yylval.val = x;
    m_pos     += i;

    return token::NUMBER;
  }

  // find identifier [a-zA-Z]([.a-zA-Z0-9]*)
  if (isalpha(c)) {
    int i;
    for (i = 0; isalnum(p[i]) || p[i] == '.'; ++i) {
      continue;
    }

#ifdef USE_STRDUP
    char *symbol = _strdup(p);
#else
    char *symbol = strndup(p, i + 1);
#endif
    symbol[i] = '\0';

    if (m_trace) vout.paranoiac("%s: accept: identifier: \"%s\" (%d)\n", __func__, symbol, i);

    yylval.sym = symbol;
    m_pos     += i;

    return token::IDENTIFIER;
  }

  ++m_pos;

  return c;
}


//====================================================================
double EvalExpr::parse()
{
  yy::parser parser(*this);

  parser.set_debug_level(m_trace);

  int retv = parser.parse();  // returns 0 if successful.

  if (retv == 0) {
    if (m_trace) vout.paranoiac("%s: accept, result = %f\n", __func__, m_result);
    return m_result;
  } else {
    if (m_trace) vout.paranoiac("%s: reject\n", __func__);

    vout.crucial("EvalExpr: parse failed.\n");
    exit(EXIT_FAILURE);

    // return double();  // should throw exception or abort.
  }
}


//====================================================================
void EvalExpr::set_result(double result)
{
  if (m_trace) vout.paranoiac("%s: result = %f\n", __func__, result);

  m_result = result;
}


//====================================================================
double EvalExpr::get_symbol_value(char const *name)
{
  return global_symbol_table.get_symbol_value(name);
}


//====================================================================
function_t EvalExpr::get_symbol_function(char const *name)
{
  return global_symbol_table.get_symbol_function(name);
}


//====================================================================
void EvalExpr::error(const std::string& msg)
{
  vout.general("EvalExpr: %s\n", msg.c_str());
}


//====================================================================
//============================================================END=====
#endif
