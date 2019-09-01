/*!
        @file    $Id:: evalexpr_parser.yy #$

        @brief

        @author  Tatsumi Aoyama (aoym)
	         $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2016-10-04 16:56:56 #$

        @version $LastChangedRevision: 1497 $
*/

%skeleton "lalr1.cc"
%require "2.7"
%defines "evalexpr_parser.h"
%define namespace "yy"
%define parser_class_name "parser"

%code requires {
#include <cmath>
#include <string>

class EvalExpr;
}

%parse-param { EvalExpr& driver }
%lex-param   { EvalExpr& driver }

%debug
%error-verbose

%union value
{
  double val;
  char *sym;
}

%code {
#include "evalexpr.h"
}

%token EOL
%token <sym> IDENTIFIER
%token <val> NUMBER
%type  <val> expr

%left '-' '+'
%left '*' '/'
%left PLS NEG
%right '^'

%printer { yyoutput << $$; } IDENTIFIER
%printer { yyoutput << $$; } <val>

%destructor { free($$); } IDENTIFIER

%%
%start line;

line: /* expr */
    | expr                       { driver.set_result($1); }

expr: NUMBER                     { $$ = $1; }
    | IDENTIFIER                 { $$ = driver.get_symbol_value($1); free($1); }
    | IDENTIFIER '(' expr ')'    { $$ = (driver.get_symbol_function($1))($3); free($1); }
    | expr '+' expr              { $$ = $1 + $3; }
    | expr '-' expr              { $$ = $1 - $3; }
    | expr '*' expr              { $$ = $1 * $3; }
    | expr '/' expr              { $$ = $1 / $3; }
    | '+' expr %prec PLS         { $$ = $2; }
    | '-' expr %prec NEG         { $$ = -$2; }
    | expr '^' expr              { $$ = pow($1, $3); }
    | '(' expr ')'               { $$ = $2; }

%%

void
yy::parser::error(const yy::parser::location_type& l, const std::string& m)
{
  driver.error(m);
}
