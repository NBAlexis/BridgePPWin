/*!
        @file    $Id:: evalexpr_global.h #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#ifndef EVALEXPR_GLOBAL_INCLUDED
#define EVALEXPR_GLOBAL_INCLUDED

#include <cmath>

// global symbol definitions.
namespace {
  const _init_function arithmetic_functions[] =
  {
    { "cos",  cos,  },
    { "sin",  sin,  },
    { "tan",  tan,  },
    { "atan", atan, },
    { "log",  log,  },
    { "exp",  exp,  },
    { "sqrt", sqrt, },
    {      0,    0, }, // sentinel
  };

  const _init_variable predefined_constants[] =
  {
    { "pi", atan(1.0) * 4, },
    {    0,             0, }, // sentinel
  };
}
#endif
