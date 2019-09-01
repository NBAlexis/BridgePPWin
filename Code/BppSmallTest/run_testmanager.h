/*!
         @file    run_testmanager.h

         @brief

         @author  Tatsumi Aoyama (aoym)
                  $LastChangedBy: aoyama $

         @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

         @version $LastChangedRevision: 1929 $
*/

#ifndef RUN_TESTMANAGER_INCLUDED
#define RUN_TESTMANAGER_INCLUDED

#ifdef USE_TESTMANAGER

#include "testManager.h"

int run_testmanager(int argc, char **argv);

// show list of registered tests
// to run prior to initializations of manager classes such as communicators
int preprocess_testmanager(int argc, char** argv);

#endif
#endif  /* RUN_TESTMANAGER_INCLUDED */
