/*!
         @file    $Id:: run_testmanager.h #$

         @brief

         @author  Tatsumi Aoyama (aoym)
                  $LastChangedBy: aoym $

         @date    $LastChangedDate:: 2017-03-14 06:41:34 #$

         @version $LastChangedRevision: 1593 $
*/

#ifndef RUN_TESTMANAGER_INCLUDED
#define RUN_TESTMANAGER_INCLUDED

#ifdef USE_TESTMANAGER

#include "Tests/testManager.h"

int run_testmanager(int argc, char **argv);

// show list of registered tests
// to run prior to initializations of manager classes such as communicators
int preprocess_testmanager(int argc, char **argv);
#endif
#endif  /* RUN_TESTMANAGER_INCLUDED */
