/*!
         @file    main.cpp

         @brief

         @author  Hideo Matsufuru (matsufuru)
                  $LastChangedBy: aoyama $

         @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

         @version $LastChangedRevision: 1929 $
*/

#include "bridge_setup.h"
#ifdef DEBUG
#include "bridge_init_factory.h"
#endif
#include "Parameters/commonParameters.h"
#include "Parameters/parameters.h"
#include "Parameters/parameterManager_YAML.h"
#include "Tools/timer.h"

#ifdef USE_TESTMANAGER
#include "run_testmanager.h"
#endif

#include "IO/bridgeIO.h"
using Bridge::vout;


const std::string filename_main_input = "main.yaml";
// const std::string filename_main_input = "stdin";

//- prototype declarations
int run_test();


//====================================================================
int main(int argc, char *argv[])
{
  // ###  initial setup  ###
  Bridge::VerboseLevel vl = Bridge::GENERAL;

  //- initialization step one
  bridge_initialize(&argc, &argv);

  std::string filename_input = filename_main_input;
  if (filename_input == "stdin") {
    vout.general(vl, "input filename : ");
    std::cin >> filename_input;
    vout.general(vl, "%s\n", filename_input.c_str());
  } else {
    vout.general(vl, "input filename : %s\n", filename_input.c_str());
  }
  vout.general(vl, "\n");

  //- load input parameters
  Parameters params_all  = ParameterManager::read(filename_input);
  Parameters params_main = params_all.lookup("Main");

  //- initialization step two: setup using parameter values,
  bridge_setup(params_main);

#ifdef USE_FACTORY
#ifdef DEBUG
//  bridge_report_factory();
#endif
#endif

  //- start timer
  unique_ptr<Timer> timer(new Timer("Main"));
  timer->start();

#ifdef USE_TESTMANAGER
  run_testmanager(argc, argv);
#else
  run_test();
#endif

  //- find total elapsed time
  timer->report();

  //- finalization step
  bridge_finalize();

  return EXIT_SUCCESS;
}


//====================================================================
//============================================================END=====
