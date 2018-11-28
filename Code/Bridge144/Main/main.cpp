/*!
         @file    $Id:: main.cpp #$

         @brief

         @author  Hideo Matsufuru (matsufuru)
                  $LastChangedBy: aoym $

         @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

         @version $LastChangedRevision: 1561 $
*/

#include "Communicator/communicator.h"
#include "ResourceManager/threadManager_OpenMP.h"

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

#ifdef USE_TESTMANAGER
  preprocess_testmanager(argc, argv);
#endif

  Communicator::init(&argc, &argv);

  // ####  banner  ####
  vout.general(vl, "Bridge++ %s\n\n", BRIDGE_VERSION);

  std::string filename_input = filename_main_input;
  if (filename_input == "stdin") {
    vout.general(vl, "input filename : ");
    std::cin >> filename_input;
    vout.general(vl, "%s\n", filename_input.c_str());
  } else {
    vout.general(vl, "input filename : %s\n", filename_input.c_str());
  }
  vout.general(vl, "\n");

  Parameters params_all  = ParameterManager::read(filename_input);
  Parameters params_main = params_all.lookup("Main");

  const std::vector<int> lattice_size     = params_main.get_int_vector("lattice_size");
  const std::vector<int> grid_size        = params_main.get_int_vector("grid_size");
  const int              number_of_thread = params_main.get_int("number_of_thread");
  const int              number_of_color  = params_main.get_int("number_of_color");
  const std::string      str_logfile      = params_main.get_string("log_filename");
  const std::string      str_ildg_logfile = params_main.get_string("ildg_log_filename");
  const std::string      str_vlevel       = params_main.get_string("verbose_level");


  //- initializations
  vl = vout.set_verbose_level(str_vlevel);
  CommonParameters::init_Vlevel(vl);

  if (str_logfile != "stdout") {
    vout.init(str_logfile);
  }

  if (str_ildg_logfile != "stdout") {
    vout.ildg_init(str_ildg_logfile);
  }


  // CommonParameters::init(lattice_size, grid_size);
  CommonParameters::init(lattice_size, grid_size, number_of_color);
  Communicator::setup();

  ThreadManager_OpenMP::init(number_of_thread);


  //- print input parameters
  vout.general(vl, "Main: input parameters\n");
  vout.general(vl, "  lattice_size     = %s\n", Parameters::to_string(lattice_size).c_str());
  vout.general(vl, "  grid_size        = %s\n", Parameters::to_string(grid_size).c_str());

  vout.general(vl, "  number of thread = %d\n", number_of_thread);
  vout.general(vl, "  number of color  = %d\n", number_of_color);
  vout.general(vl, "  logfile          = %s\n", str_logfile.c_str());
  vout.general(vl, "  ildg_logfile     = %s\n", str_ildg_logfile.c_str());
  vout.general(vl, "  vlevel           = %s\n", str_vlevel.c_str());
  vout.general(vl, "\n");

  //- input parameter check
  int err = 0;
  err += ParameterCheck::non_NULL(str_logfile);
  err += ParameterCheck::non_NULL(str_ildg_logfile);

  if (err) {
    vout.crucial(vl, "Error at main: input parameters have not been set.\n");
    exit(EXIT_FAILURE);
  }


  //- timestamp (starting time)
  unique_ptr<Timer> timer(new Timer("Main"));
  timer->timestamp();
  timer->start();

#ifdef USE_TESTMANAGER
  run_testmanager(argc, argv);
#else
  run_test();
#endif

  //- timestamp (end time)
  timer->report();
  timer->timestamp();

  ThreadManager_OpenMP::finalize();
  Communicator::finalize();

  return EXIT_SUCCESS;
}


//====================================================================
//============================================================END=====
