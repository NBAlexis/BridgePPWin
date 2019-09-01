/*!
        @file    bridge_setup.cpp

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"

//#include "bridge_setup.h"
//#include "Parameters/commonParameters.h"
//#include "Communicator/communicator.h"
//#include "ResourceManager/threadManager_OpenMP.h"
//#include "Tools/timer.h"
//#include "IO/bridgeIO.h"
//#include "bridge_init_factory.h"

#ifdef USE_GROUP_SU3
#define Nc    3
#else
#ifdef USE_GROUP_SU2
#define Nc    2
#else
#ifdef USE_GROUP_GENERAL
#define Nc    3
#endif
#endif
#endif

int bridge_initialize(int *pargc, char ***pargv)
{
  // initialize communicator
  Communicator::init(pargc, pargv);

  // initialize thread manager
  //-

  // setup factory
#if defined(USE_FACTORY) && !defined(USE_FACTORY_AUTOREGISTER)
  bridge_init_factory();
#endif

#if defined(DEBUG) && defined(USE_FACTORY)
  // bridge_report_factory();
#endif

  // show banner
  vout.general("Bridge++ %s\n\n", BRIDGE_VERSION);

  Timer::timestamp();

  return EXIT_SUCCESS;
}


int bridge_initialize(int *pargc, char ***pargv,
                      const Parameters& params
                      )
{
  bridge_initialize(pargc, pargv);

  bridge_setup(params);

  return EXIT_SUCCESS;
}


int bridge_initialize(int *pargc, char ***pargv,
                      const std::vector<int>& lattice_size,
                      const std::vector<int>& grid_size,
                      const int number_of_threads,
                      const int number_of_colors,
                      const std::string& logfile,
                      const std::string& ildg_logfile,
                      const std::string& verbose_level
                      )
{
  bridge_initialize(pargc, pargv);

  bridge_setup(lattice_size,
               grid_size,
               number_of_threads,
               number_of_colors,
               logfile,
               ildg_logfile,
               verbose_level);

  return EXIT_SUCCESS;
}


int bridge_finalize()
{
  Timer::timestamp();

  ThreadManager_OpenMP::finalize();

  Communicator::finalize();

  return EXIT_SUCCESS;
}


void bridge_setup(const Parameters& params)
{
  std::vector<int> lattice_size;
  std::vector<int> grid_size         = std::vector<int>(); // will be auto-shaped.
  int              number_of_threads = 0;                  // will be auto-adjusted.
  int              number_of_colors  = Nc;
  std::string      logfile           = "stdout";
  std::string      ildg_logfile      = "stdout";
  std::string      verbose_level     = "General";

  int ret = 0;

  ret += params.fetch_int_vector("lattice_size", lattice_size);
  ret += params.fetch_int_vector("grid_size", grid_size);
  ret += params.fetch_int("number_of_thread", number_of_threads);
  ret += params.fetch_int("number_of_color", number_of_colors);
  ret += params.fetch_string("log_filename", logfile);
  ret += params.fetch_string("ildg_log_filename", ildg_logfile);
  ret += params.fetch_string("verbose_level", verbose_level);

  bridge_setup(lattice_size,
               grid_size,
               number_of_threads,
               number_of_colors,
               logfile,
               ildg_logfile,
               verbose_level
               );
}


void bridge_setup(
  const std::vector<int>& lattice_size,
  const std::vector<int>& grid_size_hint,
  const int number_of_threads,
  const int number_of_colors,
  const std::string& logfile,
  const std::string& ildg_logfile,
  const std::string& verbose_level
  )
{
  // report
  vout.general("Main: input parameters\n");
  vout.general("  lattice_size     = %s\n", Parameters::to_string(lattice_size).c_str());
  vout.general("  grid_size        = %s\n", Parameters::to_string(grid_size_hint).c_str());
  vout.general("  number of thread = %d\n", number_of_threads);
  vout.general("  number of color  = %d\n", number_of_colors);
  vout.general("  logfile          = %s\n", logfile.c_str());
  vout.general("  ildg_logfile     = %s\n", ildg_logfile.c_str());
  vout.general("  vlevel           = %s\n", verbose_level.c_str());


  // setup communicator
  // grid may be auto-shaped.

  std::vector<int> grid_size(grid_size_hint);

  Communicator::setup(lattice_size, grid_size);

  // store to global parameters
  CommonParameters::init(lattice_size, grid_size, number_of_colors);

  // setup thread mananger
  ThreadManager_OpenMP::init(number_of_threads);

  // set log parameters
  Bridge::VerboseLevel vl = Bridge::BridgeIO::set_verbose_level(verbose_level);
  CommonParameters::init_Vlevel(vl);

  if ((logfile.size() > 0) && (logfile != "stdout")) {
    vout.init(logfile);
  }

  if ((ildg_logfile.size() > 0) && (ildg_logfile != "stdout")) {
    vout.ildg_init(ildg_logfile);
  }
}


#ifdef Nc
#undef Nc
#endif
