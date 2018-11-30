//=============================================================================
// FILENAME : SpectrumTest.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [11/29/2018 nbale]
//=============================================================================

#include "SpectrumTest.h"
#include <omp.h>

const std::string filename_main_input = "SpectrumTest.yaml";

int main(int argc, char * argv[])
{
    // ###  initial setup  ###
    Bridge::VerboseLevel vl = Bridge::GENERAL;
    Communicator::init(&argc, &argv);

    // ####  banner  ####
    Bridge::vout.general(vl, "Bridge++ %s\n\n", BRIDGE_VERSION);

    std::string filename_input = filename_main_input;
    if (filename_input == "stdin") 
    {
        Bridge::vout.general(vl, "input filename : ");
        std::cin >> filename_input;
        Bridge::vout.general(vl, "%s\n", filename_input.c_str());
    }
    else 
    {
        Bridge::vout.general(vl, "input filename : %s\n", filename_input.c_str());
    }
    Bridge::vout.general(vl, "\n");

    Parameters params_all = ParameterManager::read(filename_input);
    Parameters params_main = params_all.lookup("Main");

    const std::vector<int> lattice_size = params_main.get_int_vector("lattice_size");
    const std::vector<int> grid_size = params_main.get_int_vector("grid_size");
    const int              number_of_thread = params_main.get_int("number_of_thread");
    const int              number_of_color = params_main.get_int("number_of_color");
    const std::string      str_logfile = params_main.get_string("log_filename");
    const std::string      str_ildg_logfile = params_main.get_string("ildg_log_filename");
    const std::string      str_vlevel = params_main.get_string("verbose_level");


    //- initializations
    vl = Bridge::vout.set_verbose_level(str_vlevel);
    CommonParameters::init_Vlevel(vl);

    if (str_logfile != "stdout")
    {
        Bridge::vout.init(str_logfile);
    }

    if (str_ildg_logfile != "stdout")
    {
        Bridge::vout.ildg_init(str_ildg_logfile);
    }


    CommonParameters::init(lattice_size, grid_size, number_of_color);
    Communicator::setup();

    ThreadManager_OpenMP::init(number_of_thread);


    //- print input parameters
    Bridge::vout.general(vl, "Main: input parameters\n");
    Bridge::vout.general(vl, "  lattice_size     = %s\n", Parameters::to_string(lattice_size).c_str());
    Bridge::vout.general(vl, "  grid_size        = %s\n", Parameters::to_string(grid_size).c_str()); //seems always 1,1,1,1

    Bridge::vout.general(vl, "  number of thread = %d\n", number_of_thread);
    Bridge::vout.general(vl, "  number of color  = %d\n", number_of_color);
    Bridge::vout.general(vl, "  logfile          = %s\n", str_logfile.c_str());
    Bridge::vout.general(vl, "  ildg_logfile     = %s\n", str_ildg_logfile.c_str());
    Bridge::vout.general(vl, "  vlevel           = %s\n", str_vlevel.c_str());
    Bridge::vout.general(vl, "\n");

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_logfile);
    err += ParameterCheck::non_NULL(str_ildg_logfile);

    if (err)
    {
        vout.crucial(vl, "Error at main: input parameters have not been set.\n");
        std::cout << "press enter to quit";
        std::cin.get();
        exit(EXIT_FAILURE);
    }


    //Run the thing here
    RunJob1();

    ThreadManager_OpenMP::finalize();
    Communicator::finalize();
    std::cout << "press enter to quit";
    std::cin.get();
    return EXIT_SUCCESS;
}

//=============================================================================
// END OF FILE
//=============================================================================