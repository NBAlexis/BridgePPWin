//=============================================================================
// FILENAME : test_RandomNumbers_Schrage.h
// 
// DESCRIPTION:
//
// REVISION:
//  [11/30/2018 nbale]
//=============================================================================

#include "BppSmallTest.h"

//====================================================================
//! Test of random number generator.

namespace Test_RandomNumbers_Schrage
{
    const std::string test_name_uniform = "RandomNumbers.Schrage.Uniform";
    const std::string test_name_gaussian = "RandomNumbers.Schrage.Gaussian";
    const std::string test_name_gaussianfield = "RandomNumbers.Schrage.GaussianField";
    const std::string test_name_gobal = "RandomNumbers.Schrage.Global";

    //- test-private parameters
    namespace
    {
        const std::string filename_input_uniform = "test_RandomNumbers_Schrage_Uniform.yaml";
        const std::string filename_input_gaussian = "test_RandomNumbers_Schrage_Gaussian.yaml";
        const std::string filename_input_gaussianfield = "test_RandomNumbers_Schrage_GaussianField.yaml";
        const std::string filename_input_global = "test_RandomNumbers_Schrage_Global.yaml";
    }

    //- prototype declaration
    int uniform_calc_pi(void);
    int gaussian(void);
    int gaussian_field(void);
    int test_global(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
    namespace {
        static const bool is_registered1 = TestManager::RegisterTest(
            test_name_uniform,
            uniform_calc_pi
        );

        static const bool is_registered2 = TestManager::RegisterTest(
            test_name_gaussian,
            gaussian
        );

        static const bool is_registered3 = TestManager::RegisterTest(
            test_name_gaussianfield,
            gaussian_field
        );

        static const bool is_registered4 = TestManager::RegisterTest(
            test_name_gobal,
            test_global
        );
    }
#endif

    //====================================================================
    int uniform_calc_pi(void)
    {
        // ####  parameter setup  ####

        Parameters params_all = ParameterManager::read(filename_input_uniform);

        Parameters params_test = params_all.lookup("Test_RandomNumbers");

        int          Nseed = params_test.get_int("number_of_seeds");
        int          iseed_base = params_test.get_int("int_seed_base");
        int          Nrand = params_test.get_int("number_of_samples");
        const string str_vlevel = params_test.get_string("verbose_level");

        const bool   do_check = params_test.is_set("expected_result");
        const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

        Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

        //- print input parameters
        vout.general(vl, "  Nseed      = %d\n", Nseed);
        vout.general(vl, "  iseed_base = %d\n", iseed_base);
        vout.general(vl, "  Nrand      = %d\n", Nrand);
        vout.general(vl, "  vlevel     = %s\n", str_vlevel.c_str());
        vout.general(vl, "\n");


        // #### object setup #####
        unique_ptr<Timer> timer(new Timer(test_name_uniform));


        // ####  Execution main part  ####
        timer->start();

        vout.general(vl, "\n");
        vout.general(vl, "Monte Carlo estimate of pi:\n");
        vout.general(vl, "  number of samples = %10d\n", Nrand);
        vout.general(vl, "        seed    estimate of pi\n");

        double t1 = 0.0;
        double t2 = 0.0;
        for (int iseed = 0; iseed < Nseed; ++iseed) 
        {
            int iseed2 = iseed_base + iseed;

            unique_ptr<RandomNumbers> rand(new RandomNumbers_Schrage(iseed2));

            int Npi = 0;
            for (int i = 0; i < Nrand; ++i) {
                double rand1 = rand->get();
                double rand2 = rand->get();
                double r = rand1 * rand1 + rand2 * rand2;
                if (r < 1.0) { ++Npi; }
                //  vout.general(vl, "  %10.8f  %10.8f\n",rand1,rand2);
            }

            double pi_exp = (4.0 * Npi) / Nrand;

            t1 += pi_exp;
            t2 += pi_exp * pi_exp;

            //vout.general(vl, "  estimate of pi    = %10.8f\n",pi_exp);
            vout.general(vl, "  %10d    %14.10f\n", iseed2, pi_exp);
        }

        double api = t1 / (double)Nseed;
        double vpi = t2 / (double)Nseed - api * api;
        double dpi = sqrt(vpi);
        double epi = dpi / sqrt((double)Nseed - 1);

        double pi = 3.141592653589793;
        vout.general(vl, "  true value = %10.8f\n", pi);
        vout.general(vl, "  average    = %10.8f\n", api);
        vout.general(vl, "  variance   = %10.8f\n", vpi);
        vout.general(vl, "  deviation  = %10.8f\n", dpi);
        vout.general(vl, "  error      = %10.8f\n", epi);

        double result = api;
        // changed to check the obtained value of pi itseft. [25 May 2014 H.M.]

        timer->report();


        if (do_check)
        {
            return Test::verify(result, expected_result);
        }
        else
        {
            vout.detailed(vl, "check skipped: expected_result not set.\n\n");
            return EXIT_SKIP;
        }
    }

    int gaussian(void)
    {
        // ####  parameter setup  ####

        Parameters params_all = ParameterManager::read(filename_input_gaussian);

        Parameters params_test = params_all.lookup("Test_RandomNumbers");

        int          iseed = params_test.get_int("int_seed");
        int          Nrand = params_test.get_int("number_of_samples");
        const string str_vlevel = params_test.get_string("verbose_level");

        const bool   do_check = params_test.is_set("expected_result");
        const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

        Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

        //- print input parameters
        vout.general(vl, "  iseed  = %d\n", iseed);
        vout.general(vl, "  Nrand  = %d\n", Nrand);
        vout.general(vl, "  vlevel = %s\n", str_vlevel.c_str());
        vout.general(vl, "\n");


        // ####  object setup  ####
        unique_ptr<RandomNumbers> rand(new RandomNumbers_Schrage(iseed));
        unique_ptr<Timer>         timer(new Timer(test_name_gaussian));


        // ####  Execution main part  ####
        timer->start();

        double av = 0.0;
        double vr = 0.0;

        double rand1, rand2;
        for (int i = 0; i < Nrand; ++i) {
            rand->gauss(rand1, rand2);
            av += rand1 + rand2;
            vr += rand1 * rand1 + rand2 * rand2;
            // vout.general(vl, "  %10.8f  %10.8f\n",rand1,rand2);
        }
        av = av / (2.0 * Nrand);
        vr = vr / (2.0 * Nrand) - av * av;
        vr = sqrt(vr);

        vout.general(vl, "\n");
        vout.general(vl, "Gaussian distribution:\n");
        vout.general(vl, "  number of samples = %10d\n", Nrand);
        vout.general(vl, "  average           = %10.8f\n", av);
        vout.general(vl, "  variance          = %10.8f\n", vr);
        vout.general(vl, "  variance(expect)  = %10.8f\n", 1.0 / sqrt(2.0));

        double result = vr;

        timer->report();


        if (do_check) {
            return Test::verify(result, expected_result);
        }
        else {
            vout.detailed(vl, "check skipped: expected_result not set.\n\n");
            return EXIT_SKIP;
        }
    }

    int test_global(void)
    {
        // ####  parameter setup  ####

        Parameters params_all = ParameterManager::read(filename_input_global);

        Parameters params_test = params_all.lookup("Test_RandomNumbers");

        int          iseed = params_test.get_int("int_seed");
        const string str_vlevel = params_test.get_string("verbose_level");

        const bool   do_check = params_test.is_set("expected_result");
        const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

        Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

        //- print input parameters
        vout.general(vl, "  iseed     = %d\n", iseed);
        vout.general(vl, "  vlevel    = %s\n", str_vlevel.c_str());
        vout.general(vl, "\n");


        // #### object setup #####
        unique_ptr<Timer> timer(new Timer(test_name_gobal));


        // ####  Execution main part  ####
        timer->start();

        vout.general(vl, "\n");
        vout.general(vl, "Serial and Node-parallel test of Random Number Generator:\n");

        // sample field size
        int nin = 8;
        int nex = 4;

        int nvol = CommonParameters::Nvol();
        int lvol = CommonParameters::Lvol();

        vout.general(vl, "field size: nin = %d, nex = %d, nvol = %d, lvol = %d\n", nin, nex, nvol, lvol);

        // buffer for checks
        const size_t nsample = 1024;
        double       data[nsample];
        for (size_t i = 0; i < nsample; ++i) {
            data[i] = 0.0;
        }


        // 1. generate field in parallel
        Field field1(nin, nvol, nex);

        if (true) {
            unique_ptr<RandomNumbers> rand(new RandomNumbers_Schrage(iseed));
            rand->set_parameter_verboselevel(Bridge::DETAILED);

            // fill field with uniform random numbers
            rand->uniform_lex_global(field1);

            // generate additional random numbers to check rng state.
            for (size_t i = 0; i < nsample; ++i) {
                data[i] = rand->get();
            }
        }

        // 2. generate field at rank 0 with the same seed
        Field field2(0, 0, 0);

        if (Communicator::is_primary()) {
            unique_ptr<RandomNumbers> rand(new RandomNumbers_Schrage(iseed));

            field2.reset(nin, lvol, nex);

            double *p = field2.ptr(0);

            for (size_t i = 0, n = field2.size(); i < n; ++i) {
                *p++ = rand->get();
            }
        }

        // 3. gather parallel field to rank 0
        Field field1b(0, 0, 0);

        if (Communicator::is_primary()) {
            field1b.reset(nin, lvol, nex);
        }

        FieldIO_Text fieldio(IO_Format::Gauge::ILDG);

        fieldio.gather(&field1b, &field1);

        // 4. compare
        int err1 = 0;

        if (Communicator::is_primary()) {
            if (field1b.size() != field2.size()) {
                vout.crucial(vl, "Error at %s: field size mismatch.\n", test_name_gobal.c_str());
                exit(EXIT_FAILURE);
            }

            double *p1 = field1b.ptr(0);
            double *p2 = field2.ptr(0);

            for (size_t i = 0, n = field2.size(); i < n; ++i) {
                if (*p1++ != *p2++) ++err1;
            }
        }

        Communicator::broadcast(1, &err1, 0);
        vout.general(vl, "%s: serial and parallel: err = %d\n", test_name_gobal.c_str(), err1);

        // 5. check rng state
        int    err2 = 0;
        double buf[nsample];

        for (int ipe = 1, npe = Communicator::size(); ipe < npe; ++ipe) {
            Communicator::send_1to1(nsample, buf, data, 0, ipe, ipe);

            if (Communicator::is_primary()) {
                for (size_t i = 0; i < nsample; ++i) {
                    if (data[i] != buf[i]) ++err2;
                }
            }

            Communicator::broadcast(1, &err2, 0);
            vout.general(vl, "%s: check local state at rank %d, err = %d\n", test_name_gobal.c_str(), ipe, err2);
        }

        // 6. check save and restore
        int err3 = 0;

        if (true) {
            unique_ptr<RandomNumbers_Schrage> rand(new RandomNumbers_Schrage(iseed));
            rand->set_parameter_verboselevel(Bridge::DETAILED);

            // save current state to file
            rand->write_file("RNG_Schrage.state");

            vout.detailed("%s: number of samples = %lu\n", test_name_gobal.c_str(), nsample);

            for (size_t i = 0; i < nsample; ++i) {
                data[i] = rand->get();
            }

            // restore state from file
            rand->read_file("RNG_Schrage.state");

            // check if same series are generated.
            int err3_part = 0;

            for (size_t i = 0; i < nsample; ++i) {
                if (data[i] != rand->get()) ++err3_part;
            }

            Communicator::reduce_sum(1, &err3, &err3_part);

            vout.general(vl, "%s: save/restore test: err = %d\n", test_name_gobal.c_str(), err3);
        }

        // 7. check gaussian field
        field1.reset(nin, nvol, nex);

        if (true) {
            unique_ptr<RandomNumbers> rand(new RandomNumbers_Schrage(iseed));
            rand->set_parameter_verboselevel(Bridge::DETAILED);

            // fill field with gaussian random numbers
            rand->gauss_lex_global(field1);
        }

        if (Communicator::is_primary()) {
            unique_ptr<RandomNumbers> rand(new RandomNumbers_Schrage(iseed));

            field2.reset(nin, lvol, nex);

            double *p = field2.ptr(0);

            if (field2.nin() % 2 == 0) {
                double r1, r2;

                for (int j = 0, Nex = field2.nex(); j < Nex; ++j) {
                    for (int isite = 0, Nvol = field2.nvol(); isite < Nvol; ++isite) {
                        for (int i = 0, Nin = field2.nin(); i < Nin; i += 2) {
                            rand->gauss(r1, r2);
                            *p++ = r1;
                            *p++ = r2;
                        }
                    }
                }
            }
            else {
                double r1, r2;

                for (int j = 0, Nex = field2.nex(); j < Nex; ++j) {
                    for (int isite = 0, Nvol = field2.nvol(); isite < Nvol; ++isite) {
                        for (int i = 0, Nin = field2.nin(); i < Nin; ++i) {
                            rand->gauss(r1, r2);
                            *p++ = r1;
                        }
                    }
                }
            }
        }

        if (Communicator::is_primary()) {
            field1b.reset(nin, lvol, nex);
        }

        fieldio.gather(&field1b, &field1);

        int err4 = 0;

        if (Communicator::is_primary()) {
            if (field1b.size() != field2.size()) {
                vout.crucial(vl, "Error at %s: field size mismatch.\n", test_name_gobal.c_str());
                exit(EXIT_FAILURE);
            }

            double *p1 = field1b.ptr(0);
            double *p2 = field2.ptr(0);

            for (size_t i = 0, n = field2.size(); i < n; ++i) {
                if (*p1++ != *p2++) ++err4;
            }
        }

        Communicator::broadcast(1, &err4, 0);
        vout.general(vl, "%s: serial and parallel for gaussian: err = %d\n", test_name_gobal.c_str(), err4);


        // 8. summary
        double result = err1 + err2 + err3 + err4;

        timer->report();


        if (do_check) {
            return Test::verify(result, expected_result);
        }
        else {
            vout.detailed(vl, "check skipped: expected_result not set.\n\n");
            return EXIT_SKIP;
        }
    }

    int gaussian_field(void)
    {
        // ####  parameter setup  ####

        Parameters params_all = ParameterManager::read(filename_input_gaussianfield);

        Parameters params_test = params_all.lookup("Test_RandomNumbers");

        int          iseed = params_test.get_int("int_seed");
        const string str_vlevel = params_test.get_string("verbose_level");

        const bool   do_check = params_test.is_set("expected_result");
        const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

        Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

        //- print input parameters
        vout.general(vl, "  iseed  = %d\n", iseed);
        vout.general(vl, "  vlevel = %s\n", str_vlevel.c_str());
        vout.general(vl, "\n");


        // ####  object setup  ####
        unique_ptr<RandomNumbers> rand(new RandomNumbers_Schrage(iseed));
        unique_ptr<Timer>         timer(new Timer(test_name_gaussianfield));


        // ####  Execution main part  ####
        timer->start();

        int   Nin = 29;
        int   Nvol = CommonParameters::Nvol();
        int   Nex = 33;
        Field v(Nin, Nvol, Nex);

        double av = 0.0;
        double vr = 0.0;

        rand->gauss_lex_global(v);

        int size = v.size();
        for (int i = 0; i < size; ++i) {
            av += v.cmp(i);
            vr += v.cmp(i) * v.cmp(i);
            // vout.general(vl, "  %10.8f\n",v.cmp(i));
        }

        double av_all = Communicator::reduce_sum(av);
        double vr_all = Communicator::reduce_sum(vr);

        int global_size = CommonParameters::Lvol() * Nin * Nex;

        av = av_all / global_size;
        vr = vr_all / global_size - av * av;
        vr = sqrt(vr);

        vout.general(vl, "\n");
        vout.general(vl, "Gaussian distribution (Field):\n");
        vout.general(vl, "  number of samples = %10d\n", size);
        vout.general(vl, "  average           = %10.8f\n", av);
        vout.general(vl, "  variance          = %10.8f\n", vr);
        vout.general(vl, "  variance(expect)  = %10.8f\n", 1.0 / sqrt(2.0));

        double result = vr;

        timer->report();


        if (do_check) {
            return Test::verify(result, expected_result);
        }
        else {
            vout.detailed(vl, "check skipped: expected_result not set.\n\n");
            return EXIT_SKIP;
        }
    }

} // namespace Test_RandomNumbers
