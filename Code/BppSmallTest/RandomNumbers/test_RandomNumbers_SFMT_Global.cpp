/*!
        @file    test_RandomNumbers_SFMT_Global.cpp

        @brief

        @author  T.Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BppSmallTest.h"
#ifdef USE_NTLLIB

#include "test.h"

#include "Tools/randomNumbers_SFMT.h"

#include "Field/field_G.h"

#include "IO/fieldIO_Text.h"  // for gather field
#include "IO/io_format_gauge.h"

//====================================================================
//! Test of random number generator.

/*!
                                [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
 */

namespace Test_RandomNumbers_SFMT {
  const std::string test_name = "RandomNumbers.SFMT.Global";

  //- test-private parameters
  namespace {
    const std::string filename_input = "test_RandomNumbers_SFMT_Global.yaml";
  }

  //- prototype declaration
  int test_global(void);

#ifdef USE_TESTMANAGER_AUTOREGISTER
  namespace {
    static const bool is_registered = TestManager::RegisterTest(
      test_name,
      test_global
      );
  }
#endif

  //====================================================================
  int test_global(void)
  {
    // ####  parameter setup  ####

    const Parameters params_all = ParameterManager::read(filename_input);

    const Parameters params_test = params_all.lookup("Test_RandomNumbers");

    const int    iseed      = params_test.get_int("int_seed");
    const string str_vlevel = params_test.get_string("verbose_level");

    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  iseed     = %d\n", iseed);
    vout.general(vl, "  vlevel    = %s\n", str_vlevel.c_str());
    vout.general(vl, "\n");


    // #### object setup #####
    const unique_ptr<Timer> timer(new Timer(test_name));


    // ####  Execution main part  ####
    timer->start();

    vout.general(vl, "\n");
    vout.general(vl, "Serial and Node-parallel test of Random Number Generator:\n");

    // sample field size
    const int nin = 8;
    const int nex = 4;

    const int    nvol = CommonParameters::Nvol();
    const long_t lvol = CommonParameters::Lvol();

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
      const unique_ptr<RandomNumbers> rand(new RandomNumbers_SFMT(iseed));
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
      const unique_ptr<RandomNumbers> rand(new RandomNumbers_SFMT(iseed));

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
        vout.crucial(vl, "Error at %s: field size mismatch.\n", test_name.c_str());
        exit(EXIT_FAILURE);
      }

      double *p1 = field1b.ptr(0);
      double *p2 = field2.ptr(0);

      for (size_t i = 0, n = field2.size(); i < n; ++i) {
        if (*p1++ != *p2++) ++err1;
      }
    }

    Communicator::broadcast(1, &err1, 0);
    vout.general(vl, "%s: serial and parallel: err = %d\n", test_name.c_str(), err1);

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
      vout.general(vl, "%s: check local state at rank %d, err = %d\n", test_name.c_str(), ipe, err2);
    }

    // 6. check save and restore
    int err3 = 0;

    if (true) {
      const unique_ptr<RandomNumbers_SFMT> rand(new RandomNumbers_SFMT(iseed));
      rand->set_parameter_verboselevel(Bridge::DETAILED);

      // save current state to file
      rand->write_file("RNG_SFMT.state");

      vout.detailed("%s: number of samples = %lu\n", test_name.c_str(), nsample);

      for (size_t i = 0; i < nsample; ++i) {
        data[i] = rand->get();
      }

      // restore state from file
      rand->read_file("RNG_SFMT.state");

      // check if same series are generated.
      int err3_part = 0;

      for (size_t i = 0; i < nsample; ++i) {
        if (data[i] != rand->get()) ++err3_part;
      }

      Communicator::reduce_sum(1, &err3, &err3_part);

      vout.general(vl, "%s: save/restore test: err = %d\n", test_name.c_str(), err3);
    }

    // 7. check gaussian field
    field1.reset(nin, nvol, nex);

    if (true) {
      const unique_ptr<RandomNumbers> rand(new RandomNumbers_SFMT(iseed));
      rand->set_parameter_verboselevel(Bridge::DETAILED);

      // fill field with gaussian random numbers
      rand->gauss_lex_global(field1);
    }

    if (Communicator::is_primary()) {
      const unique_ptr<RandomNumbers> rand(new RandomNumbers_SFMT(iseed));

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
      } else {
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
        vout.crucial(vl, "Error at %s: field size mismatch.\n", test_name.c_str());
        exit(EXIT_FAILURE);
      }

      double *p1 = field1b.ptr(0);
      double *p2 = field2.ptr(0);

      for (size_t i = 0, n = field2.size(); i < n; ++i) {
        if (*p1++ != *p2++) ++err4;
      }
    }

    Communicator::broadcast(1, &err4, 0);
    vout.general(vl, "%s: serial and parallel for gaussian: err = %d\n", test_name.c_str(), err4);


    // 8. summary
    const double result = err1 + err2 + err3 + err4;

    timer->report();


    if (do_check) {
      return Test::verify(result, expected_result);
    } else {
      vout.detailed(vl, "check skipped: expected_result not set.\n\n");
      return EXIT_SKIP;
    }
  }
} // namespace Test_RandomNumbers
#endif // USE_NTLLIB
