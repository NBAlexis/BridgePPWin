//=============================================================================
// FILENAME : SpectrumTest.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [11/29/2018 nbale]
//=============================================================================

#include "SpectrumTest.h"

bool RunJob1()
{
    const std::string filename_input = "SpectrumTestJob1.yaml";

    // ####  parameter setup  ####
    int Nc = CommonParameters::Nc();
    int Nd = CommonParameters::Nd();
    int Ndim = CommonParameters::Ndim();
    int Nvol = CommonParameters::Nvol();

    Parameters params_all = ParameterManager::read(filename_input);
    Parameters params_test = params_all.lookup("Test_Spectrum");

    const int N_quark = params_test.get_int("number_of_valence_quarks");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read = params_test.get_string("gauge_config_type_input");
    const string        str_rand_type = params_test.get_string("random_number_type");
    const unsigned long seed = params_test.get_unsigned_long("seed_for_random_number");
    const string        str_vlevel = params_test.get_string("verbose_level");

    const bool   do_check = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;


    const string str_gfix_type = params_all.lookup("GaugeFixing").get_string("gauge_fixing_type");
    const string str_proj_type = params_all.lookup("Projection").get_string("projection_type");
    const string str_smear_type = params_all.lookup("Smear").get_string("smear_type");
    const string str_solver_type = params_all.lookup("Solver").get_string("solver_type");

    Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    //- print input parameters
    vout.general(vl, "  gconf_status = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read   = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  rand_type    = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed         = %lu\n", seed);
    vout.general(vl, "  vlevel       = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gfix_type    = %s\n", str_gfix_type.c_str());
    vout.general(vl, "  proj_type    = %s\n", str_proj_type.c_str());
    vout.general(vl, "  smear_type   = %s\n", str_smear_type.c_str());
    vout.general(vl, "  solver_type  = %s\n", str_solver_type.c_str());

    vector<Parameters> params_quark;

    char tmp[128];
    for (int iq = 0; iq < N_quark; ++iq) {
        memset(tmp, sizeof(char) * 128, 0);
        _itoa_s(iq + 1, tmp, 10);
        string qlabel = "Quark_" + string(tmp);
        params_quark.push_back(params_all.lookup(qlabel));
    }

    // NB. all str_gmset_type are supposed to be the same.
    string str_gmset_type = params_quark[0].lookup("Fopr").get_string("gamma_matrix_type");
    vout.general(vl, "  gmset_type  = %s\n", str_gmset_type.c_str());


    for (int iq = 0; iq < N_quark; ++iq) {
        vout.general(vl, "  Quark_%d:\n", iq + 1);
        vout.general(vl, "    source_type = %s\n",
            params_quark[iq].lookup("Source").get_string("source_type").c_str());

        vector<int> pos = params_quark[iq].lookup("Source").get_int_vector("source_position");

        for (int mu = 0; mu < Ndim; ++mu) {
            vout.general(vl, "    source_position[%d] = %d\n", mu, pos[mu]);
        }
    }
    vout.general(vl, "\n");


    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);

    if (err) {
        vout.crucial(vl, "Error !!!\n");
        exit(EXIT_FAILURE);
    }

    if (str_solver_type == "CG") {
        vout.crucial(vl, "Error !!! CG Cannot adopted \n");
        exit(EXIT_FAILURE);
    }

    RandomNumberManager::initialize(str_rand_type, seed);


    // ####  Set up a gauge configuration  ####
    unique_ptr<Field_G> U(new Field_G(Nvol, Ndim));

    if (str_gconf_status == "Hot_start") 
    {
        GaugeConfig("Random").read(U);
    }
    else 
    {
        vout.crucial(vl, "Error !!!: unsupported gconf status \"%s\"\n", str_gconf_status.c_str());
        exit(EXIT_FAILURE);
    }

    unique_ptr<Projection> proj(Projection::New(str_proj_type));
    proj->set_parameters(params_all.lookup("Projection"));

    unique_ptr<Smear> smear(Smear::New(str_smear_type, proj));
    smear->set_parameters(params_all.lookup("Smear"));

    unique_ptr<Director_Smear> dr_smear(new Director_Smear(smear));
    dr_smear->set_parameters(params_all.lookup("Director_Smear"));
    dr_smear->set_config(U);

    int     Nsmear = dr_smear->get_Nsmear();
    Field_G *Usmear = (Field_G *)dr_smear->getptr_smearedConfig(Nsmear);


    // ####  Gauge fixing  ####
    unique_ptr<Field_G> Ufix(new Field_G(Nvol, Ndim));

    unique_ptr<GaugeFixing> gfix(GaugeFixing::New(str_gfix_type));
    gfix->set_parameters(params_all.lookup("GaugeFixing"));

    gfix->fix(*Ufix, *Usmear);

    copy(*U, *Ufix);


    // ####  object setup  #####
    unique_ptr<GammaMatrixSet> gmset(GammaMatrixSet::New(str_gmset_type));

    std::vector<unique_ptr<Fopr> >   fopr(N_quark);
    std::vector<unique_ptr<Solver> > solver(N_quark);
    std::vector<unique_ptr<Fprop> >  fprop_lex(N_quark);
    std::vector<unique_ptr<Source> > source(N_quark);

    for (int iq = 0; iq < N_quark; ++iq) {
        string str_fopr_type = params_quark[iq].lookup("Fopr").get_string("fermion_type");

        fopr[iq].reset(Fopr::New(str_fopr_type, str_gmset_type));
        fopr[iq]->set_parameters(params_quark[iq].lookup("Fopr"));
        fopr[iq]->set_config(U);

        solver[iq].reset(Solver::New(str_solver_type, fopr[iq]));
        solver[iq]->set_parameters(params_all.lookup("Solver"));

        fprop_lex[iq].reset(new Fprop_Standard_lex(solver[iq]));

        string str_source_type = params_quark[iq].lookup("Source").get_string("source_type");
        source[iq].reset(Source::New(str_source_type));
        source[iq]->set_parameters(params_quark[iq].lookup("Source"));
    }

    Corr2pt_4spinor corr(gmset);
    corr.set_parameters(params_all.lookup("Corr2pt_4spinor"));


    // ####  Execution main part  ####
    typedef std::vector<Field_F> PropagatorSet;

    std::vector<PropagatorSet> sq(N_quark);
    for (int iq = 0; iq < N_quark; ++iq) {
        sq[iq].resize(Nc * Nd);

        for (int i = 0; i < Nc * Nd; ++i) {
            sq[iq][i].set(0.0);
        }
    }

    Field_F b;
    b.set(0.0);

    int    Nconv;
    double diff, diff2;

    for (int iq = 0; iq < N_quark; ++iq) {
        vout.general(vl, "Solving quark propagator, flavor = %d:\n", iq + 1);
        vout.general(vl, "  color spin   Nconv      diff           diff2\n");

        for (int ispin = 0; ispin < Nd; ++ispin) {
            for (int icolor = 0; icolor < Nc; ++icolor) {
                int idx = icolor + Nc * ispin;
                source[iq]->set(b, idx);

                fprop_lex[iq]->invert_D(sq[iq][idx], b, Nconv, diff);

                Field_F y(b);
                fopr[iq]->set_mode("D");
                fopr[iq]->mult(y, sq[iq][idx]);  // y  = fopr[iq]->mult(sq[iq][idx]);
                axpy(y, -1.0, b);                // y -= b;
                diff2 = y.norm2() / b.norm2();

                vout.general(vl, "   %2d   %2d   %6d   %12.4e   %12.4e\n",
                    icolor, ispin, Nconv, diff, diff2);
            }
        }

        vout.general(vl, "\n");
    }


    //- meson correlators
    std::vector<double> result(N_quark);

    vout.general(vl, "2-point correlator:\n");

    //- case(iq_1 == iq_2)
    for (int iq = 0; iq < N_quark; ++iq) {
        vout.general(vl, "Flavor combination = %d, %d\n", iq + 1, iq + 1);
        result[iq] = corr.meson_all(sq[iq], sq[iq]);
        vout.general(vl, "\n");
    }

    //- case(iq_1 < iq_2)
    for (int iq = 0; iq < N_quark; ++iq) {
        for (int jq = iq + 1; jq < N_quark; ++jq) {
            vout.general(vl, "Flavor combination = %d, %d\n", iq + 1, jq + 1);
            double result_2 = corr.meson_all(sq[iq], sq[jq]);
            vout.general(vl, "\n");
        }
    }

    RandomNumberManager::finalize();
}

//=============================================================================
// END OF FILE
//=============================================================================