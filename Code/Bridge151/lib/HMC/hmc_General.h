/*!
        @file    hmc_General.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef HMC_GENERAL_INCLUDED
#define HMC_GENERAL_INCLUDED

#include "action_list.h"
#include "integrator.h"
#include "langevin_Momentum.h"

#include "Measurements/Gauge/staple_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! General HMC update class BAPI.

/*
    This class BAPI defines Hybrid Monte Carlo update algorithm
    with a given molecular dynamics integrator.
    To be improved:
    - at present setting conjugate momenta explicitly assumes
      SU(3) gauge group. This should be generalized, and for
      SU(3) case assert is necessary.
                                   [25 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.           [14 Nov 2012 Y.Namekawa]
    Langevin step was moved to separate file.
                                   [07 May 2014 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                   [21 Mar 2015 Y.Namekawa]
 */


class BAPI HMC_General
{
 public:
  static const std::string class_name;

 private:
  //int m_Metropolis_test;                 //!< Metropolis test: Metropolis_test=0: no test, !=0: test
  bool m_Metropolis_test;                //!< Metropolis test: enabled if true.
  std::vector<Action *> m_action;        //!< actions
  std::vector<Director *> m_director;    //!< directors
  Integrator *m_integrator;              //!< MD integrator
  RandomNumbers *m_rand;                 //!< random number generator
  Staple_lex *m_staple;
  Langevin_Momentum *m_Langevin_P;
  double m_trajectory_length;

  Bridge::VerboseLevel m_vl;

 public:
  //! constructor with action_list, directors, and random number generator
  HMC_General(const ActionList& action_list,
              std::vector<Director *> director,
              Integrator *integrator,
              RandomNumbers *rand);

  HMC_General(const ActionList& action_list,
              std::vector<Director *> director,
              Integrator *integrator,
              unique_ptr<RandomNumbers>& rand);

  //! constructor with action_list when no director is necessary
  HMC_General(const ActionList& action_list,
              Integrator *integrator,
              RandomNumbers *rand);

  HMC_General(const ActionList& action_list,
              Integrator *integrator,
              unique_ptr<RandomNumbers>& rand);

  //! destructor
  ~HMC_General();

 private:
  // non-copyable
  HMC_General(const HMC_General&);
  HMC_General& operator=(const HMC_General&);

 public:
  void set_parameters(const Parameters& params);
  void set_parameters(const double trajectory_length, const bool Metropolis_test);
  void set_parameters(const double trajectory_length, const int Metropolis_test); // backward compatibility

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  double update(Field_G&);

  double langevin(Field_G& iP, const Field_G& U);

  //  double langevin_P(Field_G& iP);  removed [07 May 2014]

  double calc_Hamiltonian(const Field_G& iP, const Field_G& U);
  double calcH_P(const Field_G& iP);
};
#endif
