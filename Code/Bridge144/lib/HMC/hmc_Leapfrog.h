/*!
        @file    $Id:: hmc_Leapfrog.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef HMC_LEAPFROG_INCLUDED
#define HMC_LEAPFROG_INCLUDED

#include "action_list.h"
#include "integrator.h"
#include "langevin_Momentum.h"

#include "Measurements/Gauge/staple_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HMC with single level leapfrog intetgrator.

/*!
    This class implements standartd HMC with simple leapfrog
    molecular dynamics integrator.
    While more general integrator is now available, this class
    is easy to understand and to convenient for first test,
    and thus kept as it is.
                                     [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [03 Mar 2013 Y.Namekawa]
    Langevin step was moved to separate class.
                                     [07 May 2014 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                     [21 Mar 2015 Y.Namekawa]
*/


class HMC_Leapfrog
{
 public:
  static const std::string class_name;

 private:
  int    m_Nmdc;
  int    m_Nprec;
  int    m_Metropolis_test;
  double m_Estep;

  Langevin_Momentum *m_Langevin_P;
  Staple_lex        *m_staple;

  std::vector<Action *>   m_action;
  std::vector<Director *> m_director;
  RandomNumbers           *m_rand;
  Bridge::VerboseLevel    m_vl;

 public:

  //! constructor: with array of actions
  HMC_Leapfrog(std::vector<Action *> action,
               RandomNumbers *rand);

  HMC_Leapfrog(std::vector<Action *> action,
               unique_ptr<RandomNumbers>& rand);

  //! constructor: with array of actions and directors
  HMC_Leapfrog(std::vector<Action *> action,
               std::vector<Director *> director,
               RandomNumbers *rand);

  HMC_Leapfrog(std::vector<Action *> action,
               std::vector<Director *> director,
               unique_ptr<RandomNumbers>& rand);

  //! constructor: with action_list
  HMC_Leapfrog(const ActionList& action_list,
               RandomNumbers *rand);

  HMC_Leapfrog(const ActionList& action_list,
               unique_ptr<RandomNumbers>& rand);

  //! constructor: with action_list and array of directors
  HMC_Leapfrog(const ActionList& action_list,
               std::vector<Director *> director,
               RandomNumbers *rand);

  HMC_Leapfrog(const ActionList& action_list,
               std::vector<Director *> director,
               unique_ptr<RandomNumbers>& rand);

  //! destructor
  ~HMC_Leapfrog();

 private:
  // non-copyable
  HMC_Leapfrog(const HMC_Leapfrog&);
  HMC_Leapfrog& operator=(const HMC_Leapfrog&);

 public:
  void set_parameters(const Parameters& params);
  void set_parameters(double Estep, int Nmdc, int Nprec, int Metropolis_test);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  //  Metropolis_test=0: no test, !=0: test
  double update(Field_G&);

  double langevin(Field_G& iP, Field_G& U);

  //  double langevin_P(Field_G& iP);  removed. [07 May 2014]

  double calc_Hamiltonian(Field_G& iP, Field_G& U);
  double calcH_P(Field_G& iP);

  void integrate(Field_G& iP, Field_G& U);

  void update_U(double estep, Field_G& iP, Field_G& U);
  void update_P(double estep, Field_G& iP, Field_G& U);
};
#endif
