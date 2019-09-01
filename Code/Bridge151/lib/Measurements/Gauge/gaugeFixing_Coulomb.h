/*!
        @file    gaugeFixing_Coulomb.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef GAUGEFIXING_COULOMB_INCLUDED
#define GAUGEFIXING_COULOMB_INCLUDED

#include "gaugeFixing.h"

#include "Field/shiftField_eo.h"
#include "Tools/randomNumbers.h"
#include "Tools/randomNumberManager.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Coulomb gauge fixing.

/*
    This class BAPI fixes the gauge of configuration to the Coulomb gauge.
    The implementation assumes that the dimension is 4 and the
    Coulomb gauge fixing is performed within each time slice.
    The algorithm is that developed by the Los Alamos group [see the
    implementation note].
    Overrelaxation is incorporated.
    To escape the Gribov copy, if convergence is not reached on some
    timeslices within the iterations specified by Nreset, random
    gauge transformation is performed to reset the configuration on
    that timeslice.
    This is the reason that random number generator is needed at the
    construction of this class BAPI.

    The implementation is not complete:
    - only applies to SU(3) case: because of specific implementation
      of maxTr function (Cabibbo-Marinari maximization).
    - unnecessary arithmetic operations exist for the timeslices
      on which the gauge is already fixed to good precision.
    These should be improved in the version beyond test phase.
                                        [16 Feb 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.                [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                        [21 Mar 2015 Y.Namekawa]
    Staple and RandomNumbers are moved into gaugeFixing
                                        [30 Mar 2016 Y.Namekawa]
*/

class BAPI GaugeFixing_Coulomb : public GaugeFixing
{
 public:
  static const std::string class_name;

 private:
  int m_Niter;             // max iteration number
  int m_Nnaive;            // number of naive iterations
  int m_Nmeas;             // interval of measurements
  int m_Nreset;            // Number of iteration to reset the config.
  double m_Enorm;          // convergence criterion
  double m_wp;             // overrelaxation parameter

  RandomNumbers *m_rand;
  Index_eo m_index;

 public:
  GaugeFixing_Coulomb()
    : GaugeFixing(),
    m_rand(RandomNumberManager::getInstance())
  {
  }

  GaugeFixing_Coulomb(RandomNumbers *rand)
    : GaugeFixing(),
    m_rand(rand)
  {
  }

  GaugeFixing_Coulomb(unique_ptr<RandomNumbers>& rand)
    : GaugeFixing(),
    m_rand(rand.get())
  {
  }

  ~GaugeFixing_Coulomb() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Niter, const int Nnaive,
                      const int Nmeas, const int Nreset,
                      const double Enorm, const double wp);

  void fix(Field_G& Ufix, const Field_G& Uorg);

 private:
  //! one step of gauge fixing with overrelaxation parameter wp.
  void gfix_step(Field_G& Ue, Field_G& Uo, const double wp);

  void set_randomGaugeTrans(const std::valarray<double>& sg, Field_G& Geo);
  void gauge_trans_eo(Field_G& Ue, Field_G& Uo,
                      const Field_G& Geo, const int Ieo);

  void calc_SG(std::valarray<double>& sg, std::valarray<double>& Fval,
               const Field_G& Ue, const Field_G& Uo);
  void calc_DLT(Field_G& Weo,
                const Field_G& Ue, const Field_G& Uo, const int Ieo);
  void calc_W(Field_G& Weo,
              const Field_G& Ue, const Field_G& Uo, const int Ieo);

  void maxTr(Field_G&, Field_G&);
  void maxTr1(Field_G&, Field_G&);
  void maxTr2(Field_G&, Field_G&);
  void maxTr3(Field_G&, Field_G&);

  void sum_global_t(std::valarray<double>& val_global,
                    const std::valarray<double>& val_local);

#ifdef USE_FACTORY
 private:
  static GaugeFixing *create_object()
  {
    return new GaugeFixing_Coulomb();
  }

 public:
  static bool register_factory()
  {
    return GaugeFixing::Factory::Register("Coulomb", create_object);
  }
#endif
};
#endif
