/*!
        @file    energyMomentumTensor.h

        @brief

        @author  Yusuke Taniguchi  (tanigchi)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef ENERGY_MOMENTUM_TENSOR_INCLUDED
#define ENERGY_MOMENTUM_TENSOR_INCLUDED

#include "fieldStrength.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Energy Momentum Tensor measurement.

/*!
  This class BAPI measures a gauge contribution to the energy momentum tensor
  \f$O_{1\mu\nu}(x)=F_{\mu\rho}(x)F_{\nu\rho}(x)\f$.
  (See eq.(2.17) of https://arxiv.org/abs/1404.2758.)
  The renormalized energy momentum tensor is intended to be given according to
  the Suzuki method using the gradient flow.
  (See eq.(2.24) of https://arxiv.org/abs/1404.2758.)
  The field strength \f$F_{\mu\nu}\f$ is defined by a clover leaf
  or by an imaginary part of a plaquette on the lattice.
  [24 May 2017 Y.Taniguchi]

  Change measure_EMT from void to double for Test::verify
  [9 Sep 2017 Y.Namekawa]
 */

class BAPI EnergyMomentumTensor
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  std::string m_filename_output;

  double m_c_plaq;
  double m_c_rect;
  //! maximum of momentum for Fourier transformation: p_x=[0,max_mom], p_y=[-max_mom,max_mom], p_z=[-max_mom,max_mom]
  int m_max_mom;

  FieldStrength m_field_strength;

  std::vector<Field_G> m_Fmunu_plaq;
  std::vector<Field_G> m_Fmunu_1x1;
  std::vector<Field_G> m_Fmunu_1x2;
  int m_flag_field_set;

 public:
  EnergyMomentumTensor()
    : m_vl(CommonParameters::Vlevel())
  {
    m_filename_output = "stdout";
    m_flag_field_set  = 0;
  }

  ~EnergyMomentumTensor() {}

 public:
  //! setting parameters.
  virtual void set_parameters(const Parameters& params);
  void set_parameters(const double c_plaq, const double c_rect, const int max_mom);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  /*!
    Measure averaged energy momentum tensor
    \f$O_{1\mu\nu}(t_{\rm flow})=1/V_4\sum_{x}F_{\mu\rho}(t_{\rm flow},x)F_{\nu\rho}(t_{\rm flow},x)\f$ and
    \f$O_{2}(t_{\rm flow})=1/V_4\sum_{x}F_{\mu\nu}(t_{\rm flow},x)F_{\mu\nu}(t_{\rm flow},x)\f$.
    The field strength should be constructed beforehand by calling "set_field_strength()" with the flowed link U.
    The flow time t_flow in the argument is used just for print out.
   */
  double measure_EMT(const double t_flow);

  //! Measure energy momentum tensor density \f$O_{1\mu\nu}(t_{\rm flow},t)=1/V_{xyz}\sum_{x,y,z}F_{\mu\rho}(t_{\rm flow},x,y,z,t)F_{\nu\rho}(t_{\rm flow},x,y,z,t)\f$ as a function of time and print out the result using an argument t_{\rm flow}.
  double measure_EMT_at_t(const double t_flow);

  /*!
    Measure Fourier transformation of energy momentum tensor density \f$O_{1\mu\nu}(t_{\rm flow},t)=1/V_{xyz}\sum_{x,y,z}\cos(p_xx+p_yy+p_zz)F_{\mu\rho}(t_{\rm flow},x,y,z,t)F_{\nu\rho}(t_{\rm flow},x,y,z,t)\f$ as a function of time and print out the result using an argument t_flow.
    The momentum range is set by max_momentum in yaml.
    <ul>
    <li> \f$0\le p_x\le\f$ max_momentum
    <li> -max_momentum\f$\le p_y\le\f$ max_momentum
    <li> -max_momentum\f$\le p_z\le\f$ max_momentum
    </ul>
  */
  double measure_EMT_at_t_FT(const double t_flow);

  //! Measure energy momentum tensor density \f$O_{1\mu\nu}(t_{\rm flow},x)=1/V_{yzt}\sum_{y,z,t}F_{\mu\rho}(t_{\rm flow},x,y,z,t)F_{\nu\rho}(t_{\rm flow},x,y,z,t)\f$ in x direction and print out the result using an argument t_flow.
  double measure_EMT_at_x(const double t_flow);

  /*!
    Measure Fourier transformation of energy momentum tensor density \f$O_{1\mu\nu}(t_{\rm flow},x)=1/V_{yzt}\sum_{y,z,t}\cos(p_yy+p_zz+p_tt)F_{\mu\rho}(t_{\rm flow},x,y,z,t)F_{\nu\rho}(t_{\rm flow},x,y,z,t)\f$ in x direction and print out the result using an argument t_flow.
    The momentum range is set by max_momentum in yaml.
    <ul>
    <li> \f$0\le p_x\le\f$ max_momentum
    <li> -max_momentum\f$\le p_y\le\f$ max_momentum
    <li> -max_momentum\f$\le p_z\le\f$ max_momentum
    </ul>
  */
  double measure_EMT_at_x_FT(const double t_flow);

  //! Measure energy momentum tensor density \f$O_{1\mu\nu}(t_{\rm flow},y)=1/V_{xzt}\sum_{x,z,t}F_{\mu\rho}(t_{\rm flow},x,y,z,t)F_{\nu\rho}(t_{\rm flow},x,y,z,t)\f$ in y direction and print out the result using an argument t_flow.
  double measure_EMT_at_y(const double t_flow);

  /*!
    Measure Fourier transformation of energy momentum tensor density \f$O_{1\mu\nu}(t_{\rm flow},y)=1/V_{xzt}\sum_{x,z,t}\cos(p_xx+p_zz+p_tt)F_{\mu\rho}(t_{\rm flow},x,y,z,t)F_{\nu\rho}(t_{\rm flow},x,y,z,t)\f$ in y direction and print out the result using an argument t_flow.
    The momentum range is set by max_momentum in yaml.
    <ul>
    <li> \f$0\le p_x\le\f$ max_momentum
    <li> -max_momentum\f$\le p_y\le\f$ max_momentum
    <li> -max_momentum\f$\le p_z\le\f$ max_momentum
    </ul>
  */
  double measure_EMT_at_y_FT(const double t_flow);

  //! Measure energy momentum tensor density \f$O_{1\mu\nu}(t_{\rm flow},z)=1/V_{xyt}\sum_{x,y,t}F_{\mu\rho}(t_{\rm flow},x,y,z,t)F_{\nu\rho}(t_{\rm flow},x,y,z,t)\f$ in z direction and print out the result using an argument t_flow.
  double measure_EMT_at_z(const double t_flow);

  /*!
    Measure Fourier transformation of energy momentum tensor density \f$O_{1\mu\nu}(t_{\rm flow},z)=1/V_{xyt}\sum_{x,y,t}\cos(p_xx+p_yy+p_tt)F_{\mu\rho}(t_{\rm flow},x,y,z,t)F_{\nu\rho}(t_{\rm flow},x,y,z,t)\f$ in z direction and print out the result using an argument t_flow.
    The momentum range is set by max_momentum in yaml.
    <ul>
    <li> \f$0\le p_x\le\f$ max_momentum
    <li> -max_momentum\f$\le p_y\le\f$ max_momentum
    <li> -max_momentum\f$\le p_z\le\f$ max_momentum
    </ul>
  */
  double measure_EMT_at_z_FT(const double t_flow);

  //! Construct the anti-Hermitian traceless field strength \f$F_{\mu\nu}\f$ by the flowed link U. Should be called before measuring the energy-momentum tensor.
  void set_field_strength(const Field_G& U);

 private:
  //! Returns +1 for mu<nu, -1 for mu>nu and 0 otherwise.
  int factor(const int mu, const int nu);

  //! Returns array number [0-5] of vector m_Fmunu correponding to (mu,nu).
  int index_munu2i(int mu, int nu);
};
#endif
