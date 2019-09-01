/*!
        @file    topologicalCharge.h

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef TOPOLOGICALCHARGE_INCLUDED
#define TOPOLOGICALCHARGE_INCLUDED

#include "fieldStrength.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Topological Charge measurement.

/*!
  This class BAPI measures a topological charge
  defined by a clover leaf on the lattice.
                     [01 Jan 2014 Y.Namekawa]
  Add three members m_Fmunu_1x1, m_Fmunu_1x2, m_Fmunu_plaq to store
  a field strength.
  Whether they have already been stored or not is specified by a member
  m_flag_field_set.
  Add a function set_field_strength() to store a field strength into the member
  m_Fmunu.
  Move a function contract_epsilon_tensor() to fieldStrength.h
  renaming contract_2F().
  Add a function measure_topological_charge.
  Add four functions measure_topological_density_t(),
  measure_topological_density_x(), measure_topological_density_y(),
  measure_topological_density_z().
                     [25 May 2017 Y.Taniguchi]
*/

class BAPI TopologicalCharge
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
  bool m_flag_field_set;

 public:
  TopologicalCharge()
    : m_vl(CommonParameters::Vlevel())
  {
    m_filename_output = "stdout";
    m_flag_field_set  = false;
  }

  virtual ~TopologicalCharge() {}

 private:
  // non-copyable
  TopologicalCharge(const TopologicalCharge&);
  TopologicalCharge& operator=(const TopologicalCharge&);

 public:
  //! setting parameters.
  virtual void set_parameters(const Parameters& params);
  void set_parameters(const double c_plaq, const double c_rect, const int max_mom);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  //! main function to measure Topological Charge. The field strength is constructed inside the function.
  double measure(const Field_G& U);

  /*!
    Measure topological charge using the stored m_Fmunu.
    The field strength should be constructed by the link U beforehand.
    This is an upper compatible of measure() including a measurement with
    Fmunu_plaq.
    The flow time tt in the argument is used just for print out.
    The function measure() is better if you are not interested in
    the topological charge density correlation functions.
  */
  void measure_topological_charge(const double tt);

  //! Measure topological charge density corr[t]=\f$1/V_{xyz}\sum_{x,y,z}q(x,y,z,t)\f$ in temporal direction using the stored m_Fmunu and print out the result using an argument tt.
  void measure_topological_density_t(const double tt);

  //! Measure Fourier transformation of topological charge density corr[t]=\f$1/V_{xyz}\sum_{x,y,z}\cos(p_xx+p_yy+p_zz)q(x,y,z,t)\f$ using the stored m_Fmunu and print out the result using an argument tt.
  void measure_topological_density_t_FT(const double tt);

  //! Measure topological charge density corr[x]=\f$1/V_{yzt}\sum_{y,z,t}q(x,y,z,t)\f$ in x direction using the stored m_Fmunu and print out the result using an argument tt.
  void measure_topological_density_x(const double tt);

  //! Measure Fourier transformation of topological charge density corr[x]=\f$1/V_{yzt}\sum_{y,z,t}\cos(p_yy+p_zz+p_tt)q(x,y,z,t)\f$ using the stored m_Fmunu and print out the result using an argument tt.
  void measure_topological_density_x_FT(const double tt);

  //! Measure topological charge density corr[y]=\f$1/V_{xzt}\sum_{x,z,t}q(x,y,z,t)\f$ in y direction using the stored m_Fmunu and print out the result using an argument tt.
  void measure_topological_density_y(const double tt);

  //! Measure Fourier transformation of topological charge density corr[y]=\f$1/V_{xzt}\sum_{x,z,t}\cos(p_xx+p_zz+p_tt)q(x,y,z,t)\f$ using the stored m_Fmunu and print out the result using an argument tt.
  void measure_topological_density_y_FT(const double tt);

  //! Measure topological charge density corr[z]=\f$1/V_{xyt}\sum_{x,y,t}q(x,y,z,t)\f$ in z direction using the stored m_Fmunu and print out the result using an argument tt.
  void measure_topological_density_z(const double tt);

  //! Measure Fourier transformation of topological charge density corr[z]=\f$1/V_{xyt}\sum_{x,y,t}\cos(p_xx+p_yy+p_tt)q(x,y,z,t)\f$ using the stored m_Fmunu and print out the result using an argument tt.
  void measure_topological_density_z_FT(const double tt);

  //  double contract_epsilon_tensor(const Field_G& Fmunu_1, const Field_G& Fmunu_2);

  //! Construct the anti-Hermitian traceless field strength \f$F_{\mu\nu}\f$ by the link U. Should be called before measuring the topological charge and density.
  void set_field_strength(const Field_G& U);
};
#endif
