/*!
        @file    fieldStrength.h

        @brief

        @author  Yusuke Namekawa  (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FIELDSTRENGTH_INCLUDED
#define FIELDSTRENGTH_INCLUDED

#include "staple_lex.h"

#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! field strength construction.

/*!
    This class BAPI constructs a field strength, Fmunu,
    defined by a clover leaf on the lattice.
                     [03 Mar 2016 Y.Namekawa]
    Add two functions construct_Fmunu_1x1_traceless() and construct_Fmunu_1x2_traceless() to give anti-Hermitian traceless field strength.
    Add a function construct_Fmunu_plaq_traceless() which defines the anti-Hermitian traceless field strength with
    an imaginary part of the plaquette.
    Add five functions contract_2F(), contract_2F_at_t(), contract_2F_at_x(),
    contract_2F_at_y(), contract_2F_at_z().
    Add four functions contract_2F_at_t(), contract_2F_at_x(),
    contract_2F_at_y(), contract_2F_at_z() to calculate Fourier transformation.
    Add four functions global_corr_x(), global_corr_y(), global_corr_z(),
    global_corr_t() which is used in contract_2F_at_?().
                     [25 May 2017 Y.Taniguchi]
*/

class BAPI FieldStrength
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  ShiftField_lex m_shift;
  Staple_lex m_staple;

 public:
  FieldStrength()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~FieldStrength() {}

 private:
  // non-copyable
  FieldStrength(const FieldStrength&);
  FieldStrength& operator=(const FieldStrength&);

 public:
  //! Constructs the anti-Hermitian field strength with four 1x1 plquette clover leaves
  void construct_Fmunu_1x1(Field_G& Fmunu,
                           const int mu, const int nu, const Field_G& U);

  //! Constructs the anti-Hermitian field strength with eight 1x2 rectangular clover leaves
  void construct_Fmunu_1x2(Field_G& Fmunu,
                           const int mu, const int nu, const Field_G& U);

  //! Constructs the anti-Hermitian traceless field strength with four 1x1 plquette clover leaves
  void construct_Fmunu_1x1_traceless(Field_G& Fmunu,
                                     const int mu, const int nu, const Field_G& U);

  //! Constructs the anti-Hermitian traceless field strength with eight 1x2 rectangular clover leaves
  void construct_Fmunu_1x2_traceless(Field_G& Fmunu,
                                     const int mu, const int nu, const Field_G& U);

  //! Constructs the anti-Hermitian traceless field strength with an imaginary part of the plaquette.
  void construct_Fmunu_plaq_traceless(Field_G& Fmunu,
                                      const int mu, const int nu, const Field_G& U);

  //! Calculate \f$\sum_x {\rm Re}{\rm tr}(F_{\mu\nu}F_{\rho\sigma})\f$ and returns its value. Intended to be used for the topological charge and energy momentum tensor.
  double contract(const Field_G& Fmunu_1, const Field_G& Fmunu_2);

  //! Calculate corr[t]=\f$\sum_{x,y,z} {\rm Re}{\rm tr}(F_{\mu\nu}F_{\rho\sigma})\f$. Intended to be used to calculate correlations function of the topological charge and energy momentum tensor.
  void contract_at_t(std::vector<double>& corr_global,
                     const Field_G& Fmunu_1, const Field_G& Fmunu_2);

  //! Calculate corr[t]=\f$\sum_{x,y,z}cos(p_xx+p_yy+p_zz){\rm Re}{\rm tr}(F_{\mu\nu}F_{\rho\sigma})\f$. Intended to be used to calculate Fourier transformed correlations function of the topological charge and energy momentum tensor. (p_x,p_y,p_z) is given by momentum_sink.
  void contract_at_t(std::vector<double>& corr_global,
                     const Field_G& Fmunu_1, const Field_G& Fmunu_2,
                     const std::vector<int>& momentum_sink,
                     const std::vector<int>& source_position);

  //! Calculate corr[x]=\f$\sum_{y,z,t}{\rm Re}{\rm tr}(F_{\mu\nu}F_{\rho\sigma})\f$. Intended to be used to calculate correlations function of the topological charge and energy momentum tensor.
  void contract_at_x(std::vector<double>& corr_global,
                     const Field_G& Fmunu_1, const Field_G& Fmunu_2);

  //! Calculate corr[x]=\f$\sum_{y,z,t}cos(p_yy+p_zz+p_tt){\rm Re}{\rm tr}(F_{\mu\nu}F_{\rho\sigma})\f$. Intended to be used to calculate correlations function of the topological charge and energy momentum tensor.
  void contract_at_x(std::vector<double>& corr_global,
                     const Field_G& Fmunu_1, const Field_G& Fmunu_2,
                     const std::vector<int>& momentum_sink,
                     const std::vector<int>& source_position);

  //! Calculate corr[y]=\f$\sum_{x,z,t} {\rm Re}{\rm tr}(F_{\mu\nu}F_{\rho\sigma})\f$. Intended to be used to calculate correlations function of the topological charge and energy momentum tensor.
  void contract_at_y(std::vector<double>& corr_global,
                     const Field_G& Fmunu_1, const Field_G& Fmunu_2);

  //! Calculate corr[y]=\f$\sum_{x,z,t}cos(p_xx+p_zz+p_tt){\rm Re}{\rm tr}(F_{\mu\nu}F_{\rho\sigma})\f$. Intended to be used to calculate correlations function of the topological charge and energy momentum tensor.
  void contract_at_y(std::vector<double>& corr_global,
                     const Field_G& Fmunu_1, const Field_G& Fmunu_2,
                     const std::vector<int>& momentum_sink,
                     const std::vector<int>& source_position);

  //! Calculate corr[z]=\f$\sum_{x,y,t} {\rm Re}{\rm tr}(F_{\mu\nu}F_{\rho\sigma})\f$. Intended to be used to calculate correlations function of the topological charge and energy momentum tensor.
  void contract_at_z(std::vector<double>& corr_global,
                     const Field_G& Fmunu_1, const Field_G& Fmunu_2);

  //! Calculate corr[z]=\f$\sum_{x,y,t}cos(p_xx+p_yy+p_tt){\rm Re}{\rm tr}(F_{\mu\nu}F_{\rho\sigma})\f$. Intended to be used to calculate correlations function of the topological charge and energy momentum tensor.
  void contract_at_z(std::vector<double>& corr_global,
                     const Field_G& Fmunu_1, const Field_G& Fmunu_2,
                     const std::vector<int>& momentum_sink,
                     const std::vector<int>& source_position);

  //! transform node-local correlator in x into global.
  void global_corr_x(std::vector<double>& corr_global,
                     const std::vector<double>& corr_local);

  //! transform node-local correlator in y into global.
  void global_corr_y(std::vector<double>& corr_global,
                     const std::vector<double>& corr_local);

  //! transform node-local correlator in z into global.
  void global_corr_z(std::vector<double>& corr_global,
                     const std::vector<double>& corr_local);

  //! transform node-local correlator in t into global.
  void global_corr_t(std::vector<double>& corr_global,
                     const std::vector<double>& corr_local);
};

#endif
