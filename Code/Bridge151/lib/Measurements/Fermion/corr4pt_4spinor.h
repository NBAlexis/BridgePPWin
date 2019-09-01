/*!
        @file    corr4pt_4spinor.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef CORR4PT_4SPINOR_INCLUDED
#define CORR4PT_4SPINOR_INCLUDED

#include "contract_4spinor.h"

#include "Parameters/parameters.h"
#include "Tools/gammaMatrixSet.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Four-point correlator for Wilson-element_type fermions.

/*!
    Meson four-point functions were implemented.
                                  [06 Jan 2017 Y.Namekawa]
 */

class BAPI Corr4pt_4spinor
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 private:
  std::string m_filename_output;

  GammaMatrixSet *m_gmset;
  std::vector<int> m_epsilon_index;  //!< index of totally antisymmetric tensor

 public:
  Corr4pt_4spinor(GammaMatrixSet *gmset)
    : m_vl(CommonParameters::Vlevel()), m_gmset(gmset)
  {
  }

  Corr4pt_4spinor(unique_ptr<GammaMatrixSet>& gmset)
    : m_vl(CommonParameters::Vlevel()), m_gmset(gmset.get())
  {
  }

 private:
  // non-copyable
  Corr4pt_4spinor(const Corr4pt_4spinor&);
  Corr4pt_4spinor& operator=(const Corr4pt_4spinor&);

 public:
  virtual void set_parameters(const Parameters& params);

  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  double meson_all(
    const std::vector<Field_F>& sq1,
    const std::vector<Field_F>& sq2,
    const std::vector<Field_F>& sq3,
    const std::vector<Field_F>& sq4);

  void meson_correlator(
    std::vector<dcomplex>& corr_global,
    const GammaMatrix& gm_sink_12,
    const GammaMatrix& gm_sink_34,
    const GammaMatrix& gm_src_21,
    const GammaMatrix& gm_src_43,
    const std::vector<Field_F>& sq1,
    const std::vector<Field_F>& sq2,
    const std::vector<Field_F>& sq3,
    const std::vector<Field_F>& sq4);

  double meson_momentum_all(
    const std::vector<Field_F>& sq1,
    const std::vector<Field_F>& sq2,
    const std::vector<Field_F>& sq3,
    const std::vector<Field_F>& sq4,
    const std::vector<int>& source_position);

  void meson_momentum_correlator(
    std::vector<dcomplex>& corr_global,
    const std::vector<int>& momentum_sink,
    const GammaMatrix& gm_sink_12,
    const GammaMatrix& gm_sink_34,
    const GammaMatrix& gm_src_21,
    const GammaMatrix& gm_src_43,
    const std::vector<Field_F>& sq1,
    const std::vector<Field_F>& sq2,
    const std::vector<Field_F>& sq3,
    const std::vector<Field_F>& sq4,
    const std::vector<int>& source_position);

  // double proton_test(
  //   const std::vector<Field_F>& sq_u,
  //   const std::vector<Field_F>& sq_d);

  // void proton_correlator(
  //   std::vector<dcomplex>& corr_global,
  //   const GammaMatrix& gm,
  //   const std::vector<Field_F>& sq_u,
  //   const std::vector<Field_F>& sq_d);

 private:
  // void init();

  //! totally antisymmetric tensor: index.
  // int epsilon_index(int i, int n)
  // {
  //   return m_epsilon_index[i + 3 * n];
  // }

  //! totally antisymmetric tensor: value.
  // double epsilon_value(int n)
  // {
  //   return 1.0 - 2.0 * (n / 3);
  // }

  void corr_direct(std::vector<dcomplex>& corr_direct_global,
                   const GammaMatrix& gm5_gm_sink,
                   const GammaMatrix& gm_gm5_src,
                   const std::vector<Field_F>& sq1,
                   const std::vector<Field_F>& sq2);

  void corr_cross_sink(std::vector<std::vector<dcomplex> >& corr_cross_global,
                       const GammaMatrix& gm5_gm_sink,
                       const GammaMatrix& gm_gm5_src,
                       const std::vector<Field_F>& sq1,
                       const std::vector<Field_F>& sq2);

  //! transform node-local correlator in t to global.
  void global_corr_t(std::vector<dcomplex>& corr_global,
                     const std::vector<dcomplex>& corr_local);
};
#endif
