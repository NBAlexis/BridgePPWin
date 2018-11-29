/*!
@file    $Id:: corr2pt_4spinor.h #$

@brief

@author  Hideo Matsufuru (matsufuru)
$LastChangedBy: namekawa $

@date    $LastChangedDate:: 2017-03-02 16:46:34 #$

@version $LastChangedRevision: 1580 $
*/

#ifndef CORR2PT_4SPINOR_INCLUDED
#define CORR2PT_4SPINOR_INCLUDED

#include "contract_4spinor.h"

#include "Parameters/parameters.h"
#include "Tools/gammaMatrixSet.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Two-point correlator for Wilson-type fermions.

/*!
Meson correlators were implemented.
[4 Feb 2012 H.Matsufuru]
Baryon (proton) correlator was implemented by K.Nemuta.
This implementation assumes Nc=3, and some of parameters
are replaced by explicit numbers.
Better performance version:   [28 Jul 2012 H.Matsufuru].
unique_ptr is introduced to avoid memory leaks.
Add momentum of sink.         [21 Mar 2015 Y.Namekawa]
Add parameters for output.    [27 Jun 2016 Y.Namekawa]
*/

class BAPI Corr2pt_4spinor
{
public:
    static const std::string class_name;

protected:
    Bridge::VerboseLevel m_vl;

private:
    std::string m_filename_output;

    GammaMatrixSet   *m_gmset;
    std::vector<int> m_epsilon_index;  //!< index of totally antisymmetric tensor

public:
    Corr2pt_4spinor(GammaMatrixSet *gmset)
        : m_vl(CommonParameters::Vlevel()), m_gmset(gmset)
    {
        init();
    }

    Corr2pt_4spinor(unique_ptr<GammaMatrixSet>& gmset)
        : m_vl(CommonParameters::Vlevel()), m_gmset(gmset.get())
    {
        init();
    }

private:
    // non-copyable
    Corr2pt_4spinor(const Corr2pt_4spinor&);
    Corr2pt_4spinor& operator=(const Corr2pt_4spinor&);

public:
    virtual void set_parameters(const Parameters& params);

    void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

    double meson_all(
        const std::vector<Field_F>& sq1,
        const std::vector<Field_F>& sq2);

    void meson_correlator(
        std::vector<dcomplex>& corr_global,
        const GammaMatrix& gm_sink,
        const GammaMatrix& gm_src,
        const std::vector<Field_F>& sq1,
        const std::vector<Field_F>& sq2);

    double meson_momentum_all(
        const std::vector<Field_F>& sq1,
        const std::vector<Field_F>& sq2,
        const std::vector<int>& source_position);

    void meson_momentum_correlator(
        std::vector<dcomplex>& corr_global,
        const std::vector<int>& momentum_sink,
        const GammaMatrix& gm_sink,
        const GammaMatrix& gm_src,
        const std::vector<Field_F>& sq1,
        const std::vector<Field_F>& sq2,
        const std::vector<int>& source_position);

    double proton_test(
        const std::vector<Field_F>& sq_u,
        const std::vector<Field_F>& sq_d);

    void proton_correlator(
        std::vector<dcomplex>& corr_global,
        const GammaMatrix& gm,
        const std::vector<Field_F>& sq_u,
        const std::vector<Field_F>& sq_d);

private:
    void init();

    //! totally antisymmetric tensor: index.
    int epsilon_index(int i, int n)
    {
        return m_epsilon_index[i + 3 * n];
    }

    //! totally antisymmetric tensor: value.
    double epsilon_value(int n)
    {
        return 1.0 - 2.0 * (n / 3);
    }

    //! transform node-local correlator in t to global.
    void global_corr_t(std::vector<dcomplex>& corr_global,
        std::vector<dcomplex>& corr_local);
};
#endif
