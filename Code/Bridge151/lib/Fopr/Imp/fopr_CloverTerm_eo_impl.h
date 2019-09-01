/*!
        @file    fopr_CloverTerm_eo_impl.h

        @brief

        @author  UEDA, Satoru (sueda)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2019-01-22 15:20:26 #$

        @version $LastChangedRevision: 1930 $
*/

#ifndef FOPR_CLOVERTERM_EO_IMPL_IMP_INCLUDED
#define FOPR_CLOVERTERM_EO_IMPL_IMP_INCLUDED

#include "Fopr/fopr.h"

#include "Field/shiftField_eo.h"

#include "Tools/gammaMatrixSet.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Clover term operator.

/*!
    This class BAPI implements the clover term for the clover (improved
    Wilson) fermion operator.
    This part was separated from the Fopr_Clover class BAPI.
    The field strength is calculate when the function
    set_config() is called.
    The `mode' for setting fermion operator mode is now only
    defined to the case 'D'.
                [30 Sep 2012 H.Matsufuru,
                 original clover operator: 24 Dec 2011 H.M.]
    (Coding history will be recovered from trac.)
    Modify this code to work.      [03 Mar 2013 Y.Namekawa]
    Multi-threading was applied to D() and mult_csw_inv().
    Previous version explicitly implements the Dirac representaion
    of gamma-matrices and thus failed for other representation.
    In the cpp files in Imp/ and other improved performance versions,
    now chiral representation is available in addition to the Dirac.
    For these implementation, performance tuning was also applied.
                                   [31 Jul 2014 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                   [21 Mar 2015 Y.Namekawa]

    Note: mult with mode 'even' or 'odd' multiplies

      1 - csw kappa sigma_{mu nu} F_{mu nu}.

    (this is different from that of fopr_CloverTerm with mode 'D',
     which multiplies csw kappa sigma_{mu nu} F_{mu nu}. )
                                        [22 Jan 2019 I.Kanamori]
 */

namespace Imp {
  class BAPI Fopr_CloverTerm_eo : public Fopr
  {
   public:
    static const std::string class_name;

// This class BAPI returns D_ee = 1-f_ee or D_oo = 1-f_oo
   private:
    double m_kappa;
    double m_cSW;
    std::vector<int> m_boundary;
    std::string m_repr;
    std::string m_mode;

    int m_Nvol, m_Nvol2;
    int m_Ndim;
    int m_Nc, m_Nd, m_Ndm2;
    int m_NinF;

    //! Gamma Matrix and Sigma_{mu,nu} = -i [Gamma_mu, Gamma_nu] /2
    std::vector<GammaMatrix> m_GM, m_SG;

    void (Fopr_CloverTerm_eo::*m_mult)(Field&, const Field&);

    const Field_G *m_Ueo;

    Index_eo m_idx;
    ShiftField_eo m_shift_eo;

    Field_F *m_fee_inv;
    Field_F *m_foo_inv;
    Field_F m_vf, m_ff;
    Vec_SU_N v1, v2;

    //! m_T = 1 - kappa c_SW sigma F / 2
    Field_G m_T;

    //! m_T2 is used in Org-version.
    std::vector<Field_G> m_T2;

   public:
    Fopr_CloverTerm_eo(std::string repr)
      : m_Nvol(CommonParameters::Nvol()),
      m_Nvol2(m_Nvol / 2),
      m_Ndim(CommonParameters::Ndim()),
      m_Nc(CommonParameters::Nc()),
      m_Nd(CommonParameters::Nd()),
      m_Ndm2(m_Nd * m_Nd / 2),
      m_T(m_Nvol, m_Ndm2)
    {
      init(repr);
    }

    ~Fopr_CloverTerm_eo()
    {
      delete m_fee_inv;
      delete m_foo_inv;
    }

    void set_parameters(const Parameters& params);
    void set_parameters(const double kappa, const double cSW,
                        const std::vector<int> bc);

    void set_config(Field *Ueo);

    void set_config(unique_ptr<Field_G>& Ueo)
    {
      set_config(Ueo.get());
    }

    void set_mode(const std::string mode)
    {
      m_mode = mode;
    }

    std::string get_mode() const
    {
      return m_mode;
    }

    //! return D = D^dag = 1-f_ee or 1-f_oo
    const Field mult(const Field& f)
    {
      Field v(f.nin(), f.nvol(), f.nex());

      mult(v, f);
      return v;
    }

    const Field mult_dag(const Field& f)
    {
      Field v(f.nin(), f.nvol(), f.nex());

      mult_dag(v, f);
      return v;
    }

    //! return D = D^dag = 1-f_ee or 1-f_oo
    void mult(Field& v, const Field& f)
    {
      // multiplies 1-csw kappa sigma_{mu nu} F_{mu nu}
      if (m_mode == "even") {
        D(v, f, 0);
      } else if (m_mode == "odd") {
        D(v, f, 1);
      } else {
        vout.crucial("Error at %s: undefined mode = %s\n", class_name.c_str(), m_mode.c_str());
        exit(EXIT_FAILURE);
      }
    }

    void mult_dag(Field& v, const Field& f)
    {
      // multiplies 1-csw kappa sigma_{mu nu} F_{mu nu}
      mult(v, f);
    }

    void mult_isigma(Field_F&, const Field_F&,
                     const int mu, const int nu);

    // multiplies 1-csw kappa sigma_{mu nu} F_{mu nu}
    void D(Field& v, const Field& f, const int ieo);

    //! explicit implementation for Dirac representation (for Imp-version).
    void D_dirac(Field& v, const Field& f, const int ieo);

    //! explicit implementation for Chiral representation (for Imp-version).
    void D_chiral(Field& v, const Field& f, const int ieo);

    //const Field_F mult_csw_inv(const Field_F&, const int ieo);
    //const Field_G trSigmaInv(const int mu, const int nu);
    void trSigmaInv(Field_G&, const int mu, const int nu);

    // multiplies [ 1-csw kappa sigma_{mu nu} F_{mu nu} ]^{-1}
    void mult_csw_inv(Field&, const Field&, const int ieo);

    void mult_csw_inv_dirac(Field&, const Field&, const int ieo);

    void mult_csw_inv_chiral(Field&, const Field&, const int ieo);

    std::vector<double> csmatrix(const int&);

    int field_nvol() { return m_Nvol2; }
    int field_nin() { return 2 * m_Nc * m_Nd; }
    int field_nex() { return 1; }

    //! returns number of floating point operations.
    double flop_count();

   private:
    void init(const std::string repr);
    void tidyup();

    void solve_csw_inv();

    void set_csw();

    //! explicit implementation for Dirac representation (for Imp-version).
    void set_csw_dirac();

    //! explicit implementation for Chiral representation (for Imp-version).
    void set_csw_chiral();

    void mult_csw(Field_F&, const Field_F&, const int ieo);
    void set_fieldstrength(Field_G&, const int, const int);

    int sg_index(const int mu, const int nu) { return mu * m_Ndim + nu; }
  };
}
#endif /* FOPR_CLOVERTERM_EO_IMPL_IMP_INCLUDED */
