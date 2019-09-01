/*!
        @file    fopr_CloverTerm_impl.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2019-01-22 15:20:26 #$

        @version $LastChangedRevision: 1930 $
*/

#ifndef FOPR_CLOVERTERM_IMPL_IMP_INCLUDED
#define FOPR_CLOVERTERM_IMPL_IMP_INCLUDED

#include "Fopr/fopr_Wilson.h"

#include "Measurements/Gauge/staple_lex.h"

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
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    Selector is implemented.        [03 Mar 2013 Y.Namekawa]
    (Selectors are replaced with factories by Aoyama-san)
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
    A mode 'F' is added.
    Note: mult with mode 'D' or 'F' multiplies

      csw kappa sigma_{mu nu} F_{mu nu}.

    (this is different from that in fopr_CloverTerm_eo,
     which multiplies 1 - csw kappa sigma_{mu nu} F_{mu nu}. )
                                    [12 Jan, 22 Jan 2019 I.Kanamori]
 */

namespace Imp {
  class BAPI Fopr_CloverTerm : public Fopr
  {
   public:
    static const std::string class_name;

   private:
    double m_kappa;
    double m_cSW;
    std::vector<int> m_boundary;
    std::string m_repr;
    std::string m_mode;
    void (Fopr_CloverTerm::*m_csw)(Field&, const Field&);
    void (Fopr_CloverTerm::*m_gm5)(Field&, const Field&);

    int m_Nc, m_Nd, m_NinF, m_Ndim;
    int m_Nvol;

    const Field_G *m_U;  //!< pointer to gauge configuration.

    ShiftField_lex m_shift;
    Staple_lex m_staple;
    Field_G m_Cup, m_Cdn, m_v1, m_v2;                  //!< for calculation of field strength.
    Field_G m_Bx, m_By, m_Bz, m_Ex, m_Ey, m_Ez;        //!< field strength.
    // Bx = -iF(1,2), By = -iF(2,0), Bz = -iF(0,1)
    // Ex = -iF(3,0), Ey = -iF(3,1), Ez = -iF(3,2)

    std::vector<GammaMatrix> m_SG;
    GammaMatrix m_GM5;

   public:
    Fopr_CloverTerm()
      : Fopr()
    {
      init("Dirac");
    }

    Fopr_CloverTerm(const std::string repr)
      : Fopr()
    {
      init(repr);
    }

    ~Fopr_CloverTerm()
    {
      tidyup();
    }

    void set_parameters(const Parameters& params);
    void set_parameters(const double kappa, const double cSW, const std::vector<int> bc);

    void set_config(Field *U);

    void set_config(unique_ptr<Field_G>& U)
    {
      set_config(U.get());
    }

    void set_mode(const std::string mode)
    {
      m_mode = mode;
    }

    std::string get_mode() const
    {
      return m_mode;
    }

    void mult(Field& v, const Field& f)
    {
      // csw kappa sigma_{mu nu} F_{mu nu}
      if ((m_mode == "D") || (m_mode == "F")) {
        mult_sigmaF(v, f);
//  } else if(m_mode=="H"){
//    H(v,f);
      } else {
        vout.crucial("Error at %s: undefined mode = %s\n", class_name.c_str(), m_mode.c_str());
        exit(EXIT_FAILURE);
      }
    }

    void mult_dag(Field& v, const Field& f)
    {
      mult(v, f);
    }

    void mult_sigmaF(Field&, const Field&);

    void mult_gm5(Field& v, const Field& w);

    void mult_isigma(Field_F&, const Field_F&,
                     const int mu, const int nu);

    int field_nvol() { return m_Nvol; }
    int field_nin() { return 2 * m_Nc * m_Nd; }
    int field_nex() { return 1; }

    //! this returns the number of floating point operations.
    double flop_count();

   private:
    void init(const std::string repr);
    void tidyup();

    void set_csw();
    void set_fieldstrength(Field_G&, const int, const int);

    void gm5_dirac(Field&, const Field&);
    void gm5_chiral(Field&, const Field&);

    // multiplies csw kappa sigma_{mu nu} F_{mu nu}
    // NOTE: this is NOT 1 - csw kappa sigma_{mu nu} F_{mu nu}
    void mult_csw(Field&, const Field&);
    void mult_csw_dirac(Field&, const Field&);
    void mult_csw_chiral(Field&, const Field&);

    void mult_csw_dirac(Field_F&, const Field_F&);
    void mult_csw_chiral(Field_F&, const Field_F&);

    int sg_index(const int mu, const int nu) { return mu * m_Ndim + nu; }
  };
}
#endif /* FOPR_CLOVERTERM_IMPL_IMP_INCLUDED */
