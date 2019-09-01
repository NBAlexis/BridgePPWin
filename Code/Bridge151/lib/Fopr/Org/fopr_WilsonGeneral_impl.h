/*!
        @file    fopr_WilsonGeneral_impl.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2019-01-25 20:04:20 #$

        @version $LastChangedRevision: 1934 $
*/

#ifndef FOPR_WILSON_GENERAL_IMPL_ORG_INCLUDED
#define FOPR_WILSON_GENERAL_IMPL_ORG_INCLUDED

#include "Fopr/fopr.h"

#include "Field/shiftField_lex.h"

#include "Tools/gammaMatrixSet.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wilson General fermion operator.

/*!
    This class BAPI implements the Wilson General fermion operator,
    including Wilson fermion on anisotropic lattice, and
    relativistic heavy quarks.
                                    [21 Mar 2015 Y.Namekawa]
 */

namespace Org {
  class BAPI Fopr_WilsonGeneral : public Fopr
  {
   public:
    static const std::string class_name;

   private:
    // lattice parameters
    int m_Nvol;
    int m_Ndim;
    int m_Nc;
    int m_Nd;

    double m_kappa_s, m_kappa_t;
    double m_nu_s, m_r_s;
    std::vector<int> m_boundary;

    std::string m_mode;
    std::string m_repr;

    void (Fopr_WilsonGeneral::*m_mult)(Field&, const Field&);
    void (Fopr_WilsonGeneral::*m_mult_dag)(Field&, const Field&);

    const Field_G *m_U;

    std::vector<GammaMatrix> m_GM;

    ShiftField_lex m_shift;
    Field_F m_trf, m_trf2;

    Bridge::VerboseLevel m_vl;


   public:
    Fopr_WilsonGeneral() { init("Dirac"); }
    Fopr_WilsonGeneral(std::string repr) { init(repr); }
    ~Fopr_WilsonGeneral() {}

    void init(const std::string repr);

    void set_parameters(const Parameters& params);
    void set_parameters(const double kappa_s, const double kappa_t,
                        const double nu_s, const double r_s);
    void set_parameters(const double kappa_s, const double kappa_t,
                        const double nu_s, const double r_s,
                        const std::vector<int> bc);

    void set_config(Field *U)
    { m_U = (Field_G *)U; }

    void set_config(unique_ptr<Field_G>& U)
    { return set_config(U.get()); }

    void set_mode(const std::string mode);

    std::string get_mode() const;

    inline void mult(Field& v, const Field& f)
    { (this->*m_mult)(v, f); }

    inline void mult_dag(Field& v, const Field& f)
    { (this->*m_mult_dag)(v, f); }

    void mult_gm5(Field& v, const Field& f);

    void proj_chiral(Field& w, const int ex1, const Field& v, const int ex2, const int ipm);

    void D(Field& v, const Field& f);

    void D_ex(Field& v, const int ex1, const Field& f, const int ex2);

    inline void Ddag(Field& w, const Field& f)
    {
      Field w2(f.nin(), f.nvol(), f.nex());

      mult_gm5(w, f);
      D(w2, w);
      mult_gm5(w, w2);
    }

    inline void DdagD(Field& w, const Field& f)
    {
      Field w2(f.nin(), f.nvol(), f.nex());

      D(w2, f);
      mult_gm5(w, w2);
      D(w2, w);
      mult_gm5(w, w2);
    }

    inline void DDdag(Field& w, const Field& f)
    {
      Field w2(f.nin(), f.nvol(), f.nex());

      mult_gm5(w2, f);
      D(w, w2);
      mult_gm5(w2, w);
      D(w, w2);
    }

    inline void H(Field& w, const Field& f)
    {
      Field w2(f.nin(), f.nvol(), f.nex());

      D(w2, f);
      mult_gm5(w, w2);
    }

    inline void mult_undef(Field&, const Field& f)
    {
      vout.crucial(m_vl, "Error at %s: mode undefined.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    const Field_F mult_gm5p(const int mu, const Field_F& w);

    void mult_gm5p(const int mu, Field_F& v, const Field_F& w);

    // void gm5p(Field_F& v, const int mu, const Field_F& w);

    //! this returns the number of floating point operations.
    double flop_count();

    int field_nvol()
    { return CommonParameters::Nvol(); }
    int field_nin()
    { return 2 * CommonParameters::Nc() * CommonParameters::Nd(); }
    int field_nex()
    { return 1; }

    void mult_up(const int mu, Field& w, const Field& f);
    void mult_dn(const int mu, Field& w, const Field& f);

   private:
    //- prohibit copy
    Fopr_WilsonGeneral(const Fopr_WilsonGeneral&) {}
    Fopr_WilsonGeneral& operator=(const Fopr_WilsonGeneral&);

    //   void mult_p (int mu, Field_F&, const Field_F&);
    //   void mult_m (int mu, Field_F&, const Field_F&);

#ifdef USE_FACTORY
   private:
    static Fopr *create_object()
    {
      return new Fopr_WilsonGeneral();
    }

    static Fopr *create_object_with_arg(const std::string& repr)
    {
      return new Fopr_WilsonGeneral(repr);
    }

   public:
    static bool register_factory()
    {
      bool init1 = Fopr::Factory_noarg::Register("WilsonGeneral/Org", create_object);
      bool init2 = Fopr::Factory_string::Register("WilsonGeneral/Org", create_object_with_arg);

      return init1 && init2;
    }
#endif
  };
}
#endif /* FOPR_WILSON_GENERAL_IMPL_ORG_INCLUDED */
