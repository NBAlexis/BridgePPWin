/*!
        @file    $Id:: fopr_Wilson_impl.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef FOPR_WILSON_IMPL_ORG_INCLUDED
#define FOPR_WILSON_IMPL_ORG_INCLUDED

#include "Fopr/fopr.h"

#include "Tools/gammaMatrixSet.h"
#include "Field/shiftField_lex.h"

#include "Tools/mat_SU_N.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wilson fermion operator.

/*!
    This fermion operator defines the standard Wilson fermion.
    The gamma matrix representation is given as control string
    "Dirac"(default) or "Chiral" at the construction, which is
    used to construct the Fopr_Wilson instance.
    The `mode', which of D, Ddag, H, DdagD are multiplied, is
    controlled by setting the pointers to these functions,
    m_mult and m_mult_dag.
    At the beginning, they are set to point mult_undef() which
    just represent the mode has not been set.
    set_mode(string) must be called before mult() is called.
                                    [24 Dec 2011 H,Matsufuru]
 */

namespace Org {
  class Fopr_Wilson : public Fopr
  {
   public:
    static const std::string class_name;

   public:

    Fopr_Wilson() { init("Dirac"); }
    Fopr_Wilson(std::string repr) { init(repr); }
    ~Fopr_Wilson() {}

    void init(std::string repr);

    void set_parameters(const Parameters& params);
    void set_parameters(const double kappa);
    void set_parameters(const double kappa, const std::vector<int> bc);

    void set_config(Field *U)
    { m_U = (Field_G *)U; }

    void set_config(unique_ptr<Field_G>& U)
    { return set_config(U.get()); }

    void set_mode(std::string mode);

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

    const Field_F mult_gm5p(int mu, const Field_F& w);

    void mult_gm5p(int mu, Field_F& v, const Field_F& w);

//  void gm5p(Field_F& v, const int mu, const Field_F& w);

    //! this returns the number of floating point operations.
    double flop_count();

    int field_nvol()
    { return CommonParameters::Nvol(); }
    int field_nin()
    { return 2 * CommonParameters::Nc() * CommonParameters::Nd(); }
    int field_nex()
    { return 1; }


    void mult_up(int mu, Field& w, const Field& f);
    void mult_dn(int mu, Field& w, const Field& f);

    inline void fprop_normalize(Field& v)
    { scal(v, (2.0 * m_kappa)); }

    inline void fopr_normalize(Field& v)
    { scal(v, 1.0 / (2.0 * m_kappa)); }


   private:
    //- prohibit copy
    Fopr_Wilson(const Fopr_Wilson&) {}
    Fopr_Wilson& operator=(const Fopr_Wilson&);


    Bridge::VerboseLevel m_vl;

    // lattice parameters
    int m_Nvol;
    int m_Ndim;
    int m_Nc;
    int m_Nd;

    double           m_kappa;
    std::vector<int> m_boundary;

    std::string m_mode;
    std::string m_repr;

    void (Fopr_Wilson::*m_mult)(Field&, const Field&);
    void (Fopr_Wilson::*m_mult_dag)(Field&, const Field&);

    const Field_G *m_U;

    std::vector<GammaMatrix> m_GM;

    ShiftField_lex shift;
    Field_F        trf, trf2;

//   void mult_p (int mu, Field_F&, const Field_F&);
//   void mult_m (int mu, Field_F&, const Field_F&);
  };
}
#endif
