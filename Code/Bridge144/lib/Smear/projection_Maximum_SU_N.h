/*!
        @file    $Id:: projection_Maximum_SU_N.h #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/


#ifndef PROJECTION_MAXIMUM_SU_N_INCLUDED
#define PROJECTION_MAXIMUM_SU_N_INCLUDED

#include "projection.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Maximum projection to SU(N) gauge group.

/*!
    Maximum projection for SU(N) matrix by Cabibbo-Marinari
    method, SU(2) subgroup transformation.
    The code was originally written by Takashi Umeda (1997) in
    Fortran by explicitely assuming SU(3) group.
    Genralization to SU((N) was done by H.M.
                                    [09 Aug 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
 */


class Projection_Maximum_SU_N : public Projection
{
 public:
  static const std::string class_name;

 private:
  int    m_Niter;  //!< maximum iteration of maximization steps
  double m_Enorm;  //!< convergence criterion of maximization

 public:
  Projection_Maximum_SU_N()
  {
    //- defaults
    m_Niter = 100;
    m_Enorm = 1.0e-12;
  }

  ~Projection_Maximum_SU_N() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int Niter, const double Enorm);

  //! projection U = P[alpha, C, Uorg]
  void project(Field_G& U,
               double alpha,
               const Field_G& C, const Field_G& Uorg);

  //! force calculation: invalid in this class.
  void force_recursive(Field_G& Xi, Field_G& iTheta,
                       double alpha, const Field_G& Sigmap,
                       const Field_G& C, const Field_G& U);

 private:
  //- maximization of ReTr[U^\dag V].
  void maxTr(Field_G& U, Field_G& V);

  //- maximization by SU(2) subgroup.
  void maxTr_SU2(int, int, Field_G&, Field_G&, Field_G&);

  //- matrix index for convenience.
  int mindex(int i, int j, int Nc)
  {
    return i + j * Nc;
  }
};
#endif
