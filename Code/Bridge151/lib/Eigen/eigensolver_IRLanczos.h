/*!
        @file    eigensolver_IRLanczos.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef EIGENSOLVER_IRLANCZOS_INCLUDED
#define EIGENSOLVER_IRLANCZOS_INCLUDED

#include "eigensolver.h"
#include "Fopr/fopr.h"
#include "Tools/sorter.h"
#include "bridge_complex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


//! Eigenvalue solver with Implicitly Restarted Lanczos algorithm.

/*!
    This class BAPI determines eigenvalues and eigenvectors for a given
    fermion operator.
    Low- or high-lying eigenmodes are determined by chaning
    SortField class BAPI object.
                                 [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.         [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                 [21 Mar 2015 Y.Namekawa]
 */

class BAPI Eigensolver_IRLanczos : public Eigensolver
{
 public:
  static const std::string class_name;

 private:

  int m_Nk;
  int m_Np;
  int m_Niter_eigen;
  double m_Enorm_eigen;
  double m_Vthreshold;

  Fopr *m_fopr;
  Sorter *m_sorter;

 public:
  Eigensolver_IRLanczos(Fopr *fopr)
    : Eigensolver(), m_fopr(fopr), m_sorter(0) {}

  Eigensolver_IRLanczos(unique_ptr<Fopr>& fopr)
    : Eigensolver(), m_fopr(fopr.get()), m_sorter(0) {}

  ~Eigensolver_IRLanczos();

  void set_parameters(const Parameters& params);
  void set_parameters(const int Nk, const int Np, const int Niter_eigen, const double Enorm_eigen,
                      const double Vthreshold);
  void set_parameters(const std::string& sort_type,
                      const int Nk, const int Np, const int Niter_eigen, const double Enorm_eigen,
                      const double Vthreshold);

  void solve(std::vector<double>& TDa, std::vector<Field>& vk,
             int& Nsbt, int& Nconv, const Field& b);

 private:
  void step(const int Nm, const int k,
            std::vector<double>& TDa,
            std::vector<double>& TDb,
            std::vector<Field>& vk, Field& f);

  void qrtrf(std::vector<double>& TDa, std::vector<double>& TDb,
             const int Nk, const int Nm, std::vector<double>& Qt,
             const double Dsh, const int kmin, const int kmax);

  void tqri(std::vector<double>& TDa, std::vector<double>& TDb,
            const int Nk, const int Nm, std::vector<double>& Qt);

  void setUnit_Qt(const int Nm, std::vector<double>& Qt);

  void schmidt_orthogonalization(Field& w, const std::vector<Field>& vk, const int k);
};
#endif
