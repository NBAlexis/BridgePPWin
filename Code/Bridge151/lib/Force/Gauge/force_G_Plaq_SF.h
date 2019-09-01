/*!
        @file    force_G_Plaq_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef FORCE_G_PLAQ_SF_INCLUDED
#define FORCE_G_PLAQ_SF_INCLUDED

#include "force_G.h"

#include "Measurements/Gauge/staple_SF.h"

#include "IO/bridgeIO.h"
using Bridge::vout;


/*!
  \brief HMC force class BAPI for plaquette gauge force with SF BC.

  SF boundary condition is intrduced accrdong to the following policy.
  <ul>
    <li>The total size of temporal links is Nt with the Dirichlet BC.
    <li>We need t=0, ..., Nt sites for Nt links.
    <li>The boundary spatial link at t=0 is dummy.
    <ul>
      <li>Calculation of the Hamiltonian and force is overridden to given the correct result without using the spatial link at t=0 stored in the Field_G object.
      <li>The spatial link at t=0 is updated by a randomly given conjugate momentum but is not used for the Hamiltoniand nor the force.
      <li>The corresponding force for the spatial conjugate momentum at t=0 is set to zero. This conjugate momentum is not updated at the boundary.
    </ul>
    <li>The boundary spatial sites at t=Nt does not exist in the Field_G object.
    <ul>
      <li>Calculation of the Hamiltonian and force is overridden to given the correct result using an approriate matrix for the spatial link at t=Nt.
    </ul>
  </ul>
  The plaquette gauge force is given by
\f[
S_g[U]=\frac{\beta}{N_c}\sum_{x}\sum_{\mu<\nu}\omega_{\mu\nu}(x)
{\rm Re}{\rm Tr}\left(-U_{\mu\nu}(x)\right)
\f]
\f$\omega_{0i}(t=0)=\omega_{0i}(t=T)=c_t\f$,
\f$\omega_{\mu\nu}(n)=1\f$: otherwise.
  <ul>
  <li>A major difference from Force_G_Plaq class BAPI is a use of Staple_SF instead of Staple.
  <ul>
  <li>For correct Hamltonian with SF BC.
  <li>For correct force with SF BC.
  </ul>
  <li>Boundary condition can be accessed with m_phi and m_phipr.
  <li>Boundary improvement factor is stored in m_ct.
  <li> [25 Jan. 2012 Y.Taniguchi]
  </ul>
    (Coding history will be recovered from trac.)
    YAML is implemented.     [14 Nov 2012 Y.Namekawa]
*/

class BAPI Force_G_Plaq_SF : public Force_G
{
 public:
  static const std::string class_name;

 private:
  //- NB. m_U has been defined in force_G.h
  // Field_G *m_U;

  double m_beta;

  //! SF boundary condition at t=0
  double *m_phi;
  //! SF boundary condition at t=Nt
  double *m_phipr;
  //! SF boundary improvement coefficient for the plaquatte force
  double m_ct;

  Staple_SF m_staple;

 public:
  Force_G_Plaq_SF()
    : Force_G() {}

  ~Force_G_Plaq_SF() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const double beta,
                      double *phi, double *phipr, const double ct);

  //- NB. set_config has been defined in force_G.h
  // void set_config(Field *U)

  void force_core(Field&);

  /*!
    \brief Print out the gauge force for the boundary spatial links (for debugging).
  */
  void print_force(const Field_G *);

#ifdef USE_FACTORY
 private:
  static Force_G *create_object()
  {
    return new Force_G_Plaq_SF();
  }

 public:
  static bool register_factory()
  {
    return Force_G::Factory::Register("Force_G_Plaq_SF", create_object);
  }
#endif
};
#endif
