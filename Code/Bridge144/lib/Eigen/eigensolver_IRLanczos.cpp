#include "BridgeLib_Private.h"

/*!
        @file    $Id:: eigensolver_IRLanczos.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: namekawa $

        @date    $LastChangedDate:: 2017-03-03 15:12:38 #$

        @version $LastChangedRevision: 1582 $
*/

#include "eigensolver_IRLanczos.h"

using std::string;



const std::string Eigensolver_IRLanczos::class_name = "Eigensolver_IRLanczos";

//====================================================================
Eigensolver_IRLanczos::~Eigensolver_IRLanczos()
{
  if (m_sorter) delete m_sorter;
}


//====================================================================
void Eigensolver_IRLanczos::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  string str_sortfield_type;
  int    Nk, Np;
  int    Niter_eigen;
  double Enorm_eigen, Vthreshold;

  int err = 0;
  err += params.fetch_string("eigensolver_mode", str_sortfield_type);
  err += params.fetch_int("number_of_wanted_eigenvectors", Nk);
  err += params.fetch_int("number_of_working_eigenvectors", Np);
  err += params.fetch_int("maximum_number_of_iteration", Niter_eigen);
  err += params.fetch_double("convergence_criterion_squared", Enorm_eigen);
  err += params.fetch_double("threshold_value", Vthreshold);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(str_sortfield_type, Nk, Np, Niter_eigen, Enorm_eigen, Vthreshold);
}


//====================================================================
void Eigensolver_IRLanczos::set_parameters(const std::string& sort_type,
                                           int Nk, int Np,
                                           int Niter_eigen, double Enorm_eigen,
                                           double Vthreshold)
{
  if (m_sorter) delete m_sorter;

  m_sorter = new Sorter(sort_type);

  set_parameters(Nk, Np, Niter_eigen, Enorm_eigen, Vthreshold);
}


//====================================================================
void Eigensolver_IRLanczos::set_parameters(int Nk, int Np,
                                           int Niter_eigen, double Enorm_eigen,
                                           double Vthreshold)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Nk          = %d\n", Nk);
  vout.general(m_vl, "  Np          = %d\n", Np);
  vout.general(m_vl, "  Niter_eigen = %d\n", Niter_eigen);
  vout.general(m_vl, "  Enorm_eigen = %16.8e\n", Enorm_eigen);
  vout.general(m_vl, "  Vthreshold  = %16.8e\n", Vthreshold);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Nk);
  err += ParameterCheck::non_negative(Np);
  err += ParameterCheck::non_negative(Niter_eigen);
  err += ParameterCheck::square_non_zero(Enorm_eigen);
  // NB. Vthreshold == 0 is allowed.

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Nk          = Nk;
  m_Np          = Np;
  m_Niter_eigen = Niter_eigen;
  m_Enorm_eigen = Enorm_eigen;
  m_Vthreshold  = Vthreshold;
}


//====================================================================
void Eigensolver_IRLanczos::solve(std::vector<double>& TDa, std::vector<Field>& vk,
                                  int& Nsbt, int& Nconv, const Field& )
{
  int    Nk          = m_Nk;
  int    Np          = m_Np;
  int    Niter_eigen = m_Niter_eigen;
  double Enorm_eigen = m_Enorm_eigen;
  double Vthreshold  = m_Vthreshold;


  Nconv = -1;
  Nsbt  = 0;

  int Nm = Nk + Np;

  if (Nk + Np > TDa.size()) {
    vout.crucial(m_vl, "Error at %s: std::vector TDa is too small.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  } else if (Nk + Np > vk.size()) {
    vout.crucial(m_vl, "Error at %s: std::vector vk is too small.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  std::vector<double> TDb(Nm);
  std::vector<double> TDa2(Nm);
  std::vector<double> TDb2(Nm);
  std::vector<double> Qt(Nm * Nm);
  std::vector<int>    Iconv(Nm);

  int                Nin  = vk[0].nin();
  int                Nvol = vk[0].nvol();
  int                Nex  = vk[0].nex();
  std::vector<Field> B(Nm);
  for (int k = 0; k < Nm; ++k) {
    B[k].reset(Nin, Nvol, Nex);
  }

  Field f(vk[0]);
  Field v(vk[0]);

  vout.general(m_vl, "  Nk = %d  Np = %d\n", Nk, Np);
  vout.general(m_vl, "  Nm = %d\n", Nm);
  vout.general(m_vl, "  size of TDa = %d\n", TDa.size());
  vout.general(m_vl, "  size of vk  = %d\n", vk.size());

  int k1    = 1;
  int k2    = Nk;
  int kconv = 0;

  int    Kdis       = 0;
  int    Kthreshold = 0;
  double beta_k;

  //- Set initial vector
  vk[0].set(1.0);
  double vnorm = dot(vk[0], vk[0]);
  vk[0].set(1.0 / sqrt(vnorm));
  // (uniform vector)

  //- Initial Nk steps
  for (int k = 0; k < k2; ++k) {
    step(Nm, k, TDa, TDb, vk, f);
  }

  //- Restarting loop begins
  for (int iter = 0; iter < Niter_eigen; ++iter) {
    vout.detailed(m_vl, "\n iter=%d\n", iter);

    int Nm2 = Nm - kconv;

    for (int k = k2; k < Nm; ++k) {
      step(Nm, k, TDa, TDb, vk, f);
    }

    //f *= TDb[Nm - 1];
    scal(f, TDb[Nm - 1]);

    //- getting eigenvalues
    for (int k = 0; k < Nm2; ++k) {
      TDa2[k] = TDa[k + k1 - 1];
      TDb2[k] = TDb[k + k1 - 1];
    }
    setUnit_Qt(Nm, Qt);
    tqri(TDa2, TDb2, Nm2, Nm, Qt);

    //- sorting
    m_sorter->sort(TDa2, Nm);

    //- Implicitly shifted QR transformations
    setUnit_Qt(Nm, Qt);
    for (int ip = k2; ip < Nm; ++ip) {
      double Dsh  = TDa2[ip - kconv];
      int    kmin = k1;
      int    kmax = Nm;
      qrtrf(TDa, TDb, Nm, Nm, Qt, Dsh, kmin, kmax);
    }

    for (int i = 0; i < (Nk + 1); ++i) {
      B[i].set(0.0);
    }

    for (int j = k1 - 1; j < k2 + 1; ++j) {
      for (int k = 0; k < Nm; ++k) {
        axpy(B[j], Qt[k + Nm * j], vk[k]);
      }
    }

    for (int j = k1 - 1; j < k2 + 1; ++j) {
      vk[j] = B[j];
    }

    //- Compressed vector f and beta(k2)
    scal(f, Qt[Nm - 1 + Nm * (k2 - 1)]);
    axpy(f, TDb[k2 - 1], vk[k2]);

    beta_k = dot(f, f);
    beta_k = sqrt(beta_k);

    vout.detailed(m_vl, " beta(k) = %20.14f\n", beta_k);


    double beta_r = 1.0 / beta_k;
    vk[k2] = f;
    scal(vk[k2], beta_r);
    TDb[k2 - 1] = beta_k;

    //- Convergence test
    TDa2 = TDa;
    TDb2 = TDb;
    setUnit_Qt(Nm, Qt);

    tqri(TDa2, TDb2, Nk, Nm, Qt);
    for (int k = 0; k < Nk; ++k) {
      B[k].set(0.0);
    }

    for (int j = 0; j < Nk; ++j) {
      for (int k = 0; k < Nk; ++k) {
        axpy(B[j], Qt[k + j * Nm], vk[k]);
      }
    }

    Kdis       = 0;
    Kthreshold = 0;

    for (int i = 0; i < Nk; ++i) {
      m_fopr->mult(v, B[i]);
      double vnum = dot(B[i], v);
      double vden = dot(B[i], B[i]);

      //      vout.paranoiac(m_vl, " vden = %20.14e\n",vden);

      TDa2[i] = vnum / vden;
      axpy(v, -TDa2[i], B[i]);

      double vv = dot(v, v);

      vout.detailed(m_vl, "  %4d  %18.14f  %18.14e\n", i, TDa2[i], vv);

      if (vv < Enorm_eigen) {
        Iconv[Kdis] = i;
        ++Kdis;

        if (!m_sorter->comp(TDa2[i], Vthreshold)) {
          ++Kthreshold;
        }
      }
    }  // i-loop end


    vout.detailed(m_vl, " #modes converged: %d\n", Kdis);


    if (Kthreshold > 0) {
      // (there is a converged eigenvalue larger than Vthreshold.)
      Nconv = iter;

      //- Sorting
      for (int i = 0; i < Kdis; ++i) {
        TDa[i] = TDa2[Iconv[i]];
      }

      std::vector<int> idx = m_sorter->sort_index(TDa, Kdis);

      for (int i = 0; i < Kdis; ++i) {
        vk[i] = B[Iconv[idx[i]]];
      }

      Nsbt = Kdis - Kthreshold;

      vout.general(m_vl, "\n Converged:\n");
      vout.general(m_vl, "  Nconv   = %d\n", Nconv);
      vout.general(m_vl, "  beta(k) = %20.14e\n", beta_k);
      vout.general(m_vl, "  Kdis    = %d\n", Kdis);
      vout.general(m_vl, "  Nsbt    = %d\n", Nsbt);

      return;
    }
  } // end of iter loop


  if (Nconv == -1) {
    vout.crucial(m_vl, "Error at %s: NOT converged.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Eigensolver_IRLanczos::step(int Nm, int k, std::vector<double>& TDa,
                                 std::vector<double>& TDb, std::vector<Field>& vk,
                                 Field& w)
{
  if (k >= Nm) {
    vout.crucial(m_vl, "Error at %s: k is larger than Nm.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  } else if (k == 0) {  // Initial step
    m_fopr->mult(w, vk[k]);
    double alph = dot(vk[k], w);

    axpy(w, -alph, vk[k]);

    double beta = dot(w, w);
    beta = sqrt(beta);
    double beta_r = 1.0 / beta;
    vk[k + 1] = w;
    scal(vk[k + 1], beta_r);

    TDa[k] = alph;
    TDb[k] = beta;
  } else {   // Iteration step
    m_fopr->mult(w, vk[k]);
    axpy(w, -TDb[k - 1], vk[k - 1]);

    double alph = dot(vk[k], w);

    axpy(w, -alph, vk[k]);

    double beta = dot(w, w);
    beta = sqrt(beta);
    double beta_r = 1.0 / beta;
    scal(w, beta_r);

    TDa[k] = alph;
    TDb[k] = beta;

    schmidt_orthogonalization(w, vk, k);

    if (k < Nm - 1) vk[k + 1] = w;
  }
}


//====================================================================
void Eigensolver_IRLanczos::schmidt_orthogonalization(Field& w,
                                                      std::vector<Field>& vk, int k)
{
  for (int j = 0; j < k; ++j) {
    dcomplex prod = dotc(vk[j], w);
    prod *= cmplx(-1.0, 0.0);
    axpy(w, prod, vk[j]);
  }
}


//====================================================================
void Eigensolver_IRLanczos::setUnit_Qt(int Nm, std::vector<double>& Qt)
{
  for (int i = 0; i < Qt.size(); ++i) {
    Qt[i] = 0.0;
  }

  for (int k = 0; k < Nm; ++k) {
    Qt[k + k * Nm] = 1.0;
  }
}


//====================================================================
void Eigensolver_IRLanczos::tqri(std::vector<double>& TDa,
                                 std::vector<double>& TDb,
                                 int Nk, int Nm, std::vector<double>& Qt)
{
  int Niter = 100 * Nm;
  int kmin  = 1;
  int kmax  = Nk;
  // (these parameters should be tuned)


  int Nconv = -1;

  for (int iter = 0; iter < Niter; ++iter) {
    //- determination of 2x2 leading submatrix
    double dsub = TDa[kmax - 1] - TDa[kmax - 2];
    double dd   = sqrt(dsub * dsub + 4.0 * TDb[kmax - 2] * TDb[kmax - 2]);
    double Dsh  = 0.5 * (TDa[kmax - 2] + TDa[kmax - 1]
                         + fabs(dd) * (dsub / fabs(dsub)));
    // (Dsh: shift)

    //- transformation
    qrtrf(TDa, TDb, Nk, Nm, Qt, Dsh, kmin, kmax);

    //- Convergence criterion (redef of kmin and kmax)
    for (int j = kmax - 1; j >= kmin; --j) {
      double dds = fabs(TDa[j - 1]) + fabs(TDa[j]);
      if (fabs(TDb[j - 1]) + dds > dds) {
        kmax = j + 1;

        for (int jj = 0; jj < kmax - 1; ++jj) {
          double ddsjj = fabs(TDa[jj]) + fabs(TDa[jj + 1]);

          if (fabs(TDb[jj]) + ddsjj > ddsjj) {
            kmin = jj + 1;

            break;
          }
        }

        break;
      }

      if (j == kmin) {
        Nconv = iter;
        vout.paranoiac(m_vl, "  converged at iter = %d\n", Nconv);

        return;
      }
    } // end of j loop
  }   // end of iter loop

  if (Nconv == -1) {
    vout.crucial(m_vl, "Error at %s: QL method NOT converged, Niter = %d.\n", class_name.c_str(), Niter);
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Eigensolver_IRLanczos::qrtrf(std::vector<double>& TDa,
                                  std::vector<double>& TDb,
                                  int Nk, int Nm, std::vector<double>& Qt,
                                  double Dsh, int kmin, int kmax)
{
  int    k = kmin - 1;
  double x;

  double Fden = 1.0 / sqrt((TDa[k] - Dsh) * (TDa[k] - Dsh)
                           + TDb[k] * TDb[k]);
  double c = (TDa[k] - Dsh) * Fden;
  double s = -TDb[k] * Fden;

  double tmpa1 = TDa[k];
  double tmpa2 = TDa[k + 1];
  double tmpb  = TDb[k];

  TDa[k]     = c * c * tmpa1 + s * s * tmpa2 - 2.0 * c * s * tmpb;
  TDa[k + 1] = s * s * tmpa1 + c * c * tmpa2 + 2.0 * c * s * tmpb;
  TDb[k]     = c * s * (tmpa1 - tmpa2) + (c * c - s * s) * tmpb;
  x          = -s * TDb[k + 1];
  TDb[k + 1] = c * TDb[k + 1];

  for (int i = 0; i < Nk; ++i) {
    double Qtmp1 = Qt[i + Nm * k];
    double Qtmp2 = Qt[i + Nm * (k + 1)];
    Qt[i + Nm * k]       = c * Qtmp1 - s * Qtmp2;
    Qt[i + Nm * (k + 1)] = s * Qtmp1 + c * Qtmp2;
  }

  //- Givens transformations
  for (int k1 = kmin; k1 < kmax - 1; ++k1) {
    double Fden1 = 1.0 / sqrt(x * x + TDb[k1 - 1] * TDb[k1 - 1]);
    double c1    = TDb[k1 - 1] * Fden1;
    double s1    = -x * Fden1;

    double tmpa1_inner = TDa[k1];
    double tmpa2_inner = TDa[k1 + 1];
    double tmpb_inner  = TDb[k1];
    TDa[k1]     = c1 * c1 * tmpa1_inner + s1 * s1 * tmpa2_inner - 2.0 * c1 * s1 * tmpb_inner;
    TDa[k1 + 1] = s1 * s1 * tmpa1_inner + c1 * c1 * tmpa2_inner + 2.0 * c1 * s1 * tmpb_inner;
    TDb[k1]     = c1 * s1 * (tmpa1_inner - tmpa2_inner) + (c1 * c1 - s1 * s1) * tmpb_inner;
    TDb[k1 - 1] = c1 * TDb[k1 - 1] - s1 * x;
    if (k1 != kmax - 2) {
      x          = -s1 * TDb[k1 + 1];
      TDb[k1 + 1] = c1 * TDb[k1 + 1];
    }

    for (int i = 0; i < Nk; ++i) {
      double Qtmp1 = Qt[i + Nm * k1];
      double Qtmp2 = Qt[i + Nm * (k1 + 1)];
      Qt[i + Nm * k1]       = c1 * Qtmp1 - s1 * Qtmp2;
      Qt[i + Nm * (k1 + 1)] = s1 * Qtmp1 + c1 * Qtmp2;
    }
  }
}


//====================================================================
//============================================================END=====
