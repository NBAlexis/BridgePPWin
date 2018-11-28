#include "BridgeLib_Private.h"

/*!
        @file    $Id:: source_Exponential.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "source_Exponential.h"



#ifdef USE_FACTORY
namespace {
  Source *create_object()
  {
    return new Source_Exponential();
  }


  bool init = Source::Factory::Register("Exponential", create_object);
}
#endif



const std::string Source_Exponential::class_name = "Source_Exponential";

//====================================================================
void Source_Exponential::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  std::vector<int> source_position;
  double           slope, power;

  int err = 0;
  err += params.fetch_int_vector("source_position", source_position);
  err += params.fetch_double("slope", slope);
  err += params.fetch_double("power", power);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(source_position, slope, power);
}


//====================================================================
void Source_Exponential::set_parameters(const std::vector<int>& source_position,
                                        const double slope, const double power)
{
  // ####  parameter setup  ####
  int Ndim = CommonParameters::Ndim();

  //- global lattice size
  std::vector<int> Lsize(Ndim);

  Lsize[0] = CommonParameters::Lx();
  Lsize[1] = CommonParameters::Ly();
  Lsize[2] = CommonParameters::Lz();
  Lsize[3] = CommonParameters::Lt();

  //- local size
  std::vector<int> Nsize(Ndim);
  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();
  Nsize[3] = CommonParameters::Nt();

  //- print input parameters
  vout.general(m_vl, "Source for spinor field - exponential smeared:\n");
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  source_position[%d] = %d\n",
                 mu, source_position[mu]);
  }
  vout.general(m_vl, "  slope = %12.6f\n", slope);
  vout.general(m_vl, "  power = %12.6f\n", power);

  //- range check
  int err = 0;
  for (int mu = 0; mu < Ndim; ++mu) {
    // NB. Lsize[mu] > abs(source_position[mu])
    err += ParameterCheck::non_negative(Lsize[mu] - abs(source_position[mu]));
  }
  // NB. slope,power == 0 is allowed.

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  assert(source_position.size() == Ndim);

  //- store values
  m_source_position.resize(Ndim);
  for (int mu = 0; mu < Ndim; ++mu) {
    m_source_position[mu] = (source_position[mu] + Lsize[mu]) % Lsize[mu];
  }

  m_slope = slope;
  m_power = power;


  //- post-process
  const int Lvol3 = Lsize[0] * Lsize[1] * Lsize[2];
  const int Nvol3 = Nsize[0] * Nsize[1] * Nsize[2];

  m_src_func.reset(1, Nvol3, 1);
  //m_src_func = 0.0;
  m_src_func.set(0.0);

  //- PE location in t-direction.
  int tpe = m_source_position[3] / Nsize[3];

  m_in_node = false;

  if (tpe == Communicator::ipe(3)) {
    m_in_node = true;
  }

  //- fill exponential.
  //   center at source_position,
  //   range -L+1 < (x - x0) < L-1,
  //   tail folded.
  for (int z = -Lsize[2] + 1; z < Lsize[2]; ++z) {
    for (int y = -Lsize[1] + 1; y < Lsize[1]; ++y) {
      for (int x = -Lsize[0] + 1; x < Lsize[0]; ++x) {
        //- global position
        int z2 = (m_source_position[2] + z + Lsize[2]) % Lsize[2];
        int y2 = (m_source_position[1] + y + Lsize[1]) % Lsize[1];
        int x2 = (m_source_position[0] + x + Lsize[0]) % Lsize[0];

        //- PE location
        int xpe = x2 / Nsize[0];
        int ype = y2 / Nsize[1];
        int zpe = z2 / Nsize[2];

        //- local position
        int xl = x2 % Nsize[0];
        int yl = y2 % Nsize[1];
        int zl = z2 % Nsize[2];

        if (
          (xpe == Communicator::ipe(0)) &&
          (ype == Communicator::ipe(1)) &&
          (zpe == Communicator::ipe(2)) &&
          (tpe == Communicator::ipe(3)))
        {
          double r    = sqrt((double)(x * x + y * y + z * z));
          double ex   = pow(r, m_power);
          double expf = exp(-m_slope * ex);

          int lsite = xl + Nsize[0] * (yl + Nsize[1] * zl);

          m_src_func.add(0, lsite, 0, expf);
        }
      }
    }
  }

  //- normalize
  double Fnorm = 0.0;
  for (int i = 0; i < Nvol3; ++i) {
    Fnorm += m_src_func.cmp(i) * m_src_func.cmp(i);
  }
  double Fnorm_global = Communicator::reduce_sum(Fnorm);

  //m_src_func *= 1.0 / sqrt(Fnorm_global);
  scal(m_src_func, 1.0 / sqrt(Fnorm_global));

  //- check normalization
  double epsilon_criterion = CommonParameters::epsilon_criterion();

  Fnorm = 0.0;
  for (int i = 0; i < Nvol3; i++) {
    Fnorm += m_src_func.cmp(i) * m_src_func.cmp(i);
  }
  Fnorm_global = Communicator::reduce_sum(Fnorm);

  // assert(abs(sqrt(Fnorm_global) - 1.0) < epsilon_criterion);
  //- Fluctuation of order of sqrt(NPE) is acceptable. [21 May 2014]
  double epsilon_criterion2 = epsilon_criterion *
                              sqrt((double)CommonParameters::NPE());
  //  assert(abs(sqrt(Fnorm_global) - 1.0) < epsilon_criterion2);
  if (abs(sqrt(Fnorm_global) - 1.0) > epsilon_criterion2) {
    vout.crucial(m_vl, "%s: norm criterion is not satisfied.\n",
                 class_name.c_str());
    vout.crucial(m_vl, "  |sqrt(Fnorm_global) -1| = %e.\n",
                 abs(sqrt(Fnorm_global) - 1.0));
    vout.crucial(m_vl, "  epsilon_criterion2      = %e.\n",
                 epsilon_criterion2);
  }
}


//====================================================================
void Source_Exponential::set(Field& src, int j)
{
  int Ndim = CommonParameters::Ndim();

  std::vector<int> Nsize(Ndim);

  Nsize[0] = CommonParameters::Nx();
  Nsize[1] = CommonParameters::Ny();
  Nsize[2] = CommonParameters::Nz();
  Nsize[3] = CommonParameters::Nt();

  //- clear field
  src.set(0.0);

  if (m_in_node) {
    int t = m_source_position[3] % Nsize[3];

    for (int z = 0; z < Nsize[2]; ++z) {
      for (int y = 0; y < Nsize[1]; ++y) {
        for (int x = 0; x < Nsize[0]; ++x) {
          int lsite = x + Nsize[0] * (y + Nsize[1] * z);

          int isite = m_index.site(x, y, z, t);

          //XXX field layout: complex as two doubles
          src.set(2 * j, isite, 0, m_src_func.cmp(0, lsite, 0));
        }
      }
    }
  }
}


//====================================================================
//============================================================END=====
