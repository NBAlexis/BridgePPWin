/*!
        @file    quarkNumberSusceptibility_Wilson.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/
#include "BridgeLib_Private.h"
#include "quarkNumberSusceptibility_Wilson.h"

const std::string QuarkNumberSusceptibility_Wilson::class_name = "QuarkNumberSusceptibility_Wilson";

//====================================================================
void QuarkNumberSusceptibility_Wilson::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  int Nnoise;

  int err = 0;
  err += params.fetch_int("number_of_noises", Nnoise);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Nnoise);
}


//====================================================================
void QuarkNumberSusceptibility_Wilson::set_parameters(const int Nnoise)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Nnoise = %d\n", Nnoise);

  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Nnoise);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  //- store values
  m_Nnoise = Nnoise;
}


//====================================================================
double QuarkNumberSusceptibility_Wilson::measure()
{
  // measurement of quark number susceptibility:
  //  tr1 = Tr[D1*Sq]
  //  tr2 = Tr[D2*Sq]
  //  tr3 = Tr[D1*Sq*D1*Sq]
  // For definition, see the implementation note.

  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Number of noise vector = %d\n", m_Nnoise);

  const int dir_t = CommonParameters::Ndim() - 1;

  const int Nex  = m_fopr->field_nex();
  const int Nvol = m_fopr->field_nvol();
  const int Nin  = m_fopr->field_nin();

  int    Nconv = 0;
  double diff  = 1.0;

  Field xi(Nin, Nvol, Nex);
  Field v1(Nin, Nvol, Nex), v2(Nin, Nvol, Nex), v3(Nin, Nvol, Nex);

  dcomplex tr1 = cmplx(0.0, 0.0);
  dcomplex tr2 = cmplx(0.0, 0.0);
  dcomplex tr3 = cmplx(0.0, 0.0);

  for (int i_noise = 0; i_noise < m_Nnoise; ++i_noise) {
    vout.general(m_vl, "  noise vector = %d\n", i_noise);

    m_nv->set(xi);
    // noise vector was set.

    m_fprop->invert_D(v3, xi, Nconv, diff);
    vout.general(m_vl, "    Nconv = %d  diff  = %.8e\n", Nconv, diff);
    // now v3 is  M^-1 * xi.

    // mult D_1 and D_2
    v1.set(0.0);
    m_fopr->mult_up(dir_t, v1, v3);
    v2.set(0.0);
    m_fopr->mult_dn(dir_t, v2, v3);

    dcomplex tr_c1 = dotc(xi, v1);
    dcomplex tr_c2 = dotc(xi, v2);

    dcomplex tr1_c = tr_c1 - tr_c2;
    dcomplex tr2_c = tr_c1 + tr_c2;

    v3 = v1;
    axpy(v3, -1.0, v2);
    // now v3 is  D_1 * M^-1 * xi.

    v1 = v3;
    m_fprop->invert_D(v3, v1, Nconv, diff);
    vout.general(m_vl, "    Nconv = %d  diff  = %.8e\n", Nconv, diff);
    // now v3 is  M^-1 * D_1 * M^-1 * xi.

    // mult D_1
    v1.set(0.0);
    m_fopr->mult_up(dir_t, v1, v3);
    v2.set(0.0);
    m_fopr->mult_dn(dir_t, v2, v3);
    axpy(v1, -1.0, v2);
    // now v1 is  D_1 * M^-1 * D_1 * M^-1 * xi.

    dcomplex tr3_c = dotc(xi, v1);

    vout.general(m_vl, "    tr1 = (%f,%f)\n", real(tr1_c), imag(tr1_c));
    vout.general(m_vl, "    tr2 = (%f,%f)\n", real(tr2_c), imag(tr2_c));
    vout.general(m_vl, "    tr3 = (%f,%f)\n", real(tr3_c), imag(tr3_c));

    tr1 += tr1_c;
    tr2 += tr2_c;
    tr3 += tr3_c;
  }

  tr1 = tr1 / cmplx(double(m_Nnoise), 0.0);
  tr2 = tr2 / cmplx(double(m_Nnoise), 0.0);
  tr3 = tr3 / cmplx(double(m_Nnoise), 0.0);

  vout.general(m_vl, "  averaged over noise vector:\n");
  vout.general(m_vl, "    tr1 = (%f,%f)\n", real(tr1), imag(tr1));
  vout.general(m_vl, "    tr2 = (%f,%f)\n", real(tr2), imag(tr2));
  vout.general(m_vl, "    tr3 = (%f,%f)\n", real(tr3), imag(tr3));


  return real(tr1);
}


//====================================================================
//============================================================END=====
