#include "BridgeLib_Private.h"

/*!
        @file    $Id:: randomNumbers_Mseries.cpp #$

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "randomNumbers_Mseries.h"

const std::string RandomNumbers_Mseries::class_name = "RandomNumbers_Mseries";
const double      RandomNumbers_Mseries::Fnorm      = 4.656612870908988e-10;

//====================================================================
void RandomNumbers_Mseries::reset(unsigned long seed)
{
  initset((int)seed);
}


//====================================================================
void RandomNumbers_Mseries::initset(const int ndelay)
{
  vout.general(m_vl, " Random number init: ndelay =%10d\n", ndelay);

  int IP = Np;
  int IQ = Nq;

  int      MACRM = 40;
  int      MACRI = 1;
  unsigned IB[IP];

  unsigned IX = MACRI;
  for (int i = 0; i < IP; ++i) {
    IX    = IX * 69069;
    IB[i] = IX >> 31;
  }

  int JR = 0;
  int KR = IP - IQ;
  for (int j = 0; j < IP; ++j) {
    unsigned iwork = 0;

    for (int i = 0; i < 32; ++i) {
      iwork  = iwork * 2 + IB[JR];
      IB[JR] = IB[JR] ^ IB[KR];
      ++JR;
      if (JR == IP) JR = 0;
      ++KR;
      if (KR == IP) KR = 0;
    }

    w[j] = iwork >> 1;
  }

  jr = JR;
  kr = KR;

  delay3(ndelay);
}


//====================================================================
void RandomNumbers_Mseries::delay3(const int ndelay)
{
  int IP = Np;
  int IQ = Nq;

  int LAMBDA = ndelay;

  int MACRM = 40;
  int MACRI = 1;

  unsigned IWK[2 * IP - 1];
  unsigned C[2 * IP];
  unsigned IB[IP + 32];

  int MU = MACRM;

  for (int i = 0; i < IP; ++i) {
    IWK[i] = w[i];
  }

  for (int i = IP; i < 2 * IP - 1; ++i) {
    IWK[i] = IWK[i - IP] ^ IWK[i - IQ];
  }

  for (int i = 0; i < MU; ++i) {
    IB[i] = 0;
  }
  int M  = LAMBDA;
  int NB = MU - 1;

continued220:
  if (M <= IP - 1) goto continued300;
  ++NB;
  IB[NB] = M % 2;
  M      = M / 2;
  goto continued220;

continued300:
  for (int i = 0; i < IP; ++i) {
    C[i] = 0;
  }

  C[M] = 1;
  for (int j = NB; j >= 0; --j) {
    for (int i = IP - 1; i >= 0; --i) {
      C[2 * i + IB[j]]     = C[i];
      C[2 * i + 1 - IB[j]] = 0;
    }
    for (int i = 2 * IP - 1; i >= IP; --i) {
      C[i - IP] = C[i - IP] ^ C[i];
      C[i - IQ] = C[i - IQ] ^ C[i];
    }
  }

  for (int j = 0; j < IP; ++j) {
    unsigned iwork = 0;
    for (int i = 0; i < IP; ++i) {
      iwork = iwork ^ (C[i] * IWK[j + i]);
    }
    w[j] = iwork;
  }
}


//====================================================================
void RandomNumbers_Mseries::write_file(const std::string& filename)
{
  vout.general(m_vl, "%s: write down to file = %s\n", class_name.c_str(), filename.c_str());

  if (Communicator::is_primary()) {
    std::ofstream outfile;
    outfile.open(filename.c_str());

    if (!outfile) {
      vout.crucial(m_vl, "Error at %s: unable to open output file.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    outfile << " " << jr << " " << kr << std::endl;

    for (int i = 0; i < Np; i++) {
      outfile << " " << w[i] << std::endl;
    }

    outfile.close();
  }
}


//====================================================================
void RandomNumbers_Mseries::read_file(const std::string& filename)
{
  vout.general(m_vl, "%s: read from file = %s\n", class_name.c_str(), filename.c_str());

  if (Communicator::is_primary()) {
    std::ifstream infile;
    infile.open(filename.c_str());

    if (!infile) {
      vout.crucial(m_vl, "Error at %s: unable to open input file.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    infile >> jr >> kr;

    for (int i = 0; i < Np; i++) {
      infile >> w[i];
    }

    infile.close();
  }

  Communicator::broadcast(Np, w, 0);
  Communicator::broadcast(1, &jr, 0);
  Communicator::broadcast(1, &kr, 0);
}


//====================================================================
//============================================================END=====
