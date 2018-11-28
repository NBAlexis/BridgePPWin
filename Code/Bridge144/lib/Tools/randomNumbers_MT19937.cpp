#include "BridgeLib_Private.h"

/*!
        @file    $Id:: randomNumbers_MT19937.cpp #$

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "randomNumbers_MT19937.h"

const std::string RandomNumbers_MT19937::class_name = "RandomNumbers_MT19937";

//====================================================================
RandomNumbers_MT19937::RandomNumbers_MT19937(int s)
  : m_left(1)
{
  if (s < 0) {
    vout.crucial(m_vl, "Error at %s: seed must be positive.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }
  vout.general(m_vl, "%s: seed = %d\n", class_name.c_str(), s);
  init(s);
}


RandomNumbers_MT19937::RandomNumbers_MT19937(unsigned long s)
  : m_left(1)
{
  vout.general(m_vl, "%s: seed = %ul\n", class_name.c_str(), s);
  init(s);
}


RandomNumbers_MT19937::RandomNumbers_MT19937(vector<unsigned long>& key)
  : m_left(1)
{
  vout.general(m_vl, "%s: seeds = (", class_name.c_str());
  for (int i = 0; i < key.size(); ++i) {
    vout.general(m_vl, "%ul, ", key[i]);
  }
  vout.general(m_vl, ")\n");

  init(19650218UL, key);
}


RandomNumbers_MT19937::RandomNumbers_MT19937(const string& filename)
{
  read_file(filename);
}


//====================================================================
void RandomNumbers_MT19937::read_file(const string& filename)
{
  vout.detailed(m_vl, "%s: read from file %s\n", class_name.c_str(), filename.c_str());

  if (Communicator::is_primary()) {
    std::ifstream in_file(filename.c_str());
    if (!in_file) {
      vout.crucial(m_vl, "Error at %s: error: unable to open input file.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    unsigned long int n;
    for (int i = 0; i < N; ++i) {
      if (in_file >> n) {
        m_state[i] = n;
      } else {
        vout.crucial(m_vl, "Error at %s: seed data in input file is missing.\n", class_name.c_str());
        exit(EXIT_FAILURE);
      }
    }

    if (in_file >> n) {
      m_left = n;
    } else {
      vout.crucial(m_vl, "Error at %s: counter data in input file is missing.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    in_file.close();
  }

  Communicator::Base::broadcast(sizeof(int), &m_left, 0);
  Communicator::Base::broadcast(sizeof(unsigned long) * N, m_state, 0);

  m_next = m_state + N - m_left + 1;
}


//====================================================================
void RandomNumbers_MT19937::write_file(const string& filename)
{
  vout.detailed(m_vl, "%s: save random number state to file %s\n", class_name.c_str(), filename.c_str());

  if (Communicator::is_primary()) {
    std::ofstream out_file(filename.c_str());

    if (!out_file) {
      vout.crucial(m_vl, "Error at %s: unable to open output file.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; ++i) {
      out_file << m_state[i] << std::endl;
    }

    out_file << m_left << std::endl;

    out_file.close();
  }
}


//====================================================================
void RandomNumbers_MT19937::reset(unsigned long seed)
{
  return init(seed);
}


//====================================================================
void RandomNumbers_MT19937::init(unsigned long s)
{
  m_state[0] = s & 0xffffffffUL;
  for (int j = 1; j < N; ++j) {
    m_state[j]  = (1812433253UL * (m_state[j - 1] ^ (m_state[j - 1] >> 30)) + j);
    m_state[j] &= 0xffffffffUL;
  }
}


void RandomNumbers_MT19937::init(unsigned long s, vector<unsigned long>& key)
{
  init(s);
  int i          = 1;
  int j          = 0;
  int key_length = (int)key.size();

  for (int k = N > key_length ? N : key_length; k; --k) {
    // non linear
    m_state[i]
      = (m_state[i] ^ ((m_state[i - 1] ^ (m_state[i - 1] >> 30)) * 1664525UL)) + key[j] + j;
    // for WORDSIZE > 32 machines
    m_state[i] &= 0xffffffffUL;

    ++i;
    ++j;
    if (i >= N) {
      m_state[0] = m_state[N - 1];
      i          = 1;
    }
    if (j >= key_length) j = 0;
  }

  for (int k = N - 1; k; --k) {
    // non linear
    m_state[i]  = (m_state[i] ^ ((m_state[i - 1] ^ (m_state[i - 1] >> 30)) * 1566083941UL)) - i;
    m_state[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
    i++;
    if (i >= N) {
      m_state[0] = m_state[N - 1];
      i          = 1;
    }
  }
  // MSB is 1; assuring non-zero initial array
  m_state[0] = 0x80000000UL;
}


//====================================================================
void RandomNumbers_MT19937::nextState() const
{
  unsigned long *p = m_state;

  m_left = N;
  m_next = m_state;
  for (int j = N - M + 1; --j; ++p) {
    *p = p[M] ^ twist(p[0], p[1]);
  }
  for (int j = M; --j; ++p) {
    *p = p[M - N] ^ twist(p[0], p[1]);
  }

  *p = p[M - N] ^ twist(p[0], m_state[0]);
}


unsigned long
RandomNumbers_MT19937::twist(unsigned long u, unsigned long v) const
{
  // Period parameters
  const unsigned long mtrx_a = 0x9908b0dfUL; // constant vector a
  const unsigned long umask  = 0x80000000UL; // most significant w-r bits
  const unsigned long lmask  = 0x7fffffffUL; // least significant r bits

  unsigned long maxbits = (u & umask) | (v & lmask);

  return (maxbits >> 1) ^ (v & 1UL ? mtrx_a : 0UL);
}


//====================================================================
// generates a random number on [0,0xffffffff]
unsigned long RandomNumbers_MT19937::randInt32() const
{
  if (--m_left == 0) nextState();
  unsigned long y = *m_next++;

  // Tempering
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);
  return y;
}


//generates a random number on [0,0x7fffffff]
long RandomNumbers_MT19937::randInt31() const
{
  return (long)(randInt32() >> 1);
}


//====================================================================
//generates a random number on [0,1] double number
double RandomNumbers_MT19937::randDouble1() const
{
  static const double factor = 1.0 / 4294967295.0;

  return static_cast<double>(randInt32() * factor);
}


//generates a random number on [0,1) double number
double RandomNumbers_MT19937::randDouble2() const
{
  static const double factor = 1.0 / 4294967296.0;

  return static_cast<double>(randInt32() * factor);
}


// generates a random number on (0,1) double number
double RandomNumbers_MT19937::randDouble3() const
{
  static const double factor = 1.0 / 4294967296.0;

  return static_cast<double>((randInt32() + 0.5) * factor);
}


// generates a random number on [0,1) with 53-bit resolution
double RandomNumbers_MT19937::randRes53() const
{
  unsigned long       a      = randInt32() >> 5;
  unsigned long       b      = randInt32() >> 6;
  static const double factor = 1.0 / 9007199254740992.0;

  return (a * 67108864.0 + b) * factor;
}


//====================================================================
//============================================================END=====
