#include "BridgeLib_Private.h"

/*!
        @file    $Id:: randomNumberManager.cpp #$

        @brief

        @author  Tatsumi Aoyama (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/

#include "randomNumberManager.h"
#include "randomNumbers.h"
#include "randomNumbers_Mseries.h"
#include "randomNumbers_MT19937.h"

#ifdef USE_SFMTLIB
#include "randomNumbers_SFMT.h"
#endif

const std::string RandomNumberManager::class_name = "RandomNumberManager";

//================================================================
RandomNumbers *RandomNumberManager::s_rand = NULL;

//================================================================
RandomNumbers *RandomNumberManager::factory(const std::string& rng_type, unsigned long seed)
{
  if ((rng_type == "Mseries") || (rng_type == "mseries")) {
    return new RandomNumbers_Mseries(seed);
  } else if ((rng_type == "MT19937") || (rng_type == "mt19937")) {
    return new RandomNumbers_MT19937(seed);

#ifdef USE_SFMTLIB
  } else if (rng_type == "SFMT") {
    return new RandomNumbers_SFMT(seed);
#endif
  } else {
    vout.crucial("%s: unsupported random number generator: %s\n", class_name.c_str(), rng_type.c_str());
    return NULL;
  }
}


//================================================================
RandomNumbers *RandomNumberManager::getInstance()
{
  if (!s_rand) {
    vout.crucial("Error: %s: uninitialized.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  return s_rand;
}


//================================================================
bool RandomNumberManager::initialize(const std::string& rng_type, unsigned long seed)
{
  vout.detailed("%s: initialize: type = \"%s\", seed = %lu\n", class_name.c_str(), rng_type.c_str(), seed);

  if (s_rand) {
    vout.crucial("Error: %s: already initialized.\n", class_name.c_str());

    exit(EXIT_FAILURE);
    //return false;
  }

  s_rand = factory(rng_type, seed);

  if (!s_rand) {
    vout.crucial("ERROR: %s: initialize failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  return s_rand != NULL;
}


//================================================================
void RandomNumberManager::finalize()
{
  if (s_rand) {
    delete s_rand;
    s_rand = NULL;
  }
}


//================================================================
void RandomNumberManager::reset(unsigned long seed)
{
  if (!s_rand) {
    vout.crucial("Warning: %s: uninitialized.\n", class_name.c_str());
  }

  if (s_rand) {
    s_rand->reset(seed);
  }
}


//================================================================
void RandomNumberManager::restore_state(const std::string& filename)
{
  if (!s_rand) {
    vout.crucial("Warning: %s: uninitialized.\n", class_name.c_str());
  }

  if (s_rand) {
    s_rand->read_file(filename);
  }
}


//================================================================
void RandomNumberManager::save_state(const std::string& filename)
{
  if (!s_rand) {
    vout.crucial("Warning: %s: uninitialized.\n", class_name.c_str());
  }

  if (s_rand) {
    s_rand->write_file(filename);
  }
}


//================================================================
RandomNumbers *RandomNumberManager::New(const std::string& rng_type, unsigned long seed)
{
  return factory(rng_type, seed);
}


//================================================================
//================================================================
