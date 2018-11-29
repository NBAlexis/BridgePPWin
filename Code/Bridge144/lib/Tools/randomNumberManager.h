/*!
@file    $Id:: randomNumberManager.h #$

@brief

@author  Tatsumi Aoyama (aoym)
$LastChangedBy: aoym $

@date    $LastChangedDate:: 2017-02-24 00:49:46 #$

@version $LastChangedRevision: 1561 $
*/

#ifndef RANDOMNUMBER_MANAGER_INCLUDED
#define RANDOMNUMBER_MANAGER_INCLUDED

#include "randomNumbers.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Manager class for RandomNumbers

/*!
This class manages a global instance of RandomNumbers class.
[7 July 2016 T.Aoyama]

Usage:

RandomNumberManager::initialize(rng_type, seed);
-- initialize global RNG by creating an instance of rng_type with
the specified seed.

RandomNumbers* rand = RandomNumberManager::getInstance();
-- fetch the global RNG. the program terminates unless initialized.

RandomNumberManager::reset(seed);
-- reset the state of the global RNG.

RandomNumberManager::save_state(filename);
RandomNumberManager::restore_state(filename);
-- save/restore the state of the global RNG to/from the file.

RandomNumberManager::finalize();
-- finalize the global RNG by releasing the instance.


RandomNumbers* rand = RandomNumberManager::New(rng_type, seed);
-- factory for RandomNumbers class. an instance of rng_type
is new'ed and returned.
*/

class BAPI RandomNumberManager
{
public:
    static const std::string class_name;

protected:
    Bridge::VerboseLevel m_vl;

public:

    RandomNumberManager()
        : m_vl(CommonParameters::Vlevel()) {}

    virtual ~RandomNumberManager() {}

    // class methods
    static RandomNumbers *getInstance();

    static bool initialize(const std::string& rng_type, unsigned long seed);
    static void finalize();

    static void reset(unsigned long seed);

    static void restore_state(const std::string& filename);
    static void save_state(const std::string& filename);

    // factory
    static RandomNumbers *New(const std::string& rng_type, unsigned long seed);

private:
    // non-copyable
    RandomNumberManager(const RandomNumberManager&);
    RandomNumberManager& operator=(const RandomNumberManager&);

public:
    void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

private:
    static RandomNumbers *s_rand;

    static RandomNumbers *factory(const std::string& rng_type, unsigned long seed);
};
#endif
