//=============================================================================
// FILENAME : RandomNumbers_GoldenratioNoise.h
// 
// DESCRIPTION:
// DO NOT USE THIS!!!!!! TOO SLOW AND TOO BAD!!!!!!
// REVISION:
//  [11/30/2018 nbale]
//=============================================================================

#ifndef RANDOMNUMBERS_SCHRAGE_INCLUDED
#define RANDOMNUMBERS_SCHRAGE_INCLUDED

#include <string>
#include <cassert>

#include "randomNumbers.h"


#include "IO/bridgeIO.h"
using Bridge::vout;

class BAPI RandomNumbers_Schrage : public RandomNumbers
{
    static const std::string class_name;
    const double AM = (1.0 / 4294967296UL);
    unsigned int m_ulSeed;   

public:
    RandomNumbers_Schrage(const int s) : m_ulSeed(s) {;}
    RandomNumbers_Schrage(const std::string& filename) { read_file(filename); }

    ~RandomNumbers_Schrage() {}

    double get() 
    {
        m_ulSeed = (1664525UL * m_ulSeed + 1013904223UL) & 0xffffffff;
        return AM * m_ulSeed;
    }

    void write_file(const std::string& filename);
    void read_file(const std::string& filename);

    void reset(unsigned long seed) { m_ulSeed = seed; }
};

#endif /* RANDOMNUMBERS_SCHRAGE_INCLUDED */
