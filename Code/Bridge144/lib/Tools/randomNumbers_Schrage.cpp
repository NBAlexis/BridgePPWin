//=============================================================================
// FILENAME : RandomNumbers_GoldenratioNoise.h
// 
// DESCRIPTION:
//
// REVISION:
//  [11/30/2018 nbale]
//=============================================================================

#include "BridgeLib_Private.h"

const std::string RandomNumbers_Schrage::class_name = "RandomNumbers_Schrage";

//====================================================================
void RandomNumbers_Schrage::write_file(const std::string& filename)
{
    vout.general(m_vl, "%s: write down to file = %s\n", class_name.c_str(), filename.c_str());

    if (Communicator::is_primary()) 
    {
        std::ofstream outfile;
        outfile.open(filename.c_str());

        if (!outfile) 
        {
            vout.crucial(m_vl, "Error at %s: unable to open output file.\n", class_name.c_str());
            exit(EXIT_FAILURE);
        }

        outfile << " " << m_ulSeed << std::endl;
        outfile.close();
    }
}


//====================================================================
void RandomNumbers_Schrage::read_file(const std::string& filename)
{
    vout.general(m_vl, "%s: read from file = %s\n", class_name.c_str(), filename.c_str());

    if (Communicator::is_primary()) 
    {
        std::ifstream infile;
        infile.open(filename.c_str());

        if (!infile) 
        {
            vout.crucial(m_vl, "Error at %s: unable to open input file.\n", class_name.c_str());
            exit(EXIT_FAILURE);
        }

        infile >> m_ulSeed;
        infile.close();
    }
    int tobebroadcast = (int)m_ulSeed;
    Communicator::broadcast(1, &tobebroadcast, 0);
}

//====================================================================
//============================================================END=====
