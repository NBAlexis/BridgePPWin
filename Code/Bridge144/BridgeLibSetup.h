//=============================================================================
// FILENAME : BridgeLibSetup.h
// 
// DESCRIPTION:
//
// REVISION:
//  [11/28/2018 nbale]
//=============================================================================


#pragma once

#ifndef _BRIDGELIBSETUP_H_
#define _BRIDGELIBSETUP_H_

//define many things here
#define BRIDGE_WIN 1

//choose complex
#define USE_STD_COMPLEX 1

//choose communicator
//COMMUNICATOR_USE_BGNET
//COMMUNICATOR_USE_MPI
//COMMUNICATOR_USE_SINGLE
//COMMUNICATOR_USE_GPU
#define COMMUNICATOR_USE_MPI 1

//Improvements
//USE_ORG
//USE_IMP
//USE_IMP_BGQ
//USE_GROUP_SU2
//USE_GROUP_SU3
//USE_GROUP_SU_N
#define USE_GROUP_SU3 1
#define USE_IMP 1

//USE_RANDOM_MSERIES
//USE_RANDOM_MT19937
//USE_RANDOM_SFMT
#define USE_RANDOM_MT19937 1

#define USE_FACTORY


#endif //#ifndef _BRIDGELIBSETUP_H_
//=============================================================================
// END OF FILE
//=============================================================================