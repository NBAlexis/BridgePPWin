#include "BridgeLib_Private.h"

#if COMMUNICATOR_USE_BGNET

/*!
        @file    $Id: counter.cpp #$

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate: 2012-11-15 16:42:21 #$

        @version $LastChangedRevision: 1561 $
*/

#include "counter.h"

#include "libkek.h"

#include "Parameters/commonParameters.h"
#include "IO/bridgeIO.h"
using Bridge::vout;

int Counter::cid_static = 10;
//====================================================================
Counter::Counter()
{
  m_id = cid_static;
  ++cid_static;

  Bridge::VerboseLevel vl = CommonParameters::Vlevel();
  vout.general(vl, "Counter_BGNET constructed with id = %d\n", m_id);
}


//====================================================================
void Counter::start()
{
  KEK_FopCountStart(m_id);
}


//====================================================================
void Counter::finish()
{
  double time, gflops;

  finish(time, gflops);

  Bridge::VerboseLevel vl = CommonParameters::Vlevel();
  vout.general(vl, "  time = %f  GFlops = %f\n", time, gflops);
}


//====================================================================
void Counter::finish(double& time, double& gflops)
{
  unsigned long count_op;
  double        time_op;

  KEK_FopCountFinish(m_id, &count_op, &time_op);

  time   = time_op;
  gflops = 1.e-9 * double(count_op) / time_op;
}


//====================================================================
//============================================================END=====

#endif //COMMUNICATOR_USE_BGNET
