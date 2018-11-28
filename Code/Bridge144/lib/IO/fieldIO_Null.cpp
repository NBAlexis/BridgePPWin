#include "BridgeLib_Private.h"

/*!
        @file    $Id: fieldIO_Null.cpp #$

        @brief

        @author  $Tatsumi Aoyama (aoym)
                 LastChangedBy: aoym $

        @date    $LastChangedDate: 2014-10-09 17:20:36 #$

        @version $LastChangedRevision: 1561 $
*/

#include "fieldIO_Null.h"

#include "bridgeIO.h"
using Bridge::vout;

const std::string FieldIO_Null::class_name = "FieldIO_Null";

//====================================================================
void FieldIO_Null::read_file(Field *v, string filename)
{
  if (v == 0) {
    vout.crucial(m_vl, "%s: field empty.\n", __func__);
    return;
  }

  vout.detailed(m_vl, "read successful: nothing performed.\n");
}


//====================================================================
void FieldIO_Null::write_file(Field *v, std::string filename)
{
  if (v == 0) {
    vout.crucial(m_vl, "%s: field empty.\n", __func__);
    return;
  }

  vout.detailed(m_vl, "write successful: data discarded.\n");
}


//====================================================================
//============================================================END=====
