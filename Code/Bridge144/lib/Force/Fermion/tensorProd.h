/*!
        @file    $Id:: tensorProd.h #$

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 18:35:38 #$

        @version $LastChangedRevision: 1571 $
*/

#ifndef TENSORPROD_INCLUDED
#define TENSORPROD_INCLUDED

#include "Field/field_F.h"

extern void BAPI tensorProd_Field_F(Field_G&, const Field_F&, const Field_F&);
extern void BAPI tensorProd_Field_F(Field_G&, const int mu, const Field_F&, const Field_F&);
#endif
