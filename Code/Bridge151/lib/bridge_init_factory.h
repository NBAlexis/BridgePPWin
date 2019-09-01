/*!
        @file    bridge_init_factory.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef BRIDGE_INIT_FACTORY_INCLUDED
#define BRIDGE_INIT_FACTORY_INCLUDED


#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

bool bridge_init_factory();

#endif /* USE_FACTORY_AUTOREGISTER */

#ifdef DEBUG
void bridge_report_factory();

#endif

#endif /* USE_FACTORY */


#endif /* BRIDGE_INIT_FACTORY_INCLUDED */
