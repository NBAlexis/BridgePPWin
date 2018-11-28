//=============================================================================
// FILENAME : BridgeLib.h
// 
// DESCRIPTION:
//
// REVISION:
//  [11/28/2018 nbale]
//=============================================================================


#pragma once

#ifndef _BRIDGELIB_H_
#define _BRIDGELIB_H_

#if defined(_VCWIN)
#   if !defined(Z2API)
#	define __LIB_TITLE__	"BridgeLib"
#	define Z2API __DLL_IMPORT
#	ifdef DEBUG
#		define __LIB_FILE__	__LIB_TITLE__ "_d.lib"
#	else
#		define __LIB_FILE__ __LIB_TITLE__ ".lib"
#	endif
#	pragma __IMPORT_LIB(__LIB_FILE__)
#	pragma message("linking with " __LIB_FILE__ "...")
#	undef __LIB_FILE__
#	undef __LIB_TITLE__
#   endif
#else
#	define BRIDGELIBAPI  
#endif

#define _BRIDGE_LIB_PUBLIC 1
#include "BridgeLib_Private.h"


#endif //#ifndef _BRIDGELIB_H_
//=============================================================================
// END OF FILE
//=============================================================================