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

#if _VCWIN
# define __DLL_IMPORT			__declspec(dllimport)
# define Z2PRIVATE
# define __DLL_EXPORT			__declspec(dllexport)
# define __IMPORT_LIB(libname)	comment(lib, libname)
# undef FORCEINLINE
# undef CDECL
# define FORCEINLINE             __forceinline
# define CDECL                   __cdecl

# define SUPPORTS_PRAGMA_PACK 1
# define __PACK_PUSH				pack(push, 8)
# define __PACK_POP				pack(pop)
#else
# define __DLL_IMPORT			
# define CCPRIVATE
# define __DLL_EXPORT			
# define __IMPORT_LIB(libname)	
# undef FORCEINLINE
# undef CDECL
# define FORCEINLINE             __forceinline
# define CDECL                   __cdecl

# define SUPPORTS_PRAGMA_PACK 0
# define __PACK_PUSH			
# define __PACK_POP				
#endif

#if defined(_VCWIN)
#   if !defined(BAPI)
#	define __LIB_TITLE__	"BridgeLib"
#	define BAPI __DLL_IMPORT
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
#	define BAPI  
#endif

#define _BRIDGE_LIB_PUBLIC 1
#include "BridgeLib_Private.h"


#endif //#ifndef _BRIDGELIB_H_
//=============================================================================
// END OF FILE
//=============================================================================