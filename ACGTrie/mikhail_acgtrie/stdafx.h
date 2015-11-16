// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

//#include <Windows.h>
//#undef min
//#undef max
#include <assert.h>
//#include <conio.h>
//#include <iostream>
//#include <map>
//#include <tchar.h>
//#include <time.h>
//#include <stdio.h>
//#include <functional>
//#include <list>
#include <string>
#include <vector>

//#if defined(_DEBUG) || defined(CHECK_ASSERTS)
//  #undef assert
//  #define assert(exp) (void)( (exp) || (THROW_EXCEPTION("Assertion failed: "	#exp), 0) )
//  //#define assert(exp) (void)( (exp) || (THROW_EXCEPTION(FormatStr("Assertion failed: %s (file %s, line %d)",\
//		//	#exp, __FILE__, __LINE__)), 0) )
//#endif

#define THROW_EXCEPTION(mess) throw std::runtime_error(mess)

