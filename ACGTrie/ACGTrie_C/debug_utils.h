#pragma once

#include <map>

#ifndef TRACE_DETAIL_LEVEL
#  ifdef _DEBUG
#    define TRACE_DETAIL_LEVEL 4
#  else
//#    define FUNC_TIME_MEASUREMENT //d_
#    ifdef FUNC_TIME_MEASUREMENT
#      define TRACE_DETAIL_LEVEL 3
#    else
#      define TRACE_DETAIL_LEVEL 2
#    endif
#  endif
#endif

#ifndef DBG_UTILS_TRACE
// #  define DBG_UTILS_TRACE dbg_utils::trace_debug_output
// #  define DBG_UTILS_TRACE dbg_utils::trace_file
#  define DBG_UTILS_TRACE dbg_utils::trace_console
	//d_ Writing the log to file
#endif

#define TRACE_ERROR_RELEASE(X) DBG_UTILS_TRACE(dbg_utils::tdl_Error, __FUNCTION__, dbg_utils::format_msg X)

#if TRACE_DETAIL_LEVEL >= 1
#  define TRACE_ERROR(X) DBG_UTILS_TRACE(dbg_utils::tdl_Error, __FUNCTION__, dbg_utils::format_msg X)
#else
#  define TRACE_ERROR(X)
#endif

#if TRACE_DETAIL_LEVEL >= 2
#  define TRACE_NORMAL(X) DBG_UTILS_TRACE(dbg_utils::tdl_Normal, __FUNCTION__, dbg_utils::format_msg X)
#else
#  define TRACE_NORMAL(X)
#endif

#if TRACE_DETAIL_LEVEL >= 3
#  define TRACE_VERBOSE(X) DBG_UTILS_TRACE(dbg_utils::tdl_Verbose, __FUNCTION__, dbg_utils::format_msg X)
#else
#  define TRACE_VERBOSE(X)
#endif

#if TRACE_DETAIL_LEVEL >= 4
#  define TRACE_DEBUG(X) DBG_UTILS_TRACE(dbg_utils::tdl_Debug, __FUNCTION__, dbg_utils::format_msg X)
#else
#  define TRACE_DEBUG(X)
#endif

typedef const char *LPCSTR;


namespace dbg_utils
{
  enum TraceDetailLevel
  {
    tdl_Error,
    tdl_Normal,
    tdl_Verbose,
    tdl_Debug
  };

  std::string format_msg(const char* format, ...);
  std::string format_msg(const wchar_t* format, ...);
  std::string format_msg(const std::string& str);

	std::string get_trace_str_header(TraceDetailLevel level, const char* fun_name);
  void trace_debug_output(TraceDetailLevel level, const char* fun_name, const std::string& message);
  void trace_file(TraceDetailLevel level, const char* fun_name, const std::string& message);
  void trace_console(TraceDetailLevel level, const char* fun_name, const std::string& message);

	inline long long GetTimeStamp()
	{
#if (_MSC_VER >= 1400)  // VS2005 or later
		return __rdtsc();          // RDTSC is much more fast and accurate
	// QueryPerformanceCounter(&startTick_);
#else
		long long i64;

		__asm 
		{
			rdtsc;
			mov dword ptr i64, eax;
			mov dword ptr i64 + 4, edx;
		} 
		return i64;
#endif
	}

} // end of namespace dbg_utils


class func_guard 
{
public:
  func_guard(const char *funcName)
  {
		TRACE_VERBOSE((" %*s>>>>>> ENTER FUNCTION %s", 
			level_ * 2, "", 
			std::string(funcName).c_str()));
		//if (!_CrtCheckMemory())
		//{
		//	::OutputDebugString(std::string().format("########################  MEMORY FAILURE before %s\n", std::string(funcName).c_str()).ctstr());
		//}

    funcName_ = funcName;
    level_ += 1;
  }

  ~func_guard()
  {
    level_ -= 1;

    TRACE_VERBOSE(("%*s<<<<<< EXIT FUNCTION %s", 
			level_ * 2, "", 
			std::string(funcName_).c_str()));
		//if (!_CrtCheckMemory())
		//{
		//	TRACE_VERBOSE(("########################  MEMORY FAILURE after %s", std::string(funcName_).c_str()));
		//}
  }

	static void outputFunctionsStatistics();

private:
  std::string funcName_;
  static int level_;
};


class time_func_guard 
{
public:
  time_func_guard(LPCSTR funcName);
  ~time_func_guard();
	static void resetFunctionsStatistics();
	static void outputFunctionsStatistics();

private:
	struct TFunctionInfo
	{
		LARGE_INTEGER totalTicks, minTicks, maxTicks;
		int callCount;

		TFunctionInfo();
	};

	typedef std::map<LPCSTR, TFunctionInfo> TFunctionsInfo;

	static TFunctionsInfo functionsInfo_;

 // std::string funcName_;
	LPCSTR funcName_;
  static int level_;
	LARGE_INTEGER startTick_;
};

#if defined(FUNC_TIME_MEASUREMENT) 
	#define FUNC_GUARD time_func_guard fg(__FUNCTION__);
#elif defined(_DEBUG)
//	#define FUNC_GUARD func_guard fg(__FUNCTION__);
  #define FUNC_GUARD   
#else
	#define FUNC_GUARD 
#endif
