#include "stdafx.h"
#include "debug_utils.h"

#include <sstream>
#include <sys/timeb.h>
#include <time.h>
#include <windows.h>

#pragma warning (disable:4996)
                
const size_t TRACE_MAXLOGMSG = 16*1024;

namespace
{
  
std::string TraceDetailLevelToString(dbg_utils::TraceDetailLevel level)
{
  switch (level)
  {
    case dbg_utils::tdl_Error:	 return "E";
    case dbg_utils::tdl_Normal:	 return "N";
    case dbg_utils::tdl_Verbose: return "V";
    case dbg_utils::tdl_Debug:	 return "D";
  }
  return "";
}

void GetTimeStr (char *timest, short bufsize)
{
  #ifndef __BORLANDC__
    struct __timeb64 timeb;
    struct tm *today;
    size_t len;
    char st[10];

    _ftime64 (&timeb);
    today = _localtime64 ( &timeb.time );
    len=strftime (timest, bufsize>50? 50: bufsize, "%H:%M:%S", today);
    timest[len]=0;
    sprintf (st,".%03d ",timeb.millitm);
    strcat (timest,st);
  #else  // #ifndef __BORLANDC__
   // TDateTime::FormatString("dd.mm.yy hh:mm:aa");
    LongTimeFormat = "HH:MM:SS.zzz ";
    strncpy (timest, TDateTime::CurrentDateTime().TimeString().c_str(), bufsize);
  #endif  // #ifndef __BORLANDC__
}

void GetLongTimeStr (char *timest, short bufsize)
{
  #ifndef __BORLANDC__
    struct __timeb64 timeb;
    struct tm *today;
    size_t len;
    char st[10];

    _ftime64 (&timeb);
   // time( &ltime );
    today = _localtime64 ( &timeb.time );
    len=strftime (timest, bufsize>50? 50: bufsize, "%d.%m.%Y %H:%M:%S", today);
    timest[len]=0;
    sprintf (st,".%03d ",timeb.millitm);
    strcat (timest,st);
  #else  // #ifndef __BORLANDC__
   // TDateTime::FormatString("dd.mm.yy hh:mm:aa");
    strncpy (timest, TDateTime::CurrentDateTime().DateTimeString().c_str(), bufsize);
  #endif  // #ifndef __BORLANDC__
}

int GetCurrThreadId()
{
  int threadId = (int)::GetCurrentThreadId();

  if (threadId < 100000)
    return threadId;
  else  return threadId%(1<<12);
}

void GetTimeAThreadId (char *timest, short bufsize)
{
  GetTimeStr (timest, bufsize);
  _snprintf(timest + strlen(timest), bufsize - strlen(timest) - 1, "[%4d] ", GetCurrThreadId ());
}

} //end of anonymous namespace

namespace dbg_utils {

std::string format_msg(const char* format, ...)
{
  va_list argp;
  va_start(argp, format);
  char msg_buffer[TRACE_MAXLOGMSG];
  //int n = _vsnprintf(msg_buffer, sizeof(msg_buffer), format, argp);
#if (_MSC_VER >= 1400)  // VS2005 or later
  int n = _vsnprintf_s(msg_buffer, TRACE_MAXLOGMSG, TRACE_MAXLOGMSG, format, argp);
#else
  int n = _vsnprintf(msg_buffer,TRACE_MAXLOGMSG - 1, format, argp);
#endif
  va_end(argp);
  if (n < 0 || n >= TRACE_MAXLOGMSG)
    return "Trace message is too long (skipped)!";

  return msg_buffer;
}

//std::string format_msg(const wchar_t* format, ...)
//{
//  va_list argp;
//  va_start(argp, format);
//  wchar_t msg_buffer[TRACE_MAXLOGMSG];
//  //int n = _vsnwprintf(msg_buffer, sizeof(msg_buffer), format, argp);
//#if (_MSC_VER >= 1400)  // VS2005 or later
//  int n = _vsnwprintf_s(msg_buffer, TRACE_MAXLOGMSG, _TRUNCATE, format, argp);
//#else
//  int n = _vsnwprintf(msg_buffer,TRACE_MAXLOGMSG - 1, format, argp);
//#endif
//  va_end(argp);
//  if (n < 0 || n >= TRACE_MAXLOGMSG)
//    return "Trace message is too long (skipped)!";
//
//  return msg_buffer;
//}

std::string format_msg(const std::string& str)
{
  return str;
}


std::string get_trace_str_header(TraceDetailLevel level, const char* fun_name)
{
	char st[100];
	GetTimeAThreadId (st, sizeof (st));

	return std::string(st) + TraceDetailLevelToString(level) + "|" +
	  std::string(fun_name) + "|";
}

void trace_debug_output(TraceDetailLevel level, const char* fun_name, const std::string& message)
{
  if (level < tdl_Error || level > tdl_Debug)
    level = tdl_Error;

  std::string str = get_trace_str_header(level, fun_name) + message + "\n";

  ::OutputDebugStringA(str.c_str());
}

void trace_file(TraceDetailLevel level, const char* fun_name, const std::string& message)
{
//	return; //d_
  if (level < tdl_Error || level > tdl_Debug)
    level = tdl_Error;

  std::string str = get_trace_str_header(level, fun_name) + message + "\n";

  static FILE *f = NULL;

	if (!f)
	{
		if ((GetKeyState(VK_SCROLL) & 1) != 0)
			f = fopen("C:\\15_2.log", "a");           //d_ 
		else
			f = fopen("C:\\15.log", "a");

		char st[100];

		GetLongTimeStr (st, sizeof (st));
    fputs("-----------------------------------------------------------\n", f);
		fprintf(f, "Logging started at %s\n", st);
	}
  fputs(str.c_str(), f);
	fflush(f);
}

void trace_console(TraceDetailLevel level, const char* fun_name, const std::string& message)
{
  if (level < tdl_Error || level > tdl_Debug)
    level = tdl_Error;

	char st[100];
	GetTimeStr(st, sizeof (st));

	printf("%s%s\n", st, message.c_str());
	fflush(stdout);
}

} // end of namespace dbg_utils



//#define rdtsc(low,high) \
//     __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high))

int func_guard::level_ = 0;
int time_func_guard::level_ = 0;
time_func_guard::TFunctionsInfo time_func_guard::functionsInfo_;

time_func_guard::TFunctionInfo::TFunctionInfo()
{
	totalTicks.QuadPart = 0;
	minTicks.QuadPart = -1;
	maxTicks.QuadPart = -1;
	callCount = 0;
}

time_func_guard::time_func_guard(LPCSTR funcName)
{
  funcName_ = funcName;

	if (level_ <= 1)
	{
		startTick_.QuadPart = dbg_utils::GetTimeStamp();
	}
  level_ += 1;
}

extern std::string g_curDna;

time_func_guard::~time_func_guard()
{
	level_ -= 1;

	LARGE_INTEGER currTick;

	if (level_ <= 1)
	{
		currTick.QuadPart = dbg_utils::GetTimeStamp() - startTick_.QuadPart;

		TFunctionInfo& functionInfo = functionsInfo_[funcName_];

		functionInfo.callCount++;
		functionInfo.totalTicks.QuadPart += currTick.QuadPart;
		if (functionInfo.minTicks.QuadPart < 0 || functionInfo.minTicks.QuadPart > currTick.QuadPart)
			functionInfo.minTicks.QuadPart = currTick.QuadPart;									 
		if (functionInfo.maxTicks.QuadPart < 0 || functionInfo.maxTicks.QuadPart < currTick.QuadPart)
		{
			functionInfo.maxTicks.QuadPart = currTick.QuadPart;
			//if (!g_curDna.empty())
			//	printf("Maximum AddDna time for %s (%I64d cycles)\n", 
			//				 g_curDna.c_str(), currTick.QuadPart);
		}

		if (functionsInfo_.size() > 100000)
			THROW_EXCEPTION("Functions time statistics is too large");
		//if (functionInfo.callCount == 2)   //d_
		//	TRACE_VERBOSE(("FUNCTION %50s call N2", funcName_));
	}
	//else      //d_
	//{
	//	QueryPerformanceCounter(&currTick);
	//	currTick.QuadPart -= startTick_.QuadPart;
	//}

	////{ //d_
	////	LARGE_INTEGER freq; 
	////	static int callNum = 0;

	////	callNum++;
	////	QueryPerformanceFrequency(&freq);

	////	if ((funcName_.find("::QueryJunction") >= 0 ||
	////		   funcName_.find("::QueryAdjacencies") >= 0) &&
	////	  	double(currTick.QuadPart) / freq.QuadPart > 0.002)    
	////	{
	////		TRACE_VERBOSE(("FUNCTION %50s: call %8d, %12.5f s", 
	////			funcName_.c_str(),
	////			callNum, 
	////			double(currTick.QuadPart) / freq.QuadPart));
	////	}

	////	if (funcName_.find(" (eid") >= 0 || 
	////		  funcName_.find(", eid") >= 0)
	////		functionsInfo_.erase(std::string(funcName_));
	////}


	//	if (level_ <= 1)       //d_
	//{
	//	QueryPerformanceCounter(&currTick);
	//	currTick.QuadPart -= startTick_.QuadPart;

	//	TFunctionInfo& functionInfo = functionsInfo_[std::string(funcName_)];

	//	if (functionsInfo_.size() > 100000)
	//		THROW_EXCEPTION("Functions time statistics is too large");
	//	if (functionInfo.callCount == 2)   //d_
	//		TRACE_VERBOSE(("FUNCTION %50s call N2", funcName_.c_str()));
	//}

}

void time_func_guard::resetFunctionsStatistics()
{
	functionsInfo_.clear();
}

void time_func_guard::outputFunctionsStatistics()
{
	LARGE_INTEGER freq;

	QueryPerformanceFrequency(&freq);

	const double mult = 0.001;
	//const double mult = double(1000000) / freq.QuadPart;

	for (TFunctionsInfo::const_iterator it = functionsInfo_.begin(); it != functionsInfo_.end(); it++)
	{
		const TFunctionInfo& functionInfo = it->second;

		TRACE_VERBOSE(("FUNCTION %40s stat:\n  %6d calls, %12.7f Gcycles, min %10.7g Kc., max %10.7g Kc.", 
									 it->first,
									 functionInfo.callCount, 
									 double(functionInfo.totalTicks.QuadPart) / 1e9,
									 double(functionInfo.minTicks.QuadPart) * mult,
									 double(functionInfo.maxTicks.QuadPart) * mult));
		//TRACE_VERBOSE(("FUNCTION %40s stat:\n  %8d calls, %12.7f s, min %10.7g mcs, max %10.7g mcs", 
		//							 it->first,
		//							 functionInfo.callCount, 
		//							 double(functionInfo.totalTicks.QuadPart) / freq.QuadPart,
		//							 double(functionInfo.minTicks.QuadPart) * mult,
		//							 double(functionInfo.maxTicks.QuadPart) * mult));
	}
}


