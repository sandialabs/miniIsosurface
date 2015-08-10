/*
 * Log.h
 *
 *  Created on: Jul 14, 2015
 *      Author: sjmunn
 */

#ifndef LOG_H_
#define LOG_H_

#include"../includes.h"

inline std::string NowTime();

enum TLogLevel {
	logERROR,
	logWARNING,
	logYAML, 			// Essential performance info
	logINFO,			// Info (progress info for impatient users mostly)
	logDEBUG,			// If you want to add debug info of your own
	logDEBUG1,			// Built-in debugging messages
	logDEBUG_Step, 		// Step by Step debugging, only activate for very small files
};

class Log {
public:
	Log();
	virtual ~Log();
	std::ostringstream& Get(TLogLevel level = logINFO);
public:
	static TLogLevel& ReportingLevel();
	static std::string ToString(TLogLevel level);
	static TLogLevel FromString(const std::string& level);
protected:
	std::ostringstream os;
private:
	Log(const Log&);
	Log& operator =(const Log&);
};

inline Log::Log() {
}

inline std::ostringstream& Log::Get(TLogLevel level) {
	os << "- " << NowTime();
	os << " " << ToString(level) << ": ";
	os << std::string(level > logDEBUG ? level - logDEBUG : 0, '\t');
	return os;
}

inline Log::~Log() {
	os << std::endl;
	fprintf(stderr, "%s", os.str().c_str());
	fflush(stderr);
}

inline TLogLevel& Log::ReportingLevel() {
	static TLogLevel reportingLevel = logINFO;
	return reportingLevel;
}

inline std::string Log::ToString(TLogLevel level) {
	static const char* const buffer[] = { "ERROR", "WARNING", "YAML", "INFO",
			"DEBUG", "DEBUG1", "STEP_BY_STEP"};
	return buffer[level];
}

inline TLogLevel Log::FromString(const std::string& level) {
	if (level == "STEP_BY_STEP")
		return logDEBUG_Step;
	if (level == "DEBUG1")
		return logDEBUG1;
	if (level == "DEBUG")
		return logDEBUG;
	if (level == "INFO")
		return logINFO;
	if (level == "YAML")
		return logYAML;
	if (level == "WARNING")
		return logWARNING;
	if (level == "ERROR")
		return logERROR;
	Log().Get(logWARNING) << "Unknown logging level '" << level
			<< "'. Using INFO level as default.";
	return logINFO;
}

typedef Log LOG;

#define CLOG(level) \
    if (level > LOG::ReportingLevel()) ; \
    else Log().Get(level)

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)

#include <windows.h>

inline std::string NowTime()
{
	const int MAX_LEN = 200;
	char buffer[MAX_LEN];
	if (GetTimeFormatA(LOCALE_USER_DEFAULT, 0, 0,
					"HH':'mm':'ss", buffer, MAX_LEN) == 0)
	return "Error in NowTime()";

	char result[100] = {0};
	static DWORD first = GetTickCount();
	std::sprintf(result, "%s.%03ld", buffer, (long)(GetTickCount() - first) % 1000);
	return result;
}

#else

#include <sys/time.h>

inline std::string NowTime() {
	char buffer[11];
	time_t t;
	time(&t);
	tm r = { 0 };
	strftime(buffer, sizeof(buffer), "%X", localtime_r(&t, &r));
	struct timeval tv;
	gettimeofday(&tv, 0);
	char result[100] = { 0 };
	std::sprintf(result, "%s.%03ld", buffer, (long) tv.tv_usec / 1000);
	return result;
}

#endif //WIN32
#endif /* LOG_H_ */
