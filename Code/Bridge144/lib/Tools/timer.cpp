#include "BridgeLib_Private.h"

/*!
        @file    $Id:: timer.cpp #$

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
*/
#include "timer.h"

const std::string Timer::class_name = "Timer";

//====================================================================
void Timer::timestamp()
{
  const size_t buf_size = 1024;
  static char  buf[buf_size];

  time_t    current_time;
  struct tm *timep;

  current_time = time(NULL);
  timep        = localtime(&current_time);

  strftime(buf, buf_size, "%Y/%m/%d %H:%M:%S %z", timep);

  vout.general("%s: timestamp: %s\n", class_name.c_str(), buf);
}


//====================================================================
Timer::~Timer()
{
  if (m_report_on_exit) report();
}

#if BRIDGE_WIN
#include <windows.h> 


/* FILETIME of Jan 1 1970 00:00:00. */
static const unsigned __int64 epoch = ((unsigned __int64)116444736000000000ULL);

/*
* timezone information is stored outside the kernel so tzp isn't used anymore.
*
* Note: this function is not for Win32 high precision timing purpose. See
* elapsed_time().
*
* Note :  this function breaks after 2038
*
*/
int
gettimeofday(struct timeval * tp, struct timezone * )
{
    FILETIME    file_time;
    SYSTEMTIME  system_time;
    ULARGE_INTEGER ularge;

    GetSystemTime(&system_time);
    SystemTimeToFileTime(&system_time, &file_time);
    ularge.LowPart = file_time.dwLowDateTime;
    ularge.HighPart = file_time.dwHighDateTime;

    tp->tv_sec = (long)((ularge.QuadPart - epoch) / 10000000L);
    tp->tv_usec = (long)(system_time.wMilliseconds * 1000);

    return 0;
}
#endif

//====================================================================
void Timer::start()
{
  struct timeval t_start;

#ifdef USE_RUSAGE
  struct rusage ru;
  int           result = getrusage(RUSAGE_SELF, &ru);
  t_start = ru.ru_utime;
#else
  int result = gettimeofday(&t_start, 0);
#endif

  if (result) {
    vout.general("%s: warning, aquiring system clock failed.\n", class_name.c_str());
    return;
  }

  m_start = (double)t_start.tv_sec + t_start.tv_usec * 1.0e-6;

  is_started = true;
  ++m_counter;
}


//====================================================================
void Timer::stop()
{
  if (!is_started) return;

  struct timeval t_end;

#ifdef USE_RUSAGE
  struct rusage ru;
  int           result = getrusage(RUSAGE_SELF, &ru);
  t_end = ru.ru_utime;
#else
  int result = gettimeofday(&t_end, 0);
#endif

  if (result) {
    vout.general("%s: warning, aquiring system clock failed.\n", class_name.c_str());
    return;
  }

  double m_end = (double)t_end.tv_sec + t_end.tv_usec * 1.0e-6;

  m_elapsed += (m_end - m_start);

  is_started = false;
}


//====================================================================
void Timer::reset()
{
  is_started = false;
  m_elapsed  = double(0);
  m_start    = double(0);
  m_counter  = 0;
}


//====================================================================
double Timer::elapsed_sec() const
{
  return m_elapsed;
}


//====================================================================
double Timer::elapsed_msec() const
{
  return m_elapsed * 1.0e+3;
}


//====================================================================
unsigned long Timer::get_counter() const
{
  return m_counter;
}


//====================================================================
void Timer::report(const Bridge::VerboseLevel vl)
{
  stop();

  unsigned long count   = get_counter();
  double        elapsed = elapsed_sec();
  double        average = count ? elapsed / count : 0.0;

  vout.general(vl, "Elapsed time: %s: total %12.2f sec, count %4d, average %12.2f sec\n", m_id.c_str(), elapsed, count, average);
}


//==========================================================
//==================================================END=====
