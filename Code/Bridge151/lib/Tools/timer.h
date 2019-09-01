/*!
        @file    timer.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:06:23 #$

        @version $LastChangedRevision: 1929 $
*/

#ifndef TIMER_INCLUDED
#define TIMER_INCLUDED

//#define USE_RUSAGE
#undef USE_RUSAGE

#ifdef USE_RUSAGE
#include <sys/resource.h>
#endif

#include <cstdio>
#include <time.h>
#if !BRIDGE_WIN
#include <sys/time.h>
#endif

#include "IO/bridgeIO.h"
using Bridge::vout;

class BAPI Timer
{
 public:
  static const std::string class_name;

 public:
  Timer(const std::string& id = "", const bool report = false)
    : is_started(false), m_start(0), m_elapsed(0), m_counter(0),
    m_id(id),
    m_report_on_exit(report)
  {}

  ~Timer();

 private:
  // non-copyable
  Timer(const Timer&);
  Timer& operator=(const Timer&);

 public:

  void start();
  void stop();
  void reset();

  static void timestamp();

  double elapsed_sec() const;
  double elapsed_msec() const;
  unsigned long get_counter() const;

  void report(const Bridge::VerboseLevel vl = Bridge::GENERAL);

 private:
  bool is_started;
  double m_start;

  double m_elapsed;        // sec.
  unsigned long m_counter;

  std::string m_id;
  bool m_report_on_exit;
};
#endif /* _TIMER_H_ */
