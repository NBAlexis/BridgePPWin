#include "BridgeLib_Private.h"

/*!
        @file    $Id:: bridgeIO.cpp #$

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: aoym $

        @date    $LastChangedDate:: 2017-02-24 00:49:46 #$

        @version $LastChangedRevision: 1561 $
 */

#include "bridgeIO.h"
#include "Parameters/commonParameters.h"
#include "ResourceManager/threadManager_OpenMP.h"

//====================================================================
// verbose output for c style
// default verbose level, node 0

namespace Bridge {
  //====================================================================
  const std::string BridgeIO::class_name = "BridgeIO";

  //====================================================================
  BridgeIO::BridgeIO(const std::string& filename)
  {
    os_ = NULL;

#ifdef ENABLE_ILDG_TAG
    ildg_os_ = NULL;
#endif

    init(filename);

#ifdef ENABLE_ILDG_TAG
    ildg_init(filename);
#endif
  }


  //====================================================================
  BridgeIO::~BridgeIO()
  {
    tidyup_();
  }


  //====================================================================
  void BridgeIO::init(const std::string& filename)
  {
    if (os_) {
      stack_.push(os_);
    }

    if (filename == "stdout") {
      os_ = new std::ostream(std::cout.rdbuf());
    } else {
      os_ = new std::ofstream(filename.c_str());
    }

    if (!os_) {
      fprintf(stderr, "%s: init: unable to open log file \"%s\".\n", class_name.c_str(), filename.c_str());
      exit(EXIT_FAILURE);

      rewind_();
    }
  }


  //====================================================================
  void BridgeIO::init(const std::ostream& ost)
  {
    if (os_) {
      stack_.push(os_);
    }

    os_ = new std::ostream(ost.rdbuf());

    if (!os_) {
      fprintf(stderr, "%s: init: unable to open stream.\n", class_name.c_str());
      exit(EXIT_FAILURE);

      rewind_();
    }
  }


  //====================================================================
  void BridgeIO::unset()
  {
    if (os_) delete os_;

    rewind_();
  }


  //====================================================================
  void BridgeIO::rewind_()
  {
    if (stack_.size() > 0) {
      os_ = stack_.top();
      stack_.pop();
    } else {
      os_ = NULL;
    }
  }


  //====================================================================
  void BridgeIO::tidyup_()
  {
    if (os_) delete os_;

    while (stack_.size() > 0)
    {
      std::ostream *otmp = stack_.top();
      if (otmp) delete otmp;
      stack_.pop();
    }

#ifdef ENABLE_ILDG_TAG
    if (ildg_os_) delete ildg_os_;
#endif
  }


  //====================================================================
  VerboseLevel
  BridgeIO::set_verbose_level(const std::string& str)
  {
    ThreadManager_OpenMP::assert_single_thread(class_name);

    if ((str == "Crucial") || (str == "crucial") || (str == "CRUCIAL")) return Bridge::CRUCIAL;

    if ((str == "General") || (str == "general") || (str == "GENERAL")) return Bridge::GENERAL;

    if ((str == "Detailed") || (str == "detailed") || (str == "DETAILED")) return Bridge::DETAILED;

    if ((str == "Paranoiac") || (str == "paranoiac") || (str == "PARANOIAC")) return Bridge::PARANOIAC;

    if ((str == "NULL") || (str == "null")) return CommonParameters::Vlevel();

    // safe default
    return Bridge::GENERAL;
  }


  //====================================================================
  std::string
  BridgeIO::get_verbose_level(const VerboseLevel vl)
  {
    ThreadManager_OpenMP::assert_single_thread(class_name);

    switch (vl)
    {
    case Bridge::CRUCIAL:
      return "Crucial";

    case Bridge::GENERAL:
      return "General";

    case Bridge::DETAILED:
      return "Detailed";

    case Bridge::PARANOIAC:
      return "Paranoiac";

    default:
      return "NULL";
    }
  }


  //====================================================================
  void
  BridgeIO::crucial(const char *format, ...)
  {
    va_list arg;

#pragma omp barrier
#pragma omp master
    {
      va_start(arg, format);
      print(CommonParameters::Vlevel(), Bridge::CRUCIAL, 0, format, arg);
      va_end(arg);
    }
#pragma omp barrier
  }


  //====================================================================
  void
  BridgeIO::general(const char *format, ...)
  {
    va_list arg;

#pragma omp barrier
#pragma omp master
    {
      va_start(arg, format);
      print(CommonParameters::Vlevel(), Bridge::GENERAL, 0, format, arg);
      va_end(arg);
    }
#pragma omp barrier
  }


  //====================================================================
  void
  BridgeIO::detailed(const char *format, ...)
  {
    va_list arg;

#pragma omp barrier
#pragma omp master
    {
      va_start(arg, format);
      print(CommonParameters::Vlevel(), Bridge::DETAILED, 0, format, arg);
      va_end(arg);
    }
#pragma omp barrier
  }


  //====================================================================
  void
  BridgeIO::paranoiac(const char *format, ...)
  {
    va_list arg;

#pragma omp barrier
#pragma omp master
    {
      va_start(arg, format);
      print(CommonParameters::Vlevel(), Bridge::PARANOIAC, 0, format, arg);
      va_end(arg);
    }
#pragma omp barrier
  }


  //====================================================================
  // input verbose level, node 0
  void
  BridgeIO::crucial(VerboseLevel vl, const char *format, ...)
  {
    va_list arg;

#pragma omp barrier
#pragma omp master
    {
      va_start(arg, format);
      print(vl, Bridge::CRUCIAL, 0, format, arg);
      va_end(arg);
    }
#pragma omp barrier
  }


  //====================================================================
  void
  BridgeIO::general(VerboseLevel vl, const char *format, ...)
  {
    va_list arg;

#pragma omp barrier
#pragma omp master
    {
      va_start(arg, format);
      print(vl, Bridge::GENERAL, 0, format, arg);
      va_end(arg);
    }
#pragma omp barrier
  }


  //====================================================================
  void
  BridgeIO::detailed(VerboseLevel vl, const char *format, ...)
  {
#pragma omp barrier
#pragma omp master
    {
      va_list arg;

      va_start(arg, format);
      print(vl, Bridge::DETAILED, 0, format, arg);
      va_end(arg);
    }
#pragma omp barrier
  }


  //====================================================================
  void
  BridgeIO::paranoiac(VerboseLevel vl, const char *format, ...)
  {
    va_list arg;

#pragma omp barrier
#pragma omp master
    {
      va_start(arg, format);
      print(vl, Bridge::PARANOIAC, 0, format, arg);
      va_end(arg);
    }
#pragma omp barrier
  }


  //====================================================================
  // input verbose level, input node
  void
  BridgeIO::crucial(VerboseLevel vl, int node, const char *format, ...)
  {
    va_list arg;

#pragma omp barrier
#pragma omp master
    {
      va_start(arg, format);
      print(vl, Bridge::CRUCIAL, node, format, arg);
      va_end(arg);
    }
#pragma omp barrier
  }


  //====================================================================
  void
  BridgeIO::general(VerboseLevel vl, int node, const char *format, ...)
  {
    va_list arg;

#pragma omp barrier
#pragma omp master
    {
      va_start(arg, format);
      print(vl, Bridge::GENERAL, node, format, arg);
      va_end(arg);
    }
#pragma omp barrier
  }


  //====================================================================
  void
  BridgeIO::detailed(VerboseLevel vl, int node, const char *format, ...)
  {
    va_list arg;

#pragma omp barrier
#pragma omp master
    {
      va_start(arg, format);
      print(vl, Bridge::DETAILED, node, format, arg);
      va_end(arg);
    }
#pragma omp barrier
  }


  //====================================================================
  void
  BridgeIO::paranoiac(VerboseLevel vl, int node, const char *format, ...)
  {
    va_list arg;

#pragma omp barrier
#pragma omp master
    {
      va_start(arg, format);
      print(vl, Bridge::PARANOIAC, node, format, arg);
      va_end(arg);
    }
#pragma omp barrier
  }


  //====================================================================
  std::ostream& BridgeIO::getStream()
  {
    return *os_;
  }


  //====================================================================
  bool BridgeIO::isOpen()
  {
    return os_ && *os_;
  }


  //====================================================================
  inline
  void
  BridgeIO::print(VerboseLevel level, VerboseLevel write_level,
                  int node, const char *format, va_list& arg)
  {
    if ((write_level <= level) && (Communicator::nodeid() == node)) {
      if (!os_) {
        std::cerr << "ERROR: BridgeIO: no output stream." << std::endl;
        exit(EXIT_FAILURE);
      }

      vsprintf(buff_, format, arg);

      *os_ << buff_;
#ifdef DEBUG
      *os_ << std::flush;
#endif

      if (!os_->good()) {
        std::cerr << "ERROR: BridgeIO: output failed." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }


#ifdef ENABLE_ILDG_TAG
  //====================================================================
  void
  BridgeIO::ildg_init(const std::ostream& ost)
  {
    if (ildg_os_) delete ildg_os_;
    ildg_os_ = new std::ostream(ost.rdbuf());

    if (!ildg_os_) {
      fprintf(stderr, "%s: init: unable to ildg open stream.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }


  //====================================================================
  void
  BridgeIO::ildg_init(const std::string& filename)
  {
    if (ildg_os_) delete ildg_os_;

    if (filename == "stdout") {
      ildg_os_ = new std::ostream(std::cout.rdbuf());
    } else {
      ildg_os_ = new std::ofstream(filename.c_str());
    }

    if (!ildg_os_) {
      fprintf(stderr, "%s: init: unable to open ildg log file \"%s\".\n", class_name.c_str(), filename.c_str());
      exit(EXIT_FAILURE);
    }
  }


  //====================================================================
  void
  BridgeIO::ildg(const char *format, ...)
  {
    if (!ildg_os_) {
      std::cerr << "ERROR: BridgeIO: no ildg output stream." << std::endl;
      exit(EXIT_FAILURE);
    }

    va_list arg;

    va_start(arg, format);
    vsprintf(buff_, format, arg);
    va_end(arg);

    *ildg_os_ << "@ILDG:" << buff_;
#ifdef DEBUG
    *ildg_os_ << std::flush;
#endif

    if (!ildg_os_->good()) {
      std::cerr << "ERROR: BridgeIO: ildg output failed." << std::endl;
      exit(EXIT_FAILURE);
    }
  }


  //====================================================================
  std::ostream&
  BridgeIO::getILDGStream()
  {
    return *ildg_os_;
  }
#endif
  //====================================================================


  //====================================================================
  BridgeIO vout;
}

//====================================================================
//====================================================================
