#ifndef LOGMANAGER_H
/*
   Copyright (C) 2014, Dominik Gront
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   The source code is a part of BioShell library. Please cite BioShell as:
   1. D.Gront and A.Kolinski "Utility library for structural bioinformatics" Bioinformatics 2008 24(4):584-585
   2. D.Gront and A.Kolinski "BioShell - a package of tools for structural biology computations" Bioinformatics 2006 22(5):621-622 


   Any feedback is very welcome.
   http://www.bioshell.pl/   
   email: dgront @ gmail.com or dgront @ chem.uw.edu.pl (remove space)
*/

#define LOGMANAGER_H

#include <unordered_set>

#include "utils/Logger.h"

namespace utils {

class LogManager {
public:
  static LogManager & get() {
    static LogManager manager;
    return manager;
  }

  static void mute(const std::string & logger_name) {
    logger <<LogLevel::INFO << "muted channel "<<logger_name<<"\n";
    muted.insert(logger_name);
  }

  static bool is_muted(const std::string & logger_name) {
    return (muted.find(logger_name) != muted.end());
  }

  static void set_level(LogLevel & new_level) {
    manager.log_level = new_level;
  }

  static void set_level(const std::string & new_level_name) {

    std::map<LogLevel, std::string>::const_iterator it;
    for (it = log_level_names.begin(); it != log_level_names.end(); it++) {
      if (it->second.compare(new_level_name) == 0) {
        manager.log_level = it->first;
        logger << it->first << "Logging level set to " << new_level_name << "\n";
        return;
      }
    }
    logger << LogLevel::WARNING << "Unknown log-level: " << new_level_name << ", the value not set. Known levels are:";
    for (it = log_level_names.begin(); it != log_level_names.end(); it++)
      logger << " " << it->second;
    logger << "\n";
  }

  static void CRITICAL() {
    manager.log_level = LogLevel::CRITICAL;
  }
  static void SEVERE() {
    manager.log_level = LogLevel::SEVERE;
  }
  static void WARNING() {
    manager.log_level = LogLevel::WARNING;
  }
  static void FILE() {
    manager.log_level = LogLevel::FILE;
  }
  static void INFO() {
    manager.log_level = LogLevel::INFO;
  }
  static void FINE() {
    manager.log_level = LogLevel::FINE;
  }
  static void FINER() {
    manager.log_level = LogLevel::FINER;
  }
  static void FINEST() {
    manager.log_level = LogLevel::FINEST;
  }

private:
  static const LogManager manager;
  std::ostream &sink = std::cerr;
  static LogLevel log_level;
  static utils::Logger logger;
  static std::unordered_set<std::string> muted;

  LogManager() {
  }

  LogManager(LogManager const&);
  void operator=(LogManager const&);

  friend Logger;

  friend Logger &operator <<(Logger &logger, const LogLevel level);

  friend Logger &operator <<(Logger &logger, const char* message);

  friend Logger &operator <<(Logger &logger, const std::string & message);

  friend Logger &operator <<(Logger &logger, const double value);

  friend Logger &operator <<(Logger &logger, const int number);

  friend Logger &operator <<(Logger &logger, const size_t number);

  friend Logger &operator <<(Logger &logger, const unsigned long long number);

};
}

#endif // ~ UTILS_LogManager_HH