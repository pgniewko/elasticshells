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

#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>

#include "utils/Logger.h"
#include "utils/LogManager.h"
#include "utils/string_utils.h"

namespace utils {

bool Logger::is_logable(LogLevel level) const { return (level <= LogManager::log_level); }

Logger &operator <<(Logger &logger, const char* message) {

  if (logger.recent_level <= LogManager::log_level)
    if (!LogManager::is_muted(logger.module_name_)) logger.sink << message;

  return logger;
}

Logger &operator <<(Logger &logger, const size_t number) {

  if (logger.recent_level <= LogManager::log_level)
    if (!LogManager::is_muted(logger.module_name_)) logger.sink << number;

  return logger;
}


Logger &operator <<(Logger &logger, const unsigned long long number) {

  if (logger.recent_level <= LogManager::log_level)
    if (!LogManager::is_muted(logger.module_name_)) logger.sink << number;

  return logger;
}

Logger &operator <<(Logger &logger, const int number) {

  if (logger.recent_level <= LogManager::log_level)
    if (!LogManager::is_muted(logger.module_name_)) logger.sink << number;

  return logger;
}

Logger &operator <<(Logger &logger, const double value) {

  if (logger.recent_level <= LogManager::log_level) if (!LogManager::is_muted(logger.module_name_)) logger.sink << value;

  return logger;
}

Logger &operator <<(Logger &logger, const std::string & message) {

  if (logger.recent_level <= LogManager::log_level)
    if (!LogManager::is_muted(logger.module_name_)) logger.sink << message;

  return logger;
}

Logger &operator <<(Logger &logger, const LogLevel level) {

  if (LogManager::is_muted(logger.module_name_)) return logger;

  logger.recent_level = level;
  std::string color = "";
  if (logger.recent_level <= LogManager::log_level) {
    switch (level) {
      case LogLevel::CRITICAL:
      case LogLevel::SEVERE:
        color = TEXT_RED;
        logger<<"\n";
        break;
      case LogLevel::WARNING:
      case LogLevel::INFO:
        color = TEXT_YELLOW;
        break;
      case LogLevel::FILE:
      case LogLevel::FINE:
      case LogLevel::FINER:
      case LogLevel::FINEST:
        color = TEXT_WHITE;
        break;
    }
    logger.sink << color << "[" << log_level_names.at(level) << "] " << TEXT_RESET << TEXT_BOLD << logger.module_name_
        << TEXT_RESET << " ";
  }
  return logger;
}

}