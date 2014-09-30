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

#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <map>

namespace utils
{

    enum class LogLevel
    {
        CRITICAL = 0,
        SEVERE = 1,
        WARNING = 2,
        FILE = 3,
        INFO = 4,
        FINE = 5,
        FINER = 6,
        FINEST = 7,
    };

    static const std::pair<LogLevel, std::string> pairs[] =   //
    {
        std::pair<LogLevel, std::string>(LogLevel::CRITICAL, "CRITICAL"), //
        std::pair<LogLevel, std::string>(LogLevel::SEVERE, "SEVERE"), //
        std::pair<LogLevel, std::string>(LogLevel::WARNING, "WARNING"), //
        std::pair<LogLevel, std::string>(LogLevel::FILE, "FILE"), //
        std::pair<LogLevel, std::string>(LogLevel::INFO, "INFO"), //
        std::pair<LogLevel, std::string>(LogLevel::FINE, "FINE"), //
        std::pair<LogLevel, std::string>(LogLevel::FINER, "FINER"), //
        std::pair<LogLevel, std::string>(LogLevel::FINEST, "FINEST") //
    };

    static const std::map<LogLevel, std::string> log_level_names(pairs,
            pairs + sizeof(pairs) / sizeof(pairs[0]));

    class Logger
    {
        public:

            Logger(const std::string& module_name) :
                module_name_(module_name), sink(std::cerr), recent_level(LogLevel::INFO)
            {
            }

            inline const std::string& module_name() const
            {
                return module_name_;
            }

            virtual ~Logger() {}

            bool is_logable(LogLevel level) const;

            friend Logger& operator <<(Logger& logger, const LogLevel level);

            friend Logger& operator <<(Logger& logger, const char* message);

            friend Logger& operator <<(Logger& logger, const std::string& message);

            friend Logger& operator <<(Logger& logger, const double value);

            friend Logger& operator <<(Logger& logger, const int number);

            friend Logger& operator <<(Logger& logger, const size_t number);

            friend Logger& operator <<(Logger& logger, const unsigned long long number);

        private:
            const std::string  module_name_;
            std::ostream& sink;// = std::cerr;
            LogLevel recent_level;// = LogLevel::INFO;
    };

}
#endif /* UTILS_LOGGER_HH */