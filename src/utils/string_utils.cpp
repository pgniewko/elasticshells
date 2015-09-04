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

#include <cmath>

#include "utils/index.h"
#include "utils/string_utils.h"

namespace utils
{

#define BUFF_SIZE 4096
    char sprintf_buffer[BUFF_SIZE];

    std::string string_format(const std::string& fmt_str, ...)
    {
        va_list ap;
        va_start(ap, fmt_str);
        vsnprintf(sprintf_buffer, BUFF_SIZE, fmt_str.c_str(), ap);
        va_end(ap);
        return std::string(sprintf_buffer);
    }

    std::string& trim(std::string& str, const std::string delim)
    {
        std::string::size_type pos = str.find_last_not_of(delim);

        if (pos != std::string::npos)
        {
            str.erase(pos + 1);
            pos = str.find_first_not_of(delim);

            if (pos != std::string::npos)
            {
                str.erase(0, pos);
            }
        }
        else
        {
            str.erase(str.begin(), str.end());
        }

        return str;
    }

    std::vector<std::string> split(const std::string& s, const char delim)
    {
        std::vector<std::string> tokens;
        split(s, tokens, delim);
        return tokens;
    }

    template<>
    std::vector<std::string>& split<std::string>(const std::string& s, std::vector<std::string>& tokens,
            const char delim)
    {
        std::string s_copy(s);
        trim(s_copy);
        std::stringstream ss(s_copy);
        std::string item;

        while (std::getline(ss, item, delim))
        {
            if (trim(item).size() == 0)
            {
                continue;
            }

            tokens.push_back(item);
        }

        return tokens;
    }

    std::string& replace_substring(std::string& subject, const std::string& search, const std::string& replace)
    {
        if ((subject.length() == 0) || (search.length() == 0))
        {
            return subject;
        }

        size_t pos = 0;

        while ((pos = subject.find(search, pos)) != std::string::npos)
        {
            subject.replace(pos, search.length(), replace);
            pos += replace.length();
        }

        return subject;
    }

    void to_lower(std::string& s)
    {
        std::string::iterator i = s.begin();
        std::string::iterator end = s.end();

        while (i != end)
        {
            *i = std::tolower((unsigned char) * i);
            ++i;
        }
    }

    void to_upper(std::string& s)
    {
        std::string::iterator i = s.begin();
        std::string::iterator end = s.end();

        while (i != end)
        {
            *i = std::toupper((unsigned char) * i);
            ++i;
        }
    }

    double to_double(const char* p)
    {
        double r = 0.0;
        bool neg = false;

        while (*p == ' ')
        {
            ++p;
        }

        if (*p == '-')
        {
            neg = true;
            ++p;
        }

        while (*p >= '0' && *p <= '9')
        {
            r = (r * 10.0) + (*p - '0');
            ++p;
        }

        if (*p == '.')
        {
            double f = 0.0;
            int n = 0;
            ++p;

            while (*p >= '0' && *p <= '9')
            {
                f = (f * 10.0) + (*p - '0');
                ++p;
                ++n;
            }

            r += f / std::pow(10.0, n);
        }

        if (neg)
        {
            r = -r;
        }

        return r;
    }

    int to_int(const char* p)
    {
        int r = 0.0;
        bool neg = false;

        while (*p == ' ')
        {
            ++p;
        }

        if (*p == '-')
        {
            neg = true;
            ++p;
        }

        while (*p >= '0' && *p <= '9')
        {
            r = (r * 10.0) + (*p - '0');
            ++p;
        }

        if (neg)
        {
            r = -r;
        }

        return r;
    }

    std::string format_paragraph(std::vector<std::string>& words,
                                 const std::string& paragraph_pad, const std::string& line_pad, const core::index2 max_line_width)
    {
        std::string out = paragraph_pad;
        core::index2 lineLength = paragraph_pad.length();
        size_t pos = 0;

        for (std::string & w : words)
        {
            w = utils::trim(w);

            if ((pos = w.find("%N")) != w.npos)
            {
                w[pos] = ' ';
                w[pos + 1] = '\n';
                lineLength = 0;
            }

            if ((pos = w.find("%T")) != w.npos)
            {
                w[pos] = ' ';
                w[pos + 1] = '\t';
                lineLength += 4;
            }

            if (lineLength + w.length() + 1 > max_line_width)
            {
                // won't fit. Start a new line.
                if (lineLength != 0)
                {
                    out += '\n';
                    out += line_pad;
                    lineLength = line_pad.length();
                }

                // no lead space
            }
            else
            {
                /* will fit */
                if (lineLength != 0)
                {
                    // add lead space
                    out += ' ';
                    lineLength++;
                }
            }

            out += w;
            lineLength += w.length();
        } // end for

        return out;
    }

    bool has_suffix(const std::string& str, const std::string& suffix)
    {
        return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

    std::string format_paragraph(std::vector<std::string>& words,
                                 const std::string& paragraph_pad, const std::string& line_pad, const int max_line_width, const int max_first_line_width)
    {
        std::string out = paragraph_pad;
        int lineLength = paragraph_pad.length();
        size_t pos = 0;
        bool is_first_line = true;

        for (std::string & w : words)
        {
            w = utils::trim(w);

            if ((pos = w.find("%N")) != w.npos)
            {
                w[pos] = ' ';
                w[pos + 1] = '\n';
                lineLength = 0;
            }

            if ((pos = w.find("%T")) != w.npos)
            {
                w[pos] = ' ';
                w[pos + 1] = '\t';
                lineLength += 4;
            }

            const size_t mx = (is_first_line) ? max_first_line_width : max_line_width;

            if (lineLength + w.length() + 1 > mx)
            {
                // won't fit. Start a new line.
                if (lineLength != 0)
                {
                    out += '\n';
                    out += line_pad;
                    lineLength = line_pad.length();
                    is_first_line = false;
                }

                // no lead space
            }
            else
            {
                /* will fit */
                if (lineLength != 0)
                {
                    // add lead space
                    out += ' ';
                    lineLength++;
                }
            }

            out += w;
            lineLength += w.length();
        } // end for

        return out;
    }

} // ~utils
