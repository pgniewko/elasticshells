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

#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <algorithm> // for std::unique
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdarg>
#include <vector>
#include <sstream>

#include <iostream>

#include "utils/index.h"

namespace utils
{

#define TEXT_BOLD    "\x1b[1m"
#define TEXT_BLACK   "\x1b[30m"
#define TEXT_RED     "\x1b[31m"
#define TEXT_GREEN   "\x1b[32m"
#define TEXT_YELLOW  "\x1b[33m"
#define TEXT_BLUE    "\x1b[34m"
#define TEXT_MAGENTA "\x1b[35m"
#define TEXT_CYAN    "\x1b[36m"
#define TEXT_WHITE   "\x1b[37m"
#define TEXT_RESET   "\x1b[0m"

    std::string string_format(const std::string& fmt_str, ...);

    /** @brief Removes spaces (or other characters, if desired) from both ends of a given string
     *
     * @param str - the string to be modified
     * @param delim - characters to be trimmed off (by default = " \n\t\r")
     * @return the string representing the input data
     */
    std::string& trim(std::string& str, const std::string delim = " \n\t\r");

    /** @brief Splits a string into tokens based on a given delimiter.
     *
     * @param str - the input string with tokens to be converted
     * @param delim - character to be trimmed off (a space by default)
     * @return the newly created vector that holds the tokens extracted from the given string
     */
    std::vector<std::string> split(const std::string& str, const char delim);

    /** @brief Splits a string into tokens based on a given delimiter and converts them into a generic type T
     *
     * @tparam T - the final type of the data that will be converted from string
     * @param str - the input string with tokens to be converted
     * @param tokens - the resulting (converted) tokens will be stored here
     * @param delim - character to be trimmed off (a space by default)
     * @return the reference to the vector of converted data
     */
    template<typename T>
    std::vector<T>& split(const std::string& s, std::vector<T>& tokens, const char delim = ' ')
    {
        std::string s_copy(s);
        trim(s_copy);
        std::stringstream ss(s_copy);
        std::string item;
        T data;

        while (std::getline(ss, item, delim))
        {
            if (item.size() == 0)
            {
                continue;
            }

            std::stringstream ss2(item);
            ss2 >> data;
            tokens.push_back(data);
        }

        return tokens;
    }

    /** @brief Splits a string into tokens based on a given delimiter.
     *
     * This speciated method does not perform any conversion
     * @param str - the input string with tokens to be split
     * @param tokens - the resulting tokens will be stored here
     * @param delim - character to be trimmed off (a space by default)
     * @return the reference to the vector of tokens
     */
    template<>
    std::vector<std::string>& split<std::string>(const std::string& s, std::vector<std::string>& tokens, const char delim);

    /** @brief Converts any data to a string
     *
     * @tparam T - the generic type of the input data to be converted
     * @param data - the element to be converted
     * @return the string representing the input data
     */
    template<class T>
    const std::string to_string(const T& data)
    {
        std::stringstream ss;
        ss << data;
        return ss.str();
    }

    /** @brief Converts all the data from the given vector to a string
     *
     * @tparam T - the generic type of the input data to be converted
     * @param data - the element to be converted
     * @return the string representing the input data
     */
    template<class T>
    const std::string to_string(const std::vector<T>& data, const std::string& separator)
    {
        std::stringstream ss;

        for (const std::string d : data)
        {
            ss << d << separator;
        }

        return ss.str();
    }

    /** @brief Generic method to convert a string into other types.
     *
     * @tparam T - the generic type of data to be created
     * @param token - the element to be converted
     * @return output value extracetd from the given string
     */
    template<class T>
    T from_string(const std::string& token)
    {
        std::stringstream ss(token);
        T data;
        ss >> data;
        return data;
    }

    /** @brief Predicate that returns true if the two given chars are identical and equal to the template parameter
     *
     * The predicate is used by combine_same_chars() method.
     *
     * @tparam Remove - the character to be recognized by this method
     * @param first - the first of the two compared characters
     * @param second - the second of the two compared characters
     */
    template<char Remove> bool BothAre(char first, char second)
    {
        return first == second && first == Remove;
    }

    /** @brief Substitute multiple subsequent occurrences of a given character with a single char
     *
     * In the following examples, each stretch of spaces is substituted with a single space
     * @code
     *   std::string str = "Crazy     dog";
     *   std::cout << combine_same_chars<' '>(str);
     * @endcode
     *
     * @tparam Remove - the character to be recognized by this method
     * @param str - the string to be modified
     * @returns the modified string
     */
    template<char Remove>
    std::string& combine_same_chars(std::string& str)
    {
        str.erase(std::unique(str.begin(), str.end(), BothAre<Remove>), str.end());
        return str;
    }

    /** \brief Converts long strings into nice paragraphs.
     *
     * The method removes all existing tabs an new line characters and create a single line
     * string that in turn is nicely formatted to fit a given width. If you need to manually
     * break the text or introduce a tab, use the <code>%N</code> (new line) or <code>%T</code>
     * (tab) special command.
     *
     * @param words - a given line of text
     * @param max_line_width - length of a paragraph
     * @param paragraph_pad - string that is used as a prefix for the very first line in each paragraph
     * @param line_pad - string that is used as a prefix for every line, except the very first line in a paragraph
     * @return a string with \n characters that breaks the string into a paragraph.
     */
    std::string format_paragraph(std::vector<std::string>& words, const std::string& paragraph_pad,
                                 const std::string& line_pad, const int max_line_width);

    /** \brief Converts long strings into nice paragraphs.
     *
     * This is a variant of format_paragraph(std::vector<std::string> & , const std::string &,
     * const std::string &, const int) method where the first line might have different length than
     * all the other lines.
     *
     * @param words - a given line of text
     * @param max_line_width - length of a paragraph
     * @param max_first_line_width - length of the very first line in a paragraph
     * @param paragraph_pad - string that is used as a prefix for the very first line in each paragraph
     * @param line_pad - string that is used as a prefix for every line, except the very first line in a paragraph
     * @return a string with \n characters that breaks the string into a paragraph.
     */
    std::string format_paragraph(std::vector<std::string>& words, const std::string& paragraph_pad,
                                 const std::string& line_pad, const int max_line_width, const int max_first_line_width);

    template<class T>
    std::vector<T> from_string(const std::string& token, const core::index2 n_tokens)
    {
        std::stringstream ss(token);
        std::vector<T> tokens;

        for (core::index2 i = 0; i < n_tokens; i++)
        {
            T data;
            ss >> data;
            tokens.push_back(data);
        }

        return tokens;
    }

    template<class T>
    std::vector<T>& from_string(const std::string& token, const core::index2 n_tokens, std::vector<T>& destination)
    {
        std::stringstream ss(token);

        for (core::index2 i = 0; i < n_tokens; i++)
        {
            T data;
            ss >> data;
            destination.push_back(data);
        }

        return destination;
    }

    template<class T>
    T from_string(std::string& token, const T default_value)
    {
        trim(token);

        if (token.length() == 0)
        {
            return default_value;
        }

        std::stringstream ss(token);
        T data;
        ss >> data;
        return data;
    }

    template<class T>
    T from_string(const std::string& line, const int start, const int stop, const T default_value)
    {
        std::string token = line.substr(start, stop - start + 1);

//  trim(token);
        if (token.length() == 0)
        {
            return default_value;
        }

        std::stringstream ss(token);
        T data;
        ss >> data;
        return data;
    }

    template<typename T>
    std::string join_as_string(const T& v, const std::string& delimiter)
    {
        std::ostringstream s;

        for (const auto & i : v)
        {
            if (&i != &v[0])
            {
                s << delimiter;
            }

            s << i;
        }

        return s.str();
    }

    /** @brief Replaces all occurrences of a substring  with another string.
     *
     * @param subject - the string to be modified
     * @param search - the substring to be found and replaced
     * @param replace - the replacement string
     * @returns the modified string
     */
    std::string& replace_substring(std::string& subject, const std::string& search, const std::string& replace);

    void to_lower(std::string& s);
    void to_upper(std::string& s);

    /** @brief Fast method to convert a char buffer to a double value.
     *
     * @param p - the char* buffer
     * @returns double value
     */
    double to_double(const char* p);

    /** @brief Fast method to convert a char buffer to an integer value.
     *
     * @param p - the char* buffer
     * @returns int value
     */
    int to_int(const char* p);

    /** @brief Tests is a string ends with a particular suffix
     *
     * @param str - query string
     * @param suffix - desired suffix
     * @returns true if the string ends with the given cutoff value
     */
    bool has_suffix(const std::string& str, const std::string& suffix);

} //~utils

#endif