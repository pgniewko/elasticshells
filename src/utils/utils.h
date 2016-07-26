#ifndef COMMONS_H
#define	COMMONS_H

#include <vector>       // std::vector
#include <string>       // std::string
#include <algorithm>    // std::reverse
#include <string>       // std::string
#include <sstream>

#include "Environment.h"

extern std::string new_base_index(int, int = 26);

extern std::vector<std::string>& split(const std::string&, char , std::vector<std::string>&);
extern std::vector<std::string> split(const std::string&, char);
extern std::string trim( const std::string&);

#endif	/* COMMONS_H */

