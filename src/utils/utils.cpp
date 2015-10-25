#include "utils.h"

std::string new_base_index(int vindex, int base)
{
    std::vector<char> number;
    int i = 0, r;

    do
    {
        r = vindex % base;
        number.push_back( names[r] );
        vindex /= base;
        i++;
    }
    while (vindex != 0);

    while (number.size() < 4)
    {
        number.push_back( names[0] );
    }

    std::reverse(number.begin(), number.end());
    std::string strnumber(number.begin(), number.end());

    return strnumber;
}

std::vector<std::string>& split(const std::string& s, char delim, std::vector<std::string>& elems)
{
    std::stringstream ss(s);
    std::string item;

    while (std::getline(ss, item, delim))
    {
        if (item.length() > 0)
        {
            elems.push_back(item);
        }
    }

    return elems;
}

std::vector<std::string> split(const std::string& s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
