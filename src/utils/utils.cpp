#include "utils.h"

std::string new_base_index(int vindex, int base)
{
    std::vector<char> number;
    int i = 0, r;
    
    do {
        r = vindex % base;
        number.push_back( names[r] );
        vindex /= base;
        i++;
    } while(vindex != 0);
    
    while (number.size() < 4)
    {
        number.push_back( names[0] );
    }
    
    std::reverse(number.begin(), number.end());
    std::string strnumber(number.begin(), number.end());
    
    return strnumber;
}
