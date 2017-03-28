#include "AverageContactsNumber.h"

AverageContactsNumber::AverageContactsNumber(const char* name, const char* format) : Observer(name, format) {}

AverageContactsNumber::AverageContactsNumber(const AverageContactsNumber& orig) : Observer(orig) {}

AverageContactsNumber::~AverageContactsNumber() {}

void AverageContactsNumber::set_params(const int num, std::vector<std::string> args_)
{
    i_param = atoi(args_[ num + 0 ].c_str());
};

double AverageContactsNumber::observe(const Box& box, std::vector<Cell>& cells, const DomainList& dl)
{
    double contact_number = 0;
    std::size_t n = cells.size();
    
    for (std::size_t i = 0; i < n; i++)
    {
        for (std::size_t j = i+1; j < n; j++)
        {
            if ( i_param == 0)
            {
                if ( cells[i].isInContact(cells[j], box) )
                {
                    contact_number += 2.0;
                }
            }
            else if (i_param == 1)
            {
                if ( dl.isInContact(i, j, cells, box)  )
                {
                    contact_number += 2.0;
                }
            }
        }
    }
    
    double av_contact_number = contact_number / (double) n;
    
    return av_contact_number;
    
}

DerivedRegister<AverageContactsNumber> AverageContactsNumber::reg("AverageContactsNumber");