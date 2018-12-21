#include "AverageContactsNumber.h"

AverageContactsNumber::AverageContactsNumber(const char* name, const char* format) : Observer(name, format) {}

AverageContactsNumber::AverageContactsNumber(const AverageContactsNumber& orig) : Observer(orig) {}

AverageContactsNumber::~AverageContactsNumber() {}

void AverageContactsNumber::set_params(const int num, std::vector<std::string> args_)
{
};

double AverageContactsNumber::observe(const Box& box, const std::vector<Shell>& shells)
{
    if (image_not_created)
    {
        create_shells_image(box, shells);
    }
    copy_shells_data(box, shells);

    double contact_number = 0.0;
    std::size_t n = shells.size();

    for (std::size_t i = 0; i < n; i++)
    {
        for (std::size_t j = i + 1; j < n; j++)
        {
            if ( is_in_contact(i, j) )
            {
                contact_number += 2.0;
            }
        }
    }

    return contact_number / (double) n;
}

DerivedRegister<AverageContactsNumber> AverageContactsNumber::reg("AverageContactsNumber");