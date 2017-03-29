#ifndef CVTURGOR_H
#define CVTURGOR_H

#include "utils/observables/Observer.h"

class CvTurgor : public Observer
{
    public:
        CvTurgor(const char* name, const char* format);
        CvTurgor(const CvTurgor& orig);
        virtual ~CvTurgor();
    
        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&, const DomainList&);
    
    private:
        static DerivedRegister<CvTurgor> reg;
};

#endif /* CVTURGOR_H */

