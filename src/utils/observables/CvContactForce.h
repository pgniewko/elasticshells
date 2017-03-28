#ifndef CVCONTACTFORCE_H
#define CVCONTACTFORCE_H

#include "utils/observables/Observer.h"

class CvContactForce : public Observer
{    
    public:
        CvContactForce(const char*, const char*);
        CvContactForce(const CvContactForce& orig);
        virtual ~CvContactForce();
    
        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&, const DomainList&);
        
    private:
        static DerivedRegister<CvContactForce> reg;
};

#endif /* CVCONTACTFORCE_H */

