#ifndef MAXSTRAIN_H
#define	MAXSTRAIN_H

#include "utils/observables/Observer.h"

class MaxStrain : public Observer
{
    public:
        MaxStrain(const char*, const char*);
        MaxStrain(const MaxStrain& orig);
        virtual ~MaxStrain();
        
        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);
        
    private:
        static DerivedRegister<MaxStrain> reg;

};

#endif	/* MAXSTRAIN_H */

