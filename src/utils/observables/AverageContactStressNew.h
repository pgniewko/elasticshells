#ifndef AVERAGECONTACTSTRESSNEW_H
#define	AVERAGECONTACTSTRESSNEW_H

#include "utils/observables/Observer.h"

class AverageContactStressNew : public Observer
{
    public:
        AverageContactStressNew(const char*, const char*);
        AverageContactStressNew(const AverageContactStressNew& orig);
        virtual ~AverageContactStressNew();
        
        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);
        double observe(Box&, std::vector<Cell>&);
        
    private:
        static DerivedRegister<AverageContactStressNew> reg;

};

#endif	/* AVERAGECONTACTSTRESSNEW_H */

