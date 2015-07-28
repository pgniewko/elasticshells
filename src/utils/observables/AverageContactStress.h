#ifndef AVERAGECONTACTSTRESS_H
#define	AVERAGECONTACTSTRESS_H

#include "utils/observables/Observer.h"

class AverageContactStress : public Observer
{
    public:
        AverageContactStress(const char*, const char*);
        AverageContactStress(const AverageContactStress& orig);
        virtual ~AverageContactStress();
        
        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);
        double observe(Box&, std::vector<Cell>&);
    private:
        static DerivedRegister<AverageContactStress> reg;

};

#endif	/* AVERAGECONTACTSTRESS_H */