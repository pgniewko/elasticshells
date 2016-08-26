#ifndef MINSTRAIN_H
#define	MINSTRAIN_H

#include "utils/observables/Observer.h"

class MinStrain : public Observer
{
    public:
        explicit MinStrain(const char*, const char*);
        MinStrain(const MinStrain& orig);
        virtual ~MinStrain();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);

    private:
        static DerivedRegister<MinStrain> reg;

};

#endif	/* MINSTRAIN_H */

