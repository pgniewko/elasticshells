#ifndef AVERAGETURGOR_H
#define	AVERAGETURGOR_H

#include "utils/observables/Observer.h"

class AverageTurgor : public Observer
{
    public:
        explicit AverageTurgor(const char*, const char*);
        AverageTurgor(const AverageTurgor& orig);
        virtual ~AverageTurgor();

        double observe(const Box&, const std::vector<Shell>&);

    private:
        static DerivedRegister<AverageTurgor> reg;
};

#endif	/* AVERAGETURGOR_H */