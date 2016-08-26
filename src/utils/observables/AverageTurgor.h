#ifndef AVERAGETURGOR_H
#define	AVERAGETURGOR_H

#include "utils/observables/Observer.h"

class AverageTurgor : public Observer
{
    public:
        explicit AverageTurgor(const char*, const char*);
        AverageTurgor(const AverageTurgor& orig);
        virtual ~AverageTurgor();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);

    private:
        static DerivedRegister<AverageTurgor> reg;
};

#endif	/* AVERAGETURGOR_H */