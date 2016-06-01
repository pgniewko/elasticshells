#ifndef AVERAGEACTIVEAREA_H
#define	AVERAGEACTIVEAREA_H

#include "utils/observables/Observer.h"

class AverageActiveArea : public Observer
{
    public:
        AverageActiveArea(const char*, const char*);
        AverageActiveArea(const AverageActiveArea& orig);
        virtual ~AverageActiveArea();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);

    private:
        static DerivedRegister<AverageActiveArea> reg;

};

#endif	/* AVERAGEACTIVEAREA_H */

