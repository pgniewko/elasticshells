#ifndef AVERAGECONTACTAREA_H
#define	AVERAGECONTACTAREA_H

#include "utils/observables/Observer.h"

class AverageContactArea : public Observer
{
    public:
        explicit AverageContactArea(const char*, const char*);
        AverageContactArea(const AverageContactArea& orig);
        virtual ~AverageContactArea();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);

    private:
        static DerivedRegister<AverageContactArea> reg;

};

#endif	/* AVERAGECONTACTAREA_H */

