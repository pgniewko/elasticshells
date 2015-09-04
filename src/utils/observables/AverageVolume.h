#ifndef AVERAGEVOLUME_H
#define	AVERAGEVOLUME_H

#include "utils/observables/Observer.h"

class AverageVolume : public Observer
{
    public:
        AverageVolume (const char*, const char*);
        AverageVolume (const AverageVolume& orig);
        virtual ~AverageVolume();

        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);
        double observe(Box&, std::vector<Cell>&);

    private:
        static DerivedRegister<AverageVolume> reg;
};

#endif	/* AVERAGEVOLUME_H */

