#ifndef AVERAGEVOLUME_H
#define	AVERAGEVOLUME_H

#include "utils/observables/Observer.h"

class AverageVolume : public Observer
{
    public:
        explicit AverageVolume (const char*, const char*);
        AverageVolume (const AverageVolume& orig);
        virtual ~AverageVolume();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);

    private:
        static DerivedRegister<AverageVolume> reg;
};

#endif	/* AVERAGEVOLUME_H */

