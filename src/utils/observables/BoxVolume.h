#ifndef BOXVOLUME_H
#define	BOXVOLUME_H

#include "utils/observables/Observer.h"

class BoxVolume : public Observer
{
    public:
    public:
        explicit BoxVolume(const char*, const char*);
        BoxVolume(const BoxVolume& orig);
        virtual ~BoxVolume();

        double observe(const Box&, const std::vector<Shell>&);

    private:
        static DerivedRegister<BoxVolume> reg;
};

#endif	/* BOXVOLUME_H */

