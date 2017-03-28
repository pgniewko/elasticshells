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

        double observe(const Box&, std::vector<Cell>&, const DomainList&);
        void set_params(const int, std::vector<std::string>);

    private:
        static DerivedRegister<BoxVolume> reg;
};

#endif	/* BOXVOLUME_H */

