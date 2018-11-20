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

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Shell>&, const DomainList&);

    private:
        static DerivedRegister<BoxVolume> reg;
};

#endif	/* BOXVOLUME_H */

