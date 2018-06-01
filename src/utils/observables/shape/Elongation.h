#ifndef ELONGATION_H
#define ELONGATION_H

#include "geometry/Vector3D.h"
#include "utils/observables/Observer.h"
#include <limits>

class Elongation : public Observer
{
    public:
        explicit Elongation(const char*, const char*);
        Elongation(const Elongation& orig);
        virtual ~Elongation();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Shell>&, const DomainList&);

    private:
        static DerivedRegister<Elongation> reg;
};

#endif /* ELONGATION_H */

