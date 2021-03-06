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

        double observe(const Box&, const std::vector<Shell>&);

    private:
        static DerivedRegister<Elongation> reg;
};

#endif /* ELONGATION_H */

