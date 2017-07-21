#ifndef ROUNDNESS_H
#define ROUNDNESS_H

#include "geometry/Vector3D.h"
#include "geometry/numeric/Miniball.h"
#include "utils/observables/Observer.h"

#include <limits>

class Roundness : public Observer
{
    public:
        explicit Roundness(const char*, const char*);
        Roundness(const Roundness& orig);
        virtual ~Roundness();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&, const DomainList&);

    private:


        double miniball_r(Cell&);

        static DerivedRegister<Roundness> reg;
};

#endif /* ROUNDNESS_H */

