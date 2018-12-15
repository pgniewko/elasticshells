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
        double observe(const Box&, const std::vector<Shell>&);

    private:


        double miniball_r(const Shell&);

        static DerivedRegister<Roundness> reg;
};

#endif /* ROUNDNESS_H */

