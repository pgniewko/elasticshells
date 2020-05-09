#include "geometry/Vector3D.h"
#include "utils/observables/Observer.h"
#include "geometry/numeric/eig3.h"

#ifndef TENSORIAL_H
#define TENSORIAL_H

class Tensorial : public Observer
{
    public:
        explicit Tensorial(const char*, const char*);
        Tensorial(const Tensorial& orig);
        virtual ~Tensorial();

        double observe(const Box&, const std::vector<Shell>&);

    private:
        static DerivedRegister<Tensorial> reg;
};

#endif /* TENSORIAL_H */

