#ifndef ASPHERITY_H
#define	ASPHERITY_H

#include "geometry/Vector3D.h"
#include "utils/observables/Observer.h"

class Aspherity : public Observer
{
    public:
        explicit Aspherity(const char*, const char*);
        Aspherity(const Aspherity& orig);
        virtual ~Aspherity();

        double observe(const Box&, const std::vector<Shell>&);

    private:
        static DerivedRegister<Aspherity> reg;
};

#endif	/* ASPHERITY_H */

