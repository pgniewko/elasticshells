#ifndef ASPHERITY_H
#define	ASPHERITY_H

#include "geometry/Vector3D.h"
#include "utils/observables/Observer.h"

class Aspherity : public Observer
{
    public:
        Aspherity(const char*, const char*);
        Aspherity(const Aspherity& orig);
        virtual ~Aspherity();

        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);
        double observe(Box&, std::vector<Cell>&);

    private:
        static DerivedRegister<Aspherity> reg;
};

#endif	/* ASPHERITY_H */

