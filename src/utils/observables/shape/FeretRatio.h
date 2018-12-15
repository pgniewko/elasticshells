#ifndef FERETRATIO_H
#define FERETRATIO_H

#include "geometry/Vector3D.h"
#include "utils/observables/Observer.h"
#include <limits>

class FeretRatio : public Observer
{
    public:
        explicit FeretRatio(const char*, const char*);
        FeretRatio(const FeretRatio& orig);
        virtual ~FeretRatio();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, const std::vector<Shell>&);

    private:
        static DerivedRegister<FeretRatio> reg;
};

#endif /* FERETRATIO_H */

