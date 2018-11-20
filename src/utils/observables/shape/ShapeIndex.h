#ifndef SHAPEINDEX_H
#define SHAPEINDEX_H

#include "utils/observables/Observer.h"

class ShapeIndex  : public Observer
{
    public:
    public:
        explicit ShapeIndex(const char*, const char*);
        ShapeIndex(const ShapeIndex& orig);
        virtual ~ShapeIndex();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Shell>&, const DomainList&);

    private:
        static DerivedRegister<ShapeIndex> reg;

};

#endif /* SHAPEINDEX_H */

