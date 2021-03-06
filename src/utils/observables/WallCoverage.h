#ifndef WALLCOVERAGE_H
#define WALLCOVERAGE_H

#include "utils/observables/Observer.h"

class WallCoverage : public Observer
{
public:
    explicit WallCoverage(const char*, const char*);
    WallCoverage(const WallCoverage& orig);
    virtual ~WallCoverage();
    
    double observe(const Box&, const std::vector<Shell>&);
    
private:
    double contact_area(const Box&, const Shell&);
    bool is_in_contact(const Box&, const Shell&, const uint);
    bool is_touching_box(const Box&, const Vector3D&, const double, const double);
    static DerivedRegister<WallCoverage> reg;
};

#endif /* WALLCOVERAGE_H */