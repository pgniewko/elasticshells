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
        
    double observe(Box&, std::vector<Cell>&);        
    void set_params(int, ...) {return;};
    void set_params(int, std::vector<std::string>) {return;};
private:
    static DerivedRegister<Aspherity> reg;
};

#endif	/* ASPHERITY_H */

