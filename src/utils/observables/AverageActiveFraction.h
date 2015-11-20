#ifndef AVERAGEACTIVEFRACTION_H
#define	AVERAGEACTIVEFRACTION_H

#include "utils/observables/Observer.h"

class AverageActiveFraction : public Observer
{
public:
    AverageActiveFraction(const char*, const char*);
    AverageActiveFraction(const AverageActiveFraction& orig);
    virtual ~AverageActiveFraction();
    
    void set_params(const int, std::vector<std::string>);
    double observe(const Box&, std::vector<Cell>&);
        
private:
    static DerivedRegister<AverageActiveFraction> reg;

};

#endif	/* AVERAGEACTIVEFRACTION_H */

