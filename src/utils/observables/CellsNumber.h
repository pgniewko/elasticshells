#ifndef CELLSNUMBER_H
#define CELLSNUMBER_H

#include "utils/observables/Observer.h"

class CellsNumber : public Observer
{
    public:
        explicit CellsNumber(const char*, const char*);
        CellsNumber(const CellsNumber& orig);
        virtual ~CellsNumber();
    
        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&, const DomainList&);
        
    private:
        static DerivedRegister<CellsNumber> reg;

};

#endif /* CELLSNUMBER_H */

