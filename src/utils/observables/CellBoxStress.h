#ifndef CELLBOXSTRESS_H
#define	CELLBOXSTRESS_H

#include "utils/observables/Observer.h"

class CellBoxStress : public Observer
{
    public:
        explicit CellBoxStress(const char*, const char*);
        CellBoxStress(const CellBoxStress& orig);
        virtual ~CellBoxStress();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);

    private:
        static DerivedRegister<CellBoxStress> reg;

};

#endif	/* CELLBOXSTRESS_H */
