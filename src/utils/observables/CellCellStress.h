#ifndef AVERAGECONTACTSTRESSNEW_H
#define	AVERAGECONTACTSTRESSNEW_H

#include "utils/observables/Observer.h"

class CellCellStress : public Observer
{
    public:
        CellCellStress(const char*, const char*);
        CellCellStress(const CellCellStress& orig);
        virtual ~CellCellStress();
        
        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);
        double observe(Box&, std::vector<Cell>&);
        
    private:
        static DerivedRegister<CellCellStress> reg;

};

#endif	/* AVERAGECONTACTSTRESSNEW_H */

