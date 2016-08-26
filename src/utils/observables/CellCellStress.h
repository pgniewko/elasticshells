#ifndef AVERAGECONTACTSTRESSNEW_H
#define	AVERAGECONTACTSTRESSNEW_H

#include "utils/observables/Observer.h"

class CellCellStress : public Observer
{
    public:
        explicit CellCellStress(const char*, const char*);
        CellCellStress(const CellCellStress& orig);
        virtual ~CellCellStress();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);

    private:
        static DerivedRegister<CellCellStress> reg;

};

#endif	/* AVERAGECONTACTSTRESSNEW_H */

