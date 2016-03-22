#ifndef CELLSVOLUME_H
#define	CELLSVOLUME_H

#include "utils/observables/Observer.h"

class CellsVolume : public Observer
{
    public:
        CellsVolume(const char*, const char*);
        CellsVolume(const CellsVolume& orig);
        virtual ~CellsVolume();

        double observe(const Box&, std::vector<Cell>&);
        void set_params(const int, std::vector<std::string>);

    private:
        double calcVolumeFraction(const Box&, std::vector<Cell>&);
        double calcCellsVolume(std::vector<Cell>&);

        static DerivedRegister<CellsVolume> reg;

};

#endif	/* CELLSVOLUME_H */

