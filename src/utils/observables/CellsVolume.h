#ifndef CELLSVOLUME_H
#define	CELLSVOLUME_H

#include "utils/observables/Observer.h"
#include <algorithm>    // std::max

class CellsVolume : public Observer
{
    public:
        explicit CellsVolume(const char*, const char*);
        CellsVolume(const CellsVolume& orig);
        virtual ~CellsVolume();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Shell>&, const DomainList&);

    private:
        double calcVolumeFraction(const Box&, std::vector<Shell>&);
        double calcCellsVolume(std::vector<Shell>&);

        double sphereSphereIntersection(const Shell&, const Shell&);

        static DerivedRegister<CellsVolume> reg;
};

#endif	/* CELLSVOLUME_H */

