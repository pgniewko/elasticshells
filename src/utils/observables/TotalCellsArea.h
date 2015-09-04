#ifndef TOTALCELLSAREA_H
#define	TOTALCELLSAREA_H

#include "utils/observables/Observer.h"

class TotalCellsArea : public Observer
{
    public:
        TotalCellsArea(const char*, const char*);
        TotalCellsArea(const TotalCellsArea& orig);
        virtual ~TotalCellsArea();

        double observe(Box&, std::vector<Cell>&);
        void set_params(int, ...);
        void set_params(int, std::vector<std::string>);

    private:
        static DerivedRegister<TotalCellsArea> reg;

};

#endif	/* TOTALCELLSAREA_H */