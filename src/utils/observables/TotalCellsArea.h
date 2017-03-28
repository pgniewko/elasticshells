#ifndef TOTALCELLSAREA_H
#define	TOTALCELLSAREA_H

#include "utils/observables/Observer.h"

class TotalCellsArea : public Observer
{
    public:
        explicit TotalCellsArea(const char*, const char*);
        TotalCellsArea(const TotalCellsArea& orig);
        virtual ~TotalCellsArea();

        double observe(const Box&, std::vector<Cell>&, const DomainList&);
        void set_params(const int, std::vector<std::string>);

    private:
        static DerivedRegister<TotalCellsArea> reg;

};

#endif	/* TOTALCELLSAREA_H */