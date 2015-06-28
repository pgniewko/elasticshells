#ifndef TOTALCELLSAREA_H
#define	TOTALCELLSAREA_H

#include <vector>
#include "Cell.h"

class TotalCellsArea
{
    public:
        TotalCellsArea();
        TotalCellsArea(const TotalCellsArea& orig);
        virtual ~TotalCellsArea();
        static double totalCellArea(std::vector<Cell>&);
    private:

};

#endif	/* TOTALCELLSAREA_H */