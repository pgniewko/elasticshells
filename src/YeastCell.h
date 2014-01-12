#ifndef YEASTCELL_H
#define	YEASTCELL_H

#include "Cell.h"

class YeastCell : public Cell {
    public:
        YeastCell();
        YeastCell(double);
        YeastCell(const YeastCell& orig);
        virtual ~YeastCell();
        double calc_volume();
//    void set_coor(double, double, double);
        void set_params(double*, int);

    private:

};

#endif	/* YEASTCELL_H */

