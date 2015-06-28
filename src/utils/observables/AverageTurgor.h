#ifndef AVERAGETURGOR_H
#define	AVERAGETURGOR_H

#include <vector>
#include "Cell.h"

class AverageTurgor {
public:
    AverageTurgor();
    AverageTurgor(const AverageTurgor& orig);
    virtual ~AverageTurgor();
    static double populationAverageTurgor(std::vector<Cell>&);
private:
};

#endif	/* AVERAGETURGOR_H */