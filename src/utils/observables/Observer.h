#ifndef OBSERVER_H
#define	OBSERVER_H

#include <vector>

#include "Cell.h"
#include "simulation/Box.h"

class Observer {
public:
    Observer();
    Observer(const Observer& orig);
    virtual ~Observer();
    virtual double observe(Box&, std::vector<Cell>& cells) =0;

    const char* getFormat();
    const char* getName();
private:
    const char* observer_name;
    const char* output_format;

};

#endif	/* OBSERVER_H */