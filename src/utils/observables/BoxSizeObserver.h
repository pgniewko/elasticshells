#ifndef BOXSIZEOBSERVER_H
#define	BOXSIZEOBSERVER_H

#include "utils/observables/Observer.h"

class BoxSizeObserver : public Observer
{
    public:
        explicit BoxSizeObserver(const char*, const char*);
        BoxSizeObserver(const BoxSizeObserver& orig);
        virtual ~BoxSizeObserver();

        double observe(const Box&, const std::vector<Shell>&);

    private:
        static DerivedRegister<BoxSizeObserver> reg;
};

#endif	/* BOXSIZEOBSERVER_H */

