#ifndef BOXSIZEOBSERVER_H
#define	BOXSIZEOBSERVER_H

#include "utils/observables/Observer.h"

class BoxSizeObserver : public Observer
{
    public:
        explicit BoxSizeObserver(const char*, const char*);
        BoxSizeObserver(const BoxSizeObserver& orig);
        virtual ~BoxSizeObserver();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&, const DomainList&);

    private:
        static DerivedRegister<BoxSizeObserver> reg;
};

#endif	/* BOXSIZEOBSERVER_H */

