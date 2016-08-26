#ifndef BOXSIZEOBSERVER_H
#define	BOXSIZEOBSERVER_H

#include "utils/observables/Observer.h"

class BoxSizeObserver : public Observer
{
    public:
        explicit BoxSizeObserver(const char*, const char*);
        BoxSizeObserver(const BoxSizeObserver& orig);
        virtual ~BoxSizeObserver();

        double observe(const Box&, std::vector<Cell>&);
        void set_params(const int, std::vector<std::string>);

    private:
        static DerivedRegister<BoxSizeObserver> reg;

};

#endif	/* BOXSIZEOBSERVER_H */

