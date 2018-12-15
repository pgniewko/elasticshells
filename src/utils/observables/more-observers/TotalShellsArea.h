#ifndef TOTALSHELLSAREA_H
#define	TOTALSHELLSAREA_H

#include "utils/observables/Observer.h"

class TotalShellsArea : public Observer
{
    public:
        explicit TotalShellsArea(const char*, const char*);
        TotalShellsArea(const TotalShellsArea& orig);
        virtual ~TotalShellsArea();

        double observe(const Box&, std::vector<Shell>&, const DomainList&);
        void set_params(const int, std::vector<std::string>);

    private:
        static DerivedRegister<TotalShellsArea> reg;

};

#endif	/* TOTALSHELLSAREA_H */