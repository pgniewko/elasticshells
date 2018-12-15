#ifndef TOTALCONTACTAREA_H
#define	TOTALCONTACTAREA_H

#include "utils/observables/Observer.h"

class TotalContactArea : public Observer
{
    public:
        explicit TotalContactArea(const char*, const char*);
        TotalContactArea(const TotalContactArea& orig);
        virtual ~TotalContactArea();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Shell>&, const DomainList&);

    private:
        static DerivedRegister<TotalContactArea> reg;

};

#endif	/* TOTALCONTACTAREA_H */