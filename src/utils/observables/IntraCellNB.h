#ifndef INTRACELLNB_H
#define	INTRACELLNB_H

#include "utils/observables/Observer.h"

class IntraCellNB : public Observer
{
    public:
        IntraCellNB(const char*, const char*);
        IntraCellNB(const IntraCellNB& orig);
        virtual ~IntraCellNB();

        void set_params(const int, std::vector<std::string>);
        double observe(const Box&, std::vector<Cell>&);

    private:
        static DerivedRegister<IntraCellNB> reg;
};

#endif	/* INTRACELLNB_H */

