#ifndef FIELD_H
#define	FIELD_H

class Field
{
    public:
        Field();
        Field(const Field& orig);
        virtual ~Field();
        virtual void evolve() {};
};

#endif	/* FIELD_H */
