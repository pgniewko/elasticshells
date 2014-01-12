#ifndef FIELD_H
#define	FIELD_H

class Field {
    public:

        Field();
        Field(const Field& orig);
        virtual ~Field();
//    void set_dt(double dt) {this->dt = dt;};

        virtual bool is_static() {
            return (static_field);
        };
        virtual void evolve() = 0;

    private:
        bool static_field;
        double x, y, z;
        double dx, dy, dz;
        int nodes_no;
        int n_x, n_y, n_z;
//    double grid [][][];
//    double dt = 0.1;


};

#endif	/* FIELD_H */

