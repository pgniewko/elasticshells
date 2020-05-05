#ifndef VERTEXTRIANGLE_H
#define	VERTEXTRIANGLE_H

#include "Environment.h"
#include "geometry/Triangle.h"
#include "geometry/Vector3D.h"
#include "geometry/Vertex.h"

class Vertex;

class Element
{
        friend class Restarter;
        friend class Simulator;
        friend class Observer;
    public:
        Element();
        explicit Element(int, int, int);
        Element(const Element& orig);
        virtual ~Element();
        void set_id(int);
        int get_id() const;
        double area(const std::vector<Vertex>&) const;
        double area(const std::vector<Vertex>&, const Vector3D, double) const;
        

        Vector3D normal(const std::vector<Vertex>&) const;
        Vector3D centroid(const std::vector<Vertex>&) const;

        void set_params(const std::vector<Vertex>&, const double, const double, const double);
        void set_sign(int sign_) {sign = sign_;}
        const int get_sign() {return sign;}
        
        int ia = -1;
        int ib = -1;
        int ic = -1;
        int my_id = -1;
        int sign;

        friend std::ostream& operator<< (std::ostream&, const Element&);

    private:
        void set_l2(const std::vector<Vertex>&);
        void set_an(const std::vector<Vertex>&);
        void set_ki(const std::vector<Vertex>&, const double&, const double&, const double&);
        void set_ci(const std::vector<Vertex>&, const double&, const double&, const double&);
        

        double an[3];
        double L2[3];
        double ki[3];
        double ci[3];
};

#endif	/* VERTEXTRIANGLE_H */
