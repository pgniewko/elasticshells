#ifndef BOX_H
#define	BOX_H

#include "Vector3D.h"

class Box
{
    public:
        Box(void) : lx(0.0), ly(0.0), lz(0.0) {}
        Box(double l_x, double l_y, double l_z, Vector3D o);
//        virtual ~Box();

        void set(double l_x, double l_y, double l_z, Vector3D o);

        Vector3D delta(const Vector3D& p1) const;
        Vector3D delta(const Vector3D& p1, const Vector3D& p2) const;
        Vector3D image(const Vector3D& p1) const;


        double a() const {
            return lx;
        }

        double b() const {
            return ly;
        }

        double c() const {
            return lz;
        }

        Vector3D origin() const {
            return O;
        }
        double volume() const {
            return lx * ly * lz;
        }

    private:
        double lx, ly, lz;
        double ilx, ily, ilz;
        Vector3D O;
};

#endif	/*BOX_H */

