#ifndef BOX_H
#define	BOX_H


#include <iostream>
#include "Environment.h"
#include "utils/Logger.h"
#include "simulation/Scheduler.h"
#include "geometry/Vector3D.h"

class Box
{
    public:
        explicit Box(double, double, double);
        Box(const Box& orig);
        virtual ~Box();

        void setX(const double);
        double getX() const;
        void setY(const double);
        double getY() const;
        void setZ(const double);
        double getZ() const;

        bool resize(double = 1.0);

        double getVolume(const double = 0.0) const;
        double getArea(const double = 0.0) const;

        void setPbc(bool);
        void setEwall(double);
        void setNu(double);
        double getXEdge(const double = 0.0) const;
        double getYEdge(const double = 0.0) const;
        double getZEdge(const double = 0.0) const;
        double getNu() const;
        double getE() const;

        void configureScheduler(std::string);
        void setDefaultSchedule(int, int, double, double, double, double, double, double);

        static void getDistance(Vector3D&, const Vector3D&, const Vector3D&, const Box&);
        static Vector3D recenteringVector(const Vector3D&, const Box& box);

        void saveRemainingSchedule();

        bool nthTodo();

        bool pbc; // make it private
        
    private:
        double x;
        double y;
        double z;
        double E_box;
        double nu;

        static utils::Logger box_logger;

        Scheduler my_schedule;
};

// TODO: OPTIMIZE THIS FUNCTION - IT'S A CRUCIAL ONE
inline void Box::getDistance(Vector3D& dkj, const Vector3D& vj, const Vector3D& vk, const Box& box) // Vector pointing from vk to vj
{
    dkj = vk - vj;

    if (box.pbc)
    {
        double x, y, z;
        double bsx = 2 * box.getX();
        double bsy = 2 * box.getY();
        double bsz = 2 * box.getZ();
        x = round(dkj.x / bsx) *  bsx;
        y = round(dkj.y / bsy) *  bsy;
        z = round(dkj.z / bsz) *  bsz;
        dkj.x -= x;
        dkj.y -= y;
        dkj.z -= z;
    }
}

inline Vector3D Box::recenteringVector(const Vector3D& cm, const Box& box)
{
    if (box.pbc)
    {
        double x, y, z;
        double bsx = 2 * box.getX();
        double bsy = 2 * box.getY();
        double bsz = 2 * box.getZ();
        x = round(cm.x / bsx) *  bsx;
        y = round(cm.y / bsy) *  bsy;
        z = round(cm.z / bsz) *  bsz;
        Vector3D v_shift(-x, -y, -z);
        return v_shift;
    }
    else
    {
        Vector3D v_shift;
        return v_shift;
    }
}

#endif	/* BOX_H */
