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

        void set_x(const double);
        double get_x() const;
        void set_y(const double);
        double get_y() const;
        void set_z(const double);
        double get_z() const;

        void setXmin(const double);
        void setYmin(const double);
        void setZmin(const double);

        double getXmin() const;
        double getYmin() const;
        double getZmin() const;

        bool resize(double = 1.0);

        double get_volume(const double = 0.0) const;
        double get_area(const double = 0.0) const;

        void setPbc(bool);
        void setEwall(double);
        void setNu(double);
        bool pbc;
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

    private:
        double x;
        double y;
        double z;
        double x_min;
        double y_min;
        double z_min;
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
        double bsx = 2 * box.get_x();
        double bsy = 2 * box.get_y();
        double bsz = 2 * box.get_z();
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
        double bsx = 2 * box.get_x();
        double bsy = 2 * box.get_y();
        double bsz = 2 * box.get_z();
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
