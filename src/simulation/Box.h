#ifndef BOX_H
#define	BOX_H


#include <iostream>
#include "Environment.h"
#include "utils/Logger.h"
#include "simulation/Scheduler.h"

class Box
{
    public:
        Box(double, double, double);
//        Box(double, double, double, double);
        Box(const Box& orig);
        virtual ~Box();

        void setX(const double);
        double getX() const;
        void setY(const double);
        double getY() const;
        void setZ(const double);
        double getZ() const;

        void setXmax(const double);
        void setYmax(const double);
        void setZmax(const double);
        void setXmin(const double);
        void setYmin(const double);
        void setZmin(const double);

        double getXmax() const;
        double getYmax() const;
        double getZmax() const;
        double getXmin() const;
        double getYmin() const;
        double getZmin() const;

        bool resize();

        double getVolume(const double = 0.0) const;
        double getArea(const double = 0.0) const;

        void setPbc(bool);
        void setEwall(double);
        void setNu(double);
        bool pbc;
        double getXEdge(const double = 0.0) const;
        double getYEdge(const double = 0.0) const;
        double getZEdge(const double = 0.0) const;
        double getNu() const;
        double getE() const;

        void configureScheduler(char*);
        void setDefaultSchedule(int, int, double, double, double, double, double, double);

    private:
        double x;
        double y;
        double z;
        double x_max;
        double y_max;
        double z_max;
        double x_min;
        double y_min;
        double z_min;
        double E_box;
        double nu;

        static utils::Logger box_logger;

        Scheduler my_schedule;
};

#endif	/* BOX_H */
