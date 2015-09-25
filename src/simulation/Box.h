#ifndef BOX_H
#define	BOX_H

#include <iostream>

#include "Environment.h"

struct schedule_t
{
    double dx;
    double dy;
    double dz;
    double delta;
    int life_time;
    int boxs;
    int counter;
};


class Box
{
    public:
        Box(double, double, double);
        Box(double, double, double, double);
        Box(const Box& orig);
        virtual ~Box();

        void setX(const double);
        double getX() const;
        void setY(const double);
        double getY() const;
        void setZ(const double);
        double getZ() const;

        void setDx(const double);
        double getDx() const;
        void setDy(const double);
        double getDy() const;
        void setDz(const double);
        double getDz() const;

        void setXstart(const double);
        void setYstart(const double);
        void setZstart(const double);
        void setXend(const double);
        void setYend(const double);
        void setZend(const double);

        double getXstart() const;
        double getYstart() const;
        double getZstart() const;
        double getXend() const;
        double getYend() const;
        double getZend() const;

        void resize();

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

    private:
        double x;
        double y;
        double z;
        double xs;
        double ys;
        double zs;
        double xe;
        double ye;
        double ze;
        double dx;
        double dy;
        double dz;
        double E_box;
        double nu;
        
        class Scheduler
        {
            public:
                void readSchedule();
                void registerSchedules();
                void checkSchedule();
                void changeSchedule();
                int total_time;
                int recent_time;
                    
                
                
        };
};

#endif	/* BOX_H */
