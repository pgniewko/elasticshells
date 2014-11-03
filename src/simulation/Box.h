#ifndef BOX_H
#define	BOX_H

#include <iostream>

#include "Environment.h"

class Box
{
    public:
        Box(double, double, double);
        Box(double, double, double, double);
        Box(const Box& orig);
        virtual ~Box();

        void setX(const double);
        double getX();
        void setY(const double);
        double getY();
        void setZ(const double);
        double getZ();

        void setDx(const double);
        double getDx();
        void setDy(const double);
        double getDy();
        void setDz(const double);
        double getDz();

        void setXend(double);
        void setYend(double);
        void setZend(double);

        double getXend();
        double getYend();
        double getZend();

        void resize();

        double getVolume();
        double getVolume(double);
        double getArea();
        double getArea(double);
        
        void setPbc(bool);
        void setEcw(double);
        bool pbc;
        double ecw;

    private:
        double x;
        double y;
        double z;
        double xe;
        double ye;
        double ze;
        double dx;
        double dy;
        double dz;
};

#endif	/* BOX_H */