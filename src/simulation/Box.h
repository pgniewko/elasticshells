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

        void setXstart(double);
        void setYstart(double);
        void setZstart(double);
        void setXend(double);
        void setYend(double);
        void setZend(double);

        double getXstart();
        double getYstart();
        double getZstart();
        double getXend();
        double getYend();
        double getZend();

        void resize();

        //double getVolume();
        //double getArea();
        double getVolume(double = 0.0);
        double getArea(double = 0.0);

        void setPbc(bool);
        void setEwall(double);
        void setNu(double);
        bool pbc;
//        double E_box;
//        double r_;
        double getXEdge(double = 0.0);
        double getYEdge(double = 0.0);
        double getZEdge(double = 0.0);
        double getNu();
        double getE();

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
};

#endif	/* BOX_H */
