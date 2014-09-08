#ifndef BOX_H
#define	BOX_H

class Box {
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
    
    void resize();
    
private:
    double x;
    double y;
    double z;
    double dx;
    double dy;
    double dz;
};

#endif	/* BOX_H */