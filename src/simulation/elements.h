#ifndef ELEMENTS_H
#define ELEMENTS_H

struct element
{
    int ia, ib, ic;
    double an[3];
    double L2[3];
    double ki[3];
    double ci[3];
};

struct hinge
{
    int x1, x2, x3, x4;
    double D;
    double sinTheta0;
    double theta0;
};

struct object_map
{
    int cell_id;
    int vert_id;
    object_map(int ci, int vi) : cell_id(ci), vert_id(vi) {}
    bool operator==(const object_map &o) const {
        return cell_id == o.cell_id && vert_id == o.vert_id;
    }
    bool operator<(const object_map &o)  const {
        return cell_id < o.cell_id || (cell_id == o.cell_id && vert_id < o.vert_id);
    }
};

#endif /* ELEMENTS_H */

