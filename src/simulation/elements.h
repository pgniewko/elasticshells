#ifndef ELEMENTS_H
#define ELEMENTS_H

struct element
{
    int ia, ib, ic;
    double an[3];
    double L2[3];
    double ki[3];
    double ci[3];
    int sign;
};

struct hinge
{
    int v1, v2, v3, v4;
    double D;
    double sinTheta0;
    double theta0;
};

struct object_map
{
    int shell_id;
    int vert_id;
    object_map(int ci, int vi) : shell_id(ci), vert_id(vi) {}
    bool operator==(const object_map& o) const
    {
        return shell_id == o.shell_id && vert_id == o.vert_id;
    }
    bool operator<(const object_map& o)  const
    {
        return shell_id < o.shell_id || (shell_id == o.shell_id && vert_id < o.vert_id);
    }
};

#endif /* ELEMENTS_H */

