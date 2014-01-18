#include "Pairlist.h"

using namespace std;

Pairlist::Pairlist(double _mindist, int _steps)
{
    mindist = _mindist;
    steps = _steps;
    n = 0;
    margin = Vector3D(_mindist, _mindist, _mindist);
}

//Pairlist::~Pairlist()
//{
//}

int Pairlist::num_boxes() const
{
    return xb * yb * zb;
}

inline void find_minmax(const Vector3D* r, int n, Vector3D& min, Vector3D& max)
{
    min = max = r[0];

    for (int i = 1; i < n; i++)
    {
        const Vector3D* tmp = r + i;

        if (tmp->x < min.x)
        {
            min.x = tmp->x;
        }
        else if (tmp->x > max.x)
        {
            max.x = tmp->x;
        }

        if (tmp->y < min.y)
        {
            min.y = tmp->y;
        }
        else if (tmp->y > max.y)
        {
            max.y = tmp->y;
        }

        if (tmp->z < min.z)
        {
            min.z = tmp->z;
        }
        else if (tmp->z > max.z)
        {
            max.z = tmp->z;
        }
    }
}

int Pairlist::findbox(Vector3D p, const Box& dist) const
{
    p = dist.delta(p) - lmin;
    int axb = (int)( p.x / pairdist.x);
    int ayb = (int)( p.y / pairdist.y);
    int azb = (int)( p.z / pairdist.z);
    return  azb * xytotb + ayb * xb + axb;
}

void Pairlist::set_neighbors(const Box& domain)
{
    static  int oxb = - 1, oyb = -1, ozb = - 1;

    if (xb != oxb || yb != oyb || zb != ozb)
    {
        //boxes changed

        //IMPORTANT NOTE:
        //Resizing does not require to reallocate, but it means that
        //what is left has to be cleared with the loop

        nbrlist.resize(xb * yb * zb);
        fullnbrlist.resize(xb * yb * zb);

        for (int i = 0; i < xb * yb * zb; i++)
        {
            fullnbrlist[i].clear();
            nbrlist[i].clear();
        }

        xytotb = xb * yb;
        oxb = xb;
        oyb = yb;
        ozb = zb;

        for (int zi = 0; zi < zb; zi++)
            for (int yi = 0; yi < yb; yi++)
                for (int xi = 0; xi < xb; xi++)
                    for (int zz = -1; zz <= 1; zz++)
                        for (int yy = -1; yy <= 1; yy++)
                            for (int xx = -1; xx <= 1; xx++)
                            {
                                bool cond =  (xx  > 0) || (xx == 0 && yy  == 1) || ( zz > -1 && xx == 0 && yy == 0);
                                int axx = xx + xi;
                                int ayy = yy + yi;
                                int azz = zz + zi;
                                bool flag = true;

                                if (domain.a())
                                {
                                    axx = axx - int (xb * floor (float (axx) / float (xb)));
                                }
                                else if (axx < 0 || axx >= xb)
                                {
                                    flag = false;
                                }

                                if (domain.b())
                                {
                                    ayy = ayy - int (yb * floor (float (ayy) / float (yb)));
                                }
                                else if (ayy < 0 || ayy >= yb)
                                {
                                    flag = false;
                                }

                                if (domain.c())
                                {
                                    azz = azz - int (zb * floor (float (azz) / float (zb)));
                                }
                                else if (azz < 0 || azz >= zb)
                                {
                                    flag = false;
                                }

                                if (cond && flag)
                                {
                                    nbrlist[ zi * xytotb + yi * xb + xi ].push_back (  azz * xytotb + ayy * xb + axx );
                                }

                                if (flag)
                                {
                                    fullnbrlist[ zi * xytotb + yi * xb + xi ].push_back (  azz * xytotb + ayy * xb + axx );
                                }
                            }
    }
}


void Pairlist::compute(const Box& domain, const Vector3D* r, int np, int* aindex)
{
    double lx, ly, lz;
    lx = domain.a(); //set box size
    ly = domain.b();
    lz = domain.c();
    lmin = domain.origin() - Vector3D(0.5 * lx, 0.5 * ly, 0.5 * lz); //in this case domain is centered
    Vector3D min, max, side;

    if (!domain.a() || !domain.b() || !domain.c())    //at least one direction is NOT periodic
    {
        find_minmax(r, np, min, max);
        max += margin;// 0.1A margin
        min -= margin;
        side = max - min ;

        if (!domain.a())
        {
            lx = side.x;
            lmin.x = min.x;
        }

        if (!domain.b())
        {
            ly = side.y;
            lmin.y = min.y;
        }

        if (!domain.c())
        {
            lz = side.z;
            lmin.z = min.z;
        }
    }

    xb = (int)(lx / mindist);

    if (xb == 0)
    {
        xb = 1;
    }

    pairdist.x = lx / double(xb);

    yb = (int)(ly / mindist);

    if (yb == 0)
    {
        yb = 1;
    }

    pairdist.y = ly / double(yb);

    zb = (int)(lz / mindist);

    if (zb == 0)
    {
        zb = 1;
    }

    pairdist.z = lz / double(zb);

    set_neighbors(domain);

    boxlist.resize(xb * yb * zb);

    for (int i = 0; i < xb * yb * zb; i++)
    {
        boxlist[i].clear();
    }


    for (int k = 0; k < np; k++)
    {
        int box = findbox(r[k], domain);

//	if (box >= boxlist.size())
//        {
//	   cout<<" BOX ERROR "<<box<<" pos "<<r[k]<<" delta "<<domain.delta(r[k])<<endl;
//	}

        int anum;

        if (aindex != 0)
        {
            anum = aindex[k];
        }
        else
        {
            anum = k;
        }

        boxlist[box].push_back(anum);
    }

    n++;
}

void Pairlist::compute_pairs(Pairlist& pl, Box& domain, vector<Cell>& cells, double RCUT, vector<Pair>& pairs)
{
    Vector3D pos[cells.size()];

    for (int i = 0; i < cells.size(); i++)
    {
        pos[i] = cells[i].r;
    }

    pl.compute(domain, pos, cells.size());
    pairs.clear();

    for (int k = 0; k < pl.boxlist.size(); k++)
    { //iterate over boxes
        for (int nbr = 0; nbr < pl.nbrlist[k].size(); nbr++)
        {
            //iterate over nbr boxes
            int k1 =  pl.nbrlist[k][nbr];

            for (int l = 0; l < pl.boxlist[ k ].size(); l++)  //iterate over atoms
            {
                int i =  pl.boxlist[k][l];
                int startl1 = 0;

                if (k1 == k)
                {
                    startl1 = l + 1;
                }

                for (int l1 = startl1; l1 < pl.boxlist[k1].size(); l1++)
                {
                    //iterate over atoms
                    int j = pl.boxlist[k1][l1];
                    Vector3D r_kl = domain.delta(cells[i].r, cells[j].r);
                    double R_kl = r_kl.length();

                    if (R_kl < RCUT && i != j)
                    {
                        pairs.push_back(Pair(i, j));
                    }
                }
            }
        }
    }
}