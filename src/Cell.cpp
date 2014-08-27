#include "Cell.h"

Cell::Cell(int depth) : numberV(0), numberT(0)
{
    SimpleTriangulation sm(depth);
    list<Triangle> tris = sm.triangulate();
    constructVertices(tris);
    constructVTriangles(tris);
    constructTopology();
    calcVolume();
}

Cell::Cell(list<Triangle> tris) : numberV(0), numberT(0)
{
    constructVertices(tris);
    constructVTriangles(tris);
    constructTopology();
    calcVolume();
}

Cell::Cell(const Cell& orig) : cm(orig.cm), vertices(orig.vertices), triangles(orig.triangles),
                               numberV(orig.numberV), numberT(orig.numberT)
{}

Cell::~Cell() {}

void Cell::constructVertices(list<Triangle> tris)
{
    list<Vector3D> vectors;
    double xtmp, ytmp, ztmp;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
    {
        if ( isUnique(vectors, i->a) )
        {
          vectors.push_back(i->a);
          xtmp = i->a.x;
          ytmp = i->a.y;
          ztmp = i->a.z;
          vertices[numberV] = Vertex(xtmp, ytmp, ztmp);
          vertices[numberV].setId(numberV);
          numberV++;
        }
        if ( isUnique(vectors, i->b) )
        {
          vectors.push_back(i->b);
          xtmp = i->b.x;
          ytmp = i->b.y;
          ztmp = i->b.z;
          vertices[numberV] = Vertex(xtmp, ytmp, ztmp);
          vertices[numberV].setId(numberV);
          numberV++;
        }
        if ( isUnique(vectors, i->c) )
        {
          vectors.push_back(i->c);
          xtmp = i->c.x;
          ytmp = i->c.y;
          ztmp = i->c.z;
          vertices[numberV] = Vertex(xtmp, ytmp, ztmp);
          vertices[numberV].setId(numberV);
          numberV++;
        }
    }
}

bool Cell::isUnique(list<Vector3D>& vlist, const Vector3D& v)
{
    for(std::list<Vector3D>::iterator i = vlist.begin(); i != vlist.end(); ++i) 
    {
        if (i->x == v.x && i->y == v.y && i->z == v.z)
        {
            return false;
        }
    }
    return true;
}

void Cell::constructVTriangles(list<Triangle> tris)
{
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
    {
        int va = getVertex(i->a);
        int vb = getVertex(i->b);
        int vc = getVertex(i->c);
        VertexTriangle vrxt(va, vb, vc);
        triangles[numberT] = VertexTriangle(vrxt);
        triangles[numberT].setId(numberT);
        numberT++;
    }
}

int Cell::getVertex(const Vector3D v)
{
    for (int i = 0; i < numberV; i++)
    {
        if (vertices[i].xyz.x == v.x && vertices[i].xyz.y == v.y && vertices[i].xyz.z == v.z)
        {
            return vertices[i].getId();
        }
    }
    return -1;
 }


void Cell::constructTopology()
{
    for (int i = 0; i < numberT; i++)
    {   
        int aid = triangles[i].ia;
        int bid = triangles[i].ib;
        int cid = triangles[i].ic;
        
        Vector3D ab = vertices[aid].xyz - vertices[bid].xyz;
        Vector3D ac = vertices[aid].xyz - vertices[cid].xyz;
        Vector3D bc = vertices[bid].xyz - vertices[cid].xyz;
        
        double abl = ab.length();
        double acl = ac.length();
        double bcl = bc.length();
        
        int tid = triangles[i].id;

        vertices[aid].addNeighbor(bid, abl);
        vertices[aid].addNeighbor(cid, acl);

        vertices[bid].addNeighbor(aid, abl);
        vertices[bid].addNeighbor(cid, bcl);
               
        vertices[cid].addNeighbor(aid, acl);
        vertices[cid].addNeighbor(bid, bcl);
        
        vertices[aid].addTriangle(tid);
        vertices[bid].addTriangle(tid);
        vertices[cid].addTriangle(tid);
        
    }
}

void Cell::calcForces()
{
    double R0ij;
    int idxj;
    
    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < vertices[i].nneigh; j++)
        {
            R0ij = vertices[i].R0[j];
            idxj = vertices[i].neighbors[j];
            vertices[i].force += HookeanForce::calcForce(vertices[idxj].xyz, vertices[i].xyz, R0ij, gamma);
        }
    }
    
    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < numberV; j++)
        {
            if (i != j && !vertices[i].isNeighbor(j))
            {
                vertices[i].force += NbRepulsiveForce::calcForce(vertices[j].xyz, vertices[i].xyz, Rc, a);
            }
        }
    }

    
    calcCM();
    int iva, ivb, ivc;
    
    for (int i = 0; i < numberT; i++)
    {
        iva = triangles[i].ia;
        ivb = triangles[i].ib;
        ivc = triangles[i].ic;
        
        Vector3D fa = OsmoticForce::calcForce(vertices[iva].xyz, vertices[ivb].xyz, vertices[ivc].xyz, cm, dp);
        Vector3D fb = OsmoticForce::calcForce(vertices[ivb].xyz, vertices[ivc].xyz, vertices[iva].xyz, cm, dp);
        Vector3D fc = OsmoticForce::calcForce(vertices[ivc].xyz, vertices[iva].xyz, vertices[ivb].xyz, cm, dp);
        
        vertices[iva].force += fa;
        vertices[ivb].force += fb;
        vertices[ivc].force += fc;
    }
}

void Cell::calcForces(const Cell& other_cell)
{    
    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; i < other_cell.numberV; j++)
        {
            vertices[i].force += NbRepulsiveForce::calcForce(other_cell.vertices[j].xyz, vertices[i].xyz, Rc, a);
        }
    }
    
}

void Cell::printCell()
{
    cout << "Center of mass = " << cm << endl;
    cout << "Number of vertices = " << numberV << endl;
    cout << "Number of triangles = " << numberT << endl;
    cout << "SURFACE AREA = " << calcSurfaceArea() << endl;
    cout << "VOLUME = " << calcVolume() << endl;
    
    for (int i = 0; i < numberV; i++)
    {
        vertices[i].printVertex();
    }
    
    for (int i = 0; i < numberT; i++)
    {
        triangles[i].printVertexTriangle();
    }
    
}

double Cell::calcSurfaceArea()
{
    double surface = 0.0;
    for (int i = 0; i < numberT; i++)
    {
        surface += triangles[i].area(vertices);
    }

    return surface;
}

double Cell::calcVolume()
{
    double volume = 0.0;
    int va, vb, vc;
    for (int i = 0; i < numberT; i++)
    {
        va = triangles[i].ia;
        vb = triangles[i].ib;
        vc = triangles[i].ic;
        Tetrahedron tetr(vertices[va].xyz, vertices[vb].xyz, vertices[vc].xyz, cm);
        volume += tetr.volume();
    }

    return volume;
}

void Cell::calcCM()
{
    Vector3D tmp(0.0, 0.0, 0.0);
    
    double M = 0.0;
    double m;
    for (int i = 0; i < numberV; i++)
    {
        m = vertices[i].getMass();
        tmp += m * vertices[i].xyz;
        M += m;
    }

    tmp /= M;
    cm = tmp;
}

void Cell::addVelocity(const Vector3D& nv)
{
    for (int i = 0; i < numberV; i++)
    {
        vertices[i].velocity += nv;
    }
}

void Cell::addXYZ(const Vector3D& nxyz)
{
    for (int i = 0; i < numberV; i++)
    {
        vertices[i].xyz += nxyz;
    }    
}

int Cell::numberofFaces() 
{
    return numberT;
}

int Cell::numberofVertices() 
{
    return numberV;
}


void Cell::saveTriangulatedSurface(const char* filename)
{
    int index;
    ofstream os(filename);
    os << numberV << "\n" ;
    for (int i = 0; i < numberV; i++)
    {
        index = (vertices[i].getId()+1) ;
        os << "H" << index << " "<< vertices[i].xyz.x << " " << vertices[i].xyz.y << " " << vertices[i].xyz.z << "\n";
    }
    os.close();
}

void Cell::saveRenderingScript(const char* filename, const char* cellsfile)
{
    ofstream os(filename);
    os << "from pymol.cgo import *\n";
    os << "from pymol import cmd \n\n";
    os << "cmd.do(\"load " << cellsfile << ", cells\")\n";
    os << "cmd.do(\"hide all\")\n";
    os << "cmd.do(\"set sphere_color, tv_red\")\n";
    os << "cmd.do(\"set line_color, marine\")\n";
    os << "cmd.do(\"show spheres\")\n";
    os << "cmd.do(\"alter elem h, vdw=0.1\")\n";
    os << "cmd.do(\"rebuild\")\n";


    int iidx, jidx;
    
    for (int i = 0; i < numberV; i++)
    {
        iidx = (vertices[i].getId()+1);
        for (int j = 0; j < vertices[i].nneigh; j++)
        {
            jidx = (vertices[i].neighbors[j]+1);
            os << "cmd.do(\"bond /cells///UNK`/H"<< iidx << ", /cells///UNK`/H" << jidx << "\")\n";
        }
    }
    
    os << "cmd.do(\"show lines\")\n";
    os << "cmd.do(\"bg white\")\n";
    os.close();
}