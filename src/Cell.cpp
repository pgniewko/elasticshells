#include "Cell.h"

//Cell::Cell() : numberV(0), numberT(0)
//{
//}

Cell::Cell(list<Triangle> tris) : numberV(0), numberT(0)
{
    constructVertices(tris);
    constructVTriangles(tris);
    setTopology();
//    printTopology();
}

Cell::Cell(const Cell& orig) 
{
}

Cell::~Cell() 
{
}

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
           // vertices[numberV].
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

void Cell::constructVTriangles(list<Triangle> tris)
{
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
    {
        Vertex *ptra;
        Vertex *ptrb;
        Vertex *ptrc;
        ptra = getVertex(i->a);
        ptrb = getVertex(i->b);
        ptrc = getVertex(i->c);
        VertexTriangle vrxt(ptra, ptrb, ptrc);
        triangles[numberT] = VertexTriangle(vrxt);
        triangles[numberT].setId(numberT);
        numberT++;
    }
}

void Cell::setTopology()
{
//    double k0tmp;
    for (int i = 0; i < numberT; i++)
    {   
        
        Vector3D ab = triangles[i].a->xyz - triangles[i].b->xyz;
        Vector3D ac = triangles[i].a->xyz - triangles[i].c->xyz;
        Vector3D bc = triangles[i].b->xyz - triangles[i].c->xyz;
        
        double abl = ab.length();
        double acl = ac.length();
        double bcl = bc.length();
        
        int aid = triangles[i].a->getId();
        int bid = triangles[i].b->getId();
        int cid = triangles[i].c->getId();
        
        int tid = triangles[i].id;
        
        //cout << "triangle id :" << triangles[i].id << " ; ver ids "<< aid << " " << bid << " "<<cid << endl; 

        triangles[i].a->addNeighbor(bid, abl);
        triangles[i].a->addNeighbor(cid, acl);

        triangles[i].b->addNeighbor(aid, abl);
        triangles[i].b->addNeighbor(cid, bcl);

        triangles[i].c->addNeighbor(aid, acl);
        triangles[i].c->addNeighbor(bid, bcl);
        
        triangles[i].a->addTriangle(tid);
        triangles[i].b->addTriangle(tid);
        triangles[i].c->addTriangle(tid);
        
    }
}

void Cell::calcForces()
{
    cout << "licze harm. forces" << endl;
    double R0;
    double gamma = 1.0;
    int idxj;
    
    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < vertices[i].nneigh; j++)
        {
            R0 = vertices[i].R0[j];
            idxj = vertices[i].neighbors[j];
            vertices[i].force += HookeanForce::calcForce(vertices[idxj].xyz, vertices[i].xyz, R0, gamma);
        }
    }
    
   
    double Rc = 0.5;
    double a  = 1.0;
    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; j < numberV; j++)
        {
            if (i != j && !vertices[i].isBonded(j))
            {
                vertices[i].force += NbRepulsiveForce::calcForce(vertices[j].xyz, vertices[i].xyz, Rc, a);
            }
        }
    }

    double dp = 0.05;
    calcCM();
    for (int i = 0; i < numberT; i++)
    {
        Vector3D fa = OsmoticForce::calcForce(triangles[i].a->xyz, triangles[i].b->xyz, triangles[i].c->xyz, cm, dp);
        Vector3D fb = OsmoticForce::calcForce(triangles[i].b->xyz, triangles[i].c->xyz, triangles[i].a->xyz, cm, dp);
        Vector3D fc = OsmoticForce::calcForce(triangles[i].c->xyz, triangles[i].a->xyz, triangles[i].b->xyz, cm, dp);
        triangles[i].a->force += fa;
        triangles[i].b->force += fb;
        triangles[i].c->force += fc;
    }
}

void Cell::calcForces(const Cell& other_cell)
{    
    double Rc = 0.5;
    double a  = 1.0;
    for (int i = 0; i < numberV; i++)
    {
        for (int j = 0; i < other_cell.numberV; j++)
        {
            vertices[i].force += NbRepulsiveForce::calcForce(other_cell.vertices[j].xyz, vertices[i].xyz, Rc, a);
        }
    }
    
}

void Cell::printTopology()
{
 for (int i = 0; i < numberV; i++)
 {
     vertices[i].printVertex();
     //cout << i << " " << vertices[i].nneigh << endl;
 }
}

Vertex * Cell::getVertex(const Vector3D v)
{
    for (int i = 0; i < numberV; i++)
    {
        if (vertices[i].xyz.x == v.x && vertices[i].xyz.y == v.y && vertices[i].xyz.z == v.z)
        {
            return &vertices[i];
        }
    }
    return NULL;
 }


bool Cell::isUnique(list<Vector3D>& vlist, Vector3D& v)
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

double Cell::calcSurfaceArea()
{
    double surface = 0.0;
    for (int i = 0; i < numberT; i++)
    {
        surface += triangles[i].area();
    }

    return surface;
}

double Cell::calcVolume()
{
    double volume = 0.0;
    for (int i = 0; i < numberT; i++)
    {
        Tetrahedron tetr(triangles[i].a->xyz, triangles[i].b->xyz, triangles[i].c->xyz, cm);
        volume += tetr.volume();
    }

    return volume;
}

void Cell::calcCM()
{
    Vector3D tmp(0,0,0);
    
    double M = 0.0;
    for (int i = 0; i < numberV; i++)
    {
        tmp += vertices[i].xyz;
        M += vertices[i].getMass();
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
                //if(iidx > jidx )
                //{
                  os << "cmd.do(\"bond /cells///UNK`/H"<< iidx << ", /cells///UNK`/H" << jidx << "\")\n";
                //}
            
        }
    }
    
    os << "cmd.do(\"show lines\")\n";
    os << "cmd.do(\"bg white\")\n";
    os.close();
}