#include "Cell.h"

//Cell::Cell() {
//}

Cell::Cell(list<Triangle> tlist) : tris(tlist) {}

Cell::Cell(const Cell& orig) {
}

Cell::~Cell() {
}

double Cell::surfaceArea()
{
    double surface = 0.0;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) {
        surface += i->area();
    }
    return surface;
}

double Cell::volume()
{
    double volume = 0.0;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) {
        Tetrahedron tetr(i->a, i->b, i->c, cm);
        volume += tetr.volume();
    }
    return volume;
}

void Cell::calcCM()
{
    Vector3D tmp(0,0,0);
    int counter = 0;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
    {
        tmp += i->a;
        tmp += i->b;
        tmp += i->b;
        counter += 3.0;
    }
    tmp /= counter;
    cm = tmp;
}

int Cell::numberofFaces() 
{
    return tris.size();
}

int Cell::numberofVertices() 
{
    
    list<Vector3D> vertices;
    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
    {
      //  cout << "next triangle" <<endl;
        if ( isUnique(vertices, i->a) )
            vertices.push_back(i->a);
        
        if ( isUnique(vertices, i->b) )
            vertices.push_back(i->b);
        
        if ( isUnique(vertices, i->c) )
            vertices.push_back(i->c);
    //     cout << "end of next triangle" <<endl;
    }
//    return 0;
    return vertices.size();
}

bool Cell::isUnique(list<Vector3D>& vlist, Vector3D& v)
{
 //   cout << "inside unique" << endl;
    for(std::list<Vector3D>::iterator i = vlist.begin(); i != vlist.end(); ++i) 
    {
        if (i->x == v.x && i->y == v.y && i->z == v.z)
        {
            return false;
        }
    }
    return true;
}


//void Cell::createDataStructure()
//{
//    list<Vector3D> vectors;
//    int vertexId = 0;
    
//    double xtmp, ytmp, ztmp;
//    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
//    {
//        if ( isUnique(vectors, i->a) )
//        {
//            vectors.push_back(i->a);
//            xtmp = i->a.x;
//            ytmp = i->a.y;
//            ztmp = i->a.z;
//            Vertex v(xtmp, ytmp, ztmp);
//            v.setId(vertexId);
//            vertices.push_back( v );
//            vertexId++;            
//            cout << " xtmp = " << xtmp << " ytmp = " << ytmp << " ztmp = " << ztmp;
//            cout << " vertexId= " << vertexId << endl;
//        }
        
//        if ( isUnique(vectors, i->b) )
//        {
//            vectors.push_back(i->b);
//            xtmp = i->b.x;
//            ytmp = i->b.y;
//            ztmp = i->b.z;
//            Vertex v(xtmp, ytmp, ztmp);
//            v.setId(vertexId);
//            vertices.push_back( v );
//            vertexId++;
//            cout << " xtmp = " << xtmp << " ytmp = " << ytmp << " ztmp = " << ztmp;
//            cout << " vertexId= " << vertexId << endl;
//        }
        
//        if ( isUnique(vectors, i->c) )
//        {
//            vectors.push_back(i->c);
//            xtmp = i->c.x;
//            ytmp = i->c.y;
//            ztmp = i->c.z;
//            Vertex v(xtmp, ytmp, ztmp);
//            v.setId(vertexId);
//            vertices.push_back( v );
//            vertexId++;
//            cout << " xtmp = " << xtmp << " ytmp = " << ytmp << " ztmp = " << ztmp;
//            cout << " vertexId= " << vertexId << endl;
//        }
//    }
//    
//            for(std::list<Vertex>::iterator i = vertices.begin(); i != vertices.end(); ++i) 
//        {
//            i->printVertex();
//        }
//    
//    cout << "num of unique vertices(#1) = "<< vertices.size() << endl;
//    cout << "num of unique vertices(#1) = "<< vertexId << endl;
    //int tmpno = 0;
    //int counterx = 0;
    //for(std::list<Vertex>::iterator j = vertices.begin(); j != vertices.end(); ++j) 
    //{
    //    cout << counterx++ <<  ": j->xyz : " <<j->xyz <<endl;
    //    cout << j->getId() << " ";
    //    cout << j->xyz <<endl;
    //}
    //cout << "trojkaty" <<endl;
//    constructTriangles();
//    cout << "num of new triangles = " << tris.size() << endl;
    
//    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
//    {
//        i->printTriangle();
//    }
//}

//void Cell::constructTriangles()
//{
//    Vector3D * veca;
//    Vector3D * vecb;
//    Vector3D * vecc;
//    list<Triangle> newTris;
//    int triangleId = 0;
//    int counter = 0;
//    for(std::list<Triangle>::iterator i = tris.begin(); i != tris.end(); ++i) 
//    {
        
        //cout << "TROJKAT" <<endl;
        //cout << "i->a" << i->a << endl;
        //cout << "i->b" << i->b << endl;
        //cout << "i->c" << i->c << endl;
        //cout << "=============" <<endl;
//        for(std::list<Vertex>::iterator j = vertices.begin(); j != vertices.end(); ++j) 
//        {
            //cout << counter <<  ": j->xyz : ";
            //cout << j->xyz <<endl;
            //cout << "testyjemy if 1" << endl;
//            if (j->xyz.x == i->a.x && j->xyz.y == i->a.y && j->xyz.z == i->a.z)
//            {
//                j->addTriangle(triangleId);
//                veca = &j->xyz;
//            }
            //cout << "testyjemy if 2" << endl;
//            if (j->xyz.x == i->b.x && j->xyz.y == i->b.y && j->xyz.z == i->b.z)
//            {
//                j->addTriangle(triangleId);
//                vecb = &j->xyz;
//            }
            //cout << "testyjemy if 3" << endl;
//            if (j->xyz.x == i->c.x && j->xyz.y == i->c.y && j->xyz.z == i->c.z)
//            {
//                j->addTriangle(triangleId);
//                vecc = &j->xyz;
//            }
            //counter ++;
            //cout << "koniec testowania" << endl;
//        }
        //cout << "druku druku" << endl;
        //cout << *veca <<endl;
        //cout << *vecb <<endl;
        //cout << *vecc <<endl;
//        cout << "przed kopiando" << endl;
//        Triangle newTriangle(*veca, *vecb, *vecc);
//        cout << "po kopiando" << endl;
        
//        newTriangle.setId(triangleId);
//        newTris.push_back(newTriangle);
        
        
//        triangleId++;
//    }
//    tris.swap(newTris);
//}

