#include "VertexTests.h"

CPPUNIT_TEST_SUITE_REGISTRATION(VertexTests);

VertexTests::VertexTests() {}

VertexTests::~VertexTests() {}

void VertexTests::setUp() 
{
    v1 = new Vertex;
    v2 = new Vertex(0,0,1);
    v2->setMass(1.23);
    v3 = new Vertex(-1.0, 1.0, 2.0);
    
    v4 = new Vertex(11.0, 12.3, -9.99);
    v4->setId(18);
    v4->setMass(0.87);
    v4->addNeighbor(0, 0.76);
    v4->addNeighbor(0, 0.44);
    v4->addNeighbor(1, 0.44);
    v4->addNeighbor(2, 0.44);
    v4->addNeighbor(5, 0.33);
    
    vec1 = new Vector3D(1, 2, 3);
}

void VertexTests::tearDown() 
{
    delete v1;
    delete v2;
    delete v3;
    delete v4;
    //delete v5;
    //delete v6;
    
    delete vec1;
}

void VertexTests::testMethod() 
{
    
}

void VertexTests::testFailedMethod() 
{
    //CPPUNIT_ASSERT(true);
    CPPUNIT_ASSERT(false);
}

void VertexTests::testConstructor()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->xyz.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->xyz.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->xyz.z, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->force.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->force.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->force.z, DELTA14);
    CPPUNIT_ASSERT_EQUAL(-1, v1->getId());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, v1->getMass(), DELTA14);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2->xyz.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2->xyz.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, v2->xyz.z, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2->velocity.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2->velocity.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2->velocity.z, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.23, v2->getMass(), DELTA14);
    CPPUNIT_ASSERT_EQUAL(-1, v2->getId());
    CPPUNIT_ASSERT_EQUAL(0, v2->ntrian);
    CPPUNIT_ASSERT_EQUAL(0, v2->nneigh);
    CPPUNIT_ASSERT_EQUAL(0, v2->nbneigh);
    
    Vertex tmpv1(*vec1);
    CPPUNIT_ASSERT_EQUAL(-1, tmpv1.getId());
    tmpv1.setId(1);
    tmpv1.setMass(12.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tmpv1.xyz.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, tmpv1.xyz.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, tmpv1.xyz.z, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, tmpv1.getMass(), DELTA14);
    CPPUNIT_ASSERT_EQUAL(1, tmpv1.getId());
}

void VertexTests::testCopyConstructor()
{
    Vertex tmpv1(*v3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, tmpv1.xyz.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tmpv1.xyz.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, tmpv1.xyz.z, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpv1.velocity.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpv1.velocity.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpv1.velocity.z, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tmpv1.getMass(), DELTA14);
    CPPUNIT_ASSERT_EQUAL(-1, tmpv1.getId());
    CPPUNIT_ASSERT_EQUAL(0, tmpv1.ntrian);
    CPPUNIT_ASSERT_EQUAL(0, tmpv1.nneigh);
    CPPUNIT_ASSERT_EQUAL(0, tmpv1.nbneigh);
}


void VertexTests::testAddNeighbor()
{
    CPPUNIT_ASSERT_EQUAL(18, v4->getId());
    CPPUNIT_ASSERT_EQUAL(4, v4->getNumNeighbors());
    CPPUNIT_ASSERT_EQUAL(0, v4->getNumVTriangles());
    
    CPPUNIT_ASSERT_EQUAL(0, v4->getNeighborId(0));
    CPPUNIT_ASSERT_EQUAL(1, v4->getNeighborId(1));
    CPPUNIT_ASSERT_EQUAL(2, v4->getNeighborId(2));
    CPPUNIT_ASSERT_EQUAL(5, v4->getNeighborId(3));
}
