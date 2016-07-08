#include "VertexTests.h"

CPPUNIT_TEST_SUITE_REGISTRATION(VertexTests);

VertexTests::VertexTests() {}

VertexTests::~VertexTests() {}

void VertexTests::setUp()
{
    vec1 = new Vector3D(1, 2, 3);
    v1 = new Vertex;
    v2 = new Vertex(0, 0, 1);
    //v2->setMass(1.23);
    v3 = new Vertex(-1.0, 1.0, 2.0);
    v4 = new Vertex(11.0, 12.3, -9.99);
    v4->setId(18);
    //v4->setMass(0.87);
    v4->addNeighbor(0, 0.76);
    v4->addNeighbor(0, 0.44);
    v4->addNeighbor(1, 0.44);
    v4->addNeighbor(2, 0.44);
    v4->addNeighbor(5, 0.33);
    v5 = new Vertex(-11.0, 345.0, 666.0);
    v5->setId(0);
    //v5->setMass(0.99);
    v5->addTriangle(1);
    v5->addTriangle(0);
    v5->addTriangle(1);
    v5->addTriangle(0);
    v5->addTriangle(100);
    v5->addTriangle(13);
    v5->addTriangle(13);
}

void VertexTests::tearDown()
{
    delete v1;
    delete v2;
    delete v3;
    delete v4;
    delete v5;
    //delete v6;
    delete vec1;
}

void VertexTests::testMethod()
{
    CPPUNIT_ASSERT(true);
}

void VertexTests::testFailedMethod()
{
    //CPPUNIT_ASSERT(true);
    CPPUNIT_ASSERT(false);
}

void VertexTests::testConstructor()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->r_c.x, constants::delta14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->r_c.y, constants::delta14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->r_c.z, constants::delta14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->f_c.x, constants::delta14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->f_c.y, constants::delta14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->f_c.z, constants::delta14);
    CPPUNIT_ASSERT_EQUAL(-1, v1->getId());
    //CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, v1->getMass(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2->r_c.x, constants::delta14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2->r_c.y, constants::delta14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, v2->r_c.z, constants::delta14);
    //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2->v_c.x, constants::delta14);
    //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2->v_c.y, constants::delta14);
    //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2->v_c.z, constants::delta14);
    //CPPUNIT_ASSERT_DOUBLES_EQUAL(1.23, v2->getMass(), DELTA14);
    CPPUNIT_ASSERT_EQUAL(-1, v2->getId());
    CPPUNIT_ASSERT_EQUAL(0, v2->numTris);
    CPPUNIT_ASSERT_EQUAL(0, v2->numBonded);
    CPPUNIT_ASSERT_EQUAL(0, v2->numNbNeighs);
    Vertex tmpv1(*vec1);
    CPPUNIT_ASSERT_EQUAL(-1, tmpv1.getId());
    tmpv1.setId(1);
    //tmpv1.setMass(12.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tmpv1.r_c.x, constants::delta14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, tmpv1.r_c.y, constants::delta14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, tmpv1.r_c.z, constants::delta14);
    //CPPUNIT_ASSERT_DOUBLES_EQUAL(12.0, tmpv1.getMass(), DELTA14);
    CPPUNIT_ASSERT_EQUAL(1, tmpv1.getId());
}

void VertexTests::testCopyConstructor()
{
    Vertex tmpv1(*v3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, tmpv1.r_c.x, constants::delta14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tmpv1.r_c.y, constants::delta14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, tmpv1.r_c.z, constants::delta14);
    //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpv1.v_c.x, constants::delta14);
    //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpv1.v_c.y, constants::delta14);
    //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpv1.v_c.z, constants::delta14);
    //CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tmpv1.getMass(), DELTA14);
    CPPUNIT_ASSERT_EQUAL(-1, tmpv1.getId());
    CPPUNIT_ASSERT_EQUAL(0, tmpv1.numTris);
    CPPUNIT_ASSERT_EQUAL(0, tmpv1.numBonded);
    CPPUNIT_ASSERT_EQUAL(0, tmpv1.numNbNeighs);
    Vertex tmpv2(*v4);
    CPPUNIT_ASSERT_EQUAL(18, tmpv2.getId());
    CPPUNIT_ASSERT_EQUAL(4, tmpv2.getNumNeighs());
    CPPUNIT_ASSERT_EQUAL(0, tmpv2.getNumTris());
    CPPUNIT_ASSERT_EQUAL(0, tmpv2.getNeighborId(0));
    CPPUNIT_ASSERT_EQUAL(1, tmpv2.getNeighborId(1));
    CPPUNIT_ASSERT_EQUAL(2, tmpv2.getNeighborId(2));
    CPPUNIT_ASSERT_EQUAL(5, tmpv2.getNeighborId(3));
    CPPUNIT_ASSERT_EQUAL(0.76, tmpv2.getNeighborR0(0));
    CPPUNIT_ASSERT_EQUAL(0.44, tmpv2.getNeighborR0(1));
    CPPUNIT_ASSERT_EQUAL(0.44, tmpv2.getNeighborR0(2));
    CPPUNIT_ASSERT_EQUAL(0.33, tmpv2.getNeighborR0(3));
    Vertex tmpv3(*v5);
    CPPUNIT_ASSERT_EQUAL(4, tmpv3.getNumTris());
    CPPUNIT_ASSERT_EQUAL(1, tmpv3.getTriangleId(0));
    CPPUNIT_ASSERT_EQUAL(0, tmpv3.getTriangleId(1));
    CPPUNIT_ASSERT_EQUAL(100, tmpv3.getTriangleId(2));
    CPPUNIT_ASSERT_EQUAL(13, tmpv3.getTriangleId(3));
}


void VertexTests::testAddNeighbor()
{
    CPPUNIT_ASSERT_EQUAL(18, v4->getId());
    CPPUNIT_ASSERT_EQUAL(4, v4->getNumNeighs());
    CPPUNIT_ASSERT_EQUAL(0, v4->getNumTris());
    CPPUNIT_ASSERT_EQUAL(0, v4->getNeighborId(0));
    CPPUNIT_ASSERT_EQUAL(1, v4->getNeighborId(1));
    CPPUNIT_ASSERT_EQUAL(2, v4->getNeighborId(2));
    CPPUNIT_ASSERT_EQUAL(5, v4->getNeighborId(3));
    CPPUNIT_ASSERT_EQUAL(0.76, v4->getNeighborR0(0));
    CPPUNIT_ASSERT_EQUAL(0.44, v4->getNeighborR0(1));
    CPPUNIT_ASSERT_EQUAL(0.44, v4->getNeighborR0(2));
    CPPUNIT_ASSERT_EQUAL(0.33, v4->getNeighborR0(3));
}

void VertexTests::testAddTriangle()
{
    CPPUNIT_ASSERT_EQUAL(4, v5->getNumTris());
    CPPUNIT_ASSERT_EQUAL(1, v5->getTriangleId(0));
    CPPUNIT_ASSERT_EQUAL(0, v5->getTriangleId(1));
    CPPUNIT_ASSERT_EQUAL(100, v5->getTriangleId(2));
    CPPUNIT_ASSERT_EQUAL(13, v5->getTriangleId(3));
}