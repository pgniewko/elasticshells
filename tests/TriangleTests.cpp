#include "TriangleTests.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TriangleTests);

TriangleTests::TriangleTests() {}

TriangleTests::~TriangleTests() {}

void TriangleTests::setUp() 
{
    v1 = new Vector3D(1,0,0);
    v2 = new Vector3D(0,1,0);
    v3 = new Vector3D(0,0,1);
    v4 = new Vector3D(0,0,0);
    v5 = new Vector3D(-1,0,0);
    v6 = new Vector3D(1,1,1);
    
    t1 = new Triangle(*v1, *v2, *v3);
    t2 = new Triangle(*v1, *v2, *v5);
    t3 = new Triangle(*v4, *v4, *v4);
    
    vt1 = new VertexTriangle;
    vt2 = new VertexTriangle(0,1,2);
    vt3 = new VertexTriangle(0,1,3);
    vt3->setId(2);
    

}

void TriangleTests::tearDown() 
{
    delete v1;
    delete v2;
    delete v3;
    delete v4; 
    delete v5;
    delete v6;
    
    delete t1;
    delete t2;
    delete t3;
    
    delete vt1;
    delete vt2;
    delete vt3;
    
    delete vx1;
}

void TriangleTests::testMethod() {
    CPPUNIT_ASSERT(true);
}

void TriangleTests::testFailedMethod() {
    CPPUNIT_ASSERT(false);
}

void TriangleTests::testConstructors()
{
    
}

void TriangleTests::testTriangleArea()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5*SQRT3, t1->area(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, t2->area(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, t3->area(), DELTA14);
}

void TriangleTests::testTriangleCopy()
{
    Triangle tmpt1(*t1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5*SQRT3, tmpt1.area(), DELTA14);
    
    Triangle tmpt2(*t2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tmpt2.area(), DELTA14);
    
    Triangle tmpt3(*t3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpt3.area(), DELTA14);
}

void TriangleTests::testVConstructors()
{
    CPPUNIT_ASSERT_EQUAL(-1, vt1->ia);
    CPPUNIT_ASSERT_EQUAL(-1, vt1->ib);
    CPPUNIT_ASSERT_EQUAL(-1, vt1->ic);
    CPPUNIT_ASSERT_EQUAL(-1, vt1->getId());
    
    CPPUNIT_ASSERT_EQUAL(0, vt2->ia);
    CPPUNIT_ASSERT_EQUAL(1, vt2->ib);
    CPPUNIT_ASSERT_EQUAL(2, vt2->ic);
    CPPUNIT_ASSERT_EQUAL(-1, vt2->getId());
    
    CPPUNIT_ASSERT_EQUAL(0, vt3->ia);
    CPPUNIT_ASSERT_EQUAL(1, vt3->ib);
    CPPUNIT_ASSERT_EQUAL(3, vt3->ic);
    CPPUNIT_ASSERT_EQUAL(2, vt3->getId());
}

void TriangleTests::testVTriangleCopy()
{
    VertexTriangle tmpvt1(*vt3);
    CPPUNIT_ASSERT_EQUAL(0, tmpvt1.ia);
    CPPUNIT_ASSERT_EQUAL(1, tmpvt1.ib);
    CPPUNIT_ASSERT_EQUAL(3, tmpvt1.ic);
    CPPUNIT_ASSERT_EQUAL(2, tmpvt1.getId());
}

void TriangleTests::testVTriangleArea()
{
    Vertex vxarray[6] = {*v1, *v2, *v3, *v4, *v5, *v6};
//    cout << vxarray[0].xyz << endl;
//    cout << vxarray[1].xyz << endl;
//    cout << vxarray[2].xyz << endl;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, vt1->area(vxarray), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5*SQRT3, vt2->area(vxarray), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, vt3->area(vxarray), DELTA14);
    //vxarray[0] = 
}