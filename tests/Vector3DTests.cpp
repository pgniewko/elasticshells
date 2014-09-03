#include "Vector3DTests.h"


CPPUNIT_TEST_SUITE_REGISTRATION(Vector3DTests);

Vector3DTests::Vector3DTests() 
{
}

Vector3DTests::~Vector3DTests() 
{
}

void Vector3DTests::setUp() 
{
    v1 = new Vector3D(0, 0, 0);
    v2 = new Vector3D(-0, -0, -0);
    v3 = new Vector3D(1, -1, -1);
    v4 = new Vector3D(112, -123.000001, 0.987654321);

    
}

void Vector3DTests::tearDown() 
{
    delete v1;
    delete v2;
    delete v3;
    delete v4;
    delete v5;
    delete v6;
}

void Vector3DTests::testMethod() 
{
    CPPUNIT_ASSERT(true);
}

void Vector3DTests::testFailedMethod() 
{
    CPPUNIT_ASSERT(true);
    //CPPUNIT_ASSERT(false);
}

void Vector3DTests::testConstructor()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->length(), DELTA);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2->length(), DELTA);
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(SQRT3, v3->length(), DELTA);
    
    Vector3D v3d1(*v3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3->length(), v3d1.length(), DELTA);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3->x, v3d1.x, DELTA);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3->y, v3d1.y, DELTA);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3->z, v3d1.z, DELTA);
    
    
    Vector3D v3d2(*v4);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(112.0, v3d2.x, DELTA);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-123.000001, v3d2.y, DELTA);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.987654321, v3d2.z, DELTA);
}
