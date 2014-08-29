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
}

void Vector3DTests::tearDown() 
{
}

void Vector3DTests::testMethod() 
{
    CPPUNIT_ASSERT(true);
}

void Vector3DTests::testFailedMethod() 
{
     CPPUNIT_ASSERT(false);
}

void Vector3DTests::testConstructor()
{
    Vector3D v3d1(0, 0, 0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v3d1.length(), DELTA);
    
    Vector3D v3d2(-0, -0, -0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v3d2.length(), DELTA);
    
    Vector3D v3d3(1, -1, -1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(SQRT3, v3d3.length(), DELTA);
    
    Vector3D v3d4(v3d3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3d3.length(), v3d4.length(), DELTA);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3d3.x, v3d4.x, DELTA);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3d3.y, v3d4.y, DELTA);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3d3.z, v3d4.z, DELTA);
    
    
    Vector3D v3d5(112, -123.000001, 0.987654321);
    Vector3D v3d6(v3d5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(112.0, v3d6.x, DELTA);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-123.000001, v3d6.y, DELTA);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.987654321, v3d6.z, DELTA);
}
