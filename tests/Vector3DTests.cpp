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
    Vector3D v3d();
     CPPUNIT_ASSERT_EQUAL(0, v3d.length() ); 
}

