#include "ForcesTests.h"

CPPUNIT_TEST_SUITE_REGISTRATION(ForcesTests);

ForcesTests::ForcesTests()
{
}

ForcesTests::~ForcesTests()
{
}

void ForcesTests::setUp()
{
    v1 = new Vector3D(0, 0, 0);
    v2 = new Vector3D(1, 0, 0);
    v3 = new Vector3D(0, 1, 0);
    v4 = new Vector3D(0, 0, 1);
}

void ForcesTests::tearDown()
{
    
    delete v1;
    delete v2;
    delete v3;
    delete v4;
}

void ForcesTests::testMethod()
{
    CPPUNIT_ASSERT(true);
}

void ForcesTests::testFailedMethod()
{
    CPPUNIT_ASSERT(false);
}

void ForcesTests::testHookeanMagnitude()
{
    
}

void ForcesTests::testHookeanSign()
{
    
}