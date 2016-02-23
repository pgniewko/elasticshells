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

void ForcesTests::testHookeanForce()
{
    Vector3D f;
    f = HookeanForce::calcForce(*v1, *v2, 1.0, 1.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, f.length(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, f.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, f.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, f.z, DELTA14);
    f = HookeanForce::calcForce(*v1, *v2, 0.5, 1.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, f.length(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, f.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, f.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, f.z, DELTA14);
    f = HookeanForce::calcForce(*v2, *v1, 0.5, 1.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, f.length(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, f.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, f.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, f.z, DELTA14);
    f = HookeanForce::calcForce(*v2, *v3, constants::sqrt2, 1.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, f.length(), DELTA14);
    f = HookeanForce::calcForce(*v2, *v3, constants::sqrt2 / 2, 1.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(SQRT2 / 2, f.length(), DELTA14);
}

void ForcesTests::testHertzianForce()
{
    Vector3D f;
    f = HertzianRepulsion::calcForce(*v2 - *v1, 1.0, 100.0, 1.0, 1.0, 0.5, 0.5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, f.length(), DELTA14);
    f = HertzianRepulsion::calcForce(0.5 * (*v2 - *v1), 1.0 , 100.0, 1.0, 1.0, 0.5, 0.5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(47.140452079103173, f.length(), DELTA7);
}

void ForcesTests::testOsmoticForce()
{
}

void ForcesTests::testNbRepulsiveForce()
{
}