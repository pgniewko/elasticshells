#include "Vector3DTests.h"


CPPUNIT_TEST_SUITE_REGISTRATION(Vector3DTests);

Vector3DTests::Vector3DTests() {}

Vector3DTests::~Vector3DTests() {}

void Vector3DTests::setUp()
{
    v1 = new Vector3D(0, 0, 0);
    v2 = new Vector3D(-0, -0, -0);
    v3 = new Vector3D(1, -1, -1);
    v4 = new Vector3D(112, -123.000001, 0.987654321);
    v5 = new Vector3D(constants::delta14, constants::delta14, constants::delta14);
    v6 = new Vector3D(1, 0, 0);
    v7 = new Vector3D(0, 1, 0);
    v8 = new Vector3D(-1, 0, 0);
}

void Vector3DTests::tearDown()
{
    delete v1;
    delete v2;
    delete v3;
    delete v4;
    delete v5;
    delete v6;
    delete v7;
    delete v8;
}

void Vector3DTests::testMethod()
{
    CPPUNIT_ASSERT(true);
}

void Vector3DTests::testFailedMethod()
{
    //CPPUNIT_ASSERT(true);
    CPPUNIT_ASSERT(false);
}

void Vector3DTests::testConstructor()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->length(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v2->length(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(SQRT3, v3->length(), DELTA14);
    Vector3D v3d1(*v3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3->length(), v3d1.length(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3->x, v3d1.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3->y, v3d1.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v3->z, v3d1.z, DELTA14);
    Vector3D v3d2(*v4);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(112.0, v3d2.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-123.000001, v3d2.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.987654321, v3d2.z, DELTA14);
    Vector3D v3d3 = *v5;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(DELTA14, v3d3.x, 0.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(DELTA14, v3d3.y, 0.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(DELTA14, v3d3.z, 0.0);
}

void Vector3DTests::testSetLength()
{
    Vector3D v3d1(*v1);
    v3d1.set_length(1.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v3d1.length(), DELTA14);
    Vector3D v3d2(*v5);
    v3d2.set_length(1.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, v3d2.length(), DELTA14);
    Vector3D v3d3(*v2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v3d3.length(), DELTA14);
}

void Vector3DTests::testAngle()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v5->angle(*v5), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v6->angle(*v6), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(PI / 2, v6->angle(*v7), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(PI / 2, v7->angle(*v8), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(PI, v6->angle(*v8), DELTA14);
}

void Vector3DTests::testExceptions()
{
    v1->angle(*v5);
}

void Vector3DTests::testOperators()
{
    Vector3D v3d1;
    v3d1 += *v1;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v3d1.length(), DELTA14);
    v3d1 -= *v3;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, v3d1.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, v3d1.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, v3d1.z, DELTA14);
    v3d1 *= -1.0;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, v3d1.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, v3d1.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0, v3d1.z, DELTA14);
    v3d1 /= 2.0;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, v3d1.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, v3d1.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, v3d1.z, DELTA14);
}