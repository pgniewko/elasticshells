#include "BoxTests.h"

CPPUNIT_TEST_SUITE_REGISTRATION(BoxTests);

BoxTests::BoxTests()
{
}

BoxTests::~BoxTests()
{
}

void BoxTests::setUp()
{
    box1 = new Box(0, 0, 0);
    box2 = new Box(5, 5, 5);
    box3 = new Box(10, 10, 10);
    box4 = new Box(5, 6, 7);
    box5 = new Box(4, 4, 4, 0.1);
    box6 = new Box(4, 5, 6, 0.1);
    box7 = new Box(4, 5, 6, -0.1);
}

void BoxTests::tearDown()
{
    delete box1;
    delete box2;
    delete box3;
    delete box4;
    delete box5;
    delete box6;
    delete box7;
}

void BoxTests::testMethod()
{
    CPPUNIT_ASSERT(true);
}

void BoxTests::testFailedMethod()
{
    CPPUNIT_ASSERT(false);
}

void BoxTests::testConstructor()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, box1->getX(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, box2->getX(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0, box3->getY(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, box4->getY(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, box5->getZ(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, box6->getZ(), DELTA14);
}

void BoxTests::testCopyConstructor()
{
    Box tmpbox1(*box1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpbox1.getX(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpbox1.getY(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpbox1.getZ(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpbox1.getVolume(), DELTA14);
    Box tmpbox2(*box2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, tmpbox2.getX(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, tmpbox2.getY(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, tmpbox2.getZ(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpbox2.getDx(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpbox2.getDy(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpbox2.getDz(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1000.0, tmpbox2.getVolume(), DELTA14);
    Box tmpbox3(*box6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, tmpbox3.getX(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, tmpbox3.getY(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, tmpbox3.getZ(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1, tmpbox3.getDx(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1, tmpbox3.getDy(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1, tmpbox3.getDz(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(960.0, tmpbox3.getVolume(), DELTA14);
    Box tmpbox4(*box7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, tmpbox4.getX(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, tmpbox4.getY(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0, tmpbox4.getZ(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.1, tmpbox4.getDx(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.1, tmpbox4.getDy(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.1, tmpbox4.getDz(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(960.0, tmpbox4.getVolume(), DELTA14);
}

void BoxTests::testGetVolume()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, box1->getVolume(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1000.0, box2->getVolume(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(8000.0, box3->getVolume(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1680.0, box4->getVolume(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(512.0, box5->getVolume(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(960.0, box6->getVolume(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(960.0, box7->getVolume(), DELTA14);
}

void BoxTests::testGetArea()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, box1->getArea(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(600.0, box2->getArea(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2400.0, box3->getArea(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(856.0, box4->getArea(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(384.0, box5->getArea(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(592.0, box6->getArea(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(592.0, box7->getArea(), DELTA14);
}

void BoxTests::testResizing()
{
    for (int i = 0; i < 10; i++)
    {
        box1->resize();
    }

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, box1->getX(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, box1->getY(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, box1->getZ(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, box1->getVolume(), DELTA14);
    Box tmpbox1(*box1);
    tmpbox1.setDx(0.1);

    for (int i = 0; i < 10; i++)
    {
        tmpbox1.resize();
    }

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tmpbox1.getX(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpbox1.getY(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpbox1.getZ(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmpbox1.getVolume(), DELTA14);
    Box tmpbox2(*box1);
    tmpbox2.setDx(0.1);
    tmpbox2.setDy(0.2);
    tmpbox2.setDz(0.3);

    for (int i = 0; i < 10; i++)
    {
        tmpbox2.resize();
    }

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, tmpbox2.getX(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, tmpbox2.getY(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, tmpbox2.getZ(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(48.0, tmpbox2.getVolume(), DELTA7);
    Box tmpbox3(4, 5, 6, -0.1);
    tmpbox3.setXmin(1.0);
    tmpbox3.setYmin(1.0);
    tmpbox3.setZmin(1.0);

    for (int i = 0; i < 10; i++)
    {
        tmpbox3.resize();
    }

    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, tmpbox3.getX(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, tmpbox3.getY(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, tmpbox3.getZ(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(480, tmpbox3.getVolume(), DELTA7);

    for (int i = 0; i < 10; i++)
    {
        tmpbox3.resize();
    }

    CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, tmpbox3.getX(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, tmpbox3.getY(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0, tmpbox3.getZ(), DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(192.0, tmpbox3.getVolume(), DELTA7);
}