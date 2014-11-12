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
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box1->getX(), 0.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box2->getX(), 5.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box3->getY(), 10.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box4->getY(), 6.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box5->getZ(), 4.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box6->getZ(), 6.0, DELTA14);
}

void BoxTests::testCopyConstructor()
{
    Box tmpbox1(*box1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox1.getX(), 0.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox1.getY(), 0.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox1.getZ(), 0.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox1.getVolume(), 0.0, DELTA14);
    Box tmpbox2(*box2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox2.getX(), 5.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox2.getY(), 5.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox2.getZ(), 5.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox2.getDx(), 0.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox2.getDy(), 0.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox2.getDz(), 0.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox2.getVolume(), 1000.0, DELTA14);
    Box tmpbox3(*box6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox3.getX(), 4.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox3.getY(), 5.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox3.getZ(), 6.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox3.getDx(), 0.1, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox3.getDy(), 0.1, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox3.getDz(), 0.1, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox3.getVolume(), 960.0, DELTA14);
    Box tmpbox4(*box7);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox4.getX(), 4.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox4.getY(), 5.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox4.getZ(), 6.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox4.getDx(), -0.1, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox4.getDy(), -0.1, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox4.getDz(), -0.1, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(tmpbox4.getVolume(), 960.0, DELTA14);
}

void BoxTests::testGetVolume()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box1->getVolume(), 0.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box2->getVolume(), 1000.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box3->getVolume(), 8000.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box4->getVolume(), 1680.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box5->getVolume(), 512.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box6->getVolume(), 960.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box7->getVolume(), 960.0, DELTA14);
}

void BoxTests::testGetArea()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box1->getArea(), 0.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box2->getArea(), 600, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box3->getArea(), 2400.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box4->getArea(), 856.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box5->getArea(), 384.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box6->getArea(), 592.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box7->getArea(), 592.0, DELTA14);
}

void BoxTests::testResizing()
{
    for (int i = 0; i < 10; i++)
    {
        box1->resize();
    }

    CPPUNIT_ASSERT_DOUBLES_EQUAL(box1->getX(), 0.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box1->getY(), 0.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box1->getZ(), 0.0, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(box1->getVolume(), 0.0, DELTA14);
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
    tmpbox3.setXend(1.0);
    tmpbox3.setYend(1.0);
    tmpbox3.setZend(1.0);

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