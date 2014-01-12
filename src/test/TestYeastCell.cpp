#include "TestYeastCell.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TestYeastCell);

TestYeastCell::TestYeastCell()
{
}

TestYeastCell::~TestYeastCell()
{
}

void TestYeastCell::setUp()
{
    double p1 [1] = {1.0};
    double p2 [1] = {2.0};
    yc1 = new YeastCell();
    yc1->set_params(p1, 1);
    yc2 = new YeastCell();
    yc2->set_params(p2, 1);
    yc3 = new YeastCell(1.5);
}

void TestYeastCell::tearDown()
{
    delete yc1;
    delete yc2;
    delete yc3;
}

void TestYeastCell::testVolume()
{
    CPPUNIT_ASSERT(yc1->calc_volume() == 1.0);
    CPPUNIT_ASSERT(yc2->calc_volume() == 8.0);
    CPPUNIT_ASSERT(yc3->calc_volume() == 3.375);
}

void TestYeastCell::testMethod()
{
    CPPUNIT_ASSERT(true);
}

void TestYeastCell::testFailedMethod()
{
    CPPUNIT_ASSERT(true);
}

