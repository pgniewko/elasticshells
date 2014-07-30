#ifndef NEWTESTCLASS_H
#define	NEWTESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include "YeastCell.h"

class TestYeastCell : public CPPUNIT_NS::TestFixture
{
        CPPUNIT_TEST_SUITE(TestYeastCell);

        CPPUNIT_TEST(testMethod);
        CPPUNIT_TEST(testVolume);
        CPPUNIT_TEST(testFailedMethod);

        CPPUNIT_TEST_SUITE_END();

    public:
        TestYeastCell();
        virtual ~TestYeastCell();
        void setUp();
        void tearDown();

    private:
        YeastCell* yc1, *yc2, *yc3;
        void testMethod();
        void testFailedMethod();
        void testVolume();
};

#endif	/* NEWTESTCLASS_H */