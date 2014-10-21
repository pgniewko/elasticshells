#ifndef BOXTESTS_H
#define	BOXTESTS_H

#include <cppunit/extensions/HelperMacros.h>

#include "Environment.h"
#include "simulation/Box.h"

class BoxTests : public CPPUNIT_NS::TestFixture
{
        CPPUNIT_TEST_SUITE(BoxTests);

        CPPUNIT_TEST(testMethod);
        CPPUNIT_TEST(testConstructor);
        CPPUNIT_TEST(testCopyConstructor);
        CPPUNIT_TEST(testGetVolume);
        CPPUNIT_TEST(testResizing);

        CPPUNIT_TEST_FAIL(testFailedMethod);

        CPPUNIT_TEST_SUITE_END();

    public:
        BoxTests();
        virtual ~BoxTests();
        void setUp();
        void tearDown();

        void testConstructor();
        void testCopyConstructor();
        void testGetVolume();
        void testResizing();

    private:
        Box* box1, *box2, *box3, *box4, *box5, *box6, *box7;
        void testMethod();
        void testFailedMethod();
};

#endif	/* BOXTESTS_H */

