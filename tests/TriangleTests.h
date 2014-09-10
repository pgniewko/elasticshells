#ifndef TRIANGLETESTS_H
#define	TRIANGLETESTS_H

#include <cppunit/extensions/HelperMacros.h>

#include "geometry/Triangle.h"

class TriangleTests : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(TriangleTests);
    CPPUNIT_TEST(testMethod);
    
    CPPUNIT_TEST_FAIL(testFailedMethod);
    CPPUNIT_TEST_SUITE_END();

public:
    TriangleTests();
    virtual ~TriangleTests();
    void setUp();
    void tearDown();

private:
    int *example;
    void testMethod();
    void testFailedMethod();
};

#endif	/* TRIANGLETESTS_H */

