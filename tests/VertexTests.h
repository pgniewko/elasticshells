#ifndef VERTEXTESTS_H
#define	VERTEXTESTS_H

#include <cppunit/extensions/HelperMacros.h>

#include "Environment.h"
#include "geometry/Vertex.h"

class VertexTests : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(VertexTests);
    CPPUNIT_TEST(testMethod);
    
    CPPUNIT_TEST_FAIL(testFailedMethod);
    
    CPPUNIT_TEST_SUITE_END();

public:
    VertexTests();
    virtual ~VertexTests();
    void setUp();
    void tearDown();

private:
    int *example;
    void testMethod();
    void testFailedMethod();
};

#endif	/* VERTEXTESTS_H */

