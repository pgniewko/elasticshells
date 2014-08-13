#ifndef GEOMETRYTESTS_H
#define	GEOMETRYTESTS_H

#include <cppunit/extensions/HelperMacros.h>

#include "../src/geometry/Vector3D.h"

class GeometryTests : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(GeometryTests);

    CPPUNIT_TEST(testMethod);
    CPPUNIT_TEST(testFailedMethod);

    CPPUNIT_TEST_SUITE_END();

public:
    GeometryTests();
    virtual ~GeometryTests();
    void setUp();
    void tearDown();

private:
    void testMethod();
    void testFailedMethod();
};

#endif	/* GEOMETRYTESTS_H */
