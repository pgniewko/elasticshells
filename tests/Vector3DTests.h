#ifndef VECTOR3DTESTS_H
#define	VECTOR3DTESTS_H

#include <cppunit/extensions/HelperMacros.h>

class Vector3DTests : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(Vector3DTests);

    CPPUNIT_TEST(testMethod);
    CPPUNIT_TEST(testFailedMethod);

    CPPUNIT_TEST_SUITE_END();

public:
    Vector3DTests();
    virtual ~Vector3DTests();
    void setUp();
    void tearDown();

private:
    void testMethod();
    void testFailedMethod();
};

#endif	/* VECTOR3DTESTS)H */