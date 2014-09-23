#ifndef TRIANGLETESTS_H
#define	TRIANGLETESTS_H

#include <cppunit/extensions/HelperMacros.h>

#include "Environment.h"
#include "geometry/Vector3D.h"
#include "geometry/Triangle.h"
#include "geometry/Vertex.h"
#include "geometry/VertexTriangle.h"

class TriangleTests : public CPPUNIT_NS::TestFixture
{
        CPPUNIT_TEST_SUITE(TriangleTests);
        CPPUNIT_TEST(testMethod);

        CPPUNIT_TEST(testConstructors);
        CPPUNIT_TEST(testTriangleCopy);
        CPPUNIT_TEST(testTriangleArea);

        CPPUNIT_TEST(testVConstructors);
        CPPUNIT_TEST(testVTriangleCopy);
        CPPUNIT_TEST(testVTriangleArea);

        CPPUNIT_TEST_FAIL(testFailedMethod);
        CPPUNIT_TEST_SUITE_END();

    public:
        TriangleTests();
        virtual ~TriangleTests();
        void setUp();
        void tearDown();

        void testConstructors();
        void testTriangleCopy();
        void testTriangleArea();

        void testVConstructors();
        void testVTriangleCopy();
        void testVTriangleArea();

    private:
        Vector3D* v1, *v2, *v3, *v4, *v5, *v6;
        Triangle* t1, *t2, *t3;
        VertexTriangle* vt1, *vt2, *vt3;
        Vertex* vx1;
        void testMethod();
        void testFailedMethod();
};

#endif	/* TRIANGLETESTS_H */

