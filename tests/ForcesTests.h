#ifndef FORCESTESTS_H
#define	FORCESTESTS_H

#include <cppunit/extensions/HelperMacros.h>

#include "Environment.h"
#include "force/HertzianRepulsion.h"
#include "force/HookeanForce.h"
#include "force/OsmoticForce.h"
#include "geometry/Vector3D.h"
#include "geometry/Vertex.h"
#include "geometry/Tetrahedron.h"

class ForcesTests : public CPPUNIT_NS::TestFixture
{
        CPPUNIT_TEST_SUITE(ForcesTests);

        CPPUNIT_TEST(testMethod);
        CPPUNIT_TEST(testHookeanForce);
        CPPUNIT_TEST(testHertzianForce);
        CPPUNIT_TEST(testOsmoticForce);
        CPPUNIT_TEST(testNbRepulsiveForce);


        CPPUNIT_TEST_FAIL(testFailedMethod);

        CPPUNIT_TEST_SUITE_END();

    public:
        ForcesTests();
        virtual ~ForcesTests();
        void setUp();
        void tearDown();

    private:
        Vector3D* v1, *v2, *v3, *v4;
        void testMethod();
        void testFailedMethod();

        void testHookeanForce();
        void testHertzianForce();
        void testOsmoticForce();
        void testNbRepulsiveForce();
};

#endif	/* FORCESTESTS_H */