#ifndef DOMAINLISTTESTS_H
#define	DOMAINLISTTESTS_H

#include <cppunit/extensions/HelperMacros.h>

class DomainListTests : public CPPUNIT_NS::TestFixture
{
        CPPUNIT_TEST_SUITE(DomainListTests);

        CPPUNIT_TEST(testMethod);

        CPPUNIT_TEST_FAIL(testFailedMethod);

        CPPUNIT_TEST_SUITE_END();

    public:
        DomainListTests();
        virtual ~DomainListTests();
        void setUp();
        void tearDown();

    private:
        void testMethod();
        void testFailedMethod();
};

#endif	/* DOMAINLISTTESTS_H */

