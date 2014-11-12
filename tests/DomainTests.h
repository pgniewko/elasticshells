#ifndef DOMAINTESTS_H
#define	DOMAINTESTS_H

#include <cppunit/extensions/HelperMacros.h>

class DomainTests : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(DomainTests);
    
    CPPUNIT_TEST(testMethod);
    
    CPPUNIT_TEST_FAIL(testFailedMethod);
    
    CPPUNIT_TEST_SUITE_END();

public:
    DomainTests();
    virtual ~DomainTests();
    void setUp();
    void tearDown();

private:
    void testMethod();
    void testFailedMethod();
};

#endif	/* DOMAINTESTS_H */

