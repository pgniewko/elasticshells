#include "VertexTests.h"

CPPUNIT_TEST_SUITE_REGISTRATION(VertexTests);

VertexTests::VertexTests() {}

VertexTests::~VertexTests() {}

void VertexTests::setUp() 
{
    v1 = new Vertex;
        
}

void VertexTests::tearDown() 
{
    delete v1;
}

void VertexTests::testMethod() 
{
    
}

void VertexTests::testFailedMethod() 
{
    //CPPUNIT_ASSERT(true);
    CPPUNIT_ASSERT(false);
}

void VertexTests::testConstructor()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->xyz.x, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->xyz.y, DELTA14);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, v1->xyz.z, DELTA14);
    CPPUNIT_ASSERT_EQUAL(0, v1->getId());
}
