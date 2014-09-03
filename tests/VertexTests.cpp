#include "VertexTests.h"

CPPUNIT_TEST_SUITE_REGISTRATION(VertexTests);

VertexTests::VertexTests() {
}

VertexTests::~VertexTests() {
}

void VertexTests::setUp() {
    this->example = new int(1);
}

void VertexTests::tearDown() {
    delete this->example;
}

void VertexTests::testMethod() {
    CPPUNIT_ASSERT(*example == 1);
}

void VertexTests::testFailedMethod() {
    CPPUNIT_ASSERT(true);
    //CPPUNIT_ASSERT(++*example == 1);
}
