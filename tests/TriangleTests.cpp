#include "TriangleTests.h"

CPPUNIT_TEST_SUITE_REGISTRATION(TriangleTests);

TriangleTests::TriangleTests() {
}

TriangleTests::~TriangleTests() {
}

void TriangleTests::setUp() {
    this->example = new int(1);
}

void TriangleTests::tearDown() {
    delete this->example;
}

void TriangleTests::testMethod() {
    CPPUNIT_ASSERT(*example == 1);
}

void TriangleTests::testFailedMethod() {
    CPPUNIT_ASSERT(++*example == 1);
}
