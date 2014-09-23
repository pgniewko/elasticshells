#ifndef NOTIMPLEMENTEDEXCEPTION_H
#define	NOTIMPLEMENTEDEXCEPTION_H

#include <iostream>
#include <exception>
using namespace std;

class NotImplementedException : public std::runtime_error
{
    public:
        NotImplementedException(const string& what_arg) : runtime_error( what_arg ) {}
};

#endif	/* NOTIMPLEMENTEDEXCEPTION_H */

