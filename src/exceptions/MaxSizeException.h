#ifndef MAXSIZEEXCEPTION_H
#define	MAXSIZEEXCEPTION_H

#include <iostream>
#include <exception>
using namespace std;

class MaxSizeException : public std::runtime_error
{
    public:
        MaxSizeException(const string& what_arg) : runtime_error( what_arg ) {}
};

#endif	/* MAXSIZEEXCEPTION_H */