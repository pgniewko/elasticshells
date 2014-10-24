#ifndef MAXSIZEEXCEPTION_H
#define	MAXSIZEEXCEPTION_H

#include <string>
#include <iostream>
#include <exception>
#include <stdexcept>

class MaxSizeException : public std::runtime_error
{
    public:
        MaxSizeException(const std::string& what_arg) : std::runtime_error( what_arg ) {}
};

#endif	/* MAXSIZEEXCEPTION_H */