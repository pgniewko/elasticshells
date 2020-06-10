#ifndef NOTIMPLEMENTEDEXCEPTION_H
#define	NOTIMPLEMENTEDEXCEPTION_H

#include <string>
#include <iostream>
#include <exception>
#include <stdexcept>

class NotImplementedException : public std::runtime_error
{
    public:
        NotImplementedException(const std::string& what_arg) : std::runtime_error( what_arg ) {}
};

#endif	/* NOTIMPLEMENTEDEXCEPTION_H */

