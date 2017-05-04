#ifndef NOTALLOWEDEXCEPTION_H
#define NOTALLOWEDEXCEPTION_H

#include <string>
#include <iostream>
#include <exception>
#include <stdexcept>

class NotAllowedException : public std::runtime_error
{
    public:
        NotAllowedException(const std::string& what_arg) : std::runtime_error( what_arg ) {}

};

#endif /* NOTALLOWEDEXCEPTION_H */

