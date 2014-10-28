#ifndef RUNTIMEERROR_H
#define	RUNTIMEERROR_H

#include <string>
#include <iostream>
#include <exception>
#include <stdexcept>

class RunTimeError : public std::runtime_error
{
    public:
        RunTimeError(const std::string& what_arg) : std::runtime_error( what_arg ) {}
};


#endif	/* RUNTIMEERROR_H */

