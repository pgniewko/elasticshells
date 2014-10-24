#ifndef DATAEXCEPTION_H
#define	DATAEXCEPTION_H

#include <string>
#include <iostream>
#include <exception>
#include <stdexcept>

class DataException : public std::runtime_error
{
    public:
        DataException(const std::string& what_arg) : std::runtime_error( what_arg ) {}
};

#endif	/* DATAEXCEPTION_H */

