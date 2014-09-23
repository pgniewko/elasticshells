#ifndef DATAEXCEPTION_H
#define	DATAEXCEPTION_H

#include <string>
#include <stdexcept>
using namespace std;

class DataException : public std::runtime_error
{
    public:
        DataException(const string& what_arg) : runtime_error( what_arg ) {}
};

#endif	/* DATAEXCEPTION_H */

