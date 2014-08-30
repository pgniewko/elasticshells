#ifndef MAXSIZEEXCEPTION_H
#define	MAXSIZEEXCEPTION_H

#include <iostream>
#include <exception>
using namespace std;

struct MaxSizeException : public exception
{
  const char * what () const throw ()
  {
    return "MAX number of elements reached.\nNo action will be taken.";
  }
};

#endif	/* MAXSIZEEXCEPTION_H */

