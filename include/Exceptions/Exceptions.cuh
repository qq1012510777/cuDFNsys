///////////////////////////////////////////////////////////////////
// NAME:              Exceptions.h
//
// PURPOSE:           an exception and two inherited classes
//
// FUNCTIONS/OBJECTS: Fractures
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include <iostream>
#include <string>
using namespace std;

namespace cuDFNsys
{
class Exceptions
{
protected:
    // exception: message
    string msg;

public:
    // create
    Exceptions(const string &s);
    // output
    string what();
};
}; // namespace cuDFNsys

// the following exception class can be thrown e.g. when a resolution problem happens in Monte Carlo simulations.
namespace cuDFNsys
{
class ExceptionsIgnore : public cuDFNsys::Exceptions
{
public:
    ExceptionsIgnore(const string &s) : cuDFNsys::Exceptions("\033[1;33m" + s + "\033[0m"){};
};
}; // namespace cuDFNsys

// the following exception class can be thrown when the program should be discontinued immediately.
namespace cuDFNsys
{
class ExceptionsPause : public cuDFNsys::Exceptions
{
public:
    ExceptionsPause(const string &s) : cuDFNsys::Exceptions("\033[1;31m" + s + "\033[0m"){};
};
}; // namespace cuDFNsys
