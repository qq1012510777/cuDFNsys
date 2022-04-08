#include "Exceptions/Exceptions.cuh"

// ====================================================
// NAME:        Exceptions
// DESCRIPTION: create an exception.
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
cuDFNsys::Exceptions::Exceptions(const string &s)
{
    this->msg = s;
}; // Exceptions

// ====================================================
// NAME:        Exceptions::what()
// DESCRIPTION: get an exception.
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
string cuDFNsys::Exceptions::what()
{
    return this->msg;
}; // what