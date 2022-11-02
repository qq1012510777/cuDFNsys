/****************************************************************************
* cuDFNsys - simulating flow and transport in 3D fracture networks          *
* Copyright (C) 2022, Tingchang YIN, Sergio GALINDO-TORRES                  *
*                                                                           *
* This program is free software: you can redistribute it and/or modify      *
* it under the terms of the GNU Affero General Public License as            *
* published by the Free Software Foundation, either version 3 of the        *
* License, or (at your option) any later version.                           *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU Affero General Public License for more details.                       *
*                                                                           *
* You should have received a copy of the GNU Affero General Public License  *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.    *
*****************************************************************************/

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
