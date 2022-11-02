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