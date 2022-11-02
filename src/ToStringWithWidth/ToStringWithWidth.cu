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

#include "ToStringWithWidth/ToStringWithWidth.cuh"

// ====================================================
// NAME:        ToStringWithWidth
// DESCRIPTION: convert number to string with fixed length
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
template <class T>
string cuDFNsys::ToStringWithWidth(const T &val, const size_t &width)
{
    std::ostringstream oss;
    oss.width(width);
    oss.fill('0');
    oss << val;
    return oss.str();
}; // ToStringWithWidth
template string cuDFNsys::ToStringWithWidth(const int &val, const size_t &width);
template string cuDFNsys::ToStringWithWidth(const double &val, const size_t &width);
template string cuDFNsys::ToStringWithWidth(const float &val, const size_t &width);
template string cuDFNsys::ToStringWithWidth(const uint &val, const size_t &width);
template string cuDFNsys::ToStringWithWidth(const size_t &val, const size_t &width);