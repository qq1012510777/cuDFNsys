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
// NAME:              Graph.cuh
//
// PURPOSE:           Deep first search
//                    find all clusters
//
// FUNCTIONS/OBJECTS: Graph
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include <bits/stdc++.h>
#include <cmath>
#include <ctime>
#include <iostream>
#include <queue>
#include <vector>
using namespace std;

namespace cuDFNsys
{
template <typename T>
class Graph
{
public:
    // The number of fractures
    size_t V;

    // Adjacent list
    list<int> *Adj;

public:
    // constructor
    Graph(const size_t &NumOfFractures_1,
          std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> Intersection_map);
    // implement DFS
    void UseDFS(std::vector<std::vector<size_t>> &S);
    // destructor
    ~Graph()
    {
        delete[] this->Adj;
        this->Adj = NULL;
    };

private:
    void DFS(std::vector<std::vector<size_t>> &ListOfClusters);
    void DFSUtil(int s, vector<bool> &visited, vector<size_t> &onecluster);
    void addEdge(int v, int w);
};
}; // namespace cuDFNsys