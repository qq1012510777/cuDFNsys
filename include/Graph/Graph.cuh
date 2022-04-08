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
#include <bits/stdc++.h>
#include <cmath>
#include <ctime>
#include <iostream>
#include <queue>
#include <vector>
using namespace std;

namespace cuDFNsys
{
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
          std::map<pair<size_t, size_t>, pair<float3, float3>> Intersection_map);
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

// ====================================================
// NAME:        UseDFS
// DESCRIPTION: UseDFS
//
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
inline void Graph::UseDFS(std::vector<std::vector<size_t>> &S)
{
    S.clear();
    DFS(S);
    S.shrink_to_fit();
}; // UseDFS

// ====================================================
// NAME:        addEdge
// DESCRIPTION: create adjacent list
//
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
inline void Graph::addEdge(int v, int w)
{
    // Add w to v’s list.
    Adj[v].push_back(w);
}; // addEdge
}; // namespace cuDFNsys