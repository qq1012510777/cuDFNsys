///////////////////////////////////////////////////////////////////
// NAME:              UMapEdge.cuh
//
// PURPOSE:           a typedef of map
//
// FUNCTIONS/OBJECTS: UMapEdge
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "PairHash.cuh"

typedef std::unordered_map<pair<size_t, size_t>, int, PairHash> UMapEdge;
