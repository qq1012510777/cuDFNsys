///////////////////////////////////////////////////////////////////
// NAME:              PairHash.cuh
//
// PURPOSE:           a struct of the compare (customed) of Hash Table
//
// FUNCTIONS/OBJECTS: PairHash
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"

struct PairHash
{
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};