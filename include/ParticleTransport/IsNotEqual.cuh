#pragma once
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/remove.h>
#include <thrust/transform_reduce.h>

template <typename T>
struct IsNotEqual
{
    T value;
    __host__ __device__ IsNotEqual(T _value) : value(_value) {}

    __host__ __device__ bool operator()(const T &x) const { return !(x != value); }
};

