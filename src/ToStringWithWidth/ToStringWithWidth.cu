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