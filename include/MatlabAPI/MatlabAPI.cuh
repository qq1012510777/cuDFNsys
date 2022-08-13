///////////////////////////////////////////////////////////////////
// NAME:              Exceptions.h
//
// PURPOSE:           various Exceptions
//
// FUNCTIONS/OBJECTS: Fractures
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Exceptions/Exceptions.cuh"
#include "mat.h"
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace cuDFNsys
{
class MatlabAPI
{
public:
    // constructor
    MatlabAPI();
    // write data
    template <class T>
    void WriteMat(const string &FileKey_mat,
                  const string &mode,
                  const size_t &NUM_eles,
                  const size_t &rows,
                  const size_t &cols,
                  T *data_, /*column major data*/
                  const string &field_name);

    // write string
    void WriteMatString(const string &FileKey_mat,
                        const string &mode,
                        const string &content,
                        const string &fieldname);

    // read string
    string ReadMatString(const string &FileKey_mat,
                         const string &fieldname);

    // read data
    void ReadMat(const string &FileKey_mat,
                 const string &field_name,
                 Eigen::MatrixXd &data);
};
}; // namespace cuDFNsys
