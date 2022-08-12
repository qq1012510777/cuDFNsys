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
                        const string &fieldname)
    {
        const char *filename = FileKey_mat.c_str();
        MATFile *pMatFile;
        const char *mode_ = mode.c_str();
        pMatFile = matOpen(filename, mode_);

        if (!pMatFile)
            throw cuDFNsys::ExceptionsPause("cannot create mat file\n");

        mxArray *pa3;
        pa3 = mxCreateString(content.c_str());

        /*int status = */
        matPutVariable(pMatFile, fieldname.c_str(), pa3);

        // if (status == 0)
        // {
        //     string AS = "Error using matPutVariable\n";
        //     throw cuDFNsys::ExceptionsPause(AS);
        // }
        mxDestroyArray(pa3);

        matClose(pMatFile);
    };

    // read string
    string ReadMatString(const string &FileKey_mat, const string &fieldname)
    {
        const char *filename = FileKey_mat.c_str();
        MATFile *pMatFile;
        pMatFile = matOpen(filename, "r");

        if (!pMatFile)
            throw cuDFNsys::ExceptionsPause("cannot read mat file\n");

        mxArray *pa3;
        pa3 = matGetVariable(pMatFile, fieldname.c_str());

        //char str[BUFSIZE];
        //mxGetString(pa3, str, sizeof(str));

        size_t NUM_chars = mxGetNumberOfElements(pa3);
        char *strtt = new char[NUM_chars];

        mxGetString(pa3, strtt, sizeof(char) * (NUM_chars + 1));

        string si = strtt;

        delete[] strtt;
        strtt = NULL;
        mxDestroyArray(pa3);
        matClose(pMatFile);

        return si;
    };

    // read data
    void ReadMat(const string &FileKey_mat,
                 const string &field_name,
                 Eigen::MatrixXd &data)
    {
        const char *filename = FileKey_mat.c_str();
        MATFile *pMatFile;
        pMatFile = matOpen(filename, "r");
        if (!pMatFile)
            throw cuDFNsys::ExceptionsPause("cannot read mat file\n");

        mxArray *pa3;
        pa3 = matGetVariable(pMatFile, field_name.c_str());

        size_t rows = mxGetM(pa3); // rows

        size_t cols = mxGetN(pa3); // cols

        double *buffer = (double *)mxGetData(pa3);

        data = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(buffer, rows, cols);

        //mxFree(buffer);
        mxDestroyArray(pa3);
        matClose(pMatFile);
    };
};
}; // namespace cuDFNsys
