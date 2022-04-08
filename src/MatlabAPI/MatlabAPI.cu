#include "MatlabAPI/MatlabAPI.cuh"

// ====================================================
// NAME:        MatlabAPI
// DESCRIPTION: create a MatlabAPI class.
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
cuDFNsys::MatlabAPI::MatlabAPI(){
    //
}; // MatlabAPI

// ====================================================
// NAME:        MatlabAPI::WriteMat
// DESCRIPTION: write data.
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
template <class T>
void cuDFNsys::MatlabAPI::WriteMat(const string &FileKey_mat,
                                   const string &mode,
                                   const size_t &NUM_eles,
                                   const size_t &rows,
                                   const size_t &cols,
                                   T *data_, /*column major data*/
                                   const string &field_name)
{
    const char *filename = FileKey_mat.c_str();
    MATFile *pMatFile;
    const char *mode_ = mode.c_str();
    pMatFile = matOpen(filename, mode_);

    if (!pMatFile)
        throw cuDFNsys::ExceptionsPause("cannot create mat file\n");

    double pData1[NUM_eles] = {};

    mxArray *pMxArray1;
    pMxArray1 = mxCreateDoubleMatrix(rows, cols, mxREAL);

    if (!pMxArray1 || !pData1)
    {
        matClose(pMatFile);
        throw cuDFNsys::ExceptionsPause("cannot create pMxArray or pData\n");
    }

    for (size_t j = 0; j < NUM_eles; ++j)
        pData1[j] = (double)data_[j];

    memcpy((void *)(mxGetPr(pMxArray1)), (void *)pData1, sizeof(pData1));

    const char *field_ = field_name.c_str();

    matPutVariable(pMatFile, field_, pMxArray1);

    mxDestroyArray(pMxArray1);

    matClose(pMatFile);
}; // WriteMat
/************************************ Explicit Instantiate *****************************/
template void cuDFNsys::MatlabAPI::WriteMat(const string &FileKey_mat, const string &mode, const size_t &NUM_eles, const size_t &rows, const size_t &cols, int *data_, /*column major data*/ const string &field_name);
template void cuDFNsys::MatlabAPI::WriteMat(const string &FileKey_mat, const string &mode, const size_t &NUM_eles, const size_t &rows, const size_t &cols, double *data_, /*column major data*/ const string &field_name);
template void cuDFNsys::MatlabAPI::WriteMat(const string &FileKey_mat, const string &mode, const size_t &NUM_eles, const size_t &rows, const size_t &cols, float *data_, /*column major data*/ const string &field_name);
template void cuDFNsys::MatlabAPI::WriteMat(const string &FileKey_mat, const string &mode, const size_t &NUM_eles, const size_t &rows, const size_t &cols, size_t *data_, /*column major data*/ const string &field_name);
template void cuDFNsys::MatlabAPI::WriteMat(const string &FileKey_mat, const string &mode, const size_t &NUM_eles, const size_t &rows, const size_t &cols, uint *data_, /*column major data*/ const string &field_name);