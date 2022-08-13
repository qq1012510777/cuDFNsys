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

// ====================================================
// NAME:        MatlabAPI::WriteMatString
// DESCRIPTION: write data: string.
// AUTHOR:      Tingchang YIN
// DATE:        12/08/2022
// ====================================================
void cuDFNsys::MatlabAPI::WriteMatString(const string &FileKey_mat,
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
}; // WriteMatString

// ====================================================
// NAME:        MatlabAPI::ReadMatString
// DESCRIPTION: Read data: string.
// AUTHOR:      Tingchang YIN
// DATE:        12/08/2022
// ====================================================
string cuDFNsys::MatlabAPI::ReadMatString(const string &FileKey_mat, const string &fieldname)
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
}; // ReadMatString

// ====================================================
// NAME:        MatlabAPI::ReadMatString
// DESCRIPTION: Read data.
// AUTHOR:      Tingchang YIN
// DATE:        12/08/2022
// ====================================================
void cuDFNsys::MatlabAPI::ReadMat(const string &FileKey_mat,
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
}; // ReadMat