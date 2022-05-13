// ====================================================
// NAME:        ClearBadH5FileForPercolationTest
// DESCRIPTION: clear bad h5 file for percolation test.
// AUTHOR:      Tingchang YIN
// DATE:        13/05/2022
// ====================================================

#include "cuDFNsys.cuh"
#include <stdio.h>
#include <unistd.h>

int main(int argc, char *argv[])
{

    int MODEL_NO = atoi(argv[1]);
    int MC_NO_total = atoi(argv[2]);

    for (int i = 0; i < MC_NO_total; ++i)
    {
        cuDFNsys::HDF5API h5out;
        int MC_NO = i + 1;

        string filename = "Datafile_Model_" + cuDFNsys::ToStringWithWidth(MODEL_NO, 3) +
                          "_MC_" + cuDFNsys::ToStringWithWidth(MC_NO, 5) + ".h5";

        uint2 dim_e = make_uint2(1, 1);

        if (!h5out.IfH5FileExist(filename))
        {
            // do no thing
        }
        else
        {
            try
            {
                vector<double> Looptimes_k =
                    h5out.ReadDataset(filename, "N", "Loop_times");
            }
            catch (...)
            {
                // no dataset "Loop_times" existing means that I did not start loop
                // remove
                cout << "I am removing \"" << filename << "\".\n";
                std::remove(filename.c_str());
            }
        }
    };
    cout << "Finished clearing!\n";
    return 0;
};