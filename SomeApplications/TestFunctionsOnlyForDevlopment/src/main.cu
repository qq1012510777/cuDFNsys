#include "cuDFNsys.cuh"

int main(int argc, char *argv[])
{

    try
    {
        time_t t;
        time(&t);
        cuDFNsys::DFN<double> my_dfn;
        my_dfn.RandomSeed = (unsigned long)t;
        //my_dfn.LoadDFNFromCSV("InputParametersForStochasticDFN");
        my_dfn.LoadDFNFromCSV("InputParametersForDeterministicDFN");
        my_dfn.IdentifyIntersectionsClusters(true);
        my_dfn.Visualization("DFN_VISUAL_I", "DFN_VISUAL_I", "DFN_VISUAL_I",
                             true, true, true, true);
        my_dfn.StoreInH5("Class_DFN");
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        cout << e.what() << endl;
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << e.what() << endl;
    }
    catch (...)
    {
        throw;
    };
    return 0;
};