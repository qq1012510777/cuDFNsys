#include "cuDFNsys.cuh"
#include <string>
int main(int argc, char *argv[])
{
    try
    {
        time_t t;
        time(&t);

        cuDFNsys::DFN<double> my_dfn;

        if (argc >= 3)
            my_dfn.RandomSeed = (unsigned long)(atoi(argv[2]));
        else
            my_dfn.RandomSeed = (unsigned long)t;
        
        cout << "random seed: " << my_dfn.RandomSeed << endl;

        string CSVname = argv[1];

        my_dfn.LoadDFNFromCSV(CSVname);

        my_dfn.IdentifyIntersectionsClusters(true);

        my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true,
                             true, true, true);

        my_dfn.StoreInH5("Class_DFN");
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        cout << e.what() << endl;
        throw;
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << e.what() << endl;
        throw;
    }
    catch (H5::Exception &e)
    {
        cout << "H5::Exception\n";
        throw;
    }
    catch (...)
    {
        cout << "Unknown exceptions!\n";
        throw;
    }
    return 0;
}