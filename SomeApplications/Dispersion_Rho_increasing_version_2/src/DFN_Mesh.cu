#include "cuDFNsys.cuh"
#include <string>
int main(int argc, char *argv[])
{
    try
    {

        cuDFNsys::DFN<double> my_dfn;
        cuDFNsys::MeshDFN<double> meshGen;

        my_dfn.LoadClassFromH5("Class_DFN");
        my_dfn.IdentifyIntersectionsClusters(true);

        string CSVname = argv[1];
        meshGen.LoadParametersFromCSV(CSVname);
        meshGen.MeshGeneration(my_dfn);

        my_dfn.StoreInH5("Class_DFN");

        meshGen.StoreInH5("Class_MESH");

        meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                              "DFN_MESH_VISUAL", true, true);
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        cout << "Warning : " << e.what() << endl;
        throw;
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << "Warning : " << e.what() << endl;
        throw;
    }
    catch (H5::Exception &e)
    {
        cout << "Warning : "
             << "H5::Exception\n";
        throw;
    }
    catch (...)
    {
        cout << "Warning : "
             << "Unknown exceptions!\n";
        throw;
    }
    return 0;
}