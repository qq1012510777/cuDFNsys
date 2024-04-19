#include "cuDFNsys.cuh"
#include <string>
int main(int argc, char *argv[])
{
    try
    {

        cuDFNsys::DFN<double> my_dfn;
        cuDFNsys::MeshDFN<double> meshGen;
        cuDFNsys::FlowDFN<double> flow;

        my_dfn.LoadClassFromH5("Class_DFN");
        my_dfn.IdentifyIntersectionsClusters(true);
        meshGen.LoadClassFromH5("Class_MESH");

        string CSVname = argv[1];

        flow.LoadParametersFromCSV(CSVname);
        flow.FlowSimulation(my_dfn, meshGen);
        flow.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL",
                           "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL");
        flow.StoreInH5("Class_FLOW");
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