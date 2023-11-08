#include "cuDFNsys.cuh"

int main(int argc, char *argv[])
{

    try
    {
        time_t t;
        time(&t);
        cuDFNsys::DFN<double> my_dfn;
        my_dfn.RandomSeed = (unsigned long)t;
        my_dfn.LoadDFNFromCSV("InputParametersForStochasticDFN");
        //my_dfn.LoadDFNFromCSV("InputParametersForDeterministicDFN");
        my_dfn.IdentifyIntersectionsClusters(true);
        my_dfn.Visualization("DFN_VISUAL_I", "DFN_VISUAL_I", "DFN_VISUAL_I",
                             true, true, true, true);
        my_dfn.StoreInH5("Class_DFN");

        //my_dfn.SpatialPeriodicity();
        //my_dfn.IdentifyIntersectionsClusters(true);
        //my_dfn.StoreInH5("Class_DFN");

        cuDFNsys::MeshDFN<double> meshGen;
        meshGen.MinElementSize = 1;
        meshGen.MaxElementSize = 3;
        meshGen.MeshGeneration(my_dfn);
        meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                              "DFN_MESH_VISUAL", true, true);

        cout << "Length_Inlet: " << meshGen.MeshData.InletTraceLength
             << ", Length_Outlet: " << meshGen.MeshData.OutletTraceLength
             << endl;

        cuDFNsys::FlowDFN<double> flowDFN;

        flowDFN.MuOverRhoG = 1;
        flowDFN.InletHead = 100;
        flowDFN.OutletHead = 0;
        //flowDFN.IfPeriodic = true;
        //flowDFN.ConsTq = 1e-12;

        flowDFN.FlowSimulation(my_dfn, meshGen);
        flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL",
                              "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL");

        cout << flowDFN.FlowData.Permeability < < < < endl;
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