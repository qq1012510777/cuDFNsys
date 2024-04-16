#include "cuDFNsys.cuh"
#include <string>
int main(int argc, char *argv[])
{
    try
    {

        cuDFNsys::DFN<double> my_dfn;
        cuDFNsys::MeshDFN<double> meshGen;
        cuDFNsys::FlowDFN<double> flow;
        cuDFNsys::PTDFN<double> particleTracking;

        my_dfn.LoadClassFromH5("Class_DFN");
        my_dfn.IdentifyIntersectionsClusters(true);
        meshGen.LoadClassFromH5("Class_MESH");
        flow.LoadClassFromH5("Class_FLOW");

        particleTracking.LoadParametersFromCSV("PTPara");

        double Start_time = cuDFNsys::CPUSecond();
        particleTracking.ParticleTracking(my_dfn, meshGen, flow);
        cout << "PT total run time: " << cuDFNsys::CPUSecond() - Start_time << endl;
        particleTracking.Visualization(my_dfn, meshGen, flow,
                                       "DFN_DISPERSION", "DFN_DISPERSION",
                                       "DFN_FLOW_VISUAL");
        
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