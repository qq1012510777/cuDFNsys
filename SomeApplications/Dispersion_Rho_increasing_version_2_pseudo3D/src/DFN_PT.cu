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
        cout << "PT total run time: " << cuDFNsys::CPUSecond() - Start_time
             << endl;
        particleTracking.Visualization(my_dfn, meshGen, flow, "DFN_DISPERSION",
                                       "DFN_DISPERSION", "DFN_FLOW_VISUAL");
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        cout << "Warning : " << e.what() << endl;
        exit(0);
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << "Warning : " << e.what() << endl;
        exit(0);
    }
    catch (H5::Exception &e)
    {
        cout << "Warning : "
             << "H5::Exception\n";
        exit(0);
    }
    catch (H5::FileIException &e)
    {
        cout << "Warning : "
             << "H5::Exception\n";
        exit(0);
    }
    catch (...)
    {
        cout << "Warning : "
             << "Unknown exceptions!\n";
        exit(0);
    }
    return 0;
}