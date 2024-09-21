#include "cuDFNsys.cuh"
#include <string>
int main(int argc, char *argv[])
{
    try
    {

        int PercolationDirection = atoi(argv[1]);

        cuDFNsys::DFN<double> my_dfn;
        cuDFNsys::MeshDFN<double> meshGen;

        meshGen.LoadClassFromH5("Class_MESH");
        my_dfn.LoadClassFromH5("Class_DFN_" + std::to_string(PercolationDirection) + ".h5");
        meshGen.ChangePecolationDirectionAndRenumberingEdge(my_dfn.PercoDir,
                                                            my_dfn.DomainSizeX,
                                                            my_dfn.DomainDimensionRatio);
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
        throw;
    }
    return 0;
}