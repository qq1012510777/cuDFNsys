#include "cuDFNsys.cuh"
#include <string>
int main(int argc, char *argv[])
{
    try
    {
        string HDFNname = argv[1];

        int PercolationDirection = atoi(argv[2]);

        cuDFNsys::DFN<double> my_dfn;

        my_dfn.LoadClassFromH5(HDFNname);

        my_dfn.PercoDir = PercolationDirection;

        my_dfn.PercolationCluster.clear();
        if (my_dfn.IntersectionMap.size() > 0 && my_dfn.FracturesHost.size() > 1)
        {
            cuDFNsys::IdentifyPercolationCluster<double> IdentiClu{
                my_dfn.ListClusters,        // all clusters
                my_dfn.FracturesHost,       // host vector of fractures
                my_dfn.PercoDir,            // percolation direction / flow direction
                my_dfn.PercolationCluster}; // percolation cluster
        }
        else if (my_dfn.IntersectionMap.size() == 0)
        {
            my_dfn.ListClusters.resize(my_dfn.FracturesHost.size());
            for (int i = 0; i < my_dfn.FracturesHost.size(); ++i)
            {
                my_dfn.ListClusters[i].push_back(i);

                if (my_dfn.FracturesHost[i].ConnectModelSurf[my_dfn.PercoDir * 2] &&
                    my_dfn.FracturesHost[i].ConnectModelSurf[my_dfn.PercoDir * 2 + 1])
                    my_dfn.PercolationCluster.push_back(i);
            }
        }

        my_dfn.StoreInH5("Class_DFN_" + std::to_string(PercolationDirection));
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
