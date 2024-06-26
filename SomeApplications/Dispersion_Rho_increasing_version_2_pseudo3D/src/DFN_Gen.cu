#include "cuDFNsys.cuh"
#include <string>
int main(int argc, char *argv[])
{
    try
    {
        time_t t;
        time(&t);

        cuDFNsys::DFN<double> my_dfn;
        my_dfn.IfPseudo3D = true;
        
        if (argc <= 3)
        {
            if (argc == 3)
                my_dfn.RandomSeed = (unsigned long)(atoi(argv[2]));
            else
                my_dfn.RandomSeed = (unsigned long)t;
        }

        if (argc == 4)
        {
            if (atoi(argv[3]) == -1)
                my_dfn.RandomSeed = (unsigned long)t;
            else
                my_dfn.RandomSeed = (unsigned long)(atoi(argv[3]));
        }

        if (argc <= 3)
        {
            cout << "random seed: " << my_dfn.RandomSeed << endl;

            string CSVname = argv[1];

            my_dfn.LoadDFNFromCSV(CSVname);

            my_dfn.IdentifyIntersectionsClusters(true);

            my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true,
                                 true, true, true);

            my_dfn.StoreInH5("Class_DFN");
        }

        if (argc == 4)
        {
            cout << "random seed: " << my_dfn.RandomSeed << endl;
            cuDFNsys::DFN<double> my_dfn_2;
            string CSVname_1 = argv[1];
            string CSVname_2 = argv[2];
            my_dfn.LoadDFNFromCSV(CSVname_1);
            my_dfn_2.LoadDFNFromCSV(CSVname_2);

            int NUmfrac1 = my_dfn.NumFracturesTotal;

            if (my_dfn.DomainSizeX != my_dfn_2.DomainSizeX ||
                my_dfn.DomainDimensionRatio.x !=
                    my_dfn_2.DomainDimensionRatio.x ||
                my_dfn.DomainDimensionRatio.y !=
                    my_dfn_2.DomainDimensionRatio.y ||
                my_dfn.DomainDimensionRatio.z !=
                    my_dfn_2.DomainDimensionRatio.z)
            {
                throw cuDFNsys::ExceptionsPause(
                    "Domain dimensions are not consistent for generating both "
                    "deterministic and stochastic fractures!\n");
            }

            my_dfn.NumFracturesTotal += my_dfn_2.NumFracturesTotal;
            my_dfn.FracturesHost.resize(my_dfn.NumFracturesTotal);
            my_dfn.FracturesDevice.resize(my_dfn.NumFracturesTotal);

            thrust::copy(my_dfn_2.FracturesHost.begin(),
                         my_dfn_2.FracturesHost.end(),
                         my_dfn.FracturesHost.begin() + NUmfrac1);
            thrust::copy(my_dfn_2.FracturesDevice.begin(),
                         my_dfn_2.FracturesDevice.end(),
                         my_dfn.FracturesDevice.begin() + NUmfrac1);
            my_dfn.FracturesDevicePtr =
                thrust::raw_pointer_cast(my_dfn.FracturesDevice.data());
            my_dfn.IdentifyIntersectionsClusters(true);

            my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true,
                                 true, true, true);

            my_dfn.StoreInH5("Class_DFN");
        }
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
