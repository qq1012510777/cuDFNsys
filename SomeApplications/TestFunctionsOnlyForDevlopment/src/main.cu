#include "cuDFNsys.cuh"

int main(int argc, char *argv[])
{

    try
    {
        time_t t;
        time(&t);
        cuDFNsys::DFN<double> my_dfn;
        cuDFNsys::MeshDFN<double> meshGen;
        cuDFNsys::FlowDFN<double> flowDFN;

        std::ifstream File_d("./Class_DFN.h5");
        bool DFSW = File_d.good();

        if (!DFSW)
        {

            for (int i = 0; i < 10000; ++i)
            {
                my_dfn.RandomSeed = (unsigned long)t + i;
                my_dfn.LoadDFNFromCSV("InputParametersForStochasticDFN");
                my_dfn.IdentifyIntersectionsClusters(true);
                if (my_dfn.PercolationCluster.size() > 0)
                    break;
            }
            //my_dfn.LoadDFNFromCSV("InputParametersForDeterministicDFN");
            //

            my_dfn.SpatialPeriodicity();
            my_dfn.IdentifyIntersectionsClusters(true);
            //my_dfn.StoreInH5("Class_DFN");
            my_dfn.Visualization("DFN_VISUAL_I", "DFN_VISUAL_I", "DFN_VISUAL_I",
                                 true, true, true, true);

            meshGen.MinElementSize = 1;
            meshGen.MaxElementSize = 3;
            meshGen.MeshGeneration(my_dfn);
            meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                                  "DFN_MESH_VISUAL", true, true);
            flowDFN.MuOverRhoG = 1;
            flowDFN.InletHead = 100;
            flowDFN.OutletHead = 0;
            flowDFN.IfPeriodic = true;
            flowDFN.ConsTq = 1e-12;

            flowDFN.FlowSimulation(my_dfn, meshGen);
            flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL",
                                  "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL");
            flowDFN.StoreInH5("Class_FLOW");
            meshGen.StoreInH5("Class_MESH");
            my_dfn.StoreInH5("Class_DFN");
            cout << "meanV: " << flowDFN.FlowData.MeanVelocity << endl;
            cout << "maxV: " << flowDFN.FlowData.MaxVelocity << endl;
            cout << "meanGridsize: " << meshGen.MeshData.MeanGridSize << endl;
            cout << "recommended DeltaT: "
                 << pow(meshGen.MeshData.MeanGridSize, 0.5) /
                        flowDFN.FlowData.MaxVelocity / 2
                 << endl;
            //return 0;
        }
        else
        {
            my_dfn.LoadClassFromH5("Class_DFN");
            my_dfn.IdentifyIntersectionsClusters(true);
            meshGen.LoadClassFromH5("Class_MESH");
            flowDFN.LoadClassFromH5("Class_FLOW");
        }

        cuDFNsys::PTDFN<double> particleTracking;
        particleTracking.IfPeriodic = true;
        particleTracking.LoadParametersFromCSV("PT_parameters");
        // particleTracking.NumParticles = 20000;
        // particleTracking.NumTimeSteps = 1000;
        particleTracking.MolecularDiffusion =
            7.5 * flowDFN.FlowData.MeanVelocity / 700;
        particleTracking.DeltaT = pow(meshGen.MeshData.MeanGridSize, 0.5) /
                                  flowDFN.FlowData.MaxVelocity / 2;
        //particleTracking.InjectionMethod = "Flux-weighted";
        //particleTracking.OutputAllPTInformationOrFPTCurve = true;
        //particleTracking.SpacingOfControlPlanes = 3000;
        //particleTracking.IfOutputVarianceOfDisplacementsEachStep = true;
        //particleTracking.IfInjectAtCustomedPlane = true;
        //particleTracking.CustomedPlaneInjection = 10;
        //particleTracking.IfUseFluxWeightedOrEqualProbableMixingIntersection = true;

        particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN);
        particleTracking.Visualization(my_dfn, meshGen, flowDFN,
                                       "DFN_DISPERSION", "DFN_DISPERSION",
                                       "DFN_FLOW_VISUAL");
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