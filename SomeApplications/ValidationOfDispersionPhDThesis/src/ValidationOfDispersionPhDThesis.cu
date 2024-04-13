#include "cuDFNsys.cuh"
#include <string>

int main(int argc, char *argv[])
{
    time_t t;
    time(&t);
    cuDFNsys::DFN<double> my_dfn;
    cuDFNsys::MeshDFN<double> meshGen;
    cuDFNsys::FlowDFN<double> flowDFN;

    string csvname(argv[1]);

    my_dfn.LoadDFNFromCSV(csvname);

    if (my_dfn.NumFracturesTotal <= 2)
    {
        my_dfn.IdentifyIntersectionsClusters(true);
    }
    else
    {
        my_dfn.ListClusters.resize(my_dfn.NumFracturesTotal);
        my_dfn.PercolationCluster.resize(my_dfn.NumFracturesTotal);
        my_dfn.IntersectionMap.clear();
        for (int i = 0; i < my_dfn.NumFracturesTotal; ++i)
        {
            my_dfn.ListClusters[i].resize(1);
            my_dfn.ListClusters[i][0] = i;
            my_dfn.PercolationCluster[i] = i;
        }
    }
    my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true, false,
                         true, true);

    meshGen.MinElementSize = 1;
    meshGen.MaxElementSize = 3;
    meshGen.MeshGeneration(my_dfn);
    meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                          "DFN_MESH_VISUAL", true, true);

    flowDFN.MuOverRhoG = 1;
    flowDFN.InletHead = 30;
    flowDFN.OutletHead = 0;

    flowDFN.FlowSimulation(my_dfn, meshGen);
    flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL",
                          "DFN_FLOW_VISUAL");
    cout << "flowDFN.FlowData.MeanVelocity: " << flowDFN.FlowData.MeanVelocity
         << endl;
    cout << "Set up the PT parameter.csv or continue?\n";
    int hn = 0;
    cin >> hn;
    if (hn == 0)
        return 0;

    cuDFNsys::PTDFN<double> particleTracking;
    particleTracking.LoadParametersFromCSV(csvname + "_PT_parameters");
    particleTracking.TimeIntervalOutPTInformation = 10;
    //-----------------------------------advction time to cross a grid
    double FactorToDivideDeltaT = 10;

    double DeltaT_Advection =
        pow(meshGen.MeshData.MeanGridSize, 0.5) / flowDFN.FlowData.MaxVelocity;
    DeltaT_Advection /= FactorToDivideDeltaT;
    //-----------------------------------diffusion time to cross a grid
    double DeltaT_Diffusion =
        pow(meshGen.MeshData.MeanGridSize, 0.5) /
        (pow(2 * particleTracking.MolecularDiffusion, 0.5) * 4);
    DeltaT_Diffusion =
        DeltaT_Diffusion * DeltaT_Diffusion / FactorToDivideDeltaT;
    cout << "DeltaT_Advection: " << DeltaT_Advection
         << ", DeltaT_Diffusion: " << DeltaT_Diffusion << endl;

    if (DeltaT_Advection < DeltaT_Diffusion)
    {
        cout << "\n** ADVEction Time is smaller **\n\n";
        particleTracking.DeltaT = DeltaT_Advection;
    }
    else
    {
        cout << "\n** DIFFusion Time is smaller **\n\n";
        particleTracking.DeltaT = DeltaT_Diffusion;
    }
    //--------------------------------------------

    particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN);

    particleTracking.Visualization(my_dfn, meshGen, flowDFN, "DFN_DISPERSION",
                                   "DFN_DISPERSION", "DFN_FLOW_VISUAL");

    return 0;
}