#include "cuDFNsys.cuh"
int main()
{
    cuDFNsys::DFN<double> my_dfn;
    my_dfn.LoadClassFromH5("Class_DFN");
    my_dfn.Visualization("DFN_VISUAL_II", "DFN_VISUAL_II", "DFN_VISUAL_II",
                         true, true, true, true);
    cuDFNsys::MeshDFN<double> meshGen;
    meshGen.MinElementSize = 1;
    meshGen.MaxElementSize = 3;
    meshGen.MeshGeneration(my_dfn);
    meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL_I", "DFN_MESH_VISUAL_I",
                          "DFN_MESH_VISUAL_I", true, true);
    meshGen.StoreInH5("Class_MESH");
    my_dfn.StoreInH5("Class_DFN");
    return 0;
};