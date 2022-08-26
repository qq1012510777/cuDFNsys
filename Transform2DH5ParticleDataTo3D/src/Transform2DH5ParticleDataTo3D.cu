// ====================================================
// NAME:        A test case
// DESCRIPTION: Call cuDFNsys functions to do simulation and test.
// AUTHOR:      Tingchang YIN
// DATE:        30/06/2022
// ====================================================

#include "cuDFNsys.cuh"
#include <omp.h>
#include <unistd.h>

#ifdef USE_DOUBLES
typedef double _DataType_;
#else
typedef float _DataType_;
#endif

int main(int argc, char *argv[])
{

    try
    {
        _DataType_ L;
        //uint DSIZE;
        uint Nproc = 10;
        if (argv[1] != NULL)
            Nproc = atoi(argv[1]);
        Nproc += 0;
        string FracH5 = "FracturesForParticle.h5";
        string mshfile = "DFN_mesh_1.mat";

        thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host;

        cuDFNsys::InputObjectData<_DataType_> lk;
        lk.InputFractures(FracH5, Frac_verts_host, L);
        //DSIZE = Frac_verts_host.size();

        string ParticlePosition = "ParticlePositionResult/ParticlePosition";
        string DispersionInfo = "ParticlePositionResult/DispersionInfo";

        string DispersionInfo_h5 = DispersionInfo + ".h5";

        cuDFNsys::HDF5API h5g;
        vector<uint> Tem_p = h5g.ReadDataset<uint>(DispersionInfo_h5, "N", "NumOfSteps");
        uint ExistingNumsteps = Tem_p[0];

        Tem_p = h5g.ReadDataset<uint>(DispersionInfo_h5, "N", "BlockNOPresent");
        uint BlockNOPresent = Tem_p[0];

        Tem_p = h5g.ReadDataset<uint>(DispersionInfo_h5, "N", "SizeOfDataBlock");
        uint SizeOfDataBlock = Tem_p[0];

        cuDFNsys::MatlabAPI mu;

        Eigen::MatrixXd ElementFracTag;
        mu.ReadMat(mshfile, "element_Frac_Tag", ElementFracTag);
        //cout << ElementFracTag << endl;
        for (uint i = 0; i <= ExistingNumsteps; ++i)
        {
            cout << "step " << i << " is outputing ...\n";
            string filename = " ", outfilename = "";
            if (i == 0)
            {

                filename = ParticlePosition + "Init.h5";
                outfilename = ParticlePosition + "Init3D.h5";
                h5g.NewFile(outfilename);
            }
            else
            {
                filename = ParticlePosition + "Block" +
                           cuDFNsys::ToStringWithWidth((uint)ceil((double)i / ((double)SizeOfDataBlock)), 10) + ".h5";

                outfilename = ParticlePosition + "Block" +
                              cuDFNsys::ToStringWithWidth((uint)ceil((double)i / ((double)SizeOfDataBlock)), 10) + "3D.h5";
                if (i % SizeOfDataBlock == 1)
                    h5g.NewFile(outfilename);
            }
            vector<_DataType_> temp2Dpos = h5g.ReadDataset<_DataType_>(filename, "N", "Step_" + cuDFNsys::ToStringWithWidth(i, 10));

            uint numParticles = temp2Dpos.size() / 2;
            //cout << "numParticles: " << numParticles << endl;

            vector<uint> EleTag = h5g.ReadDataset<uint>(filename, "N", "ParticleIDAndElementTag_" + cuDFNsys::ToStringWithWidth(i, 10));

            //for (uint j = 0; j < EleTag.size(); ++j)
            //cout << EleTag[j] << endl;
            //cout << EleTag.size() << ", " << numParticles << endl;

            _DataType_ *Position3D = new _DataType_[numParticles * 3];

            //cout << 1 << endl;

#pragma omp parallel for schedule(static) num_threads(Nproc)
            for (uint j = 0; j < numParticles; ++j)
            {
                cuDFNsys::Vector3<_DataType_> Pos = cuDFNsys::MakeVector3(temp2Dpos[j], temp2Dpos[j + numParticles], (_DataType_)0.0);
                //cout << "2D:  " << temp2Dpos[j] << ", " << temp2Dpos[j + numParticles] << endl;

                uint EleTag_j = EleTag[numParticles + j];

                uint FracTag_j = ElementFracTag(EleTag_j - 1, 0) - 1;
                //cout << EleTag_j - 1 << ",  " << FracTag_j << endl;

                _DataType_ Rotate2DTo3D[3][3];
                Frac_verts_host[FracTag_j].RoationMatrix(Rotate2DTo3D, 23);

                Pos = cuDFNsys::ProductSquare3Vector3<_DataType_>(Rotate2DTo3D, Pos);
                Pos = cuDFNsys::MakeVector3(Pos.x + Frac_verts_host[FracTag_j].Center.x,
                                            Pos.y + Frac_verts_host[FracTag_j].Center.y,
                                            Pos.z + Frac_verts_host[FracTag_j].Center.z);

                Position3D[j] = Pos.x;
                Position3D[j + numParticles] = Pos.y;
                Position3D[j + 2 * numParticles] = Pos.z;
                //cout << Pos.x << ", " << Pos.y << ", " << Pos.z << endl;
            }
            //cout << 2 << endl;
            uint2 dim_pos = make_uint2(3, numParticles);

            h5g.AddDataset(outfilename, "N", "Step_" + cuDFNsys::ToStringWithWidth(i, 10), Position3D,
                           dim_pos);

            delete[] Position3D;
            Position3D = NULL;
        };
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