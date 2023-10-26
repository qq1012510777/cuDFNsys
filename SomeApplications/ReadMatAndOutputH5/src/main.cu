/****************************************************************************
* cuDFNsys - simulating flow and transport in 3D fracture networks          *
* Copyright (C) 2022, Tingchang YIN, Sergio GALINDO-TORRES                  *
*                                                                           *
* This program is free software: you can redistribute it and/or modify      *
* it under the terms of the GNU Affero General Public License as            *
* published by the Free Software Foundation, either version 3 of the        *
* License, or (at your option) any later version.                           *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU Affero General Public License for more details.                       *
*                                                                           *
* You should have received a copy of the GNU Affero General Public License  *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.    *
*****************************************************************************/

// ====================================================
// NAME:        A test case
// DESCRIPTION: Call cuDFNsys functions to do simulation and test.
// AUTHOR:      Tingchang YIN
// DATE:        30/06/2022
// ====================================================

#include "cuDFNsys.cuh"
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
        cuDFNsys::HDF5API h5g;
        cuDFNsys::MatlabAPI m1;

        string file_path = "ParticlePositionResult";

        // uint NumParticles = 0;

        Eigen::MatrixXd tmp;

        m1.ReadMat(file_path + "/ParticlePosition_step_0000000.mat",
                   "Delta_T",
                   tmp);

        _DataType_ Delta_T = tmp(0, 0);

        m1.ReadMat(file_path + "/ParticlePosition_step_0000000.mat",
                   "Dispersion_local",
                   tmp);
        _DataType_ Dispersion_local = tmp(0, 0);

        uint NUmsteps = 0;
        m1.ReadMat(file_path + "/ParticlePosition_step_0000000.mat",
                   "NumOfSteps",
                   tmp);
        NUmsteps = (uint)tmp(0, 0);

        string Injection_mode = m1.ReadMatString(file_path + "/ParticlePosition_step_0000000.mat", "Injection_mode");
        string Particle_mode = m1.ReadMatString(file_path + "/ParticlePosition_step_0000000.mat", "Particle_mode");

        string h5file = "ParticlePosition.h5";

        _DataType_ tmp_rr[1] = {Delta_T};
        uint2 dim_scalar = make_uint2(1, 1);

        h5g.NewFile(h5file);
        h5g.AddDataset(h5file, "N", "Delta_T", tmp_rr, dim_scalar);

        tmp_rr[0] = Dispersion_local;
        h5g.AddDataset(h5file, "N", "Dispersion_local", tmp_rr, dim_scalar);

        h5g.AddDatasetString(h5file, "N", "Injection_mode", Injection_mode);
        h5g.AddDatasetString(h5file, "N", "Particle_mode", Particle_mode);

        for (uint i = 0; i <= NUmsteps; ++i)
        {
            string matkey = file_path + "/ParticlePosition_step_" + cuDFNsys::ToStringWithWidth(i, 7) + ".mat";
            Eigen::MatrixXd databuffer;

            m1.ReadMat(matkey, "particle_position_3D_step_" + to_string(i), databuffer);

            uint2 dim_data = make_uint2(databuffer.cols(),databuffer.rows());
            double *data_ = new double[databuffer.size()];

            Eigen::Map<Eigen::MatrixXd>(data_, databuffer.rows(), databuffer.cols()) = databuffer;

            h5g.AddDataset(h5file, "N", "Step_" + cuDFNsys::ToStringWithWidth(i, 10), data_, dim_data);
            cout << i << " / " << NUmsteps << " finished\n";
            delete[] data_;
            data_ = NULL;
        }
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