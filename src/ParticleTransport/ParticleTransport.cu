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

#include "ParticleTransport/ParticleTransport.cuh"

// ====================================================
// NAME:        ParticleTransport
// DESCRIPTION: constructor of ParticleTransport
//              for t = 0
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
cuDFNsys::ParticleTransport<T>::ParticleTransport(const int &NumOfParticles,
                                                  const int &NumTimeStep,
                                                  T delta_T_,
                                                  T Dispersion_local,
                                                  thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                                                  cuDFNsys::Mesh<T> mesh,
                                                  const cuDFNsys::MHFEM<T> &fem,
                                                  uint Dir_flow,
                                                  T outletcoordinate,
                                                  const string &Particle_mode,
                                                  const string &Injection_mode,
                                                  bool if_cpu, int Nproc, bool record_time, string recordMode)
{
    IfRecordTime = record_time;
    this->Dir = Dir_flow;

    string filename_EdgesSharedEle = "EdgesSharedEle.h5";

    if (recordMode != "OutputAll" && recordMode != "FPTCurve")
        throw cuDFNsys::ExceptionsPause("Undefined Particle information record mode!\n");
    else
        this->RecordMode = recordMode;

    std::ifstream fileer(filename_EdgesSharedEle);
    bool pwqsc = fileer.good();

    if (!pwqsc)
        this->IdentifyEdgesSharedEle(mesh);
    else
    {
        cout << "Loading " << filename_EdgesSharedEle << " ...\n";
        cuDFNsys::HDF5API h5gdd;
        vector<uint> data_EdgesSharedEle = h5gdd.ReadDataset<uint>(filename_EdgesSharedEle, "N", "data");

        uint NUMoo = 1 + _NumOfSharedEleAtMost * 2;

        uint rows__s = data_EdgesSharedEle.size() / NUMoo;

        this->EdgesSharedEle.resize(rows__s);

        for (uint i = 0; i < rows__s; ++i)
        {
            EdgesSharedEle[i].NumSharedEle = data_EdgesSharedEle[i * NUMoo];
            for (uint j = 0; j < _NumOfSharedEleAtMost; ++j)
            {
                EdgesSharedEle[i].EleID[j] = data_EdgesSharedEle[i * NUMoo + j + 1];

                EdgesSharedEle[i].LocalEdgeNO[j] = data_EdgesSharedEle[i * NUMoo + j + 1 + _NumOfSharedEleAtMost];
            }
        }
        // cout << this->EdgesSharedEle[rows__s - 1].NumSharedEle;
        // cout << ", " << this->EdgesSharedEle[rows__s - 1].EleID[0] << ", " << this->EdgesSharedEle[rows__s - 1].EleID[1];
        // cout << ", " << this->EdgesSharedEle[rows__s - 1].LocalEdgeNO[0] << ", " << this->EdgesSharedEle[rows__s - 1].LocalEdgeNO[1] << endl;
        // cout << "Finish loading " << filename_EdgesSharedEle << " ...\n";
    }

    //this->IdentifyNeighborElements(mesh);

    this->InitilizeParticles(NumOfParticles, mesh, fem, Injection_mode);
    this->OutputMSD(0, Fracs, mesh);
    this->OutputParticleInfoStepByStep(0,
                                       delta_T_,
                                       Dispersion_local,
                                       Particle_mode,
                                       Injection_mode,
                                       Fracs, mesh);
    //cout << NumTimeStep << ", " << delta_T_ << ", " << Dispersion_local << " ______ \n";
    if (!if_cpu)
        this->ParticleMovement(1, NumTimeStep, delta_T_, Dispersion_local, Particle_mode,
                               Injection_mode, Fracs, mesh, fem, outletcoordinate);
    else
        this->ParticleMovementCPU(1, NumTimeStep, delta_T_, Dispersion_local, Particle_mode,
                                  Injection_mode, Fracs, mesh, fem, outletcoordinate, Nproc);
}; // ParticleTransport
template cuDFNsys::ParticleTransport<double>::ParticleTransport(const int &NumOfParticles,
                                                                const int &NumTimeStep,
                                                                double delta_T_,
                                                                double Dispersion_local,
                                                                thrust::host_vector<cuDFNsys::Fracture<double>> Fracs,
                                                                cuDFNsys::Mesh<double> mesh,
                                                                const cuDFNsys::MHFEM<double> &fem,
                                                                uint Dir_flow,
                                                                double outletcoordinate,
                                                                const string &Particle_mode,
                                                                const string &Injection_mode,
                                                                bool if_cpu, int Nproc, bool record_time, string recordMode);
template cuDFNsys::ParticleTransport<float>::ParticleTransport(const int &NumOfParticles,
                                                               const int &NumTimeStep,
                                                               float delta_T_,
                                                               float Dispersion_local,
                                                               thrust::host_vector<cuDFNsys::Fracture<float>> Fracs,
                                                               cuDFNsys::Mesh<float> mesh,
                                                               const cuDFNsys::MHFEM<float> &fem,
                                                               uint Dir_flow,
                                                               float outletcoordinate,
                                                               const string &Particle_mode,
                                                               const string &Injection_mode,
                                                               bool if_cpu, int Nproc, bool record_time, string recordMode);

// ====================================================
// NAME:        ParticleTransport
// DESCRIPTION: constructor of ParticleTransport
//              for t != 0; continute the transport
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
cuDFNsys::ParticleTransport<T>::ParticleTransport(const int &NumTimeStep,
                                                  thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                                                  cuDFNsys::Mesh<T> mesh,
                                                  const cuDFNsys::MHFEM<T> &fem,
                                                  uint Dir_flow,
                                                  T outletcoordinate,
                                                  int NumOfParticles_ii,
                                                  T delta_T_ii,
                                                  T Dispersion_local_ii,
                                                  string Particle_mode_ii,
                                                  string Injection_mode_ii,
                                                  string recordMode,
                                                  bool if_cpu, int Nproc, bool record_time)
{
    //cuDFNsys::MatlabAPI M1;
    // if (recordMode != "OutputAll" && recordMode != "FPTCurve")
    //     throw cuDFNsys::ExceptionsPause("Undefined Particle information record mode!\n");
    // else
    //     this->RecordMode = recordMode;

    IfRecordTime = record_time;

    this->Dir = Dir_flow;

    string filename_EdgesSharedEle = "EdgesSharedEle.h5";

    std::ifstream fileer(filename_EdgesSharedEle);
    bool pwqsc = fileer.good();

    if (!pwqsc)
        this->IdentifyEdgesSharedEle(mesh);
    else
    {
        cout << "Loading " << filename_EdgesSharedEle << " ...\n";
        cuDFNsys::HDF5API h5gdd;
        vector<uint> data_EdgesSharedEle = h5gdd.ReadDataset<uint>(filename_EdgesSharedEle, "N", "data");

        uint NUMoo = 1 + _NumOfSharedEleAtMost * 2;

        uint rows__s = data_EdgesSharedEle.size() / NUMoo;

        this->EdgesSharedEle.resize(rows__s);

        for (uint i = 0; i < rows__s; ++i)
        {
            EdgesSharedEle[i].NumSharedEle = data_EdgesSharedEle[i * NUMoo];
            for (uint j = 0; j < _NumOfSharedEleAtMost; ++j)
            {
                EdgesSharedEle[i].EleID[j] = data_EdgesSharedEle[i * NUMoo + j + 1];

                EdgesSharedEle[i].LocalEdgeNO[j] = data_EdgesSharedEle[i * NUMoo + j + 1 + _NumOfSharedEleAtMost];
            }
        }
        // cout << "Finish loading " << filename_EdgesSharedEle << " ...\n";
        // cout << this->EdgesSharedEle[rows__s - 1].NumSharedEle;
        // cout << ", " << this->EdgesSharedEle[rows__s - 1].EleID[0] << ", " << this->EdgesSharedEle[rows__s - 1].EleID[1];
        // cout << ", " << this->EdgesSharedEle[rows__s - 1].LocalEdgeNO[0] << ", " << this->EdgesSharedEle[rows__s - 1].LocalEdgeNO[1] << endl;
    }

    string matfile = DispersionInfo + ".h5";
    //this->IdentifyNeighborElements(mesh);

    std::ifstream file(matfile);
    bool pwq = file.good();

    if (pwq)
    {
        cout << "\n\e[1;32mContinue to perform particle transport\e[0m\n\n";

        cuDFNsys::HDF5API h5g;

        vector<uint> Tem_p = h5g.ReadDataset<uint>(matfile, "N", "NumOfSteps");
        uint ExistingNumsteps = Tem_p[0];

        Tem_p = h5g.ReadDataset<uint>(matfile, "N", "BlockNOPresent");
        this->BlockNOPresent = Tem_p[0];

        Tem_p = h5g.ReadDataset<uint>(matfile, "N", "SizeOfDataBlock");
        this->SizeOfDataBlock = Tem_p[0];

        Tem_p = h5g.ReadDataset<uint>(matfile, "N", "NumParticles");
        this->NumParticles = Tem_p[0];

        string Particle_mode = h5g.ReadDatasetString(matfile, "N", "Particle_mode");
        string Injection_mode = h5g.ReadDatasetString(matfile, "N", "Injection_mode");

        this->RecordMode = h5g.ReadDatasetString(matfile, "N", "RecordMode");

        if (this->RecordMode != "OutputAll" && this->RecordMode != "FPTCurve")
            throw cuDFNsys::ExceptionsPause("Undefined Particle information record mode!\n");

        string file_block_last = this->ParticlePosition + "Block" + cuDFNsys::ToStringWithWidth(this->BlockNOPresent, 10) + ".h5";
        string datasetname_last = "Step_" + cuDFNsys::ToStringWithWidth(ExistingNumsteps, 10);

        // cout << file_block_last << endl;
        // cout << datasetname_last << endl;

        vector<T> LastStepParticle;

        if (this->RecordMode == "OutputAll")
            LastStepParticle = h5g.ReadDataset<T>(file_block_last, "N",
                                                  datasetname_last);
        else if (this->RecordMode == "FPTCurve")
            LastStepParticle = h5g.ReadDataset<T>(ParticlePosition + "LastStep.h5", "N",
                                                  datasetname_last);

        uint NumParDyna = LastStepParticle.size() / 2;
        this->ParticlePlumes.resize(NumParDyna);

        vector<T> Ifreached;

        if (this->RecordMode == "OutputAll")
            Ifreached = h5g.ReadDataset<T>(file_block_last, "N",
                                           "ParticleIDAndElementTag_" + cuDFNsys::ToStringWithWidth(ExistingNumsteps, 10));
        else
            Ifreached = h5g.ReadDataset<T>(ParticlePosition + "LastStep.h5", "N",
                                           "ParticleIDAndElementTag_" + cuDFNsys::ToStringWithWidth(ExistingNumsteps, 10));

        cout << "Loading particles' positions at the last step\n";
        for (uint i = 0; i < NumParDyna; i++)
        {
            this->ParticlePlumes[i].Position2D.x = LastStepParticle[i];
            this->ParticlePlumes[i].Position2D.y = LastStepParticle[i + NumParDyna];
            this->ParticlePlumes[i].ElementID = Ifreached[i + NumParDyna];
            this->ParticlePlumes[i].ParticleID = Ifreached[i];
            // cout << this->ParticlePlumes[i].Position2D.x << ", ";
            // cout << this->ParticlePlumes[i].Position2D.y << ", ";
            // cout << this->ParticlePlumes[i].ElementID << ", ";
            // cout << this->ParticlePlumes[i].IfReachOutletPlane << endl;
        }
        // cout << "Finish loading particles' positions at the last step\n";

        vector<T> Tem_ps = h5g.ReadDataset<T>(matfile, "N", "Delta_T");
        T delta_T_ = Tem_ps[0];
        Tem_ps = h5g.ReadDataset<T>(matfile, "N", "Dispersion_local");
        T Dispersion_local = Tem_ps[0];

        if (!if_cpu)
            this->ParticleMovement(ExistingNumsteps + 1, NumTimeStep,
                                   delta_T_, Dispersion_local,
                                   Particle_mode,
                                   Injection_mode,
                                   Fracs, mesh, fem, outletcoordinate);
        else
            this->ParticleMovementCPU(ExistingNumsteps + 1, NumTimeStep, delta_T_, Dispersion_local, Particle_mode,
                                      Injection_mode, Fracs, mesh, fem, outletcoordinate, Nproc);
    }
    else
    {
        if (recordMode != "OutputAll" && recordMode != "FPTCurve")
            throw cuDFNsys::ExceptionsPause("Undefined Particle information record mode!\n");
        else
            this->RecordMode = recordMode;

        //cout << 1 << endl;
        this->InitilizeParticles(NumOfParticles_ii,
                                 mesh, fem, Injection_mode_ii);
        //cout << 2 << endl;
        this->OutputMSD(0, Fracs, mesh);
        this->OutputParticleInfoStepByStep(0,
                                           delta_T_ii,
                                           Dispersion_local_ii,
                                           Particle_mode_ii,
                                           Injection_mode_ii,
                                           Fracs, mesh);
        //cout << NumTimeStep << ", " << delta_T_ii << ", " << Dispersion_local << " ______ \n";
        //cout << 3 << endl;
        if (!if_cpu)
            this->ParticleMovement(1, NumTimeStep, delta_T_ii,
                                   Dispersion_local_ii,
                                   Particle_mode_ii,
                                   Injection_mode_ii,
                                   Fracs, mesh, fem,
                                   outletcoordinate);
        else
            this->ParticleMovementCPU(1, NumTimeStep, delta_T_ii,
                                      Dispersion_local_ii,
                                      Particle_mode_ii,
                                      Injection_mode_ii,
                                      Fracs, mesh, fem,
                                      outletcoordinate, Nproc);
    }
}; // ParticleTransport
template cuDFNsys::ParticleTransport<double>::ParticleTransport(const int &NumTimeStep,
                                                                thrust::host_vector<cuDFNsys::Fracture<double>> Fracs,
                                                                cuDFNsys::Mesh<double> mesh,
                                                                const cuDFNsys::MHFEM<double> &fem,
                                                                uint Dir_flow,
                                                                double outletcoordinate,
                                                                int NumOfParticles_ii,
                                                                double delta_T_ii,
                                                                double Dispersion_local_ii,
                                                                string Particle_mode_ii,
                                                                string Injection_mode_ii,
                                                                string recordMode,
                                                                bool if_cpu, int Nproc, bool record_time);
template cuDFNsys::ParticleTransport<float>::ParticleTransport(const int &NumTimeStep,
                                                               thrust::host_vector<cuDFNsys::Fracture<float>> Fracs,
                                                               cuDFNsys::Mesh<float> mesh,
                                                               const cuDFNsys::MHFEM<float> &fem,
                                                               uint Dir_flow,
                                                               float outletcoordinate,
                                                               int NumOfParticles_ii,
                                                               float delta_T_ii,
                                                               float Dispersion_local_ii,
                                                               string Particle_mode_ii,
                                                               string Injection_mode_ii,
                                                               string recordMode,
                                                               bool if_cpu, int Nproc, bool record_time);

// ====================================================
// NAME:        ParticleMovement
// DESCRIPTION: movement of particle during one time step
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
void cuDFNsys::ParticleTransport<T>::ParticleMovement(const int &init_NO_STEP,
                                                      const int &NumTimeStep,
                                                      T delta_T_,
                                                      T Dispersion_local,
                                                      const string &Particle_mode,
                                                      const string &Injection_mode,
                                                      thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                                                      cuDFNsys::Mesh<T> mesh,
                                                      const cuDFNsys::MHFEM<T> &fem,
                                                      T outletcoordinate)
{
    thrust::device_vector<cuDFNsys::Particle<T>> ParticlePlumes_DEV = this->ParticlePlumes;
    thrust::device_vector<cuDFNsys::Fracture<T>> Fracsvec_DEV = Fracs;
    thrust::device_vector<cuDFNsys::EdgeToEle> EdgesSharedEle_vec_dev = this->EdgesSharedEle;
    thrust::device_vector<uint> ElementFracTag_DEV = mesh.ElementFracTag;
    thrust::device_vector<T> Velocity_sep_edge;
    Velocity_sep_edge.reserve(fem.VelocityNormalScalarSepEdges.rows());
    thrust::device_vector<cuDFNsys::EleCoor<T>> Coordinate2D_Vec_dev = mesh.Coordinate2D;
    //thrust::device_vector<cuDFNsys::NeighborEle> NeighborEleOfOneEle_dev = this->NeighborEleOfOneEle;

    Eigen::MatrixXd Vg = fem.VelocityNormalScalarSepEdges;
    double *vc = Vg.data();
    Velocity_sep_edge.insert(Velocity_sep_edge.end(), &vc[0], &vc[fem.VelocityNormalScalarSepEdges.rows()]);

    cuDFNsys::Particle<T> *P_DEV = thrust::raw_pointer_cast(ParticlePlumes_DEV.data());
    cuDFNsys::Fracture<T> *Frac_DEV = thrust::raw_pointer_cast(Fracsvec_DEV.data());
    cuDFNsys::EdgeToEle *EdgesSharedEle_DEV = thrust::raw_pointer_cast(EdgesSharedEle_vec_dev.data());
    uint *EleToFracID_ptr = thrust::raw_pointer_cast(ElementFracTag_DEV.data());
    T *velocity_ptr = thrust::raw_pointer_cast(Velocity_sep_edge.data());
    cuDFNsys::EleCoor<T> *Coordinate2D_Vec_dev_ptr = thrust::raw_pointer_cast(Coordinate2D_Vec_dev.data());
    //cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr = thrust::raw_pointer_cast(NeighborEleOfOneEle_dev.data());

    uint NumPart_dynamic = this->ParticlePlumes.size();

    for (uint i = init_NO_STEP; i <= NumTimeStep + init_NO_STEP; ++i)
    {
        thrust::host_vector<uint> Particle_runtime_error(NumPart_dynamic, 0);
        thrust::device_vector<uint> Particle_runtime_error_dev = Particle_runtime_error;
        uint *Particle_runtime_error_dev_pnt = thrust::raw_pointer_cast(Particle_runtime_error_dev.data());

        time_t t;
        time(&t);

        cout << "The Step " << i << endl;
        double istart_b = cuDFNsys::CPUSecond();
        ParticleMovementOneTimeStepGPUKernel<T><<<NumPart_dynamic / 256 + 1, 256>>>((unsigned long)t + (unsigned long)(i * 100),
                                                                                    delta_T_,
                                                                                    Dispersion_local,
                                                                                    P_DEV,
                                                                                    Frac_DEV,
                                                                                    EdgesSharedEle_DEV,
                                                                                    Coordinate2D_Vec_dev_ptr,
                                                                                    //NeighborEleOfOneEle_dev_ptr,
                                                                                    EleToFracID_ptr,
                                                                                    velocity_ptr,
                                                                                    this->Dir,
                                                                                    outletcoordinate,
                                                                                    NumPart_dynamic,
                                                                                    mesh.Element3D.size(),
                                                                                    i,
                                                                                    Particle_runtime_error_dev_pnt);
        cudaDeviceSynchronize();
        Particle_runtime_error = Particle_runtime_error_dev;

        double ielaps_b = cuDFNsys::CPUSecond() - istart_b;

        if (IfRecordTime)
            RunTimeEveryStep.push_back(ielaps_b);

        uint sum = thrust::reduce(Particle_runtime_error.begin(), Particle_runtime_error.end(), (uint)0, thrust::plus<uint>());

        if (sum != 0)
        {
            if (this->RecordMode == "FPTCurve")
                this->OutputParticleInfoStepByStep(i - 1, // this step is a failure, so last step should be recorded
                                                   delta_T_,
                                                   Dispersion_local,
                                                   Particle_mode,
                                                   Injection_mode,
                                                   Fracs, mesh);

            throw cuDFNsys::ExceptionsPause("Error happens when particle is moving!\n");
        }

        this->ParticlePlumes = ParticlePlumes_DEV;

        double istart = cuDFNsys::CPUSecond();

        // if RecordMode == "FPTCurve", we have to find which particles have been reached after this step!
        // if RecordMode == "FPTCurve", we have to find which particles have been reached after this step!
        // if RecordMode == "FPTCurve", we have to find which particles have been reached after this step!
        // then we read the H5, and change the record the number of steps in the element
        if (this->RecordMode == "FPTCurve")
        {
            string matkeyd = ParticlePosition + "_WhichStepDoesTheParticleReached.h5";

            cuDFNsys::HDF5API h5gds;

            vector<int> WhichStepDoesTheParticleReached = h5gds.ReadDataset<int>(matkeyd, "N",
                                                                                 "WhichStepDoesTheParticleReached");

            bool if_changed = false;

            typename thrust::host_vector<cuDFNsys::Particle<T>>::iterator ityus;
            ityus = this->ParticlePlumes.begin();

            while ((ityus = thrust::find_if(ityus, this->ParticlePlumes.end(), cuDFNsys::PredicateNumOfReachedOutletParticles<T>())) != this->ParticlePlumes.end())
            {
                int idx_tq = ityus - this->ParticlePlumes.begin();

                if (this->ParticlePlumes[idx_tq].ParticleID < 0)
                    if (WhichStepDoesTheParticleReached[abs(this->ParticlePlumes[idx_tq].ParticleID)] == -1)
                    {
                        WhichStepDoesTheParticleReached[abs(this->ParticlePlumes[idx_tq].ParticleID)] = i;
                        //cout << "found particleID " << this->ParticlePlumes[idx_tq].ParticleID << " reached\n";
                        if (!if_changed)
                            if_changed = true;
                    }

                ityus++;
            }

            if (if_changed)
            {
                //cout << "chnage\n";
                uint2 dimu = make_uint2(1, WhichStepDoesTheParticleReached.size());
                h5gds.OverWrite<int>(matkeyd, "N", "WhichStepDoesTheParticleReached", WhichStepDoesTheParticleReached.data(), dimu);
            }
        }

        // int NumReachedParticle = thrust::count_if(this->ParticlePlumes.begin(),
        //                                           this->ParticlePlumes.end(),
        //                                           cuDFNsys::PredicateNumOfReachedOutletParticles<T>());

        this->ParticlePlumes.erase(thrust::remove_if(this->ParticlePlumes.begin(), this->ParticlePlumes.end(),
                                                     cuDFNsys::PredicateNumOfReachedOutletParticles<T>()),
                                   this->ParticlePlumes.end());

        this->ParticlePlumes.shrink_to_fit();

        ParticlePlumes_DEV = this->ParticlePlumes;

        NumPart_dynamic = this->ParticlePlumes.size();

        double ielaps = cuDFNsys::CPUSecond() - istart;

        cout << this->NumParticles - NumPart_dynamic << "/" << this->NumParticles << " reached outlet plane, running time: " << ielaps_b << "; counting time: " << ielaps << "s\n";

        this->OutputMSD(i, Fracs, mesh);
        if (this->RecordMode == "OutputAll")
            this->OutputParticleInfoStepByStep(i,
                                               delta_T_,
                                               Dispersion_local,
                                               Particle_mode,
                                               Injection_mode,
                                               Fracs, mesh);
        else if (this->RecordMode == "FPTCurve")
            if (i == NumTimeStep + init_NO_STEP) // the final step
                this->OutputParticleInfoStepByStep(i,
                                                   delta_T_,
                                                   Dispersion_local,
                                                   Particle_mode,
                                                   Injection_mode,
                                                   Fracs, mesh);

        if (NumPart_dynamic == 0 && Particle_mode == "Particle_tracking")
        {
            cout << "\e[1;32mAll particles reached outlet plane!\e[0m\n";
            break;
        }
        cout << "\n";
    }
}; // ParticleMovement
template void cuDFNsys::ParticleTransport<double>::ParticleMovement(const int &init_NO_STEP,
                                                                    const int &NumTimeStep,
                                                                    double delta_T_,
                                                                    double Dispersion_local,
                                                                    const string &Particle_mode,
                                                                    const string &Injection_mode,
                                                                    thrust::host_vector<cuDFNsys::Fracture<double>> Fracs,
                                                                    cuDFNsys::Mesh<double> mesh,
                                                                    const cuDFNsys::MHFEM<double> &fem,
                                                                    double outletcoordinate);
template void cuDFNsys::ParticleTransport<float>::ParticleMovement(const int &init_NO_STEP,
                                                                   const int &NumTimeStep,
                                                                   float delta_T_,
                                                                   float Dispersion_local,
                                                                   const string &Particle_mode,
                                                                   const string &Injection_mode,
                                                                   thrust::host_vector<cuDFNsys::Fracture<float>> Fracs,
                                                                   cuDFNsys::Mesh<float> mesh,
                                                                   const cuDFNsys::MHFEM<float> &fem,
                                                                   float outletcoordinate);

// ====================================================
// NAME:        ParticleMovementCPU
// DESCRIPTION: movement of particle during one time step
//              on the host side
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
void cuDFNsys::ParticleTransport<T>::ParticleMovementCPU(const int &init_NO_STEP,
                                                         const int &NumTimeStep,
                                                         T delta_T_,
                                                         T Dispersion_local,
                                                         const string &Particle_mode,
                                                         const string &Injection_mode,
                                                         thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                                                         cuDFNsys::Mesh<T> mesh,
                                                         const cuDFNsys::MHFEM<T> &fem,
                                                         T outletcoordinate,
                                                         int Nproc)
{
    //----------------------
    // cuDFNsys::Particle<T> *P_DEV = thrust::raw_pointer_cast(this->ParticlePlumes.data());
    // cuDFNsys::Fracture<T> *Frac_DEV = thrust::raw_pointer_cast(Fracs.data());
    // cuDFNsys::EdgeToEle *EdgesSharedEle_DEV = thrust::raw_pointer_cast(this->EdgesSharedEle.data());
    // cuDFNsys::EleCoor<T> *Coordinate2D_Vec_dev_ptr = thrust::raw_pointer_cast(mesh.Coordinate2D.data());
    // uint *EleToFracID_ptr = thrust::raw_pointer_cast(mesh.ElementFracTag.data());
    if (this->RecordMode != "OutputAll")
        throw cuDFNsys::ExceptionsPause("ParticleMovementCPU only support OutputAll mode\n");

    thrust::host_vector<T> Velocity_sep_edge;
    Velocity_sep_edge.reserve(fem.VelocityNormalScalarSepEdges.rows());
    Eigen::MatrixXd Vg = fem.VelocityNormalScalarSepEdges;
    double *vc = Vg.data();
    Velocity_sep_edge.insert(Velocity_sep_edge.end(), &vc[0], &vc[fem.VelocityNormalScalarSepEdges.rows()]);

    //T *velocity_ptr = thrust::raw_pointer_cast(Velocity_sep_edge.data());

    uint NumPart_dynamic = this->ParticlePlumes.size();

    for (uint i = init_NO_STEP; i <= NumTimeStep + init_NO_STEP; ++i)
    {
        thrust::host_vector<uint> Particle_runtime_error(NumPart_dynamic, 0);
        // uint *Particle_runtime_error_dev_pnt = thrust::raw_pointer_cast(Particle_runtime_error.data());

        time_t t;
        time(&t);

        cout << "The Step (CPU) " << i << endl;
        double istart_b = cuDFNsys::CPUSecond();
        //cout << "P_DEV[0].ElementID: " << (this->ParticlePlumes.data())[0].ElementID << endl;
        cuDFNsys::ParticleMovementOneTimeStepCPU(t, delta_T_, Dispersion_local,
                                                 this->ParticlePlumes.data(),
                                                 Fracs.data(),
                                                 this->EdgesSharedEle.data(),
                                                 mesh.Coordinate2D.data(),
                                                 mesh.ElementFracTag.data(),
                                                 Velocity_sep_edge.data(),
                                                 this->Dir,
                                                 outletcoordinate,
                                                 NumPart_dynamic,
                                                 mesh.Element3D.size(),
                                                 i,
                                                 Particle_runtime_error.data(),
                                                 Nproc);

        double ielaps_b = cuDFNsys::CPUSecond() - istart_b;

        if (IfRecordTime)
            RunTimeEveryStep.push_back(ielaps_b);

        uint sum = thrust::reduce(Particle_runtime_error.begin(), Particle_runtime_error.end(), (uint)0, thrust::plus<uint>());
        if (sum != 0)
            throw cuDFNsys::ExceptionsPause("Error happens when particle is moving!\n");

        double istart = cuDFNsys::CPUSecond();

        this->ParticlePlumes.erase(thrust::remove_if(this->ParticlePlumes.begin(), this->ParticlePlumes.end(),
                                                     cuDFNsys::PredicateNumOfReachedOutletParticles<T>()),
                                   this->ParticlePlumes.end());

        this->ParticlePlumes.shrink_to_fit();
        NumPart_dynamic = this->ParticlePlumes.size();

        double ielaps = cuDFNsys::CPUSecond() - istart;

        cout << this->NumParticles - NumPart_dynamic << "/" << this->NumParticles << " reached outlet plane, running time: " << ielaps_b << "; counting time: " << ielaps << "s\n";

        //cout << this->ParticlePlumes[0].ElementID << endl;
        if (NumPart_dynamic == 0 && Particle_mode == "Particle_tracking")
        {
            cout << "\e[1;32mAll particles reached outlet plane!\e[0m\n";
            break;
        }
        this->OutputParticleInfoStepByStep(i,
                                           delta_T_,
                                           Dispersion_local,
                                           Particle_mode,
                                           Injection_mode,
                                           Fracs, mesh);

        cout << "\n";
    }
    return;
}; // ParticleMovementCPU
template void cuDFNsys::ParticleTransport<double>::ParticleMovementCPU(const int &init_NO_STEP,
                                                                       const int &NumTimeStep,
                                                                       double delta_T_,
                                                                       double Dispersion_local,
                                                                       const string &Particle_mode,
                                                                       const string &Injection_mode,
                                                                       thrust::host_vector<cuDFNsys::Fracture<double>> Fracs,
                                                                       cuDFNsys::Mesh<double> mesh,
                                                                       const cuDFNsys::MHFEM<double> &fem,
                                                                       double outletcoordinate,
                                                                       int Nproc);
template void cuDFNsys::ParticleTransport<float>::ParticleMovementCPU(const int &init_NO_STEP,
                                                                      const int &NumTimeStep,
                                                                      float delta_T_,
                                                                      float Dispersion_local,
                                                                      const string &Particle_mode,
                                                                      const string &Injection_mode,
                                                                      thrust::host_vector<cuDFNsys::Fracture<float>> Fracs,
                                                                      cuDFNsys::Mesh<float> mesh,
                                                                      const cuDFNsys::MHFEM<float> &fem,
                                                                      float outletcoordinate,
                                                                      int Nproc);

// ====================================================
// NAME:        OutputMSD
// DESCRIPTION: OutputMSD curve
// AUTHOR:      Tingchang YIN
// DATE:        22/03/2023
// ====================================================
template <typename T>
void cuDFNsys::ParticleTransport<T>::OutputMSD(const uint &StepNO,
                                               thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                                               cuDFNsys::Mesh<T> mesh)
{
    // Output MSD
    cuDFNsys::HDF5API h5g;
    uint2 dim_scalar = make_uint2(1, 1);
    string MSD_file = "Dispersion_MeanSquareDisplacement.h5";
    //cout << "222222\n";
    if (StepNO == 0)
        h5g.NewFile(MSD_file);
    //cout << "3333\n";
    thrust::device_vector<cuDFNsys::Fracture<T>> Frac_verts_device;
    Frac_verts_device = Fracs;
    cuDFNsys::Fracture<T> *Frac_verts_device_ptr;
    Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());

    thrust::device_vector<uint> ElementFracTag_cuda_dev = mesh.ElementFracTag;
    uint *ElementFracTag_cuda_devptr = thrust::raw_pointer_cast(ElementFracTag_cuda_dev.data());

    // thrust::host_vector<uint> EleTag_host(this->ParticlePlumes.size());
    // for (size_t i = 0; i < EleTag_host.size(); ++i)
    //     EleTag_host[i] = this->ParticlePlumes[i].ElementID;
    // thrust::device_vector<uint> EleTag_device = EleTag_host;
    // uint *EleTag_device_ptr = thrust::raw_pointer_cast(EleTag_device.data());

    thrust::host_vector<T> Position3D(this->ParticlePlumes.size() * 3);
    thrust::device_vector<T> Position3D_dev = Position3D;
    T *Position3D_dev_ptr = thrust::raw_pointer_cast(Position3D_dev.data());

    thrust::device_vector<cuDFNsys::Particle<T>> ParticlePlumes_DEV = this->ParticlePlumes;
    cuDFNsys::Particle<T> *P_DEV = thrust::raw_pointer_cast(ParticlePlumes_DEV.data());

    cuDFNsys::Transform2DTo3DKernel<<<this->ParticlePlumes.size() / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                                    Position3D_dev_ptr,
                                                                                    P_DEV, //temp2DposCUDA_dev_ptr,
                                                                                    ElementFracTag_cuda_devptr,
                                                                                    //EleTag_device_ptr,
                                                                                    this->ParticlePlumes.size());
    cudaDeviceSynchronize();
    Position3D = Position3D_dev;

    std::vector<double> Position3D_zz1(Position3D.begin() + this->ParticlePlumes.size() * 2, Position3D.end());

    Eigen::Map<Eigen::VectorXd> Position3D_zz(Position3D_zz1.data(), Position3D_zz1.size());
    // cout << "size: " << Position3D_zz.rows() << endl;
    // for (uint i = 0; i < this->ParticlePlumes.size(); ++i)
    //     if (Position3D_zz[i] - Position3D[this->ParticlePlumes.size() * 2 + i] != 0)
    //         cout << i << ": " << Position3D_zz[i] - Position3D[this->ParticlePlumes.size() * 2 + i] << endl;

    Position3D_zz1.clear();
    Position3D_zz1.reserve(0);
    Position3D.clear();
    Position3D_dev.clear();
    Position3D.reserve(0);
    Position3D_dev.reserve(0);

    Position3D_zz = Position3D_zz - Eigen::VectorXd::Ones(Position3D_zz.rows()) * 50.0;
    Position3D_zz = Position3D_zz.cwiseAbs();

    double mean_z = Position3D_zz.mean();
    Position3D_zz = Position3D_zz - Eigen::VectorXd::Ones(Position3D_zz.rows()) * mean_z;
    double MSD_g = Position3D_zz.dot(Position3D_zz) / Position3D_zz.rows();

    h5g.AddDataset(MSD_file, "N", "MSD_" + cuDFNsys::ToStringWithWidth(StepNO, 10),
                   &MSD_g, dim_scalar);
    //------------
};
template void cuDFNsys::ParticleTransport<double>::OutputMSD(const uint &StepNO,
                                                             thrust::host_vector<cuDFNsys::Fracture<double>> Fracs,
                                                             cuDFNsys::Mesh<double> mesh);
template void cuDFNsys::ParticleTransport<float>::OutputMSD(const uint &StepNO,
                                                            thrust::host_vector<cuDFNsys::Fracture<float>> Fracs,
                                                            cuDFNsys::Mesh<float> mesh);

// ====================================================
// NAME:        OutputParticleInfoStepByStep
// DESCRIPTION: output the position of particles after a step
//              they are 2D though
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
void cuDFNsys::ParticleTransport<T>::OutputParticleInfoStepByStep(const uint &StepNO,
                                                                  const T delta_T,
                                                                  const T Dispersion_local,
                                                                  const string &Particle_mode,
                                                                  const string &Injection_mode,
                                                                  thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                                                                  cuDFNsys::Mesh<T> mesh)
{

    cuDFNsys::HDF5API h5g;

    string h5dispersioninfo = this->DispersionInfo + ".h5";

    uint2 dim_scalar = make_uint2(1, 1);

    if (StepNO == 0)
    {

        std::ifstream file("ParticlePositionResult");
        bool pwq = file.good();

        if (!pwq)
        {
            const int dir_err = system("mkdir -p ./ParticlePositionResult");
            if (-1 == dir_err)
                throw cuDFNsys::ExceptionsPause("Error creating directory of 'ParticlePositionResult'!\n");
        }

        h5g.NewFile(h5dispersioninfo);
        h5g.AddDatasetString(h5dispersioninfo, "N", "Particle_mode", Particle_mode);
        h5g.AddDatasetString(h5dispersioninfo, "N", "Injection_mode", Injection_mode);
        h5g.AddDatasetString(h5dispersioninfo, "N", "RecordMode", this->RecordMode);

        T po[1] = {delta_T};
        h5g.AddDataset(h5dispersioninfo, "N", "Delta_T", po, dim_scalar);
        po[0] = {Dispersion_local};
        h5g.AddDataset(h5dispersioninfo, "N", "Dispersion_local", po, dim_scalar);
        uint ouy[1] = {this->SizeOfDataBlock};
        h5g.AddDataset(h5dispersioninfo, "N", "SizeOfDataBlock", ouy, dim_scalar);

        ouy[0] = {BlockNOPresent};
        h5g.AddDataset(h5dispersioninfo, "N", "BlockNOPresent", ouy, dim_scalar);

        ouy[0] = {(uint)NumParticles};
        h5g.AddDataset(h5dispersioninfo, "N", "NumParticles", ouy, dim_scalar);

        //h5g.NewFile(DispersionInfo + "_MeanSquareDisplacement.h5");
    }

    uint Step[1] = {StepNO};
    if (StepNO != 0)
        h5g.OverWrite(h5dispersioninfo, "N", "NumOfSteps", Step, dim_scalar);
    else
        h5g.AddDataset(h5dispersioninfo, "N", "NumOfSteps", Step, dim_scalar);

    uint NumParDyna = this->ParticlePlumes.size();

    T *particle_position_3D;
    particle_position_3D = new T[NumParDyna * 2];
    if (particle_position_3D == NULL)
    {
        string AS = "Alloc error in ParticleTransport::OutputParticleInfoStepByStep\n";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    for (size_t i = 0; i < NumParDyna; ++i)
    {
        cuDFNsys::Vector3<T> tmpPos = cuDFNsys::MakeVector3(this->ParticlePlumes[i].Position2D.x,
                                                            this->ParticlePlumes[i].Position2D.y,
                                                            (T)0);

        uint elementID = this->ParticlePlumes[i].ElementID; // from 1
        uint FracID = mesh.ElementFracTag[elementID - 1];   // from 0

        particle_position_3D[i] = this->ParticlePlumes[i].Position2D.x;
        particle_position_3D[i + NumParDyna] = this->ParticlePlumes[i].Position2D.y;
    }

    uint2 dim_data = make_uint2(2, NumParDyna);
    uint2 dim_datauu = make_uint2(2, NumParDyna);
    uint *_IfReaded_and_ElementFracTag_ = new uint[NumParDyna * 2];
    if (_IfReaded_and_ElementFracTag_ == NULL)
    {
        string AS = "Alloc error in ParticleTransport::OutputParticleInfoStepByStep\n";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    for (size_t i = 0; i < NumParDyna; ++i)
    {
        _IfReaded_and_ElementFracTag_[i] = this->ParticlePlumes[i].ParticleID;
        _IfReaded_and_ElementFracTag_[i + NumParDyna] = this->ParticlePlumes[i].ElementID;
    }

    if (StepNO == 0)
    {
        string mat_key;
        if (this->RecordMode == "OutputAll")
        {
            mat_key = ParticlePosition + "Init.h5";

            h5g.NewFile(mat_key);
            h5g.AddDataset(mat_key, "N", "Step_" + cuDFNsys::ToStringWithWidth(StepNO, 10),
                           particle_position_3D, dim_data);
            h5g.AddDataset(mat_key, "N", "ParticleIDAndElementTag_" + cuDFNsys::ToStringWithWidth(StepNO, 10),
                           _IfReaded_and_ElementFracTag_, dim_datauu);
        }
        else if (this->RecordMode == "FPTCurve")
        {
            mat_key = ParticlePosition + "LastStep.h5";
            h5g.NewFile(mat_key);
            h5g.AddDataset(mat_key, "N", "Step_" + cuDFNsys::ToStringWithWidth(StepNO, 10),
                           particle_position_3D, dim_data);
            h5g.AddDataset(mat_key, "N", "ParticleIDAndElementTag_" + cuDFNsys::ToStringWithWidth(StepNO, 10),
                           _IfReaded_and_ElementFracTag_, dim_datauu);

            uint2 dimf2 = make_uint2(1, this->ParticlePlumes.size());
            std::vector<int> WhichStepDoesTheParticleReached(this->ParticlePlumes.size(), -1);

            string matkeyd = ParticlePosition + "_WhichStepDoesTheParticleReached.h5";
            h5g.NewFile(matkeyd);
            h5g.AddDataset(matkeyd, "N", "WhichStepDoesTheParticleReached", WhichStepDoesTheParticleReached.data(), dimf2);
        }
    }
    else
    {
        if (this->RecordMode == "OutputAll")
        {
            if (StepNO > (this->BlockNOPresent * this->SizeOfDataBlock))
            {
                this->BlockNOPresent++;
                uint ouy[1] = {this->BlockNOPresent};

                h5g.OverWrite(h5dispersioninfo, "N", "BlockNOPresent", ouy, dim_scalar);

                string mat_key = ParticlePosition + "Block" +
                                 cuDFNsys::ToStringWithWidth(this->BlockNOPresent, 10) + ".h5";
                h5g.NewFile(mat_key);
            }
            string mat_key = ParticlePosition + "Block" +
                             cuDFNsys::ToStringWithWidth(this->BlockNOPresent, 10) + ".h5";

            h5g.AddDataset(mat_key, "N", "Step_" + cuDFNsys::ToStringWithWidth(StepNO, 10),
                           particle_position_3D, dim_data);
            //cout << dim_data.x << ", " << dim_data.y << endl;
            h5g.AddDataset(mat_key, "N", "ParticleIDAndElementTag_" + cuDFNsys::ToStringWithWidth(StepNO, 10),
                           _IfReaded_and_ElementFracTag_, dim_datauu);
            //cout << dim_datauu.x << ", " << dim_datauu.y << endl;
        }
        else
        {
            string mat_key = ParticlePosition + "LastStep.h5";
            h5g.NewFile(mat_key);

            h5g.AddDataset(mat_key, "N", "Step_" + cuDFNsys::ToStringWithWidth(StepNO, 10),
                           particle_position_3D, dim_data);
            //cout << dim_data.x << ", " << dim_data.y << endl;
            h5g.AddDataset(mat_key, "N", "ParticleIDAndElementTag_" + cuDFNsys::ToStringWithWidth(StepNO, 10),
                           _IfReaded_and_ElementFracTag_, dim_datauu);
        }
    }
    delete[] particle_position_3D;
    particle_position_3D = NULL;

    delete[] _IfReaded_and_ElementFracTag_;
    _IfReaded_and_ElementFracTag_ = NULL;

    // if (StepNO > 0)
    //     exit(0);
}; // OutputParticleInfoStepByStep
template void cuDFNsys::ParticleTransport<double>::OutputParticleInfoStepByStep(const uint &StepNO,
                                                                                const double delta_T,
                                                                                const double Dispersion_local,
                                                                                const string &Particle_mode,
                                                                                const string &Injection_mode,
                                                                                thrust::host_vector<cuDFNsys::Fracture<double>> Fracs,
                                                                                cuDFNsys::Mesh<double> mesh);
template void cuDFNsys::ParticleTransport<float>::OutputParticleInfoStepByStep(const uint &StepNO,
                                                                               const float delta_T,
                                                                               const float Dispersion_local,
                                                                               const string &Particle_mode,
                                                                               const string &Injection_mode,
                                                                               thrust::host_vector<cuDFNsys::Fracture<float>> Fracs,
                                                                               cuDFNsys::Mesh<float> mesh);

// ====================================================
// NAME:        MatlabPlot
// DESCRIPTION: MatlabPlot the particles in the form of video
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
void cuDFNsys::ParticleTransport<T>::MatlabPlot(const string &mat_key,
                                                const string &command_key,
                                                const cuDFNsys::Mesh<T> &mesh,
                                                const cuDFNsys::MHFEM<T> &fem,
                                                const T &L,
                                                double3 DomainDimensionRatio,
                                                bool if_python_visualization,
                                                string PythonName_Without_suffix)
{
    if (if_python_visualization)
    {
        std::ofstream oss(PythonName_Without_suffix + ".py", ios::out);
        oss << "import h5py\n";
        oss << "import numpy as np\n";
        oss << "from mayavi import mlab as ML\n";
        oss << "import math\n";
        oss << "import gc            \n";
        oss << "f = h5py.File('" << mat_key << "')\n";
        oss << "coordinate_3D = np.array(f['coordinate_3D'][:])\n";
        oss << "element_3D = np.array(f['element_3D'][:], dtype=int)\n";
        oss << "InletP = np.array(f['InletP'][:])\n";
        oss << "OutletP = np.array(f['OutletP'][:])\n";
        oss << "L_m = f['L_m'][:][0]\n";
        oss << "DomainDimensionRatio = f['DomainDimensionRatio'][:]\n";
        oss << "pressure_eles = np.array(f['pressure_eles'][:])\n";
        oss << "f.close()\n";

        oss << "mesh = ML.triangular_mesh(coordinate_3D[0, :], coordinate_3D[1, :], coordinate_3D[2, :], np.transpose(element_3D - 1), representation='wireframe', color=(0, 0, 0), line_width=1.0)\n";
        oss << "mesh.mlab_source.dataset.cell_data.scalars = np.transpose(pressure_eles)\n";
        oss << "mesh.mlab_source.dataset.cell_data.scalars.name = 'Cell data'\n";
        oss << "mesh.mlab_source.update()\n";
        oss << "mesh.parent.update()\n";
        oss << "mesh2 = ML.pipeline.set_active_attribute(mesh, cell_scalars='Cell data')\n";
        oss << "s2 = ML.pipeline.surface(mesh2, colormap='rainbow', opacity=1)\n";
        oss << "ML.outline(extent=[-0.5 * DomainDimensionRatio[0] * L_m, 0.5 * DomainDimensionRatio[0] * L_m, -0.5 * DomainDimensionRatio[1] * L_m, 0.5 * DomainDimensionRatio[1] * L_m, -0.5 * DomainDimensionRatio[2] * L_m, 0.5 * DomainDimensionRatio[2] * L_m])\n";
        oss << "ML.axes()\n";
        oss << "s2.module_manager.scalar_lut_manager.data_range = np.array([OutletP[0], InletP[0]])\n";
        oss << "ML.colorbar(object=s2, orientation='vertical')\n";
        oss << "ML.xlabel('x (m)')\n";
        oss << "ML.ylabel('y (m)')\n";
        oss << "ML.zlabel('z (m)')\n";
        oss << "f_2 = h5py.File('./ParticlePositionResult/DispersionInfo.h5')\n";
        oss << "N_steps = int(np.array(f_2['NumOfSteps'][0]))            \n";
        oss << "N_particles = int(np.array(f_2['NumParticles'][0]))\n";
        oss << "BlockNOPresent = int(np.array(f_2['BlockNOPresent'][0]))\n";
        oss << "SizeOfDataBlock = int(np.array(f_2['SizeOfDataBlock'][0]))\n";
        oss << "f_2.close()\n";
        oss << "H5name = \"./ParticlePositionResult/ParticlePositionInit_3D.h5\"\n";
        oss << "H5name_2D = \"./ParticlePositionResult/ParticlePositionInit.h5\"\n";
        oss << "f2 = h5py.File(H5name)\n";
        oss << "f3 = h5py.File(H5name_2D)\n";
        oss << "S = np.array(f2[\"Step_\" + str(0).zfill(10)][:])\n";
        oss << "ParticleID = np.array(f3[\"ParticleIDAndElementTag_\" + str(0).zfill(10)][:])\n";
        oss << "Matrx3D_pso = np.zeros([3, N_particles])\n";
        oss << "Matrx3D_pso[:] = np.nan\n";
        oss << "Matrx3D_pso[:, ParticleID[0, :]] = S[:, :]\n";
        oss << "ShowSomeParticles = [i  for i in range(0, int(N_particles), int(N_particles * 1.0 / N_particles))]\n";
        oss << "SK = ML.points3d(Matrx3D_pso[0, ShowSomeParticles], Matrx3D_pso[1, ShowSomeParticles], Matrx3D_pso[2, ShowSomeParticles], scale_factor=0.7)\n";
        oss << "f2.close()\n";
        oss << "f3.close()\n";
        oss << "del Matrx3D_pso\n";
        oss << "del ParticleID\n";
        oss << "@ML.animate(delay=10)\n";
        oss << "def anim(N_steps):\n";
        oss << "    fe = ML.gcf() \n";
        oss << "    for i in range(1, N_steps + 1):\n";
        oss << "        print('step', i, '/', N_steps)          \n";
        oss << "        H5name = \"./ParticlePositionResult/ParticlePositionBlock\" + str(math.ceil(float(i) / float(SizeOfDataBlock))).zfill(10) + \"_3D.h5\"\n";
        oss << "        H5name_2D = \"./ParticlePositionResult/ParticlePositionBlock\" + str(math.ceil(float(i) / float(SizeOfDataBlock))).zfill(10) + \".h5\"\n";
        oss << "        f4 = h5py.File(H5name)\n";
        oss << "        f5 = h5py.File(H5name_2D)\n";
        oss << "        S = np.array(f4[\"Step_\" + str(i).zfill(10)][:])\n";
        oss << "        ParticleID = np.array(f5[\"ParticleIDAndElementTag_\" + str(i).zfill(10)][:])\n";
        oss << "        Matrx3D_pso = np.zeros([3, N_particles])\n";
        oss << "        Matrx3D_pso[:] = np.nan\n";
        oss << "        Matrx3D_pso[:, ParticleID[0, :]] = S[:, :]\n";
        oss << "        SK.mlab_source.set(x=Matrx3D_pso[0, ShowSomeParticles], y=Matrx3D_pso[1, ShowSomeParticles], z=Matrx3D_pso[2, ShowSomeParticles])\n";
        oss << "        f4.close()\n";
        oss << "        f5.close()\n";
        oss << "        del Matrx3D_pso\n";
        oss << "        del ParticleID\n";
        oss << "        del S\n";
        oss << "        yield\n";
        oss << "        gc.collect(generation=1)\n";
        oss << "anim(N_steps)\n";
        oss << "ML.show()\n";
        oss.close();

        std::ofstream osse(PythonName_Without_suffix + "_ShowTrajectories.py", ios::out);
        osse << "import h5py\n";
        osse << "import numpy as np\n";
        osse << "from mayavi import mlab as ML\n";
        osse << "import math           \n";
        osse << "f = h5py.File('" << mat_key << "')\n";
        osse << "coordinate_3D = np.array(f['coordinate_3D'][:])\n";
        osse << "element_3D = np.array(f['element_3D'][:], dtype=int)\n";
        osse << "InletP = np.array(f['InletP'][:])\n";
        osse << "OutletP = np.array(f['OutletP'][:])\n";
        osse << "L_m = f['L_m'][:][0]\n";
        osse << "DomainDimensionRatio = f['DomainDimensionRatio'][:]\n";
        osse << "pressure_eles = np.array(f['pressure_eles'][:])\n";
        osse << "f.close()           \n";
        osse << "mesh = ML.triangular_mesh(coordinate_3D[0, :], coordinate_3D[1, :], coordinate_3D[2, :], np.transpose(element_3D - 1), opacity=1)\n";
        osse << "mesh.mlab_source.dataset.cell_data.scalars = np.transpose(pressure_eles)\n";
        osse << "mesh.mlab_source.dataset.cell_data.scalars.name = 'Cell data'\n";
        osse << "mesh.mlab_source.update()\n";
        osse << "mesh.parent.update()\n";
        osse << "mesh2 = ML.pipeline.set_active_attribute(mesh, cell_scalars='Cell data')\n";
        osse << "s2 = ML.pipeline.surface(mesh2, colormap='rainbow', opacity=0.1)\n";
        osse << "ML.outline(extent=[-0.5 * DomainDimensionRatio[0] * L_m, 0.5 * DomainDimensionRatio[0] * L_m, -0.5 * DomainDimensionRatio[1] * L_m, 0.5 * DomainDimensionRatio[1] * L_m, -0.5 * DomainDimensionRatio[2] * L_m, 0.5 * DomainDimensionRatio[2] * L_m])\n";

        osse << "ML.axes()\n";
        osse << "s2.module_manager.scalar_lut_manager.data_range = np.array([OutletP[0], InletP[0]])\n";
        osse << "ML.colorbar(object=s2, orientation='vertical')\n";
        osse << "ML.xlabel('x (m)')\n";
        osse << "ML.ylabel('y (m)')\n";
        osse << "ML.zlabel('z (m)')\n";
        osse << "f_2 = h5py.File('./ParticlePositionResult/DispersionInfo.h5')\n";
        osse << "N_steps = int(np.array(f_2['NumOfSteps'][0]))\n";
        osse << "N_particles = int(np.array(f_2['NumParticles'][0]))\n";
        osse << "BlockNOPresent = int(np.array(f_2['BlockNOPresent'][0]))\n";
        osse << "SizeOfDataBlock = int(np.array(f_2['SizeOfDataBlock'][0]))\n";
        osse << "f_2.close()\n";
        osse << "H5name_2D = \"./ParticlePositionResult/ParticlePositionBlock\" + str(BlockNOPresent).zfill(10) + \".h5\"\n";
        osse << "f_3 = h5py.File(H5name_2D)\n";
        osse << "S = np.array(f_3['ParticleIDAndElementTag_' + str(N_steps).zfill(10)][0, :], dtype=int)\n";
        osse << "ReachedParticleNO = np.array([i for i in range(1, N_particles + 1)])\n";
        osse << "ReachedParticleNO = np.delete(ReachedParticleNO, S) - 1\n";
        osse << "\n";
        osse << "connections = list()\n";
        osse << "connections = [(connections + [l, l + N_particles]) for l in range(N_particles)]\n";
        osse << "\n";
        osse << "for i in range(0, N_steps + 1):\n";
        osse << "    print('Trajectories of step', i, 'are plotting')\n";
        osse << "    Matrx3D_pso = np.zeros([3, N_particles * 2])\n";
        osse << "    Matrx3D_pso[:] = np.nan\n";
        osse << "    yut = 0\n";
        osse << "    for j in range(i, i + 2):\n";
        osse << "        if(i + 1 > N_steps):\n";
        osse << "            break\n";
        osse << "        H5name = \"./ParticlePositionResult/ParticlePositionBlock\" + str(math.ceil(float(j) / int(SizeOfDataBlock))).zfill(10) + \"_3D.h5\"\n";
        osse << "        H5name_2D = \"./ParticlePositionResult/ParticlePositionBlock\" + str(math.ceil(float(j) / int(SizeOfDataBlock))).zfill(10) + \".h5\"\n";
        osse << "        if(j == 0):\n";
        osse << "            H5name = \"./ParticlePositionResult/ParticlePositionInit_3D.h5\"\n";
        osse << "            H5name_2D = \"./ParticlePositionResult/ParticlePositionInit.h5\"    \n";
        osse << "        f4 = h5py.File(H5name)\n";
        osse << "        f5 = h5py.File(H5name_2D)\n";
        osse << "        S = np.array(f4[\"Step_\" + str(j).zfill(10)][:])\n";
        osse << "        ParticleID = np.array(f5[\"ParticleIDAndElementTag_\" + str(j).zfill(10)][0, :], dtype = int)\n";
        osse << "        AS = np.in1d(ParticleID, ReachedParticleNO).nonzero()[0]\n";
        osse << "        update_ParticleID = ParticleID[AS]\n";
        osse << "        Matrx3D_pso[:, [int(o) + int(yut * N_particles) for o in update_ParticleID]] = S[:, AS]\n";
        osse << "        f4.close()\n";
        osse << "        f5.close()\n";
        osse << "        del S\n";
        osse << "        del f4\n";
        osse << "        del f5\n";
        osse << "        yut += 1\n";
        osse << "    src = ML.pipeline.scalar_scatter(Matrx3D_pso[0, :], Matrx3D_pso[1, :], Matrx3D_pso[2, :])\n";
        osse << "    src.mlab_source.dataset.lines = connections\n";
        osse << "    src.update()\n";
        osse << "    lines = ML.pipeline.stripper(src)\n";
        osse << "    ML.pipeline.surface(lines, line_width=1, opacity=1)\n";
        osse << "ML.show()\n";
        osse.close();
    }

    if (command_key != "N")
    {
        std::ofstream oss(command_key, ios::out);
        oss << "clc;\nclose all;\nclear all;\ncurrentPath = fileparts(mfilename('fullpath'));\n";
        //oss << "load('" << mat_key << "');\n";
        //oss << "currentPath = fileparts(mfilename('fullpath'));\n";
        oss << "coordinate_3D = h5read([currentPath, '/" << mat_key << "'], '/coordinate_3D');\n";
        oss << "element_3D = h5read([currentPath, '/" << mat_key << "'], '/element_3D');\n";
        oss << "pressure_eles = h5read([currentPath, '/" << mat_key << "'], '/pressure_eles');\n";

        oss << "P_out = " << fem.OutletP << "; P_in = " << fem.InletP << ";\n";
        oss << "Offset_colorbar_value_for_particle = " << L << ";\n";
        oss << "If_video = false;\n";
        oss << "L = h5read([currentPath, '/" << mat_key << "'], '/L_m');\n";
        oss << "DomainDimensionRatio = h5read([currentPath, '/" << mat_key << "'], '/DomainDimensionRatio');\n";
        oss << "cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];\n";
        oss << "cube_frame(:, 1) = 0.5 .* cube_frame(:, 1) .* DomainDimensionRatio(1); ";
        oss << "cube_frame(:, 2) = 0.5 .* cube_frame(:, 2) .* DomainDimensionRatio(2); ";
        oss << "cube_frame(:, 3) = 0.5 .* cube_frame(:, 3) .* DomainDimensionRatio(3);\n";

        oss << "figure(1); view(3); title('DFN flow (mhfem) and particle trajectory'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on\n";
        oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on\n";
        oss << endl;
        oss << "patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.1, 'facealpha', 0.1); colorbar; view(3); hold on\n";
        oss << "caxis([" << fem.OutletP << ", " << fem.InletP << "]);\n";
        oss << "axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * DomainDimensionRatio(1) * L, -1.1 / 2 * DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) * L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * DomainDimensionRatio(3) * L]);\n";
        oss << "pbaspect([DomainDimensionRatio]); hold on\n";

        oss << "N_steps = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], '/NumOfSteps');\n";
        // oss << "S = h5read([currentPath, '/ParticlePositionResult/ParticlePositionInit.h5'], '/Step_0000000000');\n";
        oss << "N_particles = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], '/NumParticles');\n";
        oss << "clear S;\n";

        oss << "BlockNOPresent = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], '/BlockNOPresent');\n";
        oss << "SizeOfDataBlock = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], '/SizeOfDataBlock');\n";

        oss << "S = h5read([currentPath, '/ParticlePositionResult/ParticlePositionBlock', num2str(BlockNOPresent, '%010d'), '.h5'], ['/ParticleIDAndElementTag_', num2str(N_steps, '%010d')]); \n";
        oss << "ReachedParticleNO = [1:1:N_particles];\n";

        oss << "MatchedValue = S(:, 1) + 1;\nReachedParticleNO(ismember(ReachedParticleNO, MatchedValue')) = []; clear S MatchedValue\n";

        oss << "newcolors = rand(size(ReachedParticleNO, 2), 3);\n";

        oss << "delta_t2 = 200;\ninit_ = 0;\nfinal_ = 0;\n";

        oss << "show_trajector_whole = false;\n";
        oss << "if (show_trajector_whole == true)\n";
        oss << "\tfor i = 0 : ceil(N_steps / delta_t2)\n";
        oss << "\t\tinit_ = final_; final_ = init_ + delta_t2;\n";
        oss << "\t\tif (final_ > N_steps)\n";
        oss << "\t\t\tfinal_ = N_steps;\n";
        oss << "\t\tend\n";
        oss << "\t\tfinal_\n";
        oss << "\t\tAK_1 = [];AK_2 = [];AK_3 = [];\n";
        oss << "\t\tfor j = init_:final_\n";
        // oss << "\t\t\tS1 = load([currentPath, '/" << ParticlePosition << "_step_', num2str(j, '%07d'),'.mat']);\n";
        oss << "\t\t\tH5name = []; H5name_2D = [];\n";
        oss << "\t\t\tif (j == 0); H5name = [currentPath, '/ParticlePositionResult/ParticlePositionInit_3D.h5']; else; H5name = [currentPath, '/ParticlePositionResult/ParticlePositionBlock', num2str(ceil(double(j) / double(SizeOfDataBlock)), '%010d'), '_3D.h5']; end;\n";
        oss << "\t\t\tif (j == 0); H5name_2D = [currentPath, '/ParticlePositionResult/ParticlePositionInit.h5']; else; H5name_2D = [currentPath, '/ParticlePositionResult/ParticlePositionBlock', num2str(ceil(double(j) / double(SizeOfDataBlock)), '%010d'), '.h5']; end;\n";

        oss << "\t\t\tS = h5read(H5name, ['/Step_', num2str(j, '%010d')]);\n";
        oss << "\t\t\tParticleID = h5read(H5name_2D, ['/ParticleIDAndElementTag_', num2str(j, '%010d')]);\n"; // /
        oss << "\t\t\tMatrx3D_pso = NaN(N_particles, 3);\n";
        oss << "\t\t\tMatrx3D_pso([ParticleID(:, 1) + 1], :) = S(:, [1 2 3]);\n";

        oss << "\t\t\tAK_1(j - init_ + 1, :) = Matrx3D_pso(:, 1);\n";
        oss << "\t\t\tAK_2(j - init_ + 1, :) = Matrx3D_pso(:, 2);\n";
        oss << "\t\t\tAK_3(j - init_ + 1, :) = Matrx3D_pso(:, 3); clear S Matrx3D_pso\n";
        oss << "\t\tend\n";
        oss << "\t\tAK_1 = AK_1(:, ReachedParticleNO);AK_2 = AK_2(:, ReachedParticleNO);AK_3 = AK_3(:, ReachedParticleNO);\n";
        oss << "\t\tcolororder(newcolors)\n";
        oss << "\t\tplot3(AK_1, AK_2, AK_3, 'linewidth', 2); hold on\n";
        oss << "\t\tif (final_ == N_steps)\n";
        oss << "\t\t\tbreak\n";
        oss << "\t\tend\n";
        oss << "\tend\n\n";
        oss << "clear AK_1 AK_2 AK_3\n";
        oss << "end\n";

        oss << "figure(2); view(3); title('DFN flow (mhfem) and particle tracking'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on\n";
        oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on\n";
        oss << endl;
        oss << "patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.2, 'facealpha', 0.9); view(3); hold on\n";
        oss << "colormap(jet)\n";
        oss << "caxis([P_out, P_in + Offset_colorbar_value_for_particle]);\n";

        oss << "axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * DomainDimensionRatio(1) * L, -1.1 / 2 * DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) * L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * DomainDimensionRatio(3) * L]);\n";
        oss << "pbaspect([DomainDimensionRatio]); hold on\n";

        oss << "Cb = colorbar;\n";
        oss << "Cb.Limits = [P_out, P_in];\n";
        oss << "Cb.Title.String = 'Hydraulic head';\n";

        oss << "hold on\n";

        oss << "newcolors = [];\nnewcolors = rand(N_particles, 1);\n";

        oss << "if (If_video == true)\n";
        oss << "\tparticles_video_object = VideoWriter([currentPath, '/moive_particles.avi']);\n";
        oss << "\tparticles_video_object.FrameRate = 10;\n";
        oss << "\topen(particles_video_object);\n";
        oss << "end\n";

        oss << "figure(2)\nfor i = 0:N_steps\n";
        oss << "\ttitle(['DFN flow (mhfem) and particle tracking; step NO = ', num2str(i)]);\n";
        oss << "\tH5name = []; H5name_2D = [];\n";
        oss << "\tif (i == 0); H5name = [currentPath, '/ParticlePositionResult/ParticlePositionInit_3D.h5']; else; H5name = [currentPath, '/ParticlePositionResult/ParticlePositionBlock', num2str(ceil(double(i) / double(SizeOfDataBlock)), '%010d'), '_3D.h5']; end;\n";
        oss << "\tif (i == 0); H5name_2D = [currentPath, '/ParticlePositionResult/ParticlePositionInit.h5']; else; H5name_2D = [currentPath, '/ParticlePositionResult/ParticlePositionBlock', num2str(ceil(double(i) / double(SizeOfDataBlock)), '%010d'), '.h5']; end;\n";

        oss << "\tS = h5read(H5name, ['/Step_', num2str(i, '%010d')]);\n";
        oss << "\tParticleID = h5read(H5name_2D, ['/ParticleIDAndElementTag_', num2str(i, '%010d')]);\n"; // /
        oss << "\tMatrx3D_pso = NaN(N_particles, 3);\n";
        oss << "\tMatrx3D_pso([ParticleID(:, 1) + 1], :) = S(:, [1 2 3]);\n";

        oss << "\tnewcolors = Matrx3D_pso(:, 3) ./ (2 * L) .* Offset_colorbar_value_for_particle / 2 + P_in + Offset_colorbar_value_for_particle / 2;\n";
        oss << "\tp_s = scatter3(Matrx3D_pso(:, 1), Matrx3D_pso(:, 2), Matrx3D_pso(:, 3), [], newcolors, 'filled'); clear Matrx3D_pso\n\n";

        oss << "\tif(i ~= 0 && If_video == true)\n";
        oss << "\t\tM=getframe(gcf);\n";
        oss << "\t\twriteVideo(particles_video_object,M);\n";
        oss << "\tend\n\n";

        oss << "\tif (i == 0); pause; else; pause(0.01); end;\n";
        oss << "\tif (i ~= N_steps); delete(p_s); end\n";
        oss << "end\n\n";

        oss << "if (If_video == true); close(particles_video_object); end\n";
    }
}; // MatlabPlot
template void cuDFNsys::ParticleTransport<double>::MatlabPlot(const string &mat_key,
                                                              const string &command_key,
                                                              const cuDFNsys::Mesh<double> &mesh,
                                                              const cuDFNsys::MHFEM<double> &fem,
                                                              const double &L, double3 DomainDimensionRatio,
                                                              bool if_python_visualization,
                                                              string PythonName_Without_suffix);
template void cuDFNsys::ParticleTransport<float>::MatlabPlot(const string &mat_key,
                                                             const string &command_key,
                                                             const cuDFNsys::Mesh<float> &mesh,
                                                             const cuDFNsys::MHFEM<float> &fem,
                                                             const float &L, double3 DomainDimensionRatio,
                                                             bool if_python_visualization,
                                                             string PythonName_Without_suffix);

// ====================================================
// NAME:        IdentifyEdgesSharedEle
// DESCRIPTION: IdentifyEdgesSharedEle
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
void cuDFNsys::ParticleTransport<T>::IdentifyEdgesSharedEle(cuDFNsys::Mesh<T> mesh)
{
    this->EdgesSharedEle.resize(mesh.Element3D.size() * 3);

    // create unordered_map for each edges
    typedef struct SharedElements
    {
        uint ElementID[_NumOfSharedEleAtMost] = {0};
        uint LocalEdgeNO[_NumOfSharedEleAtMost] = {0};
    } SharedEle;

    std::unordered_map<pair<size_t, size_t>, SharedEle, PairHash> UMapEdge2Ele;
    UMapEdge2Ele.reserve(mesh.Element3D.size() * 3);

    for (size_t i = 0; i < mesh.Element3D.size(); ++i)
    {
        uint *pts = &mesh.Element3D[i].x;

        for (size_t j = 0; j < 3; ++j)
        {
            size_t node1 = pts[j];
            size_t node2 = pts[(j + 1) % 3];

            pair<size_t, size_t> key_;
            key_.first = node1 < node2 ? node1 : node2;
            key_.second = node1 > node2 ? node1 : node2;

            if (UMapEdge2Ele.find(key_) == UMapEdge2Ele.end())
            {
                UMapEdge2Ele[key_].ElementID[0] = i + 1;
                UMapEdge2Ele[key_].LocalEdgeNO[0] = j;
            }
            else
            {
                for (size_t k = 0; k < _NumOfSharedEleAtMost; ++k)
                {
                    if (UMapEdge2Ele[key_].ElementID[k] == 0)
                    {
                        UMapEdge2Ele[key_].ElementID[k] = i + 1;
                        UMapEdge2Ele[key_].LocalEdgeNO[k] = j;
                        break;
                    };

                    if (k == _NumOfSharedEleAtMost - 1)
                    {
                        throw cuDFNsys::ExceptionsIgnore("Find intersections shared by more than four elements!");
                    }
                }
            };
        }
    };

    // identify
    for (size_t i = 0; i < this->EdgesSharedEle.size(); ++i)
    {
        uint *pts = &(mesh.Element3D[i / 3].x);

        size_t node1 = pts[i % 3];
        size_t node2 = pts[(i % 3 + 1) % 3];

        pair<size_t, size_t> key_;
        key_.first = node1 < node2 ? node1 : node2;
        key_.second = node1 > node2 ? node1 : node2;

        for (size_t k = 0; k < _NumOfSharedEleAtMost; ++k)
        {
            if (UMapEdge2Ele[key_].ElementID[k] != 0)
            {
                this->EdgesSharedEle[i].EleID[k] = UMapEdge2Ele[key_].ElementID[k];
                this->EdgesSharedEle[i].LocalEdgeNO[k] = UMapEdge2Ele[key_].LocalEdgeNO[k];
                this->EdgesSharedEle[i].NumSharedEle++;
            }
        }

        // for (size_t k = 0; k < this->EdgesSharedEle[i].NumSharedEle; ++k)
        // {
        //     cout << UMapEdge2Ele[key_].ElementID[k] << (k == this->EdgesSharedEle[i].NumSharedEle - 1 ? "\n" : ", ");
        //     if(UMapEdge2Ele[key_].ElementID[k] > mesh.Element3D.size())
        //     {
        //         throw cuDFNsys::ExceptionsIgnore("Out of range of element ID\n");
        //     }
        // }
    }

    // output a h5 file
    cuDFNsys::HDF5API h5g;
    string filename_ = "EdgesSharedEle.h5";

    h5g.NewFile(filename_);

    vector<uint> data_EdgesSharedEle(EdgesSharedEle.size() * (1 + _NumOfSharedEleAtMost + _NumOfSharedEleAtMost));

    uint NUMoo = (1 + _NumOfSharedEleAtMost + _NumOfSharedEleAtMost);

    for (uint i = 0; i < EdgesSharedEle.size(); ++i)
    {
        data_EdgesSharedEle[i * NUMoo] = EdgesSharedEle[i].NumSharedEle;
        for (uint j = 0; j < _NumOfSharedEleAtMost; ++j)
        {
            data_EdgesSharedEle[i * NUMoo + j + 1] = EdgesSharedEle[i].EleID[j];
            data_EdgesSharedEle[i * NUMoo + j + 1 + _NumOfSharedEleAtMost] =
                EdgesSharedEle[i].LocalEdgeNO[j];
        }
    }

    uint2 dim_u = make_uint2(EdgesSharedEle.size(), NUMoo);

    h5g.AddDataset(filename_, "N", "data", data_EdgesSharedEle.data(), dim_u);

}; // IdentifyEdgesSharedEle
template void cuDFNsys::ParticleTransport<double>::IdentifyEdgesSharedEle(cuDFNsys::Mesh<double> mesh);
template void cuDFNsys::ParticleTransport<float>::IdentifyEdgesSharedEle(cuDFNsys::Mesh<float> mesh);

// ====================================================
// NAME:        InitilizeParticles
// DESCRIPTION: InitilizeParticles
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
void cuDFNsys::ParticleTransport<T>::InitilizeParticles(const int &NumOfParticles,
                                                        cuDFNsys::Mesh<T> mesh,
                                                        const cuDFNsys::MHFEM<T> &fem,
                                                        const string &Injection_mode)
{
    thrust::host_vector<uint> NumParticlesEachEdge(mesh.InletEdgeNOLen.size(), 0);

    //double sum_speed = abs(fem.VelocityNormalScalarSepEdges.sum());
    //int sizeofinletedge = fem.VelocityNormalScalarSepEdges.size();
    //cout << "sum_speed: " << sum_speed << ", " << sizeofinletedge << endl;
    if (Injection_mode == "Flux-weighted" || Injection_mode == "Resident")
        for (size_t i = 0; i < mesh.InletEdgeNOLen.size(); ++i)
        {
            uint EdgeNO = (uint)mesh.InletEdgeNOLen[i].x;
            T length_ = mesh.InletEdgeNOLen[i].y;

            if (Injection_mode == "Flux-weighted")
            {

                T proportion_ = abs(fem.VelocityNormalScalarSepEdges(EdgeNO - 1, 0)) * length_ / fem.QIn;
                //cout << "edgeno: " << EdgeNO << ", " << proportion_ << endl;
                /// flux-weighted

                NumParticlesEachEdge[i] = (uint)round(proportion_ * NumOfParticles);

                //if (abs(fem.VelocityNormalScalarSepEdges(EdgeNO - 1, 0)) / sum_speed < 1.0 / (double)sizeofinletedge)
                //NumParticlesEachEdge[i] = 0.0;
            }
            else if (Injection_mode == "Resident")
            {
                // uniform
                T proportion_ = abs(length_ / fem.InletLength);
                NumParticlesEachEdge[i] = (uint)round(proportion_ * NumOfParticles);
            }
        }
    else if (Injection_mode == "Point")
    {
        int a;
        srand((unsigned)time(NULL));
        //cout << "NumOfParticles: " << NumOfParticles << endl;
        for (uint i = 0; i < 100; ++i)
        {
            a = rand() % (int)mesh.InletEdgeNOLen.size();
            //printf("a = %d\n", a);
            uint EdgeNO = (uint)mesh.InletEdgeNOLen[a].x;
            if (fem.VelocityNormalScalarSepEdges(EdgeNO - 1, 0) < 0)
                NumParticlesEachEdge[a] = NumOfParticles;
            break;
        }
    }
    else
        throw cuDFNsys::ExceptionsPause("Undefined particle injection mode!\n");

    this->NumParticles = thrust::reduce(thrust::host, NumParticlesEachEdge.begin(), NumParticlesEachEdge.end(), 0);
    T fraction_ = (T)NumOfParticles / (T)this->NumParticles;
    //cout << this->NumParticles << ", " << fraction_ << endl;

    for (size_t i = 0; i < NumParticlesEachEdge.size(); ++i)
        NumParticlesEachEdge[i] = round(NumParticlesEachEdge[i] * fraction_);

    this->NumParticles = thrust::reduce(thrust::host, NumParticlesEachEdge.begin(), NumParticlesEachEdge.end(), 0);
    this->ParticlePlumes.resize(this->NumParticles);

    //cout << "this->NumParticles: " << this->NumParticles << " / " << NumOfParticles << endl;
    //exit(0);

    uint TmpCountParticle = 0;
    uint pARTICLEid_T = 0;
    for (size_t i = 0; i < NumParticlesEachEdge.size(); ++i)
    {
        uint EdgeNO = (uint)mesh.InletEdgeNOLen[i].x;     // from 1
        uint elementID = (EdgeNO - 1) / 3 + 1;            // from 1
        uint FracID = mesh.ElementFracTag[elementID - 1]; // from 0

        uint LocalEdgeNO = (EdgeNO - 1) % 3;

        // cuDFNsys::Vector2<T> position_ = cuDFNsys::MakeVector2((T)0.5f * (mesh.Coordinate2D[elementID - 1].x[LocalEdgeNO] + mesh.Coordinate2D[elementID - 1].x[(LocalEdgeNO + 1) % 3]),
        //                                                        (T)0.5f * (mesh.Coordinate2D[elementID - 1].y[LocalEdgeNO] + mesh.Coordinate2D[elementID - 1].y[(LocalEdgeNO + 1) % 3]));
        // // move the position a little toward to the interior of grid
        // {
        //     cuDFNsys::Vector2<T> inwardNormal = cuDFNsys::MakeVector2(-(mesh.Coordinate2D[elementID - 1].y[(LocalEdgeNO + 1) % 3] - mesh.Coordinate2D[elementID - 1].y[LocalEdgeNO]),
        //                                                               mesh.Coordinate2D[elementID - 1].x[(LocalEdgeNO + 1) % 3] - mesh.Coordinate2D[elementID - 1].x[LocalEdgeNO]);
        //     T norm_inward = sqrt(inwardNormal.x * inwardNormal.x + inwardNormal.y * inwardNormal.y);
        //     inwardNormal.x /= norm_inward;
        //     inwardNormal.y /= norm_inward;
        //     position_.x += (inwardNormal.x * 1e-5);
        //     position_.y += (inwardNormal.y * 1e-5);
        // };
        // cuDFNsys::Particle<T> tmpP;
        // tmpP.ElementID = elementID;
        // tmpP.Position2D = position_;
        // thrust::host_vector<cuDFNsys::Particle<T>> tmpPVEC(NumParticlesEachEdge[i], tmpP);

        thrust::host_vector<cuDFNsys::Particle<T>> tmpPVEC(NumParticlesEachEdge[i]);
        cuDFNsys::Vector2<T> Vec2D;
        Vec2D.x = mesh.Coordinate2D[elementID - 1].x[(LocalEdgeNO + 1) % 3] - mesh.Coordinate2D[elementID - 1].x[LocalEdgeNO];
        Vec2D.y = mesh.Coordinate2D[elementID - 1].y[(LocalEdgeNO + 1) % 3] - mesh.Coordinate2D[elementID - 1].y[LocalEdgeNO];
        T norm_e = sqrt(Vec2D.x * Vec2D.x + Vec2D.y * Vec2D.y);
        Vec2D.x /= norm_e, Vec2D.y /= norm_e;

        norm_e /= NumParticlesEachEdge[i] + 1;

        for (uint j = 1; j <= NumParticlesEachEdge[i]; ++j)
        {
            tmpPVEC[j - 1].ElementID = elementID;
            cuDFNsys::Vector2<T> position_rr;
            position_rr.x = mesh.Coordinate2D[elementID - 1].x[LocalEdgeNO] + Vec2D.x * norm_e * j,
            position_rr.y = mesh.Coordinate2D[elementID - 1].y[LocalEdgeNO] + Vec2D.y * norm_e * j;
            tmpPVEC[j - 1].Position2D = position_rr;

            tmpPVEC[j - 1].ParticleID = pARTICLEid_T;

            pARTICLEid_T++;
        }

        thrust::copy(tmpPVEC.begin(),
                     tmpPVEC.end(),
                     this->ParticlePlumes.begin() + TmpCountParticle);

        TmpCountParticle += NumParticlesEachEdge[i];
    }
}; // InitilizeParticles
template void cuDFNsys::ParticleTransport<double>::InitilizeParticles(const int &NumOfParticles,
                                                                      cuDFNsys::Mesh<double> mesh,
                                                                      const cuDFNsys::MHFEM<double> &fem,
                                                                      const string &Injection_mode);
template void cuDFNsys::ParticleTransport<float>::InitilizeParticles(const int &NumOfParticles,
                                                                     cuDFNsys::Mesh<float> mesh,
                                                                     const cuDFNsys::MHFEM<float> &fem,
                                                                     const string &Injection_mode);

// ====================================================
// NAME:        InitilizeParticles
// DESCRIPTION: InitilizeParticles
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
void cuDFNsys::ParticleTransport<T>::IdentifyNeighborElements(cuDFNsys::Mesh<T> mesh)
{
    this->NeighborEleOfOneEle.resize(mesh.Element3D.size());

    std::unordered_map<uint, cuDFNsys::NeighborEle> SharedNodeEle;
    SharedNodeEle.reserve(mesh.Coordinate3D.size());

    //cout << "t 1\n";
    for (uint i = 0; i < mesh.Element3D.size(); ++i)
    {
        uint *Ptr = &(mesh.Element3D[i].x);

        for (uint j = 0; j < 3; ++j)
        {
            //cout << "Ptr[j]: " << Ptr[j] << ", mesh.Coordinate3D.size() = " << mesh.Coordinate3D.size() << endl;
            uint nodeID = Ptr[j] - 1;
            //cout << "\t\t~1" << endl;
            SharedNodeEle[nodeID].NumNeighborEle++;
            if (SharedNodeEle[nodeID].NumNeighborEle > _NumOfNeighborEleAtMost)
            {
                string AS = "The number of elements sharing the same node exceeds the allowed value\n";
                throw cuDFNsys::ExceptionsPause(AS);
            }
            //cout << "\t\t~2" << endl;
            SharedNodeEle[nodeID].EleID[SharedNodeEle[nodeID].NumNeighborEle - 1] = i + 1;
            //cout << "\t\t~3" << endl;
            //cout << "\t" << SharedNodeEle[nodeID].NumNeighborEle << ", " << SharedNodeEle[nodeID].EleID[SharedNodeEle[nodeID].NumNeighborEle - 1] << "\n";
        }
    }

    for (uint i = 0; i < mesh.Element3D.size(); ++i)
    {
        uint *Ptr = &(mesh.Element3D[i].x);

        std::set<uint> KL;

        for (uint j = 0; j < 3; ++j)
        {
            for (uint k = 0; k < SharedNodeEle[Ptr[j] - 1].NumNeighborEle; ++k)
            {
                if (SharedNodeEle[Ptr[j] - 1].EleID[k] != i + 1)
                    KL.insert(SharedNodeEle[Ptr[j] - 1].EleID[k]);
            }
        }

        this->NeighborEleOfOneEle[i].NumNeighborEle = KL.size();

        if (KL.size() > _NumOfNeighborEleAtMost)
        {
            string AS = "The number of neighbering elements of an element exceeds the allowed value\n";
            throw cuDFNsys::ExceptionsPause(AS);
        }
        //cout << "element ID " << i + 1 << ": ";
        uint oo = 0;
        for (std::set<uint>::iterator j = KL.begin(); j != KL.end(); ++j)
        {
            this->NeighborEleOfOneEle[i].EleID[oo] = *j;
            oo++;
            // cout << *j << (oo == KL.size() ? "\n" : ", ");
        }
    }
}; // IdentifyNeighborElements
template void cuDFNsys::ParticleTransport<double>::IdentifyNeighborElements(cuDFNsys::Mesh<double> mesh);
template void cuDFNsys::ParticleTransport<float>::IdentifyNeighborElements(cuDFNsys::Mesh<float> mesh);
