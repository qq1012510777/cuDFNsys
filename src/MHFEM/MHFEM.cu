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

#include "MHFEM/MHFEM.cuh"

// ====================================================
// NAME:        MHFEM
// DESCRIPTION: MHFEM constructor
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
template <typename T>
cuDFNsys::MHFEM<T>::MHFEM(
    const cuDFNsys::Mesh<T> &mesh,
    const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs, const T &inlet_p_,
    const T &outlet_p_, const int &dir_, const T &L,
    double3 DomainDimensionRatio, bool if_CPU, int Nproc, bool if_periodic,
    T muOverRhoG, T constant_qq)
{
    cout << "\tMHFEM ing\n";
    Eigen::setNbThreads(1);
    this->Dir = dir_;
    this->InletP = inlet_p_;
    this->OutletP = outlet_p_;
    this->IfPeriodic = if_periodic;
    this->MuOverRhoG = muOverRhoG;
    // cout << "constant_qq: " << constant_qq << endl;
    this->ConsTq = constant_qq;
    // cout << " this->ConsTq: " << this->ConsTq << endl;

    if (this->IfPeriodic &&
        (abs(mesh.InletTraceLength - mesh.OutletTraceLength) > 1e-5))
        throw cuDFNsys::ExceptionsPause("The DFN is not spatially periodic");

    this->Implementation(mesh, Fracs, if_CPU, Nproc);

    this->QError = (abs(QIn - QOut) / (QIn > QOut ? QIn : QOut)) * 100.0f;
    double *Rt = &DomainDimensionRatio.x;

    //--------------------the mean and max velocity of elements
    //--------------------the mean and max velocity of elements
    //--------------------the mean and max velocity of elements
    T *velocity_center_grid = new T[mesh.Element3D.size() * 3];
    if (velocity_center_grid == NULL)
    {
        string AS = "Alloc error in MHFEM::MHFEM\n";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    Eigen::VectorXd NormOfVelocity;
    NormOfVelocity.resize(mesh.Element3D.size(), 1);

    T volume_frac_sum = 0;
    this->MaxVelocity = 0;

    T sum_head_times_area = 0;

    for (uint i = 0; i < mesh.Element3D.size(); ++i)
    {
        uint3 EdgeNO = make_uint3((i + 1) * 3 - 3, (i + 1) * 3 - 2,
                                  (i + 1) * 3 - 1); // from 0
        cuDFNsys::Vector3<T> Velocity_ = cuDFNsys ::MakeVector3(
            (T)this->VelocityNormalScalarSepEdges(EdgeNO.x, 0),
            (T)this->VelocityNormalScalarSepEdges(EdgeNO.y, 0),
            (T)this->VelocityNormalScalarSepEdges(EdgeNO.z, 0));
        cuDFNsys::Vector2<T> Vertexes[3];
        Vertexes[0] = cuDFNsys::MakeVector2(mesh.Coordinate2D[i].x[0],
                                            mesh.Coordinate2D[i].y[0]);
        Vertexes[1] = cuDFNsys::MakeVector2(mesh.Coordinate2D[i].x[1],
                                            mesh.Coordinate2D[i].y[1]);
        Vertexes[2] = cuDFNsys::MakeVector2(mesh.Coordinate2D[i].x[2],
                                            mesh.Coordinate2D[i].y[2]);

        cuDFNsys::Vector2<T> Center_p = cuDFNsys::MakeVector2(
            1.0f / 3.0f * (Vertexes[0].x + Vertexes[1].x + Vertexes[2].x),
            1.0f / 3.0f * (Vertexes[0].y + Vertexes[1].y + Vertexes[2].y));

        cuDFNsys::Vector2<T> velocity_p =
            cuDFNsys::ReconstructVelocityGrid<T>(Center_p, Vertexes, Velocity_);

        cuDFNsys::Vector3<T> velocity_p_3D =
            cuDFNsys::MakeVector3(velocity_p.x, velocity_p.y, (T)0);

        T R_mat[3][3];
        cuDFNsys::Fracture<T> FII = Fracs[mesh.ElementFracTag[i]];
        FII.RoationMatrix(R_mat, 23);
        velocity_p_3D =
            cuDFNsys::ProductSquare3Vector3<T>(R_mat, velocity_p_3D);

        velocity_center_grid[i] = velocity_p_3D.x;
        velocity_center_grid[i + mesh.Element3D.size()] = velocity_p_3D.y;
        velocity_center_grid[i + 2 * mesh.Element3D.size()] = velocity_p_3D.z;

        T b_f = pow(FII.Conductivity * 12, 1.0 / 3.0);

        T area_this_element =
            cuDFNsys::Triangle2DArea<T>(Vertexes[0], Vertexes[1], Vertexes[2]);

        volume_frac_sum += (area_this_element * b_f);

        sum_head_times_area += (PressureEles(i, 0) * area_this_element);

        NormOfVelocity(i, 0) =
            pow(velocity_p.x * velocity_p.x + velocity_p.y * velocity_p.y,
                0.5) /
            b_f; // q / b

        this->MaxVelocity =
            (NormOfVelocity(i, 0) > this->MaxVelocity ? NormOfVelocity(i, 0)
                                                      : this->MaxVelocity);

        NormOfVelocity(i, 0) *= (area_this_element * b_f);
    };

    if (!this->IfPeriodic) // K = -Q * (mu / (rho * g)) / (Area_inlet * pressure_gradient)
    {
        this->Permeability =
            -0.5f * (QOut + QIn) * this->MuOverRhoG /
            (L * L * Rt[(this->Dir + 1) % 3] * Rt[(this->Dir + 2) % 3]) /
            ((outlet_p_ - inlet_p_) / (L * Rt[this->Dir]));
        // cout << "Original Permeability (non-periodic): " << this->Permeability
        //      << endl;
        // cout << "Q / L^2 = "
        //      << 0.5f * (QOut + QIn) /
        //             (L * L * Rt[(this->Dir + 1) % 3] * Rt[(this->Dir + 2) % 3])
        //      << endl;
        // cout << "Integrated Darcian V = "
        //      << NormOfVelocity.sum() /
        //             (L * L * L * DomainDimensionRatio.x *
        //              DomainDimensionRatio.y * DomainDimensionRatio.z)
        //      << endl;
        // this->Permeability =
        //     -NormOfVelocity.sum() /
        //     (L * L * L * DomainDimensionRatio.x * DomainDimensionRatio.y *
        //      DomainDimensionRatio.z) *
        //     this->MuOverRhoG /
        //     (sum_head_times_area /
        //      (L * L * L * DomainDimensionRatio.x * DomainDimensionRatio.y *
        //       DomainDimensionRatio.z));
    }

    if (this->IfPeriodic) // K = -Q * (mu / (rho * g)) / (Area_inlet * pressure_gradient)
    {
        this->Permeability =
            -(0.5f * (QOut + QIn)) * this->MuOverRhoG /
            (sum_head_times_area /
             (L * L * L * DomainDimensionRatio.x * DomainDimensionRatio.y *
              DomainDimensionRatio.z)) /
            (L * L * Rt[(this->Dir + 1) % 3] * Rt[(this->Dir + 2) % 3]);
        // this->Permeability =
        //     -NormOfVelocity.sum() /
        //     (L * L * L * DomainDimensionRatio.x * DomainDimensionRatio.y *
        //      DomainDimensionRatio.z) *
        //     this->MuOverRhoG /
        //     (sum_head_times_area /
        //      (L * L * L * DomainDimensionRatio.x * DomainDimensionRatio.y *
        //       DomainDimensionRatio.z));
    }

    // now it is the maximum velocity
    //this->MaxVelocity = NormOfVelocity.maxCoeff(); //NormOfVelocity.sum() / mesh.Element3D.size();
    this->MeanVelocity =
        NormOfVelocity.sum() / volume_frac_sum; // % mesh.Element3D.size();

    delete[] velocity_center_grid;
    velocity_center_grid = NULL;

}; // MHFEM
template cuDFNsys::MHFEM<double>::MHFEM(
    const cuDFNsys::Mesh<double> &mesh,
    const thrust::host_vector<cuDFNsys::Fracture<double>> &Fracs,
    const double &inlet_p_, const double &outlet_p_, const int &dir_,
    const double &L, double3 DomainDimensionRatio, bool if_CPU, int Nproc,
    bool if_periodic, double muOverRhoG, double constant_qq);
template cuDFNsys::MHFEM<float>::MHFEM(
    const cuDFNsys::Mesh<float> &mesh,
    const thrust::host_vector<cuDFNsys::Fracture<float>> &Fracs,
    const float &inlet_p_, const float &outlet_p_, const int &dir_,
    const float &L, double3 DomainDimensionRatio, bool if_CPU, int Nproc,
    bool if_periodic, float muOverRhoG, float constant_qq);

// ====================================================
// NAME:        MatlabPlot
// DESCRIPTION: plot mhfem result
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
template <typename T>
double2 cuDFNsys::MHFEM<T>::MatlabPlot(
    const string &mat_key, const string &command_key,
    thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
    const cuDFNsys::Mesh<T> &mesh, const T &L, bool if_python_visualization,
    string PythonName_Without_suffix, double3 DomainDimensionRatio)
{
    //cuDFNsys::MatlabAPI M1;

    cuDFNsys::HDF5API h5gg;

    h5gg.NewFile(mat_key);

    size_t node_num = mesh.Coordinate3D.size();

    T *ptr_coordinates_3D;
    ptr_coordinates_3D = new T[node_num * 3];

    if (ptr_coordinates_3D == NULL)
    {
        string AS = "Alloc error in MHFEM::MatlabPlot\n";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    for (size_t i = 0; i < node_num; ++i)
    {
        ptr_coordinates_3D[i] = mesh.Coordinate3D[i].x;
        ptr_coordinates_3D[i + node_num] = mesh.Coordinate3D[i].y;
        ptr_coordinates_3D[i + node_num * 2] = mesh.Coordinate3D[i].z;
    }
    // M1.WriteMat(mat_key, "w", node_num * 3,
    //             node_num, 3, ptr_coordinates_3D, "coordinate_3D");
    uint2 dim_f = make_uint2(3, node_num);
    h5gg.AddDataset(mat_key, "N", "coordinate_3D", ptr_coordinates_3D, dim_f);

    delete[] ptr_coordinates_3D;
    ptr_coordinates_3D = NULL;

    dim_f = make_uint2(1, 1);
    T yu[1] = {L};
    h5gg.AddDataset(mat_key, "N", "L_m", yu, dim_f);

    yu[0] = this->InletP;
    h5gg.AddDataset(mat_key, "N", "InletP", yu, dim_f);
    yu[0] = this->OutletP;
    h5gg.AddDataset(mat_key, "N", "OutletP", yu, dim_f);

    uint2 dim_ds = make_uint2(3, 1);
    double DomainDimensionRatio_d[3] = {
        DomainDimensionRatio.x, DomainDimensionRatio.y, DomainDimensionRatio.z};
    h5gg.AddDataset(mat_key, "N", "DomainDimensionRatio",
                    DomainDimensionRatio_d, dim_ds);

    //---------------------
    size_t ele_num = mesh.Element3D.size();
    T *ptr_element_3D = new T[ele_num * 3];
    if (ptr_element_3D == NULL)
    {
        string AS = "Alloc error in MHFEM::MatlabPlot\n";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    for (size_t i = 0; i < ele_num; ++i)
    {
        ptr_element_3D[i] = mesh.Element3D[i].x;
        ptr_element_3D[i + ele_num] = mesh.Element3D[i].y;
        ptr_element_3D[i + ele_num * 2] = mesh.Element3D[i].z;
    }
    // M1.WriteMat(mat_key, "u", ele_num * 3,
    //             ele_num, 3, ptr_element_3D, "element_3D");
    dim_f = make_uint2(3, ele_num);
    h5gg.AddDataset(mat_key, "N", "element_3D", ptr_element_3D, dim_f);

    delete[] ptr_element_3D;
    ptr_element_3D = NULL;
    //--------------------------

    T *pressure_ELEs = new T[ele_num];
    if (pressure_ELEs == NULL)
    {
        string AS = "Alloc error in MHFEM::MatlabPlot\n";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    std::copy(this->PressureEles.data(), this->PressureEles.data() + ele_num,
              pressure_ELEs);

    // M1.WriteMat(mat_key, "u", ele_num,
    //             ele_num, 1, pressure_ELEs,
    //             "pressure_eles");
    dim_f = make_uint2(1, ele_num);
    h5gg.AddDataset(mat_key, "N", "pressure_eles", pressure_ELEs, dim_f);

    delete[] pressure_ELEs;
    pressure_ELEs = NULL;

    //------------------------------
    //T *normal_veloc = new T[ele_num * 3];
    //if (normal_veloc == NULL)
    //{
    //    string AS = "Alloc error in MHFEM::MatlabPlot\n";
    //    throw cuDFNsys::ExceptionsPause(AS);
    //}
    //
    //std::copy(this->VelocityNormalScalarSepEdges.data(),
    //          this->VelocityNormalScalarSepEdges.data() + ele_num * 3,
    //          normal_veloc);
    //
    //M1.WriteMat(mat_key, "u", ele_num * 3,
    //            ele_num * 3, 1, normal_veloc,
    //            "normal_velocity");
    //
    //delete[] normal_veloc;
    //normal_veloc = NULL;

    //----------------
    T *velocity_center_grid = new T[mesh.Element3D.size() * 3];
    if (velocity_center_grid == NULL)
    {
        string AS = "Alloc error in MHFEM::MatlabPlot\n";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    // Eigen::VectorXd NormOfVelocity;
    // NormOfVelocity.resize(mesh.Element3D.size(), 1);

    // T area_sum = 0;
    // double maxV = 0;
    for (uint i = 0; i < mesh.Element3D.size(); ++i)
    {
        uint3 EdgeNO = make_uint3((i + 1) * 3 - 3, (i + 1) * 3 - 2,
                                  (i + 1) * 3 - 1); // from 0
        cuDFNsys::Vector3<T> Velocity_ = cuDFNsys ::MakeVector3(
            (T)this->VelocityNormalScalarSepEdges(EdgeNO.x, 0),
            (T)this->VelocityNormalScalarSepEdges(EdgeNO.y, 0),
            (T)this->VelocityNormalScalarSepEdges(EdgeNO.z, 0));
        cuDFNsys::Vector2<T> Vertexes[3];
        Vertexes[0] = cuDFNsys::MakeVector2(mesh.Coordinate2D[i].x[0],
                                            mesh.Coordinate2D[i].y[0]);
        Vertexes[1] = cuDFNsys::MakeVector2(mesh.Coordinate2D[i].x[1],
                                            mesh.Coordinate2D[i].y[1]);
        Vertexes[2] = cuDFNsys::MakeVector2(mesh.Coordinate2D[i].x[2],
                                            mesh.Coordinate2D[i].y[2]);

        cuDFNsys::Vector2<T> Center_p = cuDFNsys::MakeVector2(
            1.0f / 3.0f * (Vertexes[0].x + Vertexes[1].x + Vertexes[2].x),
            1.0f / 3.0f * (Vertexes[0].y + Vertexes[1].y + Vertexes[2].y));

        cuDFNsys::Vector2<T> velocity_p =
            cuDFNsys::ReconstructVelocityGrid<T>(Center_p, Vertexes, Velocity_);

        cuDFNsys::Vector3<T> velocity_p_3D =
            cuDFNsys::MakeVector3(velocity_p.x, velocity_p.y, (T)0);

        T R_mat[3][3];
        Fracs[mesh.ElementFracTag[i]].RoationMatrix(R_mat, 23);
        velocity_p_3D =
            cuDFNsys::ProductSquare3Vector3<T>(R_mat, velocity_p_3D);

        velocity_center_grid[i] = velocity_p_3D.x;
        velocity_center_grid[i + mesh.Element3D.size()] = velocity_p_3D.y;
        velocity_center_grid[i + 2 * mesh.Element3D.size()] = velocity_p_3D.z;

        // T area_this_element =
        //     cuDFNsys::Triangle2DArea<T>(Vertexes[0], Vertexes[1], Vertexes[2]);
        // area_sum += area_this_element;
        // NormOfVelocity(i, 0) =
        //     pow(velocity_p.x * velocity_p.x + velocity_p.y * velocity_p.y,
        //         0.5) /
        //     pow(Fracs[mesh.ElementFracTag[i]].Conductivity * 12, 1.0 / 3.0);
        // maxV = (NormOfVelocity(i, 0) > maxV ? NormOfVelocity(i, 0) : maxV);
        // NormOfVelocity(i, 0) *= area_this_element;
    };

    // now it is the maximum velocity
    //double maxV = NormOfVelocity.maxCoeff(); //NormOfVelocity.sum() / mesh.Element3D.size();
    // double meanV = NormOfVelocity.sum() / area_sum;

    h5gg.AddDataset<double>(mat_key, "N", "maxV", &this->MaxVelocity,
                            make_uint2(1, 0));
    h5gg.AddDataset<double>(mat_key, "N", "meanV", &this->MeanVelocity,
                            make_uint2(1, 0));

    double2 JKLDSF = make_double2(this->MeanVelocity, this->MaxVelocity);

    // M1.WriteMat(mat_key, "u", mesh.Element3D.size() * 3,
    //             mesh.Element3D.size(), 3, velocity_center_grid,
    //             "velocity_center_grid");
    dim_f = make_uint2(3, mesh.Element3D.size());
    h5gg.AddDataset(mat_key, "N", "velocity_center_grid", velocity_center_grid,
                    dim_f);

    delete[] velocity_center_grid;
    velocity_center_grid = NULL;

    //-----------
    std::vector<uint> ElementFracTag;

    ElementFracTag.insert(ElementFracTag.end(), mesh.ElementFracTag.begin(),
                          mesh.ElementFracTag.end());

    uint2 dim_fqsc = make_uint2(1, ElementFracTag.size());
    h5gg.AddDataset(mat_key, "N", "ElementFracTag", ElementFracTag.data(),
                    dim_fqsc);

    //
    std::vector<T> ElementAperture(ElementFracTag.size());
    for (uint i = 0; i < ElementAperture.size(); ++i)
        ElementAperture[i] =
            pow(Fracs[ElementFracTag[i]].Conductivity * 12.0, 1.0 / 3.0);
    h5gg.AddDataset(mat_key, "N", "ElementAperture", ElementAperture.data(),
                    dim_fqsc);

    if (1)
    {
        size_t NUM_ele = mesh.Element3D.size();

        T *edge_attri = new T[NUM_ele * 6];
        if (edge_attri == NULL)
        {
            string AS = "Alloc error in Mesh::MatlabPlot\n";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        for (size_t i = 0; i < NUM_ele; ++i)
        {
            edge_attri[i] = mesh.Element3D[i].x;
            edge_attri[i + NUM_ele] = mesh.Element3D[i].y;
            edge_attri[i + 2 * NUM_ele] = mesh.Element3D[i].z;

            edge_attri[i + 3 * NUM_ele] = mesh.EdgeAttri[i].e[0];
            edge_attri[i + 4 * NUM_ele] = mesh.EdgeAttri[i].e[1];
            edge_attri[i + 5 * NUM_ele] = mesh.EdgeAttri[i].e[2];
        }
        //M1.WriteMat(mat_key, "u", NUM_ele * 6,
        //            NUM_ele, 6, edge_attri, "edge_attri");

        uint2 dim_fwww = make_uint2(6, NUM_ele);
        h5gg.AddDataset(mat_key, "N", "edge_attri", edge_attri, dim_fwww);

        delete[] edge_attri;
        edge_attri = NULL;
    }

    //-------------------
    if (if_python_visualization)
    {
        std::ofstream oss(PythonName_Without_suffix + ".py", ios::out);
        oss << "import h5py\n";
        oss << "import numpy as np\n";
        oss << "from mayavi import mlab as ML\n";
        oss << "from numpy import linalg as LA\n";

        oss << "f = h5py.File('" << mat_key << "')\n";
        oss << "coordinate_3D = np.array(f['coordinate_3D'][:])\n";
        oss << "element_3D = np.array(f['element_3D'][:], dtype=int)\n";
        oss << "InletP = np.array(f['InletP'][:])\n";
        oss << "OutletP = np.array(f['OutletP'][:])\n";
        oss << "L_m = f['L_m'][:][0]\n";
        oss << "DomainDimensionRatio = f['DomainDimensionRatio'][:]\n";
        oss << "pressure_eles = np.array(f['pressure_eles'][:])\n";
        oss << "velocity_center_grid = "
               "np.array(f['velocity_center_grid'][:])\n";
        oss << "vec_norm = LA.norm(velocity_center_grid, axis=0)\n";

        oss << "f.close()\n";
        oss << "mesh = ML.triangular_mesh(coordinate_3D[0, :], "
               "coordinate_3D[1, :], coordinate_3D[2, :], "
               "np.transpose(element_3D-1), representation='wireframe', "
               "color=(0, 0, 0), line_width=1.0)\n";

        if (!this->IfPeriodic)
            oss << "mesh.mlab_source.dataset.cell_data.scalars = "
                   "np.transpose(pressure_eles)\n";
        else
            oss << "mesh.mlab_source.dataset.cell_data.scalars = "
                   "np.transpose(vec_norm)\n";

        oss << "mesh.mlab_source.dataset.cell_data.scalars.name = 'Cell "
               "data'\n";
        oss << "mesh.mlab_source.update()\n";
        oss << "mesh.parent.update()\n";
        oss << "mesh2 = ML.pipeline.set_active_attribute(mesh, "
               "cell_scalars='Cell data')\n";
        oss << "s2 = ML.pipeline.surface(mesh2, colormap='rainbow', "
               "opacity=0.8)\n";
        oss << "ML.outline(extent=[-0.5 * DomainDimensionRatio[0] * L_m, 0.5 * "
               "DomainDimensionRatio[0] * L_m, -0.5 * DomainDimensionRatio[1] "
               "* L_m, 0.5 * DomainDimensionRatio[1] * L_m, -0.5 * "
               "DomainDimensionRatio[2] * L_m, 0.5 * DomainDimensionRatio[2] * "
               "L_m])\n";
        oss << "ML.axes()\n";
        if (!this->IfPeriodic)
        {
            oss << "s2.module_manager.scalar_lut_manager.data_range = "
                   "np.array([OutletP[0], InletP[0]])\n";
            oss << "ML.colorbar(object=s2, orientation='vertical', "
                   "title='Hydraulic head [L]')\n";
        }
        else
        {
            oss << "s2.module_manager.scalar_lut_manager.data_range = "
                   "np.array([np.min(vec_norm), 0.01 * np.max(vec_norm)])\n";
            oss << "ML.colorbar(object=s2, orientation='vertical', "
                   "title='Norm of velocity [LL/T]')\n";
        }
        oss << "ML.xlabel('x (m)')\n";
        oss << "ML.ylabel('y (m)')\n";
        oss << "ML.zlabel('z (m)')\n";

        oss << "CenterELE = np.zeros([3, element_3D.shape[1]])\n";
        oss << "CenterELE[0, :] = 1.0 / 3.0 * (coordinate_3D[0, element_3D[0, "
               ":]-1] + coordinate_3D[0, element_3D[1, :]-1] + "
               "coordinate_3D[0, element_3D[2, :]-1])\n";
        oss << "CenterELE[1, :] = 1.0 / 3.0 * (coordinate_3D[1, element_3D[0, "
               ":]-1] + coordinate_3D[1, element_3D[1, :]-1] + "
               "coordinate_3D[1, element_3D[2, :]-1])\n";
        oss << "CenterELE[2, :] = 1.0 / 3.0 * (coordinate_3D[2, element_3D[0, "
               ":]-1] + coordinate_3D[2, element_3D[1, :]-1] + "
               "coordinate_3D[2, element_3D[2, :]-1])\n";
        oss << "ML.quiver3d(CenterELE[0, :], CenterELE[1, :], CenterELE[2, :], "
               "velocity_center_grid[0, :], velocity_center_grid[1, :], "
               "velocity_center_grid[2, :])\n";

        oss << "ML.show()\n";
        oss.close();
    }

    //-----------------
    if (command_key != "N")
    {
        std::ofstream oss(command_key, ios::out);
        oss << "clc;\nclose all;\nclear all;\n";
        //oss << "load('" << mat_key << "');\n";
        // oss << "L = 0.5 * " << L << ";\n";
        oss << "currentPath = fileparts(mfilename('fullpath'));\n";
        oss << "coordinate_3D = h5read([currentPath, '/" << mat_key
            << "'], '/coordinate_3D');\n";
        oss << "element_3D = h5read([currentPath, '/" << mat_key
            << "'], '/element_3D');\n";
        oss << "velocity_center_grid = h5read([currentPath, '/" << mat_key
            << "'], '/velocity_center_grid');\n";
        oss << "pressure_eles = h5read([currentPath, '/" << mat_key
            << "'], '/pressure_eles');\n";
        oss << "VelocityNorm = vecnorm(velocity_center_grid');\n\n";

        oss << "L = h5read([currentPath, '/" << mat_key << "'], '/L_m');\n";
        oss << "DomainDimensionRatio = h5read([currentPath, '/" << mat_key
            << "'], '/DomainDimensionRatio');\n";
        oss << "cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, "
               "-L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, "
               "-L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, "
               "-L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, "
               "L, -L; -L,L, L];\n";
        oss << "cube_frame(:, 1) = 0.5 .* cube_frame(:, 1) .* "
               "DomainDimensionRatio(1); ";
        oss << "cube_frame(:, 2) = 0.5 .* cube_frame(:, 2) .* "
               "DomainDimensionRatio(2); ";
        oss << "cube_frame(:, 3) = 0.5 .* cube_frame(:, 3) .* "
               "DomainDimensionRatio(3);\n";

        oss << "figure(1); view(3); title('DFN flow (mhfem)'); xlabel('x "
               "(m)'); ylabel('y (m)'); zlabel('z (m)'); hold on\n";
        oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 "
               "10 11 12; 13 14 15 16], 'FaceVertexCData', "
               "zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', "
               "'EdgeAlpha', 1, 'facealpha', 0); hold on\n";
        oss << endl;
        oss << "axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * "
               "DomainDimensionRatio(1) * L, -1.1 / 2 * "
               "DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) "
               "* L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * "
               "DomainDimensionRatio(3) * L]);\n";
        oss << "pbaspect([DomainDimensionRatio]); hold on\n";

        //oss << "centers_ele = zeros(size(element_3D, 1), 3);\n";
        //oss << "centers_ele(:, 1) = 0.5 * (coordinate_3D(element_3D(:, 1), 1) + coordinate_3D(element_3D(:, 2), 1) + coordinate_3D(element_3D(:, 3), 1));\n";
        //oss << "centers_ele(:, 2) = 0.5 * (coordinate_3D(element_3D(:, 1), 2) + coordinate_3D(element_3D(:, 2), 2) + coordinate_3D(element_3D(:, 3), 2));\n";
        //oss << "centers_ele(:, 3) = 0.5 * (coordinate_3D(element_3D(:, 1), 3) + coordinate_3D(element_3D(:, 2), 3) + coordinate_3D(element_3D(:, 3), 3));\n";
        //oss << endl;

        //oss << "center_s_edge = zeros(size(shared_edge_NO, 1), 3);\n";
        //oss << "center_s_edge(:, 1) = 0.5 * (coordinate_3D(shared_edge_NO(:, 1), 1) + coordinate_3D(shared_edge_NO(:, 2), 1));\n";
        //oss << "center_s_edge(:, 2) = 0.5 * (coordinate_3D(shared_edge_NO(:, 1), 2) + coordinate_3D(shared_edge_NO(:, 2), 2));\n";
        //oss << "center_s_edge(:, 3) = 0.5 * (coordinate_3D(shared_edge_NO(:, 1), 3) + coordinate_3D(shared_edge_NO(:, 2), 3));\n";

        //oss << "tool_ = [centers_ele; center_s_edge];\n";
        //oss << "pressure_vert = griddata(tool_(:, 1), tool_(:, 2), tool_(:, 3), [pressure_eles; pressure_shared_edge], coordinate_3D(:, 1), coordinate_3D(:, 2), coordinate_3D(:, 3), 'nearest');\n";
        if (!this->IfPeriodic)
        {
            oss << "patch('Vertices', coordinate_3D, 'Faces', element_3D, "
                   "'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', "
                   "'EdgeAlpha', 0.9, 'facealpha', 1); view(3); hold "
                   "on\n";
            oss << "Cb = colorbar;\n";
            // oss << "Cb.Limits = [" << this->OutletP << ", " << this->InletP << "];\n";
            oss << "Cb.Title.String = 'Hydraulic head [L]';\n";
            oss << "caxis([" << this->OutletP << ", " << this->InletP
                << "]);\n\n";
        }
        else
        {
            oss << "patch('Vertices', coordinate_3D, 'Faces', element_3D, "
                   "'FaceVertexCData', VelocityNorm', 'FaceColor', 'flat', "
                   "'EdgeAlpha', 0.9, 'facealpha', 1); view(3); hold "
                   "on\n";
            // s
            oss << "Cb = colorbar;\n";
            //oss << "Cb.Limits = [min(VelocityNorm), max(VelocityNorm)];\n";
            oss << "Cb.Title.String = 'Norm of velocity [LL/T]';\n";
            oss << "caxis([min(VelocityNorm), 0.01 * max(VelocityNorm)]);\n\n";
        }

        //oss << "element_3D_re = zeros(size(element_3D, 1) * 3, 3);\n";
        //oss << "element_3D_re([1:3:end], :) = element_3D;\n";
        //oss << "element_3D_re([2:3:end], :) = element_3D(:, [2 3 1]);\n";
        //oss << "element_3D_re([3:3:end], :) = element_3D(:, [3 1 2]);\n\n";
        //
        //oss << "Tri_plane_normal = cross([coordinate_3D(element_3D_re(:, 2), :) - coordinate_3D(element_3D_re(:, 1), :)], [coordinate_3D(element_3D_re(:, 3), :) - coordinate_3D(element_3D_re(:, 2), :)]);\n";
        //oss << "Outward_normal = cross([coordinate_3D(element_3D_re(:, 2), :) - coordinate_3D(element_3D_re(:, 1), :)], Tri_plane_normal);\n";
        //oss << "Outward_normal = Outward_normal ./ norm(Outward_normal);\n";
        //oss << "Outward_normal = Outward_normal .* normal_velocity;\n";
        //
        //oss << "center_sep_edge = 0.5 * [coordinate_3D(element_3D_re(:, 2), :) + coordinate_3D(element_3D_re(:, 1), :)];\n";
        //
        //oss << "quiver3(center_sep_edge(:, 1),center_sep_edge(:, 2),center_sep_edge(:, 3), Outward_normal(:, 1),Outward_normal(:, 2),Outward_normal(:, 3), 4, 'LineWidth', 1.5, 'color', 'r');\n";

        oss << "CenterELE = zeros(size(element_3D, 1), 3);\n";
        oss << "CenterELE(:, 1) = 1/3 * (coordinate_3D(element_3D(:, 1), 1) + "
               "coordinate_3D(element_3D(:, 2), 1) + "
               "coordinate_3D(element_3D(:, 3), 1));\n";
        oss << "CenterELE(:, 2) = 1/3 * (coordinate_3D(element_3D(:, 1), 2) + "
               "coordinate_3D(element_3D(:, 2), 2) + "
               "coordinate_3D(element_3D(:, 3), 2));\n";
        oss << "CenterELE(:, 3) = 1/3 * (coordinate_3D(element_3D(:, 1), 3) + "
               "coordinate_3D(element_3D(:, 2), 3) + "
               "coordinate_3D(element_3D(:, 3), 3));\n";
        oss << "quiver3(CenterELE(:, 1), CenterELE(:, 2), CenterELE(:, 3), "
               "velocity_center_grid(:, 1),velocity_center_grid(:, "
               "2),velocity_center_grid(:, 3), 4, 'LineWidth', 1.5, 'color', "
               "'r');\n";
        oss << "ElementAperture = h5read([currentPath, '/" << mat_key
            << "'], '/ElementAperture');\n";
        oss << "VelocityNorm = [vecnorm(velocity_center_grid')]' ./ "
               "ElementAperture;\n";
        oss << "len1=[vecnorm([coordinate_3D(element_3D(:, 1), :) - "
               "coordinate_3D(element_3D(:, 2), :)]')]';\n";
        oss << "len2=[vecnorm([coordinate_3D(element_3D(:, 3), :) - "
               "coordinate_3D(element_3D(:, 2), :)]')]';\n";
        oss << "len3=[vecnorm([coordinate_3D(element_3D(:, 1), :) - "
               "coordinate_3D(element_3D(:, 3), :)]')]';\n";
        oss << "P_ss = (len1+len2+len3)*0.5;\n";
        oss << "Area_ss=(P_ss .* (P_ss-len1) .* (P_ss-len2) .* (P_ss-len3)) .^ 0.5;\n";
        oss << "meanFractureVelocity = sum(VelocityNorm .* Area_ss .* "
               "ElementAperture) ./ "
               "(sum(Area_ss .* ElementAperture))\n";
        oss.close();
    }
    return JKLDSF;
}; // MHFEM
template double2 cuDFNsys::MHFEM<double>::MatlabPlot(
    const string &mat_key, const string &command_key,
    thrust::host_vector<cuDFNsys::Fracture<double>> Fracs,
    const cuDFNsys::Mesh<double> &mesh, const double &L,
    bool if_python_visualization, string PythonName_Without_suffix,
    double3 DomainDimensionRatio);
template double2 cuDFNsys::MHFEM<float>::MatlabPlot(
    const string &mat_key, const string &command_key,
    thrust::host_vector<cuDFNsys::Fracture<float>> Fracs,
    const cuDFNsys::Mesh<float> &mesh, const float &L,
    bool if_python_visualization, string PythonName_Without_suffix,
    double3 DomainDimensionRatio);

// ====================================================
// NAME:        Implementation
// DESCRIPTION: Implementation of mhfem
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
template <typename T>
void cuDFNsys::MHFEM<T>::Implementation(
    const cuDFNsys::Mesh<T> &mesh,
    const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs, bool if_CPU,
    int Nproc)
{
    size_t NUM_sep_edges = mesh.Element3D.size() * 3,
           NUM_eles = mesh.Element3D.size(),
           NUM_glob_interior_edges = mesh.NumInteriorEdges,
           NUM_INLET_EDGES = mesh.NumInletEdges,
           NUM_OUTLET_EDGES = mesh.NumOutletEdges,
           NUM_NEUMMANN_EDGES = mesh.NumNeumannEdges;

    size_t NumEdges = NUM_glob_interior_edges + NUM_INLET_EDGES +
                      NUM_OUTLET_EDGES + NUM_NEUMMANN_EDGES;

    //cout << "NUM_sep_edges: " << NUM_sep_edges << endl;
    //cout << "NUM_eles: " << NUM_eles << endl;
    //cout << "NUM_glob_interior_edges: " << NUM_glob_interior_edges << endl;
    pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> II_K;

    if (!if_CPU)
        II_K = this->AssembleOnGPU(mesh, Fracs, this->InletP, this->OutletP);
    else
        II_K = this->AssembleOnCPU(mesh, Fracs, this->InletP, this->OutletP,
                                   Nproc);

    int RowsC = NumEdges;

    if (!this->IfPeriodic)
        RowsC -= (NUM_INLET_EDGES + NUM_OUTLET_EDGES);

    Eigen::SparseMatrix<double> A(NUM_sep_edges, NUM_sep_edges);
    Eigen::SparseMatrix<double> B(NUM_eles, NUM_sep_edges);
    Eigen::SparseMatrix<double> C(RowsC, NUM_sep_edges);

    // cout << "K:\n"
    //      << MatrixXd(II_K.first.block(0, 0, NUM_sep_edges + NUM_eles + RowsC,
    //                                   NUM_sep_edges + NUM_eles + RowsC))
    //      << endl;
    // cout << "\n\n-----------\n\nb:\n"
    //      << MatrixXd(
    //             II_K.second.block(0, 0, NUM_sep_edges + NUM_eles + RowsC, 1))
    //      << endl;

    A = II_K.first.block(0, 0, NUM_sep_edges, NUM_sep_edges);
    B = II_K.first.block(NUM_sep_edges, 0, NUM_eles, NUM_sep_edges);
    C = II_K.first.block(NUM_sep_edges + NUM_eles, 0, RowsC, NUM_sep_edges);
    //cout << "A:\n" << MatrixXd(A) << endl;
    //cout << "B:\n" << MatrixXd(B) << endl;
    //cout << "C:\n" << MatrixXd(C) << endl;
    Eigen::SparseMatrix<double> g = II_K.second.block(0, 0, NUM_sep_edges, 1);

    Eigen::SparseMatrix<double> f =
        II_K.second.block(NUM_sep_edges, 0, NUM_eles, 1);

    Eigen::SparseMatrix<double> s =
        II_K.second.block(NUM_sep_edges + NUM_eles, 0, RowsC, 1);

    Eigen::SparseMatrix<double> A_inv(A.rows(), A.cols());
    A_inv.reserve(VectorXi::Constant(A.cols(), 3));

    double iStart_fem = cuDFNsys::CPUSecond();
    for (int i = 0; i < A.rows() / 3; ++i)
    {
        MatrixXd A_block = MatrixXd(A.block(i * 3, i * 3, 3, 3));
        A_block = A_block.inverse();
        /* 
        for (int j = i * 3, jq = 0; j < i * 3 + 3; ++j, ++jq)
        for (int k = i * 3, kq = 0; k < i * 3 + 3; ++k, ++kq)
        A_inv.insert(j, k) = A_block(jq, kq);*/

        for (int j = 0; j < 3; ++j)
        {
            SparseVector<double> SK(A.rows());
            SK.reserve(3);

            SK.insert(i * 3) = A_block.col(j)[0];
            SK.insert(i * 3 + 1) = A_block.col(j)[1];
            SK.insert(i * 3 + 2) = A_block.col(j)[2];

            A_inv.innerVector(i * 3 + j) = SK;
        }
    }
    cout << "\t\tRunning time of calculating A_inv: "
         << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";

    iStart_fem = cuDFNsys::CPUSecond();
    Eigen::SparseMatrix<double> C_tps = C.transpose();
    //cout << 4 << endl;
    Eigen::SparseMatrix<double> B_tps = B.transpose();

    Eigen::SparseMatrix<double> Wq = B * A_inv;

    Eigen::SparseMatrix<double> U = Wq * B_tps;
    cout << "\t\tRunning time of matrix multiplication 1: "
         << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";

    iStart_fem = cuDFNsys::CPUSecond();
    for (int i = 0; i < U.rows(); ++i)
        U.coeffRef(i, i) = 1.0 / U.coeffRef(i, i);
    cout << "\t\tRunning time of calculating U: "
         << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";

    iStart_fem = cuDFNsys::CPUSecond();
    Eigen::SparseMatrix<double> Re = C * A_inv;

    Eigen::SparseMatrix<double> Eq = Re * B_tps;

    Eigen::SparseMatrix<double> Sd = Wq * C_tps;
    //cout << 9 << endl;
    Eigen::SparseMatrix<double> D = Re * C_tps - Eq * U * Sd;

    Eigen::SparseMatrix<double> r = Re * g + Eq * U * (f - Wq * g) - s;
    D.makeCompressed();
    r.makeCompressed();
    cout << "\t\tRunning time of matrix multiplication 2: "
         << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";

    //cout << "\tsolving ...\n";

    iStart_fem = cuDFNsys::CPUSecond();
    UmfPackLU<Eigen::SparseMatrix<double>> solver;
    //SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(D);
    if (solver.info() != Success)
        throw cuDFNsys::ExceptionsPause("Compute matrix is wrong!\n");
    this->PressureInteriorEdge = solver.solve(r);

    // cout << "this->AllEdges:\n";
    // cout << this->PressureInteriorEdge << endl;
    cout << "\t\tRunning time of solving matrix: "
         << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";

    iStart_fem = cuDFNsys::CPUSecond();

    this->PressureEles = U * (Wq * g - Sd * this->PressureInteriorEdge - f);
    this->VelocityNormalScalarSepEdges =
        A_inv *
        (g - B_tps * this->PressureEles - C_tps * this->PressureInteriorEdge);

    // cout << this->PressureEles << endl;
    //cout << "\tcalculating flux ...\n";
    //cout << "\n\nin:\n";
    // cout << "VelocityNormalScalarSepEdges: \n";
    // cout << VelocityNormalScalarSepEdges << endl;
    for (size_t i = 0; i < mesh.InletEdgeNOLen.size(); ++i)
    {
        size_t sep_EDGE_no = mesh.InletEdgeNOLen[i].x;
        T len = mesh.InletEdgeNOLen[i].y;

        T veloc_length =
            abs((T)this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0) * len);
        this->QIn += veloc_length;
        this->InletLength += len;
        // cout << "inlet len: " << len << ";\tq: " << this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0) << "; sep_EDGE_no: " << sep_EDGE_no + 1 << "\n";
        if (this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0) > 0)
        {
            cout << "\t\t*** Note: One inlet velocity is > 0: "
                 << this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0);
            cout << ", sep_EDGE_no: " << sep_EDGE_no;
            cout << ", elementID (from 1): "
                 << ceil((1.0 * sep_EDGE_no + 1.0) / 3.0);
            cout << ", element Pressure: "
                 << this->PressureEles(
                        int(ceil((1.0 * sep_EDGE_no + 1.0) / 3.0)), 0)
                 << endl;
        }

        //opp += veloc_length;
    }
    //cout << "\n\nout:\n";
    for (size_t i = 0; i < mesh.OutletEdgeNOLen.size(); ++i)
    {
        size_t sep_EDGE_no = mesh.OutletEdgeNOLen[i].x;
        T len = mesh.OutletEdgeNOLen[i].y;

        T veloc_length =
            abs((T)this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0) * len);
        this->QOut += veloc_length;
        this->OutletLength += len;
        // cout << "outlet len: " << len << ";\tq: " << this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0) << "; sep_EDGE_no: " << sep_EDGE_no + 1 << "\n";
        if (this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0) < 0)
        {
            cout << "\t\t*** Note: One outlet velocity is < 0: "
                 << this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0);
            cout << ", sep_EDGE_no: " << sep_EDGE_no;
            cout << ", elementID (from 1): "
                 << ceil((1.0 * sep_EDGE_no + 1.0) / 3.0);
            cout << ", element Pressure: "
                 << this->PressureEles(
                        int(ceil((1.0 * sep_EDGE_no + 1.0) / 3.0)), 0)
                 << endl;
        }
    }
    cout << "\t\tRunning time of post treatments: "
         << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";
}; // Implementation
template void cuDFNsys::MHFEM<double>::Implementation(
    const cuDFNsys::Mesh<double> &mesh,
    const thrust::host_vector<cuDFNsys::Fracture<double>> &Fracs, bool if_CPU,
    int Nproc);
template void cuDFNsys::MHFEM<float>::Implementation(
    const cuDFNsys::Mesh<float> &mesh,
    const thrust::host_vector<cuDFNsys::Fracture<float>> &Fracs, bool if_CPU,
    int Nproc);

// ====================================================
// NAME:        AssembleOnGPU
// DESCRIPTION: Assemble matrix on GPU
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
template <typename T>
pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
cuDFNsys::MHFEM<T>::AssembleOnGPU(
    const cuDFNsys::Mesh<T> &mesh,
    const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs, T P_in, T P_out)
{

    int NUM_sep_edges = mesh.Element3D.size() * 3,
        NUM_eles = mesh.Element3D.size(),
        NUM_glob_interior_edges = mesh.NumInteriorEdges,
        NUM_INLET_EDGES = mesh.NumInletEdges,
        NUM_OUTLET_EDGES = mesh.NumOutletEdges,
        NUM_NEUMMANN_EDGES = mesh.NumNeumannEdges;

    int Dim = NUM_sep_edges + NUM_eles + NUM_glob_interior_edges +
              NUM_INLET_EDGES + NUM_OUTLET_EDGES + NUM_NEUMMANN_EDGES;
    // cout << "Dim: " << Dim << endl;
    //cout << "NUM_sep_edges: " << NUM_sep_edges << ", NUM_eles: " << NUM_eles << ", NUM_glob_interior_edges: " << NUM_glob_interior_edges << endl;

    int Frac_NUM = mesh.Element2D.size();

    size_t tmp_e = 0;
    int iyt = 0;
    thrust::host_vector<T> Conduc_Frac(NUM_eles);

    for (std::vector<size_t>::iterator it_fracID = mesh.FracID->begin();
         it_fracID != mesh.FracID->end(); it_fracID++)
    {
        T conduct_tmp = Fracs[*(it_fracID)].Conductivity;

        for (int j = 0; j < mesh.Element2D[iyt].size(); ++j)
        {
            Conduc_Frac[tmp_e] = conduct_tmp;
            //cout << "conduct_tmp " <<  conduct_tmp << endl;
            tmp_e++;
        };

        iyt++;
    }

    double iStart_fem = cuDFNsys::CPUSecond();

    thrust::device_vector<cuDFNsys::EleCoor<T>> Coodin_2D_dev;
    Coodin_2D_dev = mesh.Coordinate2D;
    cuDFNsys::EleCoor<T> *coord_2D_dev_ptr =
        thrust::raw_pointer_cast(Coodin_2D_dev.data());

    cuDFNsys::EleEdgeAttri *Edge_attri_dev_ptr;
    thrust::device_vector<cuDFNsys::EleEdgeAttri> Edge_attri_dev;
    Edge_attri_dev = mesh.EdgeAttri;
    Edge_attri_dev_ptr = thrust::raw_pointer_cast(Edge_attri_dev.data());

    thrust::device_vector<T> Conduc_Frac_dev = Conduc_Frac;
    T *Conduc_Frac_dev_ptr = thrust::raw_pointer_cast(Conduc_Frac_dev.data());

    thrust::host_vector<cuDFNsys::Triplet<T>> tri_h;
    thrust::device_vector<cuDFNsys::Triplet<T>> tri_dev(18 * NUM_eles);
    cuDFNsys::Triplet<T> *tri_dev_ptr =
        thrust::raw_pointer_cast(tri_dev.data());

    //cout << "this->MuOverRhoG: " << this->MuOverRhoG << endl;
    // cout << "this->ConsTq " << this->ConsTq << endl;
    cuDFNsys::AssembleOnGPUKernel<T><<<NUM_eles / 256 + 1, 256>>>(
        tri_dev_ptr, coord_2D_dev_ptr, Edge_attri_dev_ptr, Conduc_Frac_dev_ptr,
        NUM_sep_edges, NUM_eles, Dim, P_in, P_out, this->IfPeriodic,
        (T)this->MuOverRhoG, this->ConsTq);
    cudaDeviceSynchronize();
    this->TripletTime = cuDFNsys::CPUSecond() - iStart_fem;
    cout << "\t\tRunning time of GPU triplets: " << this->TripletTime
         << "sec\n";

    iStart_fem = cuDFNsys::CPUSecond();

    tri_h = tri_dev;
    Eigen::SparseMatrix<double> K(Dim, Dim);
    SparseMatrix<double> b(Dim, 1);
    K.reserve(VectorXi::Constant(Dim, 7));
    b.reserve(NUM_INLET_EDGES + NUM_OUTLET_EDGES + NUM_NEUMMANN_EDGES);

    for (size_t i = 0; i < NUM_eles * 18; ++i)
    {
        if (tri_h[i].row != -1)
        {
            //cout << "ele: " << i / 18 + 1 << endl;
            //cout << "DIM: " << Dim << "; ";
            //cout << "NUM_sep_edges: " << NUM_sep_edges << "; ";
            //cout << "NUM_eles: " << NUM_eles << "; ";
            //cout << "NUM_glob_interior_edges: " << NUM_glob_interior_edges << "; ";
            //cout << tri_h[i].row << ", " << tri_h[i].col << ", " << tri_h[i].val << endl;
            if (tri_h[i].col != Dim + 2)
                K.insert(tri_h[i].row, tri_h[i].col) = tri_h[i].val;
            else
            {
                b.insert(tri_h[i].row, 0) = tri_h[i].val;
                //cout << "b: " << tri_h[i].row << ", " << tri_h[i].val << endl;
            }
        }
    }
    cout << "\t\tRunning time of Matrix assembly and generation: "
         << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";
    return std::make_pair(K, b);
}; // AssembleOnGPU
template pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
cuDFNsys::MHFEM<double>::AssembleOnGPU(
    const cuDFNsys::Mesh<double> &mesh,
    const thrust::host_vector<cuDFNsys::Fracture<double>> &Fracs, double P_in,
    double P_out);
template pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
cuDFNsys::MHFEM<float>::AssembleOnGPU(
    const cuDFNsys::Mesh<float> &mesh,
    const thrust::host_vector<cuDFNsys::Fracture<float>> &Fracs, float P_in,
    float P_out);

// ====================================================
// NAME:        AssembleOnCPU
// DESCRIPTION: Assemble matrix on GPU
// AUTHOR:      Tingchang YIN
// DATE:        05/11/2022
// ====================================================
template <typename T>
pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
cuDFNsys::MHFEM<T>::AssembleOnCPU(
    const cuDFNsys::Mesh<T> &mesh,
    const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs, T P_in, T P_out,
    int Nproc)
{
    throw cuDFNsys::ExceptionsPause(
        "Use of `cuDFNsys::MHFEM<T>::AssembleOnCPU` is not recommended");

    int NUM_sep_edges = mesh.Element3D.size() * 3,
        NUM_eles = mesh.Element3D.size(),
        NUM_glob_interior_edges = mesh.NumInteriorEdges,
        NUM_INLET_EDGES = mesh.NumInletEdges,
        NUM_OUTLET_EDGES = mesh.NumOutletEdges;

    int Dim = NUM_sep_edges + NUM_eles + NUM_glob_interior_edges;
    // cout << "NUM_sep_edges: " << NUM_sep_edges << ", NUM_eles: " << NUM_eles << ", NUM_glob_interior_edges: " << NUM_glob_interior_edges << endl;

    int Frac_NUM = mesh.Element2D.size();

    size_t tmp_e = 0;
    int iyt = 0;
    thrust::host_vector<T> Conduc_Frac(NUM_eles);

    for (std::vector<size_t>::iterator it_fracID = mesh.FracID->begin();
         it_fracID != mesh.FracID->end(); it_fracID++)
    {
        T conduct_tmp = Fracs[*(it_fracID)].Conductivity;

        for (int j = 0; j < mesh.Element2D[iyt].size(); ++j)
        {
            Conduc_Frac[tmp_e] = conduct_tmp;
            //cout << "conduct_tmp " <<  conduct_tmp << endl;
            tmp_e++;
        };

        iyt++;
    }

    //---
    double iStart_fem = cuDFNsys::CPUSecond();

    thrust::host_vector<cuDFNsys::Triplet<T>> tri_h(21 * NUM_eles);

#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (int i = 0; i < NUM_eles; ++i)
    {
        int I[3] = {i * 3 + 1, // 2
                    i * 3 + 2, // 3
                    i * 3};    // 1

        cuDFNsys::EleCoor<T> coord = mesh.Coordinate2D[i];

        T A[3][3];
        cuDFNsys::StimaA<T>(coord, A);

        size_t j = i * 21;
        size_t j_tt = j;

        for (size_t ik = 0; ik < 3; ++ik)
            for (size_t jk = 0; jk < 3; ++jk)
            {
                tri_h[j].row = I[ik];
                tri_h[j].col = I[jk];
                tri_h[j].val = -1.0f / Conduc_Frac[i] * A[ik][jk];

                int edge_1 = tri_h[j].row % 3;
                int edge_2 = tri_h[j].col % 3;
                if (mesh.EdgeAttri[i].e[edge_1] == 2 &&
                    mesh.EdgeAttri[i].e[edge_2] == 2 && edge_1 == edge_2)
                    tri_h[j].val = 1.0;
                else if (edge_1 != edge_2 &&
                         (mesh.EdgeAttri[i].e[edge_1] == 2 ||
                          mesh.EdgeAttri[i].e[edge_2] == 2))
                    tri_h[j].row = -1;

                j++;
            }

        T B[3] = {0};
        B[0] = pow(pow(coord.x[2] - coord.x[1], 2) +
                       pow(coord.y[2] - coord.y[1], 2),
                   0.5);
        B[1] = pow(pow(coord.x[0] - coord.x[2], 2) +
                       pow(coord.y[0] - coord.y[2], 2),
                   0.5);
        B[2] = pow(pow(coord.x[1] - coord.x[0], 2) +
                       pow(coord.y[1] - coord.y[0], 2),
                   0.5);

        for (size_t ik = 0; ik < 3; ++ik)
        {
            tri_h[j].row = I[ik];
            tri_h[j].col = i + NUM_sep_edges;
            tri_h[j].val = B[ik];
            j++;
            tri_h[j].row = i + NUM_sep_edges;
            tri_h[j].col = I[ik];
            tri_h[j].val = B[ik];
            j++;

            int edge_1 = tri_h[j - 2].row % 3;
            if (mesh.EdgeAttri[i].e[edge_1] == 2)
            {
                tri_h[j - 2].row = -1;
                tri_h[j - 1].row = -1;
            }
        }

        T P_in_out[2] = {P_in, P_out};

        for (size_t ik = 0; ik < 3; ++ik)
        {
            int ek = mesh.EdgeAttri[i].e[ik];
            int NO_ = mesh.EdgeAttri[i].no[ik] - 1;

            if (ek == 3)
            {
                tri_h[j].row = NO_ + NUM_sep_edges + NUM_eles;
                tri_h[j].col = I[(ik + 2) % 3];
                tri_h[j].val = -B[(ik + 2) % 3];
                j++;
                tri_h[j].col = NO_ + NUM_sep_edges + NUM_eles;
                tri_h[j].row = I[(ik + 2) % 3];
                tri_h[j].val = -B[(ik + 2) % 3];
                j++;
            }
            else if (ek == 0 || ek == 1) // in or out
            {
                tri_h[j].row = NO_;
                tri_h[j].col =
                    NUM_sep_edges + NUM_eles + NUM_glob_interior_edges + 2;
                tri_h[j].val = P_in_out[ek] * B[(ik + 2) % 3];
                j++;
            }
            //else if (ek == 2) // neumann
            //neuman_sep_dev[NO_] = 1 + NO_;
        };

        for (size_t k = j; k < j_tt + 21; ++k)
            tri_h[k].row = -1;
    }

    //------------------------------
    //------------------------------
    //------------------------------
    //------------------------------
    this->TripletTime = cuDFNsys::CPUSecond() - iStart_fem;
    cout << "\t\tRunning time of CPU triplets with Nproc = " << Nproc << ": "
         << this->TripletTime << "sec\n";

    iStart_fem = cuDFNsys::CPUSecond();

    Eigen::SparseMatrix<double> K(Dim, Dim);
    SparseMatrix<double> b(Dim, 1);
    K.reserve(VectorXi::Constant(Dim, 5));
    b.reserve(NUM_INLET_EDGES + NUM_OUTLET_EDGES);

    for (size_t i = 0; i < NUM_eles * 21; ++i)
    {
        if (tri_h[i].row != -1)
        {
            //cout << "ele: " << i / 21 + 1 << endl;
            //cout << "DIM: " << Dim << "; ";
            //cout << "NUM_sep_edges: " << NUM_sep_edges << "; ";
            //cout << "NUM_eles: " << NUM_eles << "; ";
            //cout << "NUM_glob_interior_edges: " << NUM_glob_interior_edges << "; ";
            //cout << tri_h[i].row << ", " << tri_h[i].col << ", " << tri_h[i].val << endl;
            if (tri_h[i].col != Dim + 2)
                K.insert(tri_h[i].row, tri_h[i].col) = tri_h[i].val;
            else
            {
                b.insert(tri_h[i].row, 0) = tri_h[i].val;
                //cout << "b: " << tri_h[i].row << ", " << tri_h[i].val << endl;
            }
        }
    }
    cout << "\t\tRunning time of Matrix assembly and generation: "
         << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";
    return std::make_pair(K, b);
}; // AssembleOnCPU
template pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
cuDFNsys::MHFEM<double>::AssembleOnCPU(
    const cuDFNsys::Mesh<double> &mesh,
    const thrust::host_vector<cuDFNsys::Fracture<double>> &Fracs, double P_in,
    double P_out, int Nproc);
template pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
cuDFNsys::MHFEM<float>::AssembleOnCPU(
    const cuDFNsys::Mesh<float> &mesh,
    const thrust::host_vector<cuDFNsys::Fracture<float>> &Fracs, float P_in,
    float P_out, int Nproc);
