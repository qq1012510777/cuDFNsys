#include "MHFEM/MHFEM.cuh"

// ====================================================
// NAME:        MHFEM
// DESCRIPTION: MHFEM constructor
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
cuDFNsys::MHFEM::MHFEM(const cuDFNsys::Mesh &mesh,
                       const thrust::host_vector<cuDFNsys::Fracture> &Fracs,
                       const float &inlet_p_,
                       const float &outlet_p_,
                       const int &dir_,
                       const float &L)
{
    Eigen::setNbThreads(1);
    this->Dir = dir_;
    this->InletP = inlet_p_;
    this->OutletP = outlet_p_;
    this->Implementation(mesh, Fracs);
    this->QError = (abs(QIn - QOut) / (QIn > QOut ? QIn : QOut)) * 100.0f;
    this->Permeability = 0.5f * (QOut + QIn) / L / (inlet_p_ - outlet_p_);
}; // MHFEM

// ====================================================
// NAME:        MatlabPlot
// DESCRIPTION: plot mhfem result
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
void cuDFNsys::MHFEM::MatlabPlot(const string &mat_key,
                                 const string &command_key,
                                 thrust::host_vector<cuDFNsys::Fracture> Fracs,
                                 const cuDFNsys::Mesh &mesh,
                                 const float &L)
{
    cuDFNsys::MatlabAPI M1;

    size_t node_num = mesh.Coordinate3D.size();

    float *ptr_coordinates_3D;
    ptr_coordinates_3D = new float[node_num * 3];

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
    M1.WriteMat(mat_key, "w", node_num * 3,
                node_num, 3, ptr_coordinates_3D, "coordinate_3D");

    delete[] ptr_coordinates_3D;
    ptr_coordinates_3D = NULL;

    //---------------------
    size_t ele_num = mesh.Element3D.size();
    float *ptr_element_3D = new float[ele_num * 3];
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
    M1.WriteMat(mat_key, "u", ele_num * 3,
                ele_num, 3, ptr_element_3D, "element_3D");

    delete[] ptr_element_3D;
    ptr_element_3D = NULL;
    //--------------------------

    float *pressure_ELEs = new float[ele_num];
    if (pressure_ELEs == NULL)
    {
        string AS = "Alloc error in MHFEM::MatlabPlot\n";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    std::copy(this->PressureEles.data(),
              this->PressureEles.data() + ele_num,
              pressure_ELEs);

    M1.WriteMat(mat_key, "u", ele_num,
                ele_num, 1, pressure_ELEs,
                "pressure_eles");

    delete[] pressure_ELEs;
    pressure_ELEs = NULL;

    //------------------------------
    //float *normal_veloc = new float[ele_num * 3];
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
    float *velocity_center_grid = new float[mesh.Element3D.size() * 3];
    if (velocity_center_grid == NULL)
    {
        string AS = "Alloc error in MHFEM::MatlabPlot\n";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    for (uint i = 0; i < mesh.Element3D.size(); ++i)
    {
        uint3 EdgeNO = make_uint3((i + 1) * 3 - 3, (i + 1) * 3 - 2, (i + 1) * 3 - 1); // from 0
        float3 Velocity_ = make_float3((float)this->VelocityNormalScalarSepEdges(EdgeNO.x, 0),
                                       (float)this->VelocityNormalScalarSepEdges(EdgeNO.y, 0),
                                       (float)this->VelocityNormalScalarSepEdges(EdgeNO.z, 0));
        float2 Vertexes[3];
        Vertexes[0] = make_float2(mesh.Coordinate2D[i].x[0], mesh.Coordinate2D[i].y[0]);
        Vertexes[1] = make_float2(mesh.Coordinate2D[i].x[1], mesh.Coordinate2D[i].y[1]);
        Vertexes[2] = make_float2(mesh.Coordinate2D[i].x[2], mesh.Coordinate2D[i].y[2]);

        float2 Center_p = make_float2(1.0f / 3.0f * (Vertexes[0].x + Vertexes[1].x + Vertexes[2].x), 1.0f / 3.0f * (Vertexes[0].y + Vertexes[1].y + Vertexes[2].y));

        float2 velocity_p = cuDFNsys::ReconstructVelocityGrid(Center_p, Vertexes, Velocity_);

        float3 velocity_p_3D = make_float3(velocity_p.x, velocity_p.y, 0);

        float R_mat[3][3];
        Fracs[mesh.ElementFracTag[i]].RoationMatrix(R_mat, 23);
        velocity_p_3D = cuDFNsys::ProductSquare3Float3(R_mat, velocity_p_3D);

        velocity_center_grid[i] = velocity_p_3D.x;
        velocity_center_grid[i + mesh.Element3D.size()] = velocity_p_3D.y;
        velocity_center_grid[i + 2 * mesh.Element3D.size()] = velocity_p_3D.z;
    };

    M1.WriteMat(mat_key, "u", mesh.Element3D.size() * 3,
                mesh.Element3D.size(), 3, velocity_center_grid,
                "velocity_center_grid");
    delete[] velocity_center_grid;
    velocity_center_grid = NULL;

    //-----------------
    std::ofstream oss(command_key, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    oss << "load('" << mat_key << "');\n";
    oss << "L = 0.5 * " << L << ";\n";
    oss << "cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];\n";
    oss << "figure(1); view(3); title('DFN flow (mhfem)'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on\n";
    oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on\n";
    oss << endl;
    oss << "xlim([-1.1*L, 1.1*L]);\n";
    oss << "ylim([-1.1*L, 1.1*L]);\n";
    oss << "zlim([-1.1*L, 1.1*L]);\nhold on\n";
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
    oss << "patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1); colorbar; view(3); hold on\n";
    oss << "caxis([" << this->OutletP << ", " << this->InletP << "]);\n\n";

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
    oss << "CenterELE(:, 1) = 1/3 * (coordinate_3D(element_3D(:, 1), 1) + coordinate_3D(element_3D(:, 2), 1) + coordinate_3D(element_3D(:, 3), 1));\n";
    oss << "CenterELE(:, 2) = 1/3 * (coordinate_3D(element_3D(:, 1), 2) + coordinate_3D(element_3D(:, 2), 2) + coordinate_3D(element_3D(:, 3), 2));\n";
    oss << "CenterELE(:, 3) = 1/3 * (coordinate_3D(element_3D(:, 1), 3) + coordinate_3D(element_3D(:, 2), 3) + coordinate_3D(element_3D(:, 3), 3));\n";
    oss << "quiver3(CenterELE(:, 1), CenterELE(:, 2), CenterELE(:, 3), velocity_center_grid(:, 1),velocity_center_grid(:, 2),velocity_center_grid(:, 3), 4, 'LineWidth', 1.5, 'color', 'r');\n";
    oss.close();
}; // MHFEM

// ====================================================
// NAME:        Implementation
// DESCRIPTION: Implementation of mhfem
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
void cuDFNsys::MHFEM::Implementation(const cuDFNsys::Mesh &mesh,
                                     const thrust::host_vector<cuDFNsys::Fracture> &Fracs)
{
    size_t NUM_sep_edges = mesh.Element3D.size() * 3,
           NUM_eles = mesh.Element3D.size(),
           NUM_glob_interior_edges = mesh.NumInteriorEdges;

    pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> II_K;
    II_K = this->AssembleOnGPU(mesh, Fracs, this->InletP, this->OutletP);

    Eigen::SparseMatrix<double> A(NUM_sep_edges, NUM_sep_edges);
    Eigen::SparseMatrix<double> B(NUM_eles, NUM_sep_edges);
    Eigen::SparseMatrix<double> C(NUM_glob_interior_edges, NUM_sep_edges);

    A = II_K.first.block(0, 0, NUM_sep_edges, NUM_sep_edges);
    B = II_K.first.block(NUM_sep_edges, 0, NUM_eles, NUM_sep_edges);
    C = II_K.first.block(NUM_sep_edges + NUM_eles, 0, NUM_glob_interior_edges, NUM_sep_edges);

    Eigen::SparseMatrix<double> g = II_K.second.block(0, 0, NUM_sep_edges, 1);

    Eigen::SparseMatrix<double> f = II_K.second.block(NUM_sep_edges, 0, NUM_eles, 1);

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
    cout << "\t\tRunning time of calculating A_inv: " << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";

    iStart_fem = cuDFNsys::CPUSecond();
    Eigen::SparseMatrix<double> C_tps = C.transpose();
    //cout << 4 << endl;
    Eigen::SparseMatrix<double> B_tps = B.transpose();

    Eigen::SparseMatrix<double> Wq = B * A_inv;

    Eigen::SparseMatrix<double> U = Wq * B_tps;
    cout << "\t\tRunning time of matrix multiplication 1: " << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";

    iStart_fem = cuDFNsys::CPUSecond();
    for (int i = 0; i < U.rows(); ++i)
        U.coeffRef(i, i) = 1.0 / U.coeffRef(i, i);
    cout << "\t\tRunning time of calculating U: " << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";

    iStart_fem = cuDFNsys::CPUSecond();
    Eigen::SparseMatrix<double> Re = C * A_inv;

    Eigen::SparseMatrix<double> Eq = Re * B_tps;

    Eigen::SparseMatrix<double> Sd = Wq * C_tps;
    //cout << 9 << endl;
    Eigen::SparseMatrix<double> D = Re * C_tps - Eq * U * Sd;

    Eigen::SparseMatrix<double> r = Re * g + Eq * U * (f - Wq * g);
    D.makeCompressed();
    r.makeCompressed();
    cout << "\t\tRunning time of matrix multiplication 2: " << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";

    //cout << "\tsolving ...\n";
    iStart_fem = cuDFNsys::CPUSecond();
    UmfPackLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(D);
    this->PressureInteriorEdge = solver.solve(r);
    cout << "\t\tRunning time of solving matrix: " << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";

    iStart_fem = cuDFNsys::CPUSecond();

    this->PressureEles = U * (Wq * g - Sd * this->PressureInteriorEdge - f);
    this->VelocityNormalScalarSepEdges = A_inv * (g - B_tps * this->PressureEles - C_tps * this->PressureInteriorEdge);

    //cout << this->PressureEles << endl;
    //cout << "\tcalculating flux ...\n";
    //cout << "\n\nin:\n";
    for (size_t i = 0; i < mesh.InletEdgeNOLen.size(); ++i)
    {
        size_t sep_EDGE_no = mesh.InletEdgeNOLen[i].x - 1;
        float len = mesh.InletEdgeNOLen[i].y;

        float veloc_length = abs((float)this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0) * len);
        this->QIn += veloc_length;
        this->InletLength += len;
        //cout << "len: " << len << ";\tq: " << this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0) << "; sep_EDGE_no: " << sep_EDGE_no + 1 << "\n";
        //if (this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0) > 0)
        //opp += veloc_length;
    }
    //cout << "\n\nout:\n";
    for (size_t i = 0; i < mesh.OutletEdgeNOLen.size(); ++i)
    {
        size_t sep_EDGE_no = mesh.OutletEdgeNOLen[i].x - 1;
        float len = mesh.OutletEdgeNOLen[i].y;

        float veloc_length = abs((float)this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0) * len);
        this->QOut += veloc_length;
        this->OutletLength += len;
        //cout << "len: " << len << ";\tq: " << this->VelocityNormalScalarSepEdges(sep_EDGE_no, 0) << "; sep_EDGE_no: " << sep_EDGE_no + 1 << "\n";
    }
    cout << "\t\tRunning time of post treatments: " << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";
}; // Implementation

// ====================================================
// NAME:        AssembleOnGPU
// DESCRIPTION: Assemble matrix on GPU
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> cuDFNsys::MHFEM::AssembleOnGPU(const cuDFNsys::Mesh &mesh,
                                                                                              const thrust::host_vector<cuDFNsys::Fracture> &Fracs,
                                                                                              float P_in,
                                                                                              float P_out)
{
    double iStart_fem = cuDFNsys::CPUSecond();

    int NUM_sep_edges = mesh.Element3D.size() * 3,
        NUM_eles = mesh.Element3D.size(),
        NUM_glob_interior_edges = mesh.NumInteriorEdges,
        NUM_INLET_EDGES = mesh.NumInletEdges,
        NUM_OUTLET_EDGES = mesh.NumOutletEdges;

    int Dim = NUM_sep_edges + NUM_eles + NUM_glob_interior_edges;

    int Frac_NUM = mesh.Element2D.size();

    size_t tmp_e = 0;
    int iyt = 0;
    thrust::host_vector<float> Conduc_Frac(NUM_eles);

    for (std::vector<size_t>::iterator it_fracID = mesh.FracID->begin();
         it_fracID != mesh.FracID->end();
         it_fracID++)
    {
        float conduct_tmp = Fracs[*(it_fracID)].Conductivity;

        for (int j = 0; j < mesh.Element2D[iyt].size(); ++j)
        {
            Conduc_Frac[tmp_e] = conduct_tmp;
            //cout << "conduct_tmp " <<  conduct_tmp << endl;
            tmp_e++;
        };

        iyt++;
    }

    thrust::device_vector<cuDFNsys::EleCoor> Coodin_2D_dev;
    Coodin_2D_dev = mesh.Coordinate2D;
    cuDFNsys::EleCoor *coord_2D_dev_ptr = thrust::raw_pointer_cast(Coodin_2D_dev.data());

    cuDFNsys::EleEdgeAttri *Edge_attri_dev_ptr;
    thrust::device_vector<cuDFNsys::EleEdgeAttri> Edge_attri_dev;
    Edge_attri_dev = mesh.EdgeAttri;
    Edge_attri_dev_ptr = thrust::raw_pointer_cast(Edge_attri_dev.data());

    thrust::device_vector<float> Conduc_Frac_dev = Conduc_Frac;
    float *Conduc_Frac_dev_ptr = thrust::raw_pointer_cast(Conduc_Frac_dev.data());

    thrust::host_vector<cuDFNsys::Triplet> tri_h;
    thrust::device_vector<cuDFNsys::Triplet> tri_dev(21 * NUM_eles);
    cuDFNsys::Triplet *tri_dev_ptr = thrust::raw_pointer_cast(tri_dev.data());

    cuDFNsys::AssembleOnGPUKernel<<<NUM_eles / 256 + 1, 256>>>(tri_dev_ptr,
                                                               coord_2D_dev_ptr,
                                                               Edge_attri_dev_ptr,
                                                               Conduc_Frac_dev_ptr,
                                                               NUM_sep_edges,
                                                               NUM_eles,
                                                               NUM_glob_interior_edges,
                                                               P_in,
                                                               P_out);
    cudaDeviceSynchronize();
    cout << "\t\tRunning time of GPU assemble: " << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";

    iStart_fem = cuDFNsys::CPUSecond();

    tri_h = tri_dev;
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
    cout << "\t\tRunning time of Matrix generation: " << cuDFNsys::CPUSecond() - iStart_fem << "sec\n";
    return std::make_pair(K, b);
}; // AssembleOnGPU
