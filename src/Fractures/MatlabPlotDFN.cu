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

#include "Fractures/MatlabPlotDFN.cuh"

// ====================================================
// NAME:        cuDFNsys::MatlabPlotDFN::MatlabPlotDFN
// DESCRIPTION: Plot DFN with m and mat file
//              m is command file; mat is data file
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
template <typename T>
cuDFNsys::MatlabPlotDFN<T>::MatlabPlotDFN(
    string mat_key,     // mat file name
    string command_key, // m file name
    thrust::host_vector<cuDFNsys::Fracture<T>>
        Frac_verts_host, // Vector of Fracture
    std::map<pair<size_t, size_t>,
             pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>>
        Intersection_host, // vector of Intersection
    std::vector<std::vector<size_t>>
        ListClusters, // List of clusters: a cluster is a vector
    std::vector<size_t>
        Percolation_cluster, // Percolation cluster contains the cluster vector ID
    bool If_show_truncated_frac, // if show truncated fractures
    bool If_show_intersection, bool If_show_cluster, bool If_show_orientation,
    T L, int dir, bool if_python_visualization,
    string PythonName_Without_suffix, double3 DomainDimensionRatio_d)
{
    //cuDFNsys::MatlabAPI M1;
    cuDFNsys::HDF5API h5gg;

    h5gg.NewFile(mat_key);

    int NUM_Frac = Frac_verts_host.size();
    int sum_NUM_verts = 4 * NUM_Frac;
    T *Frac_NUM_verts = new T[NUM_Frac];
    if (Frac_NUM_verts == NULL)
    {
        string AS = "Alloc error in MatlabPlotDFN!";
        throw cuDFNsys::ExceptionsPause(AS);
    }
    std::fill_n(Frac_NUM_verts, NUM_Frac, 4);

    if (If_show_truncated_frac == true)
    {
        sum_NUM_verts = 0;
        for (int i = 0; i < NUM_Frac; ++i)
        {
            Frac_NUM_verts[i] = Frac_verts_host[i].NumVertsTruncated;
            sum_NUM_verts += Frac_NUM_verts[i];
        }
    }

    uint2 dim_f = make_uint2(1, NUM_Frac);
    h5gg.AddDataset(mat_key, "N", "Frac_NUM_verts", Frac_NUM_verts, dim_f);
    //M1.WriteMat(mat_key, "w", NUM_Frac, NUM_Frac, 1, Frac_NUM_verts, "Frac_NUM_verts");
    delete[] Frac_NUM_verts;
    Frac_NUM_verts = NULL;

    // ----------
    dim_f = make_uint2(1, 1);
    T modelsize[1] = {L};
    h5gg.AddDataset(mat_key, "N", "L_m", modelsize, dim_f);

    // int NumCLusters = ListClusters.size();
    // h5gg.AddDataset(mat_key, "N", "NumClusters", &NumCLusters, dim_f);

    //------------
    uint2 dim_ds = make_uint2(3, 1);
    double DomainDimensionRatio[3] = {DomainDimensionRatio_d.x,
                                      DomainDimensionRatio_d.y,
                                      DomainDimensionRatio_d.z};
    h5gg.AddDataset(mat_key, "N", "DomainDimensionRatio", DomainDimensionRatio,
                    dim_ds);

    //-----------------
    T *R_ = new T[NUM_Frac];
    if (R_ == NULL)
    {
        string AS = "Alloc error in MatlabPlotDFN!";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    for (int i = 0; i < NUM_Frac; ++i)
        R_[i] = Frac_verts_host[i].Radius;

    //dim_f = make_uint2(1, NUM_Frac);
    h5gg.AddDataset(mat_key, "N", "R", R_, dim_f);
    // M1.WriteMat(mat_key, "u", NUM_Frac, NUM_Frac, 1,
    //             R_, "R");
    delete[] R_;
    R_ = NULL;
    //----------------------

    T *verts = new T[sum_NUM_verts * 3];
    if (verts == NULL)
    {
        string AS = "Alloc error in MatlabPlotDFN!";
        throw cuDFNsys::ExceptionsPause(AS);
    }
    int tmp_cc = 0;

    for (int i = 0; i < NUM_Frac; ++i)
    {
        int hj = 4;

        if (If_show_truncated_frac == true)
            hj = Frac_verts_host[i].NumVertsTruncated;

        for (int j = 0; j < hj; ++j)
        {
            if (If_show_truncated_frac == false)
            {
                verts[tmp_cc] = Frac_verts_host[i].Verts3D[j].x;
                verts[tmp_cc + sum_NUM_verts] = Frac_verts_host[i].Verts3D[j].y;
                verts[tmp_cc + sum_NUM_verts * 2] =
                    Frac_verts_host[i].Verts3D[j].z;
                tmp_cc++;
            }
            else
            {
                verts[tmp_cc] = Frac_verts_host[i].Verts3DTruncated[j].x;
                verts[tmp_cc + sum_NUM_verts] =
                    Frac_verts_host[i].Verts3DTruncated[j].y;
                verts[tmp_cc + sum_NUM_verts * 2] =
                    Frac_verts_host[i].Verts3DTruncated[j].z;
                tmp_cc++;
            }
        }
    }

    dim_f = make_uint2(3, sum_NUM_verts);
    h5gg.AddDataset(mat_key, "N", "verts", verts, dim_f);
    //M1.WriteMat(mat_key, "u", sum_NUM_verts * 3, sum_NUM_verts, 3,
    //            verts, "verts");
    delete[] verts;
    verts = NULL;

    //--------------------------------------
    if (If_show_intersection == true && Intersection_host.size() >= 1)
    {
        int NUM_intersections = Intersection_host.size();

        T *Intersections_ = new T[NUM_intersections * 8];
        if (Intersections_ == NULL)
        {
            string AS = "Alloc error in MatlabPlotDFN!";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        int tmp_jj = 0;
        for (auto i = Intersection_host.begin(); i != Intersection_host.end();
             ++i)
        {
            Intersections_[tmp_jj] = i->second.first.x;
            Intersections_[tmp_jj + NUM_intersections] = i->second.first.y;
            Intersections_[tmp_jj + NUM_intersections * 2] = i->second.first.z;
            Intersections_[tmp_jj + NUM_intersections * 3] = i->second.second.x;
            Intersections_[tmp_jj + NUM_intersections * 4] = i->second.second.y;
            Intersections_[tmp_jj + NUM_intersections * 5] = i->second.second.z;
            Intersections_[tmp_jj + NUM_intersections * 6] = i->first.first;
            Intersections_[tmp_jj + NUM_intersections * 7] = i->first.second;
            tmp_jj++;
        }

        dim_f = make_uint2(8, NUM_intersections);
        h5gg.AddDataset(mat_key, "N", "intersections", Intersections_, dim_f);
        //M1.WriteMat(mat_key, "u", NUM_intersections * 6, NUM_intersections, 6,
        // Intersections_, "intersections");

        delete[] Intersections_;
        Intersections_ = NULL;
    }

    if (If_show_cluster == true)
    {
        std::vector<uint> TMP_S;
        TMP_S.reserve(Percolation_cluster.size());
        std::transform(Percolation_cluster.begin(), Percolation_cluster.end(),
                       std::back_inserter(TMP_S),
                       [](size_t value)
                       { return static_cast<uint>(value + 1); });
        //cout << "Percolation_cluster" << endl;
        h5gg.AddDataset(mat_key, "N", "PercolationClusters", TMP_S.data(),
                        make_uint2(TMP_S.size(), 0));
        //cout << "Percolation_cluster222" << endl;
        int NumberClustersyy = ListClusters.size();
        h5gg.AddDataset(mat_key, "N", "NumClusters", &NumberClustersyy,
                        make_uint2(1, 0));
        //cout << "222\n";
        for (uint i = 0; i < ListClusters.size(); ++i)
        {
            std::vector<uint> TMP_O; //(ListClusters[i].size());
            TMP_O.reserve(ListClusters[i].size());
            //cout << "i: " << i << endl;
            std::transform(ListClusters[i].begin(), ListClusters[i].end(),
                           std::back_inserter(TMP_O),
                           [](size_t value)
                           { return static_cast<uint>(value + 1); });
            //cout << "ListClusters[i]" << endl;
            h5gg.AddDataset(mat_key, "N", "Cluster_" + std::to_string(i + 1),
                            TMP_O.data(), make_uint2(TMP_O.size(), 1));
        }
        //size_t Max_cluster_size = 0;
        //for (size_t i = 0; i < ListClusters.size(); ++i)
        //    Max_cluster_size = ListClusters[i].size() > Max_cluster_size ? ListClusters[i].size() : Max_cluster_size;
        //T *clustersList = new T[ListClusters.size() * Max_cluster_size];
        //if (clustersList == NULL)
        //{
        //    string AS = "Alloc error in MatlabPlotDFN!";
        //    throw cuDFNsys::ExceptionsPause(AS);
        //}
        ////cout << "Max_cluster_size: " << Max_cluster_size << endl;
        //int tmpff = 0;
        //for (size_t j = 0; j < Max_cluster_size; ++j)
        //    for (size_t i = 0; i < ListClusters.size(); ++i)
        //    {
        //        if (j < ListClusters[i].size())
        //            clustersList[tmpff] = ListClusters[i][j] + 1;
        //        else
        //            clustersList[tmpff] = -1.0;
        //        tmpff++;
        //    }
        ////M1.WriteMat(mat_key, "u", ListClusters.size() * Max_cluster_size, ListClusters.size(), Max_cluster_size, clustersList, "ListClusters");
        //dim_f = make_uint2(Max_cluster_size, ListClusters.size());
        //h5gg.AddDataset(mat_key, "N", "ListClusters", clustersList, dim_f);
        //delete[] clustersList;
        //clustersList = NULL;
        ////-----------
        //T *percolationcluster = new T[Percolation_cluster.size()];
        //if (percolationcluster == NULL)
        //{
        //    string AS = "Alloc error in MatlabPlotDFN!";
        //    throw cuDFNsys::ExceptionsPause(AS);
        //}
        //for (size_t i = 0; i < Percolation_cluster.size(); ++i)
        //    percolationcluster[i] = Percolation_cluster[i] + 1;
        ////M1.WriteMat(mat_key, "u", Percolation_cluster.size(),
        //// Percolation_cluster.size(), 1, percolationcluster, "PercolationClusters");
        //dim_f = make_uint2(1, Percolation_cluster.size());
        //h5gg.AddDataset(mat_key, "N", "PercolationClusters", percolationcluster, dim_f);
        //delete[] percolationcluster;
        //percolationcluster = NULL;
    };

    if (If_show_orientation == true)
    {
        T *orientation = new T[NUM_Frac * 2];
        if (orientation == NULL)
        {
            string AS = "Alloc error in MatlabPlotDFN!";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        for (size_t i = 0; i < NUM_Frac; ++i)
        {
            orientation[i] = Frac_verts_host[i].Phi();
            orientation[i + NUM_Frac] = Frac_verts_host[i].Theta();
            //cout << orientation[i] << ", " << orientation[i + NUM_Frac] << endl;
        }
        // M1.WriteMat(mat_key, "u", NUM_Frac * 2, NUM_Frac, 2, orientation, "Polar_Orientation");

        dim_f = make_uint2(2, NUM_Frac);
        h5gg.AddDataset(mat_key, "N", "Polar_Orientation", orientation, dim_f);

        delete[] orientation;
        orientation = NULL;
    };
    //------------------------------------------------
    if (if_python_visualization)
    {
        std::ofstream oss(PythonName_Without_suffix + ".py", ios::out);
        oss << "import h5py\n";
        oss << "import numpy as np\n";
        oss << "from mayavi import mlab as ML\n";
        oss << "\n";
        oss << "f = h5py.File('" << mat_key << "')\n";

        oss << "\n";
        oss << "verts = np.array(f['verts'][:])\n";
        oss << "Frac_NUM_verts = np.array(f['Frac_NUM_verts'][:])\n";
        oss << "verts = np.array(f['verts'][:])\n";
        oss << "L_m = f['L_m'][0]\n";
        oss << "DomainDimensionRatio = f['DomainDimensionRatio'][:]\n";

        if (If_show_cluster)
        {
            // oss << "ListClusters = "
            //        "np.array(f['ListClusters'][:])\nnp.delete(ListClusters, 0, "
            //        "axis=0)\n";
            oss << "NumClusters=np.array(f['NumClusters'][:])\n";
            oss << "ListClusters=[]\n";
            oss << "for i in range(NumClusters[0]):\n";
            oss << "\tListClusters.append(np.array(f[\"Cluster_\" + str(i + 1)][:]))\n";

            oss << "PercolationClusters = "
                   "np.array(f['PercolationClusters'][:])\n";
        }
        if (If_show_orientation)
            oss << "Polar_Orientation = np.array(f['Polar_Orientation'][:])\n";
        oss << endl;
        if (If_show_intersection && Intersection_host.size() >= 1)
        {
            oss << "\nintersections = np.array(f['intersections'][:])\n";
            oss << "NumIntersections = intersections.shape[1] if "
                   "intersections.ndim != 1 else 1\n";
            oss << "if intersections.ndim == 1:\n";
            oss << "\tintersections=intersections[:, np.newaxis]\n";
            oss << "Intersection =  "
                   "np.concatenate((np.transpose(intersections[[0, 1, 2], :]), "
                   "np.transpose(intersections[[3, 4, 5], :])), axis=0)\n";
            oss << "Connection_ = list()\n";
            oss << "Connection_ = [(Connection_ + [i, i + NumIntersections]) "
                   "for i in range(NumIntersections)]\n";
            oss << "src = ML.pipeline.scalar_scatter(Intersection[:, 0], "
                   "Intersection[:, 1], Intersection[:, 2])\n";
            oss << "src.mlab_source.dataset.lines = Connection_\n";
            oss << "src.update()\n";
            oss << "lines = ML.pipeline.stripper(src)\n";
            oss << "ML.pipeline.surface(lines, color=(1, 0, 0), line_width=5, "
                   "opacity=1)\n\n";
        }

        oss << "f.close()\n";
        oss << "\n";
        oss << "NUM_Fracs = Frac_NUM_verts.shape[1]\n";
        oss << "\n";
        oss << "poo = 0\n";
        oss << "structure_ = list()\n";
        if (If_show_cluster)
            oss << "scalar_t = np.zeros([verts.shape[1]])\n";

        oss << "for i in range(NUM_Fracs):\n";
        oss << "\tNumVerticesSingleFrac = int(Frac_NUM_verts[0, i])\n";
        if (If_show_cluster)
        {
            //oss << "\tpos = np.where(ListClusters == i + 1)\n";
            oss << "\tpos_t = 0\n";
            oss << "\tindex_yt=0\n";
            oss << "\tfor j in range(NumClusters[0]):\n";
            oss << "\t\tindex_yt = np.where(ListClusters[j] == i + 1)\n";
            oss << "\t\tif index_yt[0].size:\n";
            oss << "\t\t\tpos_t = j\n";
            oss << "\t\t\tbreak\n";
            oss << "\tif not (pos_t + 1 in PercolationClusters):\n";
            oss << "\t\tscalar_t[[s for s in range(poo, poo + "
                   "NumVerticesSingleFrac)]] = pos_t\n";
            oss << "\telse:\n";
            oss << "\t\tscalar_t[[s for s in range(poo, poo + "
                   "NumVerticesSingleFrac)]] = 1.5 * "
                   "NumClusters\n";
        }
        oss << "\tfor j in range(NumVerticesSingleFrac - 2):\n";
        oss << "\t\tstructure_ = structure_ + [[poo, poo + (j + 1) % "
               "NumVerticesSingleFrac, poo + (j + 2) % "
               "NumVerticesSingleFrac]]\n";
        oss << "\tpoo += NumVerticesSingleFrac\n";
        oss << "\n";
        oss << "ML.triangular_mesh(verts[0, :], verts[1, :], verts[2, :], "
               "structure_, scalars=verts[2, :], opacity=1)\n";
        oss << "ML.outline(extent=[-0.5 * DomainDimensionRatio[0] * L_m, 0.5 * "
               "DomainDimensionRatio[0] * L_m, -0.5 * DomainDimensionRatio[1] "
               "* L_m, 0.5 * DomainDimensionRatio[1] * L_m, -0.5 * "
               "DomainDimensionRatio[2] * L_m, 0.5 * DomainDimensionRatio[2] * "
               "L_m])\n";
        oss << "ML.axes()\n";
        oss << "ML.colorbar(orientation='vertical')\n";
        oss << "ML.xlabel('x (m)')\n";
        oss << "ML.ylabel('y (m)')\n";
        oss << "ML.zlabel('z (m)')\n";
        oss << "ML.show()\n\n";

        if (If_show_cluster)
        {
            oss << "ML.triangular_mesh(verts[0, :], verts[1, :], verts[2, :], "
                   "structure_, scalars=scalar_t, opacity=1)\n";
            oss << "ML.outline(extent=[-0.5 * DomainDimensionRatio[0] * L_m, "
                   "0.5 * DomainDimensionRatio[0] * L_m, -0.5 * "
                   "DomainDimensionRatio[1] * L_m, 0.5 * "
                   "DomainDimensionRatio[1] * L_m, -0.5 * "
                   "DomainDimensionRatio[2] * L_m, 0.5 * "
                   "DomainDimensionRatio[2] * L_m])\n";

            oss << "ML.axes()\n";
            oss << "ML.colorbar(orientation='vertical')\n";
            oss << "ML.xlabel('x (m)')\n";
            oss << "ML.ylabel('y (m)')\n";
            oss << "ML.zlabel('z (m)')\n";
            oss << "ML.show()\n";
        }

        if (If_show_orientation)
        {
            oss << "\nimport matplotlib.pyplot as plt\n";
            oss << "fig, ax = plt.subplots(subplot_kw={'projection': "
                   "'polar'})\n";
            oss << "ax.plot(Polar_Orientation[0, :], Polar_Orientation[1, :], "
                   "'ko')\n";
            oss << "ax.set_xticks([0, 1.5708, 3.1416, 4.7124])\n";
            oss << "ax.set_xticklabels([r'0', r\"$0.5\\pi$\", r\"$\\pi$\", "
                   "r\"$1.5\\pi$\"])\n";
            oss << "ax.set_rmax(1.5708)\n";
            oss << "ax.set_rticks([0.5236, 1.0472, 1.5708])\n";
            oss << "ax.set_yticklabels([r\"$0.25\\pi$\", r\"$0.75\\pi$\", "
                   "r\"$\\pi$\"])\n";
            oss << "ax.set_rlabel_position(10)\n";
            oss << "ax.grid(True)\n";
            oss << "ax.set_title(\"Orientations\", va='bottom')\n";
            oss << "plt.show()\n";
        }
        oss.close();
    }

    //----------------------------------------------
    if (command_key != "N")
    {
        std::ofstream oss(command_key, ios::out);
        oss << "clc;\nclose all;\nclear all;\n";
        oss << "currentPath = fileparts(mfilename('fullpath'));\n";
        // oss << "load('" << mat_key << "');\n";
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
        oss << "figure(1); view(3); title('Discete fracture network'); "
               "xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on\n";

        oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 "
               "10 11 12; 13 14 15 16], 'FaceVertexCData', "
               "zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', "
               "'EdgeAlpha', 1, 'facealpha', 0); hold on\n";

        oss << "Frac_NUM_verts = h5read([currentPath, '/" << mat_key
            << "'], '/Frac_NUM_verts');\n";
        oss << "verts = h5read([currentPath, '/" << mat_key
            << "'], '/verts');\n";

        oss << "\nMAX_num_fracs = max(Frac_NUM_verts);\n";
        oss << "NUM_fracs = size(Frac_NUM_verts, 1);\n";
        oss << "element = NaN(NUM_fracs, MAX_num_fracs);\n\n";

        oss << "tmpcc = 1;\n";
        oss << "for i = 1:NUM_fracs\n";
        oss << "\ttmpkk = Frac_NUM_verts(i);\n";
        oss << "\tfor j = 1:tmpkk\n";
        oss << "\t\telement(i, j) = tmpcc; tmpcc = tmpcc + 1;\n";
        oss << "\tend\n";
        oss << "end\n\n";
        oss << "patch('Vertices', verts, 'Faces', element, 'FaceVertexCData', "
               "verts(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, "
               "'facealpha', 1); view(3); colorbar; hold on;\n";

        if (If_show_intersection == true && Intersection_host.size() >= 1)
        {
            oss << "\nintersections = h5read([currentPath, '/" << mat_key
                << "'], '/intersections');\n";
            oss << "if (size(intersections, 2) == 1); "
                   "intersections=intersections'; end;\n";
            oss << "intersections_verts = [intersections(:, [1 2 3]); "
                   "intersections(:, [4 5 6])];\n";
            oss << "intersections_structures = [[1:size(intersections, 1)]', "
                   "[1:size(intersections, 1)]' + size(intersections, 1)];\n";
            oss << "patch('Vertices', intersections_verts, 'Faces', "
                   "intersections_structures, 'FaceVertexCData', "
                   "ones(size(intersections_verts, 1), 1), 'FaceColor', "
                   "'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'linewidth', 3, "
                   "'edgecolor', 'r'); view(3); colorbar; hold on;\n";
        }

        oss << "axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * "
               "DomainDimensionRatio(1) * L, -1.1 / 2 * "
               "DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) "
               "* L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * "
               "DomainDimensionRatio(3) * L]);\n";
        oss << "pbaspect([DomainDimensionRatio]); hold on\n";
        if (If_show_cluster == true)
        {

            oss << "\nNumClusters = h5read([currentPath, '/" << mat_key
                << "'], '/NumClusters');\n";
            oss << "ListClusters = {};\n";
            oss << "for i=1:NumClusters\n";
            oss << "\tListClusters{i, 1} = h5read([currentPath, '/" << mat_key
                << "'], ['/Cluster_', num2str(i)]);\n";
            oss << "end\n";
            oss << "PercolationClusters = h5read([currentPath, '/" << mat_key
                << "'], '/PercolationClusters');\n";
            oss << "figure(2); title('DFN highlighting clusters'); xlabel('x "
                   "(m)'); ylabel('y (m)'); zlabel('z (m)'); view(3); hold "
                   "on\n";
            oss << "axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * "
                   "DomainDimensionRatio(1) * L, -1.1 / 2 * "
                   "DomainDimensionRatio(2) * L, 1.1 / 2 * "
                   "DomainDimensionRatio(2) * L, -1.1 / 2 * "
                   "DomainDimensionRatio(3) * L, 1.1 / 2 * "
                   "DomainDimensionRatio(3) * L]);\n";
            oss << "pbaspect([DomainDimensionRatio]); hold on\n";
            oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 "
                   "8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', "
                   "zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', "
                   "'EdgeAlpha', 1, 'facealpha', 0); hold on\n";
            oss << "colorValue = zeros(size(element, 1), 1);\n";
            oss << "for i = 1:size(ListClusters, 1)\n";
            oss << "\tKK = ListClusters{i, 1};\n";
            oss << "\tif(ismember(PercolationClusters, i))\n";
            oss << "\t\tcolorValue(KK) = 1.2;\n";
            oss << "\telse\n";
            oss << "\t\tcolorValue(KK) = rand;\n";
            oss << "\tend\n";
            oss << "end\n";
            oss << "patch('Vertices', verts, 'Faces', element, "
                   "'FaceVertexCData', colorValue, 'FaceColor', 'flat', "
                   "'EdgeAlpha', 0.9, 'facealpha', 1); view(3); colorbar; hold "
                   "on;\n";
            oss << "colorbar;\n";
        }

        if (If_show_orientation == true)
        {
            oss << "\nfigure(3);\n";
            oss << "Polar_Orientation = h5read([currentPath, '/" << mat_key
                << "'], '/Polar_Orientation');\n";
            oss << "if (size(Polar_Orientation, 2) == 1 && "
                   "size(Polar_Orientation, 1) == 2); "
                   "Polar_Orientation=Polar_Orientation'; end;\n";
            oss << "polarscatter(Polar_Orientation(:, 1), Polar_Orientation(:, "
                   "2), 's', 'filled'); rlim([0 0.5*pi]);\n";
            oss << "rticks([pi / 12, 2 * pi / 12, 3 * pi / 12, 4 * pi / 12, 5 "
                   "* pi / 12, 6 * pi / 12 ]);\n";
            oss << "title(['Fractures', '''',' orientations']); hold on\n";
            //oss << "set(gca,'thetaticklabel',[]);\n";
            oss << "set(gca,'rticklabel',[]);";
        }

        oss << "\n\n\n%if R values have a lognormal distribution, uncommect "
               "the following-------\n";
        oss << "%% R = h5read([currentPath, '/" << mat_key << "'], '/R');\n";
        oss << "%% figure(4);\n";
        oss << "%% nbins = 20;\n";
        oss << "%% histfit(R, nbins, 'lognormal');\n";
        oss << "%% pd=fitdist(R,'lognormal')\n";
        oss << "%% Ex = exp(pd.mu + pd.sigma^2*0.5)\n";
        oss << "%% Dx = exp(2*pd.mu+pd.sigma^2)*(exp(pd.sigma^2)-1)\n";

        oss.close();
    }
}; // MatlabPlotDFN
template cuDFNsys::MatlabPlotDFN<double>::MatlabPlotDFN(
    string mat_key,     // mat file name
    string command_key, // m file name
    thrust::host_vector<cuDFNsys::Fracture<double>>
        Frac_verts_host, // Vector of Fracture
    std::map<pair<size_t, size_t>,
             pair<cuDFNsys::Vector3<double>, cuDFNsys::Vector3<double>>>
        Intersection_host, // vector of Intersection
    std::vector<std::vector<size_t>>
        ListClusters, // List of clusters: a cluster is a vector
    std::vector<size_t>
        Percolation_cluster, // Percolation cluster contains the cluster vector ID
    bool If_show_truncated_frac, // if show truncated fractures
    bool If_show_intersection, bool If_show_cluster, bool If_show_orientation,
    double L, int dir, bool if_python_visualization,
    string PythonName_Without_suffix, double3 DomainDimensionRatio_d);
template cuDFNsys::MatlabPlotDFN<float>::MatlabPlotDFN(
    string mat_key,     // mat file name
    string command_key, // m file name
    thrust::host_vector<cuDFNsys::Fracture<float>>
        Frac_verts_host, // Vector of Fracture
    std::map<pair<size_t, size_t>,
             pair<cuDFNsys::Vector3<float>, cuDFNsys::Vector3<float>>>
        Intersection_host, // vector of Intersection
    std::vector<std::vector<size_t>>
        ListClusters, // List of clusters: a cluster is a vector
    std::vector<size_t>
        Percolation_cluster, // Percolation cluster contains the cluster vector ID
    bool If_show_truncated_frac, // if show truncated fractures
    bool If_show_intersection, bool If_show_cluster, bool If_show_orientation,
    float L, int dir, bool if_python_visualization,
    string PythonName_Without_suffix, double3 DomainDimensionRatio_d);