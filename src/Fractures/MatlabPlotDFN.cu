#include "Fractures/MatlabPlotDFN.cuh"

// ====================================================
// NAME:        cuDFNsys::MatlabPlotDFN::MatlabPlotDFN
// DESCRIPTION: Plot DFN with m and mat file
//              m is command file; mat is data file
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
cuDFNsys::MatlabPlotDFN::MatlabPlotDFN(string mat_key,
                                       string command_key,
                                       thrust::host_vector<cuDFNsys::Fracture> Frac_verts_host,
                                       MapIntersection Intersection_host,
                                       std::vector<std::vector<size_t>> ListClusters,
                                       std::vector<size_t> Percolation_cluster,
                                       bool If_show_truncated_frac,
                                       bool If_show_intersection,
                                       bool If_show_cluster,
                                       bool If_show_orientation,
                                       float L,
                                       int dir)
{
    cuDFNsys::MatlabAPI M1;

    int NUM_Frac = Frac_verts_host.size();
    int sum_NUM_verts = 4 * NUM_Frac;
    float *Frac_NUM_verts = new float[NUM_Frac];
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

    M1.WriteMat(mat_key, "w", NUM_Frac, NUM_Frac, 1, Frac_NUM_verts, "Frac_NUM_verts");
    delete[] Frac_NUM_verts;
    Frac_NUM_verts = NULL;
    //-----------------
    float *R_ = new float[NUM_Frac];
    if (R_ == NULL)
    {
        string AS = "Alloc error in MatlabPlotDFN!";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    for (int i = 0; i < NUM_Frac; ++i)
        R_[i] = Frac_verts_host[i].Radius;

    M1.WriteMat(mat_key, "u", NUM_Frac, NUM_Frac, 1,
                R_, "R");
    delete[] R_;
    R_ = NULL;
    //----------------------
    float *verts = new float[sum_NUM_verts * 3];
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
                verts[tmp_cc + sum_NUM_verts * 2] = Frac_verts_host[i].Verts3D[j].z;
                tmp_cc++;
            }
            else
            {
                verts[tmp_cc] = Frac_verts_host[i].Verts3DTruncated[j].x;
                verts[tmp_cc + sum_NUM_verts] = Frac_verts_host[i].Verts3DTruncated[j].y;
                verts[tmp_cc + sum_NUM_verts * 2] = Frac_verts_host[i].Verts3DTruncated[j].z;
                tmp_cc++;
            }
        }
    }

    M1.WriteMat(mat_key, "u", sum_NUM_verts * 3, sum_NUM_verts, 3,
                verts, "verts");
    delete[] verts;
    verts = NULL;

    //--------------------------------------
    if (If_show_intersection == true)
    {
        int NUM_intersections = Intersection_host.size();

        float *Intersections_ = new float[NUM_intersections * 6];
        if (Intersections_ == NULL)
        {
            string AS = "Alloc error in MatlabPlotDFN!";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        int tmp_jj = 0;
        for (std::map<pair<size_t, size_t>, pair<float3, float3>>::iterator i = Intersection_host.begin();
             i != Intersection_host.end(); ++i)
        {
            Intersections_[tmp_jj] = i->second.first.x;
            Intersections_[tmp_jj + NUM_intersections] = i->second.first.y;
            Intersections_[tmp_jj + NUM_intersections * 2] = i->second.first.z;
            Intersections_[tmp_jj + NUM_intersections * 3] = i->second.second.x;
            Intersections_[tmp_jj + NUM_intersections * 4] = i->second.second.y;
            Intersections_[tmp_jj + NUM_intersections * 5] = i->second.second.z;

            tmp_jj++;
        }

        M1.WriteMat(mat_key, "u", NUM_intersections * 6, NUM_intersections, 6,
                    Intersections_, "intersections");

        delete[] Intersections_;
        Intersections_ = NULL;
    }

    if (If_show_cluster == true)
    {
        size_t Max_cluster_size = 0;
        for (size_t i = 0; i < ListClusters.size(); ++i)
            Max_cluster_size = ListClusters[i].size() > Max_cluster_size ? ListClusters[i].size() : Max_cluster_size;

        float *clustersList = new float[ListClusters.size() * Max_cluster_size];
        if (clustersList == NULL)
        {
            string AS = "Alloc error in MatlabPlotDFN!";
            throw cuDFNsys::ExceptionsPause(AS);
        }
        //cout << "Max_cluster_size: " << Max_cluster_size << endl;
        int tmpff = 0;
        for (size_t j = 0; j < Max_cluster_size; ++j)
            for (size_t i = 0; i < ListClusters.size(); ++i)
            {
                if (j < ListClusters[i].size())
                    clustersList[tmpff] = ListClusters[i][j] + 1;
                else
                    clustersList[tmpff] = -1.0;

                tmpff++;
            }

        M1.WriteMat(mat_key, "u", ListClusters.size() * Max_cluster_size, ListClusters.size(), Max_cluster_size, clustersList, "ListClusters");

        delete[] clustersList;
        clustersList = NULL;
        //-----------

        float *percolationcluster = new float[Percolation_cluster.size()];
        if (percolationcluster == NULL)
        {
            string AS = "Alloc error in MatlabPlotDFN!";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        for (size_t i = 0; i < Percolation_cluster.size(); ++i)
            percolationcluster[i] = Percolation_cluster[i] + 1;

        M1.WriteMat(mat_key, "u", Percolation_cluster.size(), Percolation_cluster.size(), 1, percolationcluster, "PercolationClusters");

        delete[] percolationcluster;
        percolationcluster = NULL;
    };

    if (If_show_orientation == true)
    {
        float *orientation = new float[NUM_Frac * 2];
        if (orientation == NULL)
        {
            string AS = "Alloc error in MatlabPlotDFN!";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        for (size_t i = 0; i < NUM_Frac; ++i)
        {
            orientation[i] = Frac_verts_host[i].Phi();
            orientation[i + NUM_Frac] = Frac_verts_host[i].Theta();
        }
        M1.WriteMat(mat_key, "u", NUM_Frac * 2, NUM_Frac, 2, orientation, "Polar_Orientation");
        delete[] orientation;
        orientation = NULL;
    };

    //----------------------------------------------
    std::ofstream oss(command_key, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    oss << "load('" << mat_key << "');\n";
    oss << "L = 0.5 * " << L << ";\n";
    oss << "cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];\n";
    oss << "figure(1); view(3); title('Discete fracture network'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on\n";
    oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on\n";

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
    oss << "patch('Vertices', verts, 'Faces', element, 'FaceVertexCData', verts(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1); view(3); colorbar; hold on;\n";

    if (If_show_intersection == true)
    {
        oss << "\nintersections_verts = [intersections(:, [1 2 3]); intersections(:, [4 5 6])];\n";
        oss << "intersections_structures = [[1:size(intersections, 1)]', [1:size(intersections, 1)]' + size(intersections, 1)];\n";
        oss << "patch('Vertices', intersections_verts, 'Faces', intersections_structures, 'FaceVertexCData', ones(size(intersections_verts, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'linewidth', 3, 'edgecolor', 'r'); view(3); colorbar; hold on;\n";
    }

    oss << "axis([-1.1 * L,  1.1 * L, -1.1 * L, 1.1 * L, -1.1 * L, 1.1 * L]);\n";

    if (If_show_cluster == true)
    {
        oss << "\nfigure(2); title('DFN highlighting clusters'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); view(3); axis([-1.1 * L,  1.1 * L, -1.1 * L, 1.1 * L, -1.1 * L, 1.1 * L]); hold on\n";
        //oss << "ListClusters(ListClusters==-1) = NaN;\n";
        oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on\n";
        oss << "colorValue = zeros(size(element, 1), 1);\n";
        oss << "for i = 1:size(ListClusters, 1)\n";
        oss << "\tKK = ListClusters(i, :); KK(KK == -1) = [];\n";
        oss << "\tif(ismember(PercolationClusters, i))\n";
        oss << "\t\tcolorValue(KK) = 1.2;\n";
        oss << "\telse\n";
        oss << "\t\tcolorValue(KK) = rand;\n";
        oss << "\tend\n";
        oss << "end\n";
        oss << "patch('Vertices', verts, 'Faces', element, 'FaceVertexCData', colorValue, 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1); view(3); colorbar; hold on;\n";
        oss << "colorbar;\n";
    }

    if (If_show_orientation == true)
    {
        oss << "\nfigure(3);\n";
        oss << "polarscatter(Polar_Orientation(:, 1), Polar_Orientation(:, 2), 's', 'filled'); rlim([0 0.5*pi]);\n";
        oss << "rticks([pi / 12, 2 * pi / 12, 3 * pi / 12, 4 * pi / 12, 5 * pi / 12, 6 * pi / 12 ]);\n";
        oss << "title(['Fractures', '''',' orientations']); hold on\n";
        //oss << "set(gca,'thetaticklabel',[]);\n";
        oss << "set(gca,'rticklabel',[]);";
    }

    oss << "\n\n\n%if R values have a lognormal distribution, uncommect the following-------\n";
    oss << "%figure(4);\n";
    oss << "%nbins = 20;\n";
    oss << "%histfit(R, nbins, 'lognormal');\n";
    oss << "%pd=fitdist(R,'lognormal')\n";
    oss << "%Ex = exp(pd.mu + pd.sigma^2*0.5)\n";
    oss << "%Dx = exp(2*pd.mu+pd.sigma^2)*(exp(pd.sigma^2)-1)\n";

    oss.close();
}; // MatlabPlotDFN