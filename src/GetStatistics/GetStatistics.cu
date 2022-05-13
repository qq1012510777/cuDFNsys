#include "GetStatistics/GetStatistics.cuh"

// ====================================================
// NAME:        GetStatistics
// DESCRIPTION: Get P30, P32, permeability and so on
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
void cuDFNsys::GetStatistics(const thrust::host_vector<cuDFNsys::Fracture> &Frac_verts_host,
                             const std::map<pair<size_t, size_t>, pair<float3, float3>> &Intersection_map,
                             const std::vector<std::vector<size_t>> &ListClusters,
                             const std::vector<size_t> &Percolation_cluster,
                             const float L,
                             float &P33_total_B,
                             float &P33_connected_B,
                             float &Ratio_of_P33_B,
                             float &P33_largest_cluster_B,
                             float &P32_total_B,
                             float &P32_connected_B,
                             float &Ratio_of_P32_B,
                             float &P32_largest_cluster_B,
                             float &P30_B,
                             float &P30_connected_B,
                             float &Ratio_of_P30_B,
                             float &P30_largest_cluster_B,
                             float &Percolation_probability_B,
                             float &n_I_B)
{
    P32_total_B = 0;
    P30_B = 0;
    P32_connected_B = 0;

    float volume_DFN = L * L * L;

    float area_total = 0;
    float volume_total = 0;
    for (int i = 0; i < Frac_verts_host.size(); ++i)
    {
        area_total += pow(Frac_verts_host[i].Radius, 2) * 2;
        volume_total += pow(Frac_verts_host[i].Radius, 2) * 2 * pow(Frac_verts_host[i].Conductivity * 12.0, 1.0 / 3.0);
        //cout << pow(Frac_verts_host[i].Conductivity * 12, 1 / 3) << ", ";
        //printf("%.17f, ", pow(Frac_verts_host[i].Conductivity * 12, 1.0 / 3.0));
    }
    //cout << endl;

    vector<size_t> FRACS_perco;

    for (int i = 0; i < Percolation_cluster.size(); ++i)
    {
        FRACS_perco.insert(FRACS_perco.end(), ListClusters[Percolation_cluster[i]].begin(),
                           ListClusters[Percolation_cluster[i]].end());
    }

    float area_K = 0;
    float volume_K = 0;
    for (int i = 0; i < FRACS_perco.size(); ++i)
    {
        area_K += pow(Frac_verts_host[FRACS_perco[i]].Radius, 2) * 2;
        volume_K += pow(Frac_verts_host[FRACS_perco[i]].Radius, 2) * 2 * pow(Frac_verts_host[FRACS_perco[i]].Conductivity * 12, 1.0 / 3.0);
    }

    int Tag_largest_cluster = 0;
    int clustersize_max = 0;
    for (int i = 0; i < ListClusters.size(); ++i)
    {
        if (clustersize_max < ListClusters[i].size())
        {
            clustersize_max = ListClusters[i].size();
            Tag_largest_cluster = i;
        }
    }
    float area_s = 0;
    float volume_s = 0;
    for (int i = 0; i < clustersize_max; ++i)
    {
        area_s += pow(Frac_verts_host[ListClusters[Tag_largest_cluster][i]].Radius, 2) * 2;
        volume_s += pow(Frac_verts_host[ListClusters[Tag_largest_cluster][i]].Radius, 2) * 2 *
                    pow(Frac_verts_host[ListClusters[Tag_largest_cluster][i]].Conductivity * 12, 1.0 / 3.0);
    };

    P32_total_B = area_total / volume_DFN;
    P32_connected_B = area_K / volume_DFN;
    Ratio_of_P32_B = area_K / area_total;
    P32_largest_cluster_B = area_s / volume_DFN;

    P30_B = Frac_verts_host.size() / volume_DFN;
    P30_connected_B = FRACS_perco.size() / volume_DFN;
    Ratio_of_P30_B = FRACS_perco.size() * 1.0 / Frac_verts_host.size();
    P30_largest_cluster_B = clustersize_max * 1.0 / Frac_verts_host.size();

    P33_total_B = volume_total / volume_DFN;
    P33_connected_B = volume_K / volume_DFN;
    Ratio_of_P33_B = volume_K / volume_total;
    P33_largest_cluster_B = volume_s / volume_DFN;

    if (Percolation_cluster.size() > 0)
        Percolation_probability_B = 1;
    else
        Percolation_probability_B = 0;

    n_I_B = 1.0f * Intersection_map.size() / (Frac_verts_host.size() * 1.0f);
};