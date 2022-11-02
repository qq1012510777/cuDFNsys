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

#include "Fractures/RemoveDeadEndFrac.cuh"

// ====================================================
// NAME:        RemoveDeadEndFrac
// DESCRIPTION: Remove dead end fractures that
//              do not conduct flow
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
cuDFNsys::RemoveDeadEndFrac<T>::RemoveDeadEndFrac(std::vector<size_t> &One_cluster,
                                                  std::vector<pair<int, int>> &Intersection_pair,
                                                  const size_t &dir,
                                                  thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs,
                                                  std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> Intersection_map)
{
    size_t inlet_id = dir * 2, outlet_id = dir * 2 + 1;

    size_t frac_max = *std::max_element(One_cluster.begin(), One_cluster.end());

    Eigen::SparseMatrix<double> adjacent_matrix(frac_max + 3, frac_max + 3);
    adjacent_matrix.reserve(VectorXi::Constant(frac_max + 3, 5));

    for (size_t i = 0; i < One_cluster.size(); ++i)
    {
        size_t FracTag = One_cluster[i];

        if (Fracs[FracTag].ConnectModelSurf[inlet_id])
        {
            adjacent_matrix.insert(FracTag, frac_max + 1) = 2;
            adjacent_matrix.insert(frac_max + 1, FracTag) = 2; // the second-last
        }

        if (Fracs[FracTag].ConnectModelSurf[outlet_id])
        {
            adjacent_matrix.insert(FracTag, frac_max + 2) = 2;
            adjacent_matrix.insert(frac_max + 2, FracTag) = 2; // the last one
        }
    }

    std::sort(One_cluster.begin(), One_cluster.end(), std::greater<size_t>());

    for (size_t i = 0; i < One_cluster.size() - 1; ++i)
    {
        size_t FracTag_1 = One_cluster[i];

        for (size_t j = i + 1; j < One_cluster.size(); ++j)
        {
            size_t FracTag_2 = One_cluster[j];

            //-------

            auto ity = Intersection_map.find(std::make_pair(FracTag_1, FracTag_2));

            if (ity != Intersection_map.end())
            {
                adjacent_matrix.insert(FracTag_2, FracTag_1) = 2;
                adjacent_matrix.insert(FracTag_1, FracTag_2) = 2;
            }
        }
    }

    adjacent_matrix.makeCompressed();
    bool tees = true;

    vector<bool> if_deleted(frac_max + 1, false);

    //int Non_zeros = adjacent_matrix.nonZeros();

    while (tees == true)
    {
        for (size_t i = 0; i <= frac_max; ++i)
        {
            SparseVector<double> sub_mat = adjacent_matrix.innerVector(i);

            if (sub_mat.nonZeros() == 1 && if_deleted[i] == false) // only connect to one frac
            {
                //cout << "found i = " << i << endl;
                std::vector<size_t>::iterator itr = std::find(One_cluster.begin(), One_cluster.end(), i);
                //cout << "itr: " << *itr << endl;
                size_t local_index = std::distance(One_cluster.begin(), itr);

                //cout << "local_index: " << local_index << "; one_cluster size: " << One_cluster.size() << endl;

                One_cluster.erase(One_cluster.begin() + local_index);
                //cout << "remove frac: " << i << endl;

                if_deleted[i] = true;

                adjacent_matrix.row(i) *= 0;
                adjacent_matrix.col(i) *= 0;

                adjacent_matrix.prune(0.01);

                //cout << "Non_zeros 1: " << Non_zeros << ", Non_zeros 2: " << adjacent_matrix.nonZeros() << endl;
                //Non_zeros = adjacent_matrix.nonZeros();
                break;
            }

            if (i == frac_max)
                tees = false;
        }
    }

    int DSIZE = One_cluster.size();

    Intersection_pair.reserve(DSIZE * floor((DSIZE - 1) / 2) + (DSIZE - 1) % 2 * DSIZE * 0.5);

    for (size_t i = 0; i <= frac_max; ++i)
    {
        SparseVector<double> sub_mat = adjacent_matrix.innerVector(i);
        if (sub_mat.nonZeros() > 0)
            for (SparseVector<double>::InnerIterator it(sub_mat); it; ++it)
                if (it.row() < i)
                    Intersection_pair.push_back(std::make_pair(i, it.row()));
                else if (it.row() >= frac_max + 1)
                    break;
    }

    Intersection_pair.shrink_to_fit();

    //  delete deadend fractures in host_vector of cuDFNsys::ractures
    //  thrust::host_vector<cuDFNsys::Fracture> &Fracs
    std::map<int, int> Map_ID1_to_ID2;
    thrust::host_vector<cuDFNsys::Fracture<T>> FracsII(DSIZE);
    for (size_t i = 0; i < DSIZE; ++i)
    {
        FracsII[i] = Fracs[One_cluster[i]];
        Map_ID1_to_ID2.insert(std::make_pair(One_cluster[i], i));
    }
    Fracs.clear();
    Fracs.resize(DSIZE);
    Fracs = FracsII;
    Fracs.shrink_to_fit();
    One_cluster.shrink_to_fit();

    // std::vector<size_t> &One_cluster,
    // std::vector<pair<int, int>> &Intersection_pair

    for (size_t i = 0; i < One_cluster.size(); ++i)
        One_cluster[i] = Map_ID1_to_ID2[One_cluster[i]];

    for (size_t i = 0; i < Intersection_pair.size(); ++i)
    {
        Intersection_pair[i].first = Map_ID1_to_ID2[Intersection_pair[i].first];
        Intersection_pair[i].second = Map_ID1_to_ID2[Intersection_pair[i].second];
    }
}; // RemoveDeadEndFrac
template cuDFNsys::RemoveDeadEndFrac<double>::RemoveDeadEndFrac(std::vector<size_t> &One_cluster,
                                                                std::vector<pair<int, int>> &Intersection_pair,
                                                                const size_t &dir,
                                                                thrust::host_vector<cuDFNsys::Fracture<double>> &Fracs,
                                                                std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<double>, cuDFNsys::Vector3<double>>> Intersection_map);
template cuDFNsys::RemoveDeadEndFrac<float>::RemoveDeadEndFrac(std::vector<size_t> &One_cluster,
                                                               std::vector<pair<int, int>> &Intersection_pair,
                                                               const size_t &dir,
                                                               thrust::host_vector<cuDFNsys::Fracture<float>> &Fracs,
                                                               std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<float>, cuDFNsys::Vector3<float>>> Intersection_map);