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

#include "cuDFNsys.cuh"

// ====================================================
// NAME:        cuDFNsys::DFN<T>::FractureGeneration
// DESCRIPTION: create a DFN
// AUTHOR:      Tingchang YIN
// DATE:        20/10/2023
// ====================================================
template <typename T>
void cuDFNsys::DFN<T>::FractureGeneration()
{
    this->NumFracturesTotal = 0;
    // the total number of fractures
    for (int i = 0; i < this->NumFractures.size(); ++i)
        this->NumFracturesTotal += this->NumFractures[i];
    if (this->NumFracturesTotal == 0)
        throw cuDFNsys::ExceptionsPause(
            "The number of fractures of each group is not set!");
    this->FracturesHost.resize(NumFracturesTotal);
    this->FracturesDevice = this->FracturesHost;

    int TmpCount = 0;
    for (int i = 0; i < this->NumFractures.size(); ++i)
    {
        thrust::host_vector<cuDFNsys::Fracture<T>> Frac_host_sub(
            this->NumFractures[i]);
        thrust::device_vector<cuDFNsys::Fracture<T>> Frac_device_sub(
            this->NumFractures[i]);
        cuDFNsys::Fracture<T> *Frac_device_sub_ptr =
            thrust::raw_pointer_cast(Frac_device_sub.data());

        cuDFNsys::Fractures<T><<<this->NumFractures[i] / 256 + 1, 256>>>(
            Frac_device_sub_ptr,   // the pointer to device vector
            this->RandomSeed,      // seed
            this->NumFractures[i], // number of fracture
            this->DomainSizeX,     // domain size
            this->ModeOfSizeDistribution
                [i], // distribution pattern of fracture sizes
            this->SizeDistributionParameters
                [i],        // parameters of distribution of fracture sizes
            this->Kappa[i], // kappa value of fisher distribution
            this->Beta[i],  // beta
            this->Gamma[i], // gamma
            this->DomainDimensionRatio, // ratio of domain dimensions
            this->MeanOrientationOfFisherDistribution[i]); // mean orientations

        cudaDeviceSynchronize();
        // wait until the device function finish
        Frac_host_sub = Frac_device_sub;
        thrust::copy(thrust::host, Frac_host_sub.begin(), Frac_host_sub.end(),
                     this->FracturesHost.begin() + TmpCount);
        TmpCount += this->NumFractures[i];
    };
    this->FracturesDevice = this->FracturesHost;
    this->FracturesDevicePtr =
        thrust::raw_pointer_cast(this->FracturesDevice.data());
}; // FractureGeneration
template void cuDFNsys::DFN<double>::FractureGeneration();
template void cuDFNsys::DFN<float>::FractureGeneration();

// ====================================================
// NAME:        cuDFNsys::DFN<T>::IdentifyIntersectionsClusters
// DESCRIPTION: intersections and clusters
// AUTHOR:      Tingchang YIN
// DATE:        20/10/2023
// ====================================================
template <typename T>
void cuDFNsys::DFN<T>::IdentifyIntersectionsClusters(
    const bool &IfTruncatedFractures)
{
    this->ListClusters.clear();
    this->PercolationCluster.clear();
    this->IntersectionMap.clear();

    cuDFNsys::IdentifyIntersection<T> identifyInters{
        (size_t)this->NumFracturesTotal, // number of fractures
        this->FracturesDevicePtr, // pointer of device vector of fractures
        IfTruncatedFractures, // if you want to use truncated fractures? here is false,
        this->IntersectionMap};

    if (this->IntersectionMap.size() > 0 && this->FracturesHost.size() > 1)
    {
        cuDFNsys::Graph<T> G{(size_t)this->NumFracturesTotal,
                             this->IntersectionMap};
        G.UseDFS(this->ListClusters);
        cuDFNsys::IdentifyPercolationCluster<T> IdentiClu{
            this->ListClusters,        // all clusters
            this->FracturesHost,       // host vector of fractures
            this->PercoDir,            // percolation direction / flow direction
            this->PercolationCluster}; // percolation cluster
    }
    else if (this->IntersectionMap.size() == 0)
    {
        this->ListClusters.resize(this->FracturesHost.size());
        for (int i = 0; i < this->FracturesHost.size(); ++i)
        {
            this->ListClusters[i].push_back(i);

            if (this->FracturesHost[i].ConnectModelSurf[this->PercoDir * 2] &&
                this->FracturesHost[i].ConnectModelSurf[this->PercoDir * 2 + 1])
                this->PercolationCluster.push_back(i);
        }
    }
}; // IdentifyIntersectionsClusters
template void cuDFNsys::DFN<double>::IdentifyIntersectionsClusters(
    const bool &IfTruncatedFractures);
template void cuDFNsys::DFN<float>::IdentifyIntersectionsClusters(
    const bool &IfTruncatedFractures);

// ====================================================
// NAME:        cuDFNsys::DFN<T>::FractureGeneration
// DESCRIPTION: Visualization
// AUTHOR:      Tingchang YIN
// DATE:        20/10/2023
// ====================================================
template <typename T>
void cuDFNsys::DFN<T>::Visualization(const string &MatlabScriptName,
                                     const string &PythonScriptName,
                                     const string &HDF5FileName,
                                     const bool &IfShowTruncatedFractures,
                                     const bool &IfShowIntersections,
                                     const bool &IfHightlighAllClusters,
                                     const bool &IfShowOrientationDistribution)
{
    cuDFNsys::MatlabPlotDFN<T> As{
        HDF5FileName + ".h5", // the .h5 file, with suffix
        MatlabScriptName +
            ".m",            // matlab script to visualize the DFN, with suffix
        this->FracturesHost, // host vector of fractures
        this->IntersectionMap,         // intersection map
        this->ListClusters,            // clusters
        this->PercolationCluster,      // No. or say, ID, of percolating cluster
        IfShowTruncatedFractures,      // if show truncated fractures?
        IfShowIntersections,           // if show intersections?
        IfHightlighAllClusters,        // if show clusters?
        IfShowOrientationDistribution, // if show orientations data?
        this->DomainSizeX,             // domain size
        this->PercoDir,                // flow direction
        true, // true means I also want to see DFN with python script, a .py file will be generated
        PythonScriptName, // the name of python script, without suffix
        this->DomainDimensionRatio};
}; // Visualization
template void cuDFNsys::DFN<double>::Visualization(
    const string &MatlabScriptName, const string &PythonScriptName,
    const string &HDF5FileName, const bool &IfShowTruncatedFractures,
    const bool &IfShowIntersections, const bool &IfHightlighAllClusters,
    const bool &IfShowOrientationDistribution);
template void cuDFNsys::DFN<float>::Visualization(
    const string &MatlabScriptName, const string &PythonScriptName,
    const string &HDF5FileName, const bool &IfShowTruncatedFractures,
    const bool &IfShowIntersections, const bool &IfHightlighAllClusters,
    const bool &IfShowOrientationDistribution);

// ====================================================
// NAME:        cuDFNsys::DFN<T>::StoreInH5
// DESCRIPTION: store this class in .h5
// AUTHOR:      Tingchang YIN
// DATE:        21/10/2023
// ====================================================
template <typename T>
void cuDFNsys::DFN<T>::StoreInH5(const string &ClassNameH5)
{
    cuDFNsys::OutputObjectData<T> lk_out;
    lk_out.OutputFractures(ClassNameH5 + ".h5", this->FracturesHost,
                           this->DomainSizeX, this->DomainDimensionRatio);

    cuDFNsys::HDF5API h5gg;
    {
        int NUM_intersections = this->IntersectionMap.size();

        T *Intersections_ = new T[NUM_intersections * 8];
        if (Intersections_ == NULL)
        {
            string AS = "Alloc error in MatlabPlotDFN!";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        int tmp_jj = 0;
        for (auto i = this->IntersectionMap.begin();
             i != this->IntersectionMap.end(); ++i)
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

        uint2 dim_f = make_uint2(8, NUM_intersections);
        h5gg.AddDataset(ClassNameH5 + ".h5", "N", "Intersections",
                        Intersections_, dim_f);
        delete[] Intersections_;
        Intersections_ = NULL;

        h5gg.AddDataset(ClassNameH5 + ".h5", "N", "PercoDir", &this->PercoDir,
                        make_uint2(1, 0));
    }
    {
        size_t Max_cluster_size = 0;
        for (size_t i = 0; i < this->ListClusters.size(); ++i)
            Max_cluster_size = this->ListClusters[i].size() > Max_cluster_size
                                   ? this->ListClusters[i].size()
                                   : Max_cluster_size;

        T *clustersList = new T[this->ListClusters.size() * Max_cluster_size];
        if (clustersList == NULL)
        {
            string AS = "Alloc error in MatlabPlotDFN!";
            throw cuDFNsys::ExceptionsPause(AS);
        }
        //cout << "Max_cluster_size: " << Max_cluster_size << endl;
        int tmpff = 0;
        for (size_t j = 0; j < Max_cluster_size; ++j)
            for (size_t i = 0; i < this->ListClusters.size(); ++i)
            {
                if (j < this->ListClusters[i].size())
                    clustersList[tmpff] = this->ListClusters[i][j] + 1;
                else
                    clustersList[tmpff] = -1.0;

                tmpff++;
            }

        //M1.WriteMat(mat_key, "u", ListClusters.size() * Max_cluster_size, ListClusters.size(), Max_cluster_size, clustersList, "ListClusters");
        uint2 dim_f = make_uint2(Max_cluster_size, this->ListClusters.size());
        h5gg.AddDataset(ClassNameH5 + ".h5", "N", "ListClusters", clustersList,
                        dim_f);

        delete[] clustersList;
        clustersList = NULL;
        //-----------

        T *percolationcluster = new T[this->PercolationCluster.size()];
        if (percolationcluster == NULL)
        {
            string AS = "Alloc error in MatlabPlotDFN!";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        for (size_t i = 0; i < this->PercolationCluster.size(); ++i)
            percolationcluster[i] = this->PercolationCluster[i] + 1;

        dim_f = make_uint2(1, this->PercolationCluster.size());
        h5gg.AddDataset(ClassNameH5 + ".h5", "N", "PercolationClusters",
                        percolationcluster, dim_f);

        delete[] percolationcluster;
        percolationcluster = NULL;

        int NumberClustersyy = (int)this->ListClusters.size();
        h5gg.AddDataset(ClassNameH5 + ".h5", "N", "NumClusters",
                        &NumberClustersyy, make_uint2(1, 0));

        uint OKLDS = (this->IfPeriodic ? 1 : 0);
        h5gg.AddDataset(ClassNameH5 + ".h5", "N", "IfPeriodic", &OKLDS,
                        make_uint2(1, 0));
    }

}; //
template void cuDFNsys::DFN<double>::StoreInH5(const string &ClassNameH5);
template void cuDFNsys::DFN<float>::StoreInH5(const string &ClassNameH5);

// ====================================================
// NAME:        cuDFNsys::DFN<T>::LoadClassFromH5
// DESCRIPTION: Load class data from a h5 file
// AUTHOR:      Tingchang YIN
// DATE:        21/10/2023
// ====================================================
template <typename T>
void cuDFNsys::DFN<T>::LoadClassFromH5(const string &ClassNameH5)
{
    cuDFNsys::InputObjectData<T> lk;
    lk.InputFractures(ClassNameH5 + ".h5", this->FracturesHost,
                      this->DomainSizeX, this->DomainDimensionRatio);
    this->FracturesDevice = this->FracturesHost;
    this->FracturesDevicePtr =
        thrust::raw_pointer_cast(this->FracturesDevice.data());
    this->NumFracturesTotal = this->FracturesHost.size();
    //------
    {
        cuDFNsys::HDF5API hdf5Class;
        // std::vector<int> NumClusters =
        //     hdf5Class.ReadDataset<int>(ClassNameH5 + ".h5", "N", "NumClusters");
        //
        // this->ListClusters.resize(NumClusters[0]);
        //
        // std::vector<T> Temp_Variable =
        //     hdf5Class.ReadDataset<T>(ClassNameH5 + ".h5", "N", "ListClusters");
        //
        // for (int i = 0; i < NumClusters[0]; ++i)
        //     this->ListClusters[i].reserve(Temp_Variable.size() /
        //                                   NumClusters[0]);
        // for (int i = 0; i < Temp_Variable.size(); ++i)
        //     if (Temp_Variable[i] != -1)
        //         this->ListClusters[i % NumClusters[0]].push_back(
        //             Temp_Variable[i] - 1);
        // for (int i = 0; i < NumClusters[0]; ++i)
        //     this->ListClusters[i].shrink_to_fit();
        //
        // Temp_Variable = hdf5Class.ReadDataset<T>(ClassNameH5 + ".h5", "N",
        //                                          "PercolationClusters");
        // this->PercolationCluster.resize(Temp_Variable.size());
        // for (int i = 0; i < Temp_Variable.size(); ++i)
        //     this->PercolationCluster[i] = Temp_Variable[i] - 1;
        //
        // //------INTERSECTION PAIR---
        // Temp_Variable =
        //     hdf5Class.ReadDataset<T>(ClassNameH5 + ".h5", "N", "Intersections");
        // int NumIntersections = Temp_Variable.size() / 8;
        // for (int i = 0; i < NumIntersections; ++i)
        //     this->IntersectionMap.insert(std::make_pair(
        //         std::make_pair((uint)Temp_Variable[i + 6 * NumIntersections],
        //                        (uint)Temp_Variable[i + 7 * NumIntersections]),
        //         std::make_pair(cuDFNsys::MakeVector3(
        //                            Temp_Variable[i],
        //                            Temp_Variable[i + NumIntersections],
        //                            Temp_Variable[i + NumIntersections * 2]),
        //                        cuDFNsys::MakeVector3(
        //                            Temp_Variable[i + NumIntersections * 3],
        //                            Temp_Variable[i + NumIntersections * 4],
        //                            Temp_Variable[i + NumIntersections * 5]))));
        //------------
        std::vector<int> UO =
            hdf5Class.ReadDataset<int>(ClassNameH5 + ".h5", "N", "PercoDir");
        this->PercoDir = UO[0];

        std::vector<uint> UF =
            hdf5Class.ReadDataset<uint>(ClassNameH5 + ".h5", "N", "IfPeriodic");
        this->IfPeriodic = (UF[0] == 0 ? false : true);
    };
};
template void cuDFNsys::DFN<double>::LoadClassFromH5(const string &ClassNameH5);
template void cuDFNsys::DFN<float>::LoadClassFromH5(const string &ClassNameH5);

// ====================================================
// NAME:        cuDFNsys::DFN<T>::SpatialPeriodicity
// DESCRIPTION: make the DFN become spatially periodic
// AUTHOR:      Tingchang YIN
// DATE:        06/11/2023
// ====================================================
template <typename T>
void cuDFNsys::DFN<T>::SpatialPeriodicity()
{
    int O_size = this->FracturesHost.size();

    this->IfPeriodic = true;

    thrust::host_vector<cuDFNsys::Fracture<T>> Fractures_copy =
        this->FracturesHost;

    this->FracturesHost.insert(this->FracturesHost.end(),
                               Fractures_copy.begin(), Fractures_copy.end());
    Fractures_copy.clear();
    this->FracturesDevice = this->FracturesHost;

    this->FracturesDevicePtr =
        thrust::raw_pointer_cast(this->FracturesDevice.data());

    this->NumFracturesTotal *= 2;

    cuDFNsys::FracturesForSpatialPeriodicity<T>
        <<<this->NumFracturesTotal / 256 + 1, 256>>>(
            this->FracturesDevicePtr, (int)this->NumFracturesTotal,
            this->DomainSizeX, this->PercoDir, this->DomainDimensionRatio);
    cudaDeviceSynchronize();
    this->FracturesHost = this->FracturesDevice;
    this->FracturesDevicePtr =
        thrust::raw_pointer_cast(this->FracturesDevice.data());

    if (this->PercoDir == 0)
        this->DomainSizeX *= 2, this->DomainDimensionRatio.y /= 2,
            this->DomainDimensionRatio.z /= 2;
    else if (this->PercoDir == 1)
        this->DomainDimensionRatio.y *= 2;
    else
        this->DomainDimensionRatio.z *= 2;

    cout << "cuDFNsys::DFN<T>::SpatialPeriodicity: The DFN becomes spatially "
            "periodic now! The domain is elongated in the percolation "
            "direction. The total number of "
            "fracture is changed to "
         << this->NumFracturesTotal << " from " << O_size << endl;
};
template void cuDFNsys::DFN<double>::SpatialPeriodicity();
template void cuDFNsys::DFN<float>::SpatialPeriodicity();

// ====================================================
// NAME:        cuDFNsys::DFN<T>::LoadDFNFromCSV
// DESCRIPTION: load input parameters of DFNs from .csv file
// AUTHOR:      Tingchang YIN
// DATE:        03/11/2023
// ====================================================
template <typename T>
void cuDFNsys::DFN<T>::LoadDFNFromCSV(const string &xlsxNameWithoutSuffix)
{
    std::ifstream ifs;
    ifs.open(xlsxNameWithoutSuffix + ".csv", ios::in);

    if (!ifs.is_open())
    {
        ifs.close();
        throw cuDFNsys::ExceptionsPause("Fail to open " +
                                        xlsxNameWithoutSuffix + ".csv");
    }

    vector<string> items(12);
    for (int i = 0; i < items.size(); ++i)
    {
        string temp_k;
        getline(ifs, temp_k);
        items[i] = temp_k;
        //cout << temp_k << endl;
    }
    ifs.close();

    string temp_s;
    //--------line 2
    istringstream istrss(items[0]);
    getline(istrss, temp_s, ',');
    getline(istrss, temp_s, ',');
    bool IfStochastic = (atoi(temp_s.c_str()) == 1 ? true : false);

    if (IfStochastic)
    {
        //--------line 2
        istringstream istrss(items[1]);
        getline(istrss, temp_s, ',');
        getline(istrss, temp_s, ',');
        this->DomainSizeX = atof(temp_s.c_str());

        //--------line 3
        istrss.str(items[2]);
        getline(istrss, temp_s, ',');
        getline(istrss, temp_s, ',');
        this->DomainDimensionRatio.x = atof(temp_s.c_str());
        getline(istrss, temp_s, ',');
        this->DomainDimensionRatio.y = atof(temp_s.c_str());
        getline(istrss, temp_s, ',');
        this->DomainDimensionRatio.z = atof(temp_s.c_str());

        //--------line 4
        istrss.str(items[3]);
        getline(istrss, temp_s, ',');
        getline(istrss, temp_s, ',');
        this->PercoDir = atoi(temp_s.c_str());

        //--------line 5
        istrss.str(items[4]);
        getline(istrss, temp_s, ',');
        getline(istrss, temp_s, ',');
        int NumGroups = atoi(temp_s.c_str());
        this->NumFractures.resize(NumGroups);
        this->Kappa.resize(NumGroups);
        this->MeanOrientationOfFisherDistribution.resize(NumGroups);
        this->Beta.resize(NumGroups);
        this->Gamma.resize(NumGroups);
        this->ModeOfSizeDistribution.resize(NumGroups);
        this->SizeDistributionParameters.resize(NumGroups);

        //--------line 6
        istrss.str(items[5]);
        getline(istrss, temp_s, ',');
        for (int i = 0; i < NumGroups; ++i)
        {
            getline(istrss, temp_s, ',');
            this->NumFractures[i] = atoi(temp_s.c_str());
        }

        //--------line 7
        istrss.str(items[6]);
        getline(istrss, temp_s, ',');
        for (int i = 0; i < NumGroups; ++i)
        {
            getline(istrss, temp_s, ',');
            this->Kappa[i] = atof(temp_s.c_str());
        }

        //--------line 8
        istrss.str(items[7]);
        getline(istrss, temp_s, ',');
        for (int i = 0; i < NumGroups; ++i)
        {
            getline(istrss, temp_s, ',');
            this->MeanOrientationOfFisherDistribution[i].x =
                atof(temp_s.c_str());
            getline(istrss, temp_s, ',');
            this->MeanOrientationOfFisherDistribution[i].y =
                atof(temp_s.c_str());
            getline(istrss, temp_s, ',');
            this->MeanOrientationOfFisherDistribution[i].z =
                atof(temp_s.c_str());

            T norm_meanOri =
                pow(this->MeanOrientationOfFisherDistribution[i].x *
                            this->MeanOrientationOfFisherDistribution[i].x +
                        this->MeanOrientationOfFisherDistribution[i].y *
                            this->MeanOrientationOfFisherDistribution[i].y +
                        this->MeanOrientationOfFisherDistribution[i].z *
                            this->MeanOrientationOfFisherDistribution[i].z,
                    0.5);
            this->MeanOrientationOfFisherDistribution[i].x /= norm_meanOri,
                this->MeanOrientationOfFisherDistribution[i].y /= norm_meanOri,
                this->MeanOrientationOfFisherDistribution[i].z /= norm_meanOri;
        }

        //--------line 9
        istrss.str(items[8]);
        getline(istrss, temp_s, ',');
        for (int i = 0; i < NumGroups; ++i)
        {
            getline(istrss, temp_s, ',');
            this->ModeOfSizeDistribution[i] = atoi(temp_s.c_str());
        }

        //--------line 10
        istrss.str(items[9]);
        getline(istrss, temp_s, ',');
        for (int i = 0; i < NumGroups; ++i)
        {
            getline(istrss, temp_s, ',');
            this->SizeDistributionParameters[i].x = atof(temp_s.c_str());
            getline(istrss, temp_s, ',');
            this->SizeDistributionParameters[i].y = atof(temp_s.c_str());
            getline(istrss, temp_s, ',');
            this->SizeDistributionParameters[i].z = atof(temp_s.c_str());
            getline(istrss, temp_s, ',');
            this->SizeDistributionParameters[i].w = atof(temp_s.c_str());
        }

        //--------line 11
        istrss.str(items[10]);
        getline(istrss, temp_s, ',');
        for (int i = 0; i < NumGroups; ++i)
        {
            getline(istrss, temp_s, ',');
            this->Beta[i] = atof(temp_s.c_str());
        }

        //--------line 12
        istrss.str(items[11]);
        getline(istrss, temp_s, ',');
        for (int i = 0; i < NumGroups; ++i)
        {
            getline(istrss, temp_s, ',');
            this->Gamma[i] = atof(temp_s.c_str());
        }

        cout << "Load input parameters from "
             << xlsxNameWithoutSuffix + ".csv for a stochastic DFN\n";
        cout << "Domain size is "
             << this->DomainSizeX * this->DomainDimensionRatio.x << ", "
             << this->DomainSizeX * this->DomainDimensionRatio.y << ", "
             << this->DomainSizeX * this->DomainDimensionRatio.z << "\n";
        cout << "Percolation direction is " << this->PercoDir << endl;

        for (int i = 0; i < this->NumFractures.size(); ++i)
        {
            cout << "For group " << i << ":\n";
            cout << "\t"
                 << "NumFractures: " << NumFractures[i] << endl;
            cout << "\t"
                 << "Kappa: " << Kappa[i] << endl;
            cout << "\t"
                 << "MeanOrientationOfFisherDistribution: "
                 << MeanOrientationOfFisherDistribution[i].x << ", "
                 << MeanOrientationOfFisherDistribution[i].y << ", "
                 << MeanOrientationOfFisherDistribution[i].z << endl;
            cout << "\t"
                 << "Beta: " << Beta[i] << endl;
            cout << "\t"
                 << "Gamma: " << Gamma[i] << endl;
            cout << "\t"
                 << "ModeOfSizeDistribution: " << ModeOfSizeDistribution[i]
                 << endl;
            cout << "\t"
                 << "SizeDistributionParameters: "
                 << SizeDistributionParameters[i].x << ", "
                 << SizeDistributionParameters[i].y << ", "
                 << SizeDistributionParameters[i].z << ", "
                 << SizeDistributionParameters[i].w << endl;
        }
        this->FractureGeneration();
    }
    else
    {
        //
        //--------line 2
        istringstream istrss(items[1]);
        getline(istrss, temp_s, ',');
        getline(istrss, temp_s, ',');
        this->DomainSizeX = atof(temp_s.c_str());

        //--------line 3
        istrss.str(items[2]);
        getline(istrss, temp_s, ',');
        getline(istrss, temp_s, ',');
        this->DomainDimensionRatio.x = atof(temp_s.c_str());
        getline(istrss, temp_s, ',');
        this->DomainDimensionRatio.y = atof(temp_s.c_str());
        getline(istrss, temp_s, ',');
        this->DomainDimensionRatio.z = atof(temp_s.c_str());

        //--------line 4
        istrss.str(items[3]);
        getline(istrss, temp_s, ',');
        getline(istrss, temp_s, ',');
        this->PercoDir = atoi(temp_s.c_str());

        //--------line 5
        istrss.str(items[4]);
        getline(istrss, temp_s, ',');
        getline(istrss, temp_s, ',');
        this->NumFracturesTotal = atoi(temp_s.c_str());

        int DataNumForOneFracture = 12;
        thrust::host_vector<T> Data_deterministic(DataNumForOneFracture *
                                                  this->NumFracturesTotal, 0);

        // -------//--------line > 5
        // cout << "NumFrac: " << this->NumFracturesTotal << endl;
        // ----------------
        std::ifstream ifs_II;
        ifs_II.open(xlsxNameWithoutSuffix + ".csv", ios::in);

        if (!ifs_II.is_open())
        {
            ifs_II.close();
            throw cuDFNsys::ExceptionsPause("Fail to open " +
                                            xlsxNameWithoutSuffix + ".csv");
        }

        items.resize(5 + this->NumFracturesTotal);
        for (int i = 0; i < items.size(); ++i)
        {
            string temp_k;
            getline(ifs_II, temp_k);
            items[i] = temp_k;
            //cout << temp_k << endl;
        }
        ifs_II.close();
        //-----------------

        for (int i = 0; i < NumFracturesTotal; ++i)
        {

            istrss.str(items[5 + i]);
            getline(istrss, temp_s, ',');

            string PSLD = istrss.str();
            //cout << "Fracture " << i << ": " << PSLD << endl;

            int commaCount = std::count(PSLD.begin(), PSLD.end(), ',');

            //cout << commaCount << endl;
            int JKLS = (commaCount >= 12 ? DataNumForOneFracture
                                         : DataNumForOneFracture - 3);

            for (int j = 0; j < JKLS; ++j)
            {
                getline(istrss, temp_s, ',');
                Data_deterministic[i * DataNumForOneFracture + j] =
                    atof(temp_s.c_str());
                //cout << atof(temp_s.c_str()) << (j == JKLS - 1 ? "\n" : ", ");
            }
        }

        this->FracturesHost.resize(this->NumFracturesTotal);
        this->FracturesDevice = this->FracturesHost;
        this->FracturesDevicePtr =
            thrust::raw_pointer_cast(FracturesDevice.data());

        thrust::device_vector<T> Data_deterministic_device = Data_deterministic;
        T *Data_f = thrust::raw_pointer_cast(Data_deterministic_device.data());

        //exit(0);

        cuDFNsys::FracturesDeterministic<T>
            <<<this->NumFracturesTotal / 256 + 1, 256>>>(
                this->FracturesDevicePtr, this->RandomSeed, Data_f,
                this->NumFracturesTotal, this->DomainSizeX, this->PercoDir,
                this->DomainDimensionRatio);
        cudaDeviceSynchronize();

        this->FracturesHost = this->FracturesDevice;
        this->FracturesDevicePtr =
            thrust::raw_pointer_cast(FracturesDevice.data());

        cout << "Load input parameters from "
             << xlsxNameWithoutSuffix + ".csv for a deterministic DFN with "
             << this->NumFracturesTotal << " fractures\n";

        cout << "First fracture: ";
        for (int j = 0; j < DataNumForOneFracture; ++j)
            cout << Data_deterministic[j] << (j == DataNumForOneFracture - 1 ? "\n" : ", ");
        cout << "...\nLast fracture: ";
        for (int j = 0; j < DataNumForOneFracture; ++j)
            cout << Data_deterministic[Data_deterministic.size() -
                                       DataNumForOneFracture + j]
                 << (j == DataNumForOneFracture - 1 ? "\n" : ", ");
    }
};
template void
cuDFNsys::DFN<double>::LoadDFNFromCSV(const string &xlsxNameWithoutSuffix);
template void
cuDFNsys::DFN<float>::LoadDFNFromCSV(const string &xlsxNameWithoutSuffix);

// ====================================================
// NAME:        cuDFNsys::MeshDFN<T>::MeshGeneration
// DESCRIPTION: MeshGeneration
// AUTHOR:      Tingchang YIN
// DATE:        20/10/2023
// ====================================================
template <typename T>
void cuDFNsys::MeshDFN<T>::MeshGeneration(DFN<T> &my_dfn)
{
    if (my_dfn.PercolationCluster.size() == 0)
    {
        cout << ("No percolation cluster exists! The MeshDFN.MeshGeneration is "
                 "discontinued")
             << endl;
        return;
    }

    cuDFNsys::GetAllPercolatingFractures GetPer{
        my_dfn.PercolationCluster, my_dfn.ListClusters, this->FracsPercol};
    std::vector<pair<int, int>> IntersectionPair_percol;
    cuDFNsys::RemoveDeadEndFrac<T> RDEF{
        this->FracsPercol,    IntersectionPair_percol, (size_t)my_dfn.PercoDir,
        my_dfn.FracturesHost, my_dfn.IntersectionMap,  false};

    my_dfn.NumFracturesTotal = my_dfn.FracturesHost.size();
    my_dfn.FracturesDevice = my_dfn.FracturesHost;
    my_dfn.FracturesDevicePtr =
        thrust::raw_pointer_cast(my_dfn.FracturesDevice.data());

    this->MeshData = cuDFNsys::Mesh<T>{
        my_dfn.FracturesHost, IntersectionPair_percol,    &(this->FracsPercol),
        this->MinElementSize, this->MaxElementSize,       my_dfn.PercoDir,
        my_dfn.DomainSizeX,   my_dfn.DomainDimensionRatio};
    this->MeanGridSize = this->MeshData.MeanGridSize;
    // cout << "mesh.Element3D.size(): " << this->MeshData.Element3D.size() << endl;
};
template void cuDFNsys::MeshDFN<double>::MeshGeneration(DFN<double> &my_dfn);
template void cuDFNsys::MeshDFN<float>::MeshGeneration(DFN<float> &my_dfn);

// ====================================================
// NAME:        cuDFNsys::MeshDFN<T>::Visualization
// DESCRIPTION: Visualization of mesh
// AUTHOR:      Tingchang YIN
// DATE:        20/10/2023
// ====================================================
template <typename T>
void cuDFNsys::MeshDFN<T>::Visualization(DFN<T> my_dfn,
                                         const string &MatlabScriptName,
                                         const string &PythonScriptName,
                                         const string &HDF5FileName,
                                         const bool &IfCheck2DCoordinatesOfMesh,
                                         const bool &IfCheckEdgeAttributes)
{
    if (my_dfn.PercolationCluster.size() == 0)
    {
        cout << ("No percolation cluster exists! The MeshDFN.Visualization is "
                 "discontinued")
             << endl;
        return;
    }

    this->MeanGridSize = this->MeshData.MatlabPlot(
        HDF5FileName + ".h5", MatlabScriptName + ".m", my_dfn.FracturesHost,
        my_dfn.DomainSizeX, IfCheck2DCoordinatesOfMesh, IfCheckEdgeAttributes,
        true, PythonScriptName, my_dfn.DomainDimensionRatio);
}; //
template void cuDFNsys::MeshDFN<double>::Visualization(
    DFN<double> my_dfn, const string &MatlabScriptName,
    const string &PythonScriptName, const string &HDF5FileName,
    const bool &IfCheck2DCoordinatesOfMesh, const bool &IfCheckEdgeAttributes);
template void cuDFNsys::MeshDFN<float>::Visualization(
    DFN<float> my_dfn, const string &MatlabScriptName,
    const string &PythonScriptName, const string &HDF5FileName,
    const bool &IfCheck2DCoordinatesOfMesh, const bool &IfCheckEdgeAttributes);

// ====================================================
// NAME:        cuDFNsys::MeshDFN<T>::StoreInH5
// DESCRIPTION: store this class in .h5
// AUTHOR:      Tingchang YIN
// DATE:        23/10/2023
// ====================================================
template <typename T>
void cuDFNsys::MeshDFN<T>::StoreInH5(const string &ClassNameH5)
{
    cuDFNsys::OutputObjectData<T> lk_out;
    lk_out.OutputMesh(ClassNameH5 + ".h5", this->MeshData, this->FracsPercol);
}; // cuDFNsys::MeshDFN<T>::StoreInH5
template void cuDFNsys::MeshDFN<double>::StoreInH5(const string &ClassNameH5);
template void cuDFNsys::MeshDFN<float>::StoreInH5(const string &ClassNameH5);

// ====================================================
// NAME:        cuDFNsys::MeshDFN<T>::LoadClassFromH5
// DESCRIPTION: Load class data from a h5 file
// AUTHOR:      Tingchang YIN
// DATE:        23/10/2023
// ====================================================
template <typename T>
void cuDFNsys::MeshDFN<T>::LoadClassFromH5(const string &ClassNameH5)
{
    cuDFNsys::HDF5API hdf5Class;
    std::vector<uint> Fracs_percol_II =
        hdf5Class.ReadDataset<uint>(ClassNameH5 + ".h5", "N", "Fracs_percol");
    this->FracsPercol.resize(Fracs_percol_II.size());
    std::copy(Fracs_percol_II.begin(), Fracs_percol_II.end(),
              this->FracsPercol.data());
    cuDFNsys::InputObjectData<T> lk;
    lk.InputMesh(ClassNameH5 + ".h5", this->MeshData, &this->FracsPercol);

    this->MeanGridSize = this->MeshData.MeanGridSize;
}; // cuDFNsys::MeshDFN<T>::LoadClassFromH5
template void
cuDFNsys::MeshDFN<double>::LoadClassFromH5(const string &ClassNameH5);
template void
cuDFNsys::MeshDFN<float>::LoadClassFromH5(const string &ClassNameH5);

// ====================================================
// NAME:        cuDFNsys::FlowDFN<T>::FlowSimulation
// DESCRIPTION: Flow simulation
// AUTHOR:      Tingchang YIN
// DATE:        20/10/2023
// ====================================================
template <typename T>
void cuDFNsys::FlowDFN<T>::FlowSimulation(cuDFNsys::DFN<T> my_dfn,
                                          cuDFNsys::MeshDFN<T> my_mesh)
{
    if (my_dfn.PercolationCluster.size() == 0)
    {
        cout << ("No percolation cluster exists! The FlowDFN.FlowSimulation is "
                 "discontinued")
             << endl;
        return;
    }

    // this->IfPeriodic = my_dfn.CheckIfPeriodic();
    // cout << "my_dfn.FracturesHost.size = " << my_dfn.FracturesHost.size()
    //      << endl;
    // cout << "my_dfn.PercoDir: " << my_dfn.PercoDir << endl;
    // cout << "this->InletHead: " << this->InletHead << endl;
    // cout << "this->OutletHead: " << this->OutletHead << endl;
    // cout << "my_dfn.DomainDimensionRatio: " << my_dfn.DomainDimensionRatio.x
    //      << ", " << my_dfn.DomainDimensionRatio.y << ", "
    //      << my_dfn.DomainDimensionRatio.z << endl;
    // cout << "1 this->ConsTq" << this->ConsTq << endl;
    this->FlowData = cuDFNsys::MHFEM<T>{my_mesh.MeshData,
                                        my_dfn.FracturesHost,
                                        this->InletHead,
                                        this->OutletHead,
                                        my_dfn.PercoDir,
                                        my_dfn.DomainSizeX,
                                        my_dfn.DomainDimensionRatio,
                                        false,
                                        0,
                                        this->IfPeriodic,
                                        this->MuOverRhoG,
                                        this->ConsTq};
};
template void
cuDFNsys::FlowDFN<double>::FlowSimulation(cuDFNsys::DFN<double> my_dfn,
                                          cuDFNsys::MeshDFN<double> my_mesh);
template void
cuDFNsys::FlowDFN<float>::FlowSimulation(cuDFNsys::DFN<float> my_dfn,
                                         cuDFNsys::MeshDFN<float> my_mesh);

// ====================================================
// NAME:        cuDFNsys::FlowDFN<T>::Visualization
// DESCRIPTION: Flow Visualization
// AUTHOR:      Tingchang YIN
// DATE:        20/10/2023
// ====================================================
template <typename T>
void cuDFNsys::FlowDFN<T>::Visualization(cuDFNsys::DFN<T> my_dfn,
                                         cuDFNsys::MeshDFN<T> my_mesh,
                                         const string &MatlabScriptName,
                                         const string &PythonScriptName,
                                         const string &HDF5FileName)
{
    if (my_dfn.PercolationCluster.size() == 0)
    {
        cout << ("No percolation cluster exists! The FlowDFN.Visualization is "
                 "discontinued")
             << endl;
        return;
    }

    double2 velocities = this->FlowData.MatlabPlot(
        HDF5FileName + ".h5", MatlabScriptName + ".m", my_dfn.FracturesHost,
        my_mesh.MeshData, my_dfn.DomainSizeX, true, PythonScriptName,
        my_dfn.DomainDimensionRatio);
    this->MeanVelocity = velocities.x;
    this->MaxVelocity = velocities.y;
};
template void cuDFNsys::FlowDFN<double>::Visualization(
    cuDFNsys::DFN<double> my_dfn, cuDFNsys::MeshDFN<double> my_mesh,
    const string &MatlabScriptName, const string &PythonScriptName,
    const string &HDF5FileName);
template void cuDFNsys::FlowDFN<float>::Visualization(
    cuDFNsys::DFN<float> my_dfn, cuDFNsys::MeshDFN<float> my_mesh,
    const string &MatlabScriptName, const string &PythonScriptName,
    const string &HDF5FileName);

// ====================================================
// NAME:        cuDFNsys::FlowDFN<T>::StoreInH5
// DESCRIPTION: store this class in .h5
// AUTHOR:      Tingchang YIN
// DATE:        23/10/2023
// ====================================================
template <typename T>
void cuDFNsys::FlowDFN<T>::StoreInH5(const string &ClassNameH5)
{
    cuDFNsys::OutputObjectData<T> lk_out;
    lk_out.OutputMHFEM(ClassNameH5 + ".h5", this->FlowData);

    this->MaxVelocity = this->FlowData.MaxVelocity;
    this->MeanVelocity = this->FlowData.MeanVelocity;
    this->ConsTq = this->FlowData.ConsTq;
}; // cuDFNsys::FlowDFN<T>::StoreInH5
template void cuDFNsys::FlowDFN<double>::StoreInH5(const string &ClassNameH5);
template void cuDFNsys::FlowDFN<float>::StoreInH5(const string &ClassNameH5);

// ====================================================
// NAME:        cuDFNsys::FlowDFN<T>::LoadClassFromH5
// DESCRIPTION: Load class data from a h5 file
// AUTHOR:      Tingchang YIN
// DATE:        23/10/2023
// ====================================================
template <typename T>
void cuDFNsys::FlowDFN<T>::LoadClassFromH5(const string &ClassNameH5)
{
    cuDFNsys::InputObjectData<T> lk;
    lk.InputMHFEM(ClassNameH5 + ".h5", this->FlowData);
    this->InletHead = this->FlowData.InletP;
    this->OutletHead = this->FlowData.OutletP;
    this->IfPeriodic = this->FlowData.IfPeriodic;
    this->MuOverRhoG = this->FlowData.MuOverRhoG;
}; // cuDFNsys::FlowDFN<T>::LoadClassFromH5
template void
cuDFNsys::FlowDFN<double>::LoadClassFromH5(const string &ClassNameH5);
template void
cuDFNsys::FlowDFN<float>::LoadClassFromH5(const string &ClassNameH5);

// ====================================================
// NAME:        cuDFNsys::PTDFN<T>::ParticleTracking
// DESCRIPTION: ParticleTracking
// AUTHOR:      Tingchang YIN
// DATE:        21/10/2023
// ====================================================
template <typename T>
void cuDFNsys::PTDFN<T>::ParticleTracking(cuDFNsys::DFN<T> my_dfn,
                                          cuDFNsys::MeshDFN<T> my_mesh,
                                          cuDFNsys::FlowDFN<T> my_flow)
{
    if (my_dfn.PercolationCluster.size() == 0)
    {
        cout << ("No percolation cluster exists! The PTDFN.ParticleTracking is "
                 "discontinued")
             << endl;
        return;
    }

    cuDFNsys::OutputObjectData<T> lk;
    lk.OutputFractures("FracturesForParticle.h5", my_dfn.FracturesHost,
                       my_dfn.DomainSizeX);

    this->PTData = cuDFNsys::ParticleTransport<T>{
        this->NumTimeSteps,
        my_dfn.FracturesHost,
        my_mesh.MeshData,
        my_flow.FlowData,
        (uint)my_dfn.PercoDir,
        (T)-0.5 * (T)(&my_dfn.DomainDimensionRatio.x)[my_dfn.PercoDir] *
            my_dfn.DomainSizeX,
        this->NumParticles,
        this->DeltaT,
        this->MolecularDiffusion,
        "Particle_tracking",
        InjectionMethod,
        (this->OutputAllPTInformationOrFPTCurve ? "OutputAll" : "FPTCurve"),
        false,
        0,
        false,
        this->SpacingOfControlPlanes,
        this->IfOutputVarianceOfDisplacementsEachStep,
        this->IfInjectAtCustomedPlane,
        this->CustomedPlaneInjection,
        this->IfUseFluxWeightedOrEqualProbableMixingIntersection,
        this->IfPeriodic};
}; // ParticleTracking
template void
cuDFNsys::PTDFN<double>::ParticleTracking(cuDFNsys::DFN<double> my_dfn,
                                          cuDFNsys::MeshDFN<double> my_mesh,
                                          cuDFNsys::FlowDFN<double> my_flow);
template void
cuDFNsys::PTDFN<float>::ParticleTracking(cuDFNsys::DFN<float> my_dfn,
                                         cuDFNsys::MeshDFN<float> my_mesh,
                                         cuDFNsys::FlowDFN<float> my_flow);

// ====================================================
// NAME:        cuDFNsys::PTDFN<T>::ParticleTracking
// DESCRIPTION: Visualization of ParticleTracking
// AUTHOR:      Tingchang YIN
// DATE:        21/10/2023
// ====================================================
template <typename T>
void cuDFNsys::PTDFN<T>::Visualization(cuDFNsys::DFN<T> my_dfn,
                                       cuDFNsys::MeshDFN<T> my_mesh,
                                       cuDFNsys::FlowDFN<T> my_flow,
                                       const string &MatlabScriptName,
                                       const string &PythonScriptName,
                                       const string &HDF5FileNameOfFlowDFN)
{
    if (my_dfn.PercolationCluster.size() == 0)
    {
        cout << ("No percolation cluster exists! The PTDFN.Visualization is "
                 "discontinued")
             << endl;
        return;
    }

    this->PTData.MatlabPlot(
        HDF5FileNameOfFlowDFN + ".h5", MatlabScriptName + ".m",
        my_mesh.MeshData, my_flow.FlowData, my_dfn.DomainSizeX,
        my_dfn.DomainDimensionRatio, true, PythonScriptName);
};
template void cuDFNsys::PTDFN<double>::Visualization(
    cuDFNsys::DFN<double> my_dfn, cuDFNsys::MeshDFN<double> my_mesh,
    cuDFNsys::FlowDFN<double> my_flow, const string &MatlabScriptName,
    const string &PythonScriptName, const string &HDF5FileNameOfFlowDFN);
template void cuDFNsys::PTDFN<float>::Visualization(
    cuDFNsys::DFN<float> my_dfn, cuDFNsys::MeshDFN<float> my_mesh,
    cuDFNsys::FlowDFN<float> my_flow, const string &MatlabScriptName,
    const string &PythonScriptName, const string &HDF5FileNameOfFlowDFN);

// ====================================================
// NAME:        cuDFNsys::PTDFN<T>::LoadParametersFromCSV
// DESCRIPTION: LoadParametersFromCSV
// AUTHOR:      Tingchang YIN
// DATE:        03/11/2023
// ====================================================
template <typename T>
void cuDFNsys::PTDFN<T>::LoadParametersFromCSV(const string &CSVName)
{
    std::ifstream ifs;
    ifs.open(CSVName + ".csv", ios::in);

    if (!ifs.is_open())
    {
        ifs.close();
        throw cuDFNsys::ExceptionsPause("Fail to open " + CSVName + ".csv");
    }

    vector<string> items(11);
    for (int i = 0; i < items.size(); ++i)
    {
        string temp_k;
        getline(ifs, temp_k);
        items[i] = temp_k;
        //cout << temp_k << endl;
    }
    ifs.close();

    string temp_s;
    //--------line 1
    istringstream istrss(items[0]);
    getline(istrss, temp_s, ',');
    getline(istrss, temp_s, ',');
    this->NumParticles = atoi(temp_s.c_str());

    //--------line 2
    istrss.str(items[1]);
    getline(istrss, temp_s, ',');
    getline(istrss, temp_s, ',');
    this->NumTimeSteps = atoi(temp_s.c_str());

    //--------line 3
    istrss.str(items[2]);
    getline(istrss, temp_s, ',');
    getline(istrss, temp_s, ',');
    this->MolecularDiffusion = atof(temp_s.c_str());

    //--------line 4
    istrss.str(items[3]);
    getline(istrss, temp_s, ',');
    getline(istrss, temp_s, ',');
    this->DeltaT = atof(temp_s.c_str());

    //--------line 5
    istrss.str(items[4]);
    getline(istrss, temp_s, ',');
    getline(istrss, temp_s, ',');
    this->InjectionMethod = temp_s;

    //--------line 6
    istrss.str(items[5]);
    getline(istrss, temp_s, ',');
    getline(istrss, temp_s, ',');
    this->OutputAllPTInformationOrFPTCurve =
        (atoi(temp_s.c_str()) == 1 ? true : false);

    //--------line 7
    istrss.str(items[6]);
    getline(istrss, temp_s, ',');
    getline(istrss, temp_s, ',');
    this->SpacingOfControlPlanes = atof(temp_s.c_str());

    //--------line 8
    istrss.str(items[7]);
    getline(istrss, temp_s, ',');
    getline(istrss, temp_s, ',');

    this->IfOutputVarianceOfDisplacementsEachStep =
        (atoi(temp_s.c_str()) == 1 ? true : false);

    //--------line 9
    istrss.str(items[8]);
    getline(istrss, temp_s, ',');
    getline(istrss, temp_s, ',');
    this->IfInjectAtCustomedPlane = (atoi(temp_s.c_str()) == 1 ? true : false);

    //--------line 10
    istrss.str(items[9]);
    getline(istrss, temp_s, ',');
    getline(istrss, temp_s, ',');
    this->CustomedPlaneInjection = atof(temp_s.c_str());

    //--------line 11
    istrss.str(items[10]);
    getline(istrss, temp_s, ',');
    getline(istrss, temp_s, ',');
    this->IfUseFluxWeightedOrEqualProbableMixingIntersection =
        (atoi(temp_s.c_str()) == 1 ? true : false);

    cout << "Load PT parameters from " << CSVName << ".csv:\n";
    cout << "NumParticles: " << (NumParticles) << endl;
    cout << "NumTimeSteps: " << (NumTimeSteps) << endl;
    cout << "MolecularDiffusion: " << (MolecularDiffusion) << endl;
    cout << "DeltaT: " << (DeltaT) << endl;
    cout << "InjectionMethod: " << (InjectionMethod) << endl;
    cout << "OutputAllPTInformationOrFPTCurve: "
         << (OutputAllPTInformationOrFPTCurve ? "OutputAll" : "FPTCurve")
         << endl;
    cout << "SpacingOfControlPlanes: " << (SpacingOfControlPlanes) << endl;
    cout << "IfOutputVarianceOfDisplacementsEachStep: "
         << (IfOutputVarianceOfDisplacementsEachStep ? "yes" : "no") << endl;
    cout << "IfInjectAtCustomedPlane: "
         << (IfInjectAtCustomedPlane ? "yes" : "no") << endl;
    cout << "CustomedPlaneInjection: " << (CustomedPlaneInjection) << endl;
    cout << "IfUseFluxWeightedOrEqualProbableMixingIntersection: "
         << (IfUseFluxWeightedOrEqualProbableMixingIntersection
                 ? "outgoingfluxweighted"
                 : "equalproble")
         << endl;
};
template void
cuDFNsys::PTDFN<double>::LoadParametersFromCSV(const string &CSVName);
template void
cuDFNsys::PTDFN<float>::LoadParametersFromCSV(const string &CSVName);