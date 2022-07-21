#include "Fractures/IdentifyIntersection.cuh"

// ====================================================
// NAME:        IdentifyIntersection
// DESCRIPTION: Identify Intersections of fractures
//              in a  DFN with CPU
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
template <typename T>
cuDFNsys::IdentifyIntersection<T>::IdentifyIntersection(thrust::host_vector<cuDFNsys::Fracture<T>> verts,
                                                        const bool &if_trucncated,
                                                        std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> &Intersection_map)
{
    Intersection_map.clear();

    for (int i = 1; i < verts.size(); ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            cuDFNsys::Vector3<T> dist_two_frac = cuDFNsys::MakeVector3<T>(verts[i].Center.x - verts[j].Center.x,
                                                                          verts[i].Center.y - verts[j].Center.y,
                                                                          verts[i].Center.z - verts[j].Center.z);

            T ddis = pow(dist_two_frac.x * dist_two_frac.x +
                             dist_two_frac.y * dist_two_frac.y +
                             dist_two_frac.z * dist_two_frac.z,
                         0.5);

            if (ddis > (verts[i].Radius + verts[j].Radius))
                continue;
            //---------------

            T Mat_2_to_3[3][3];
            verts[i].RoationMatrix(Mat_2_to_3, 23); //= verts[i].Roation_Matrix_2Dto3D;
            T Mat_3_to_2[3][3];
            verts[i].RoationMatrix(Mat_3_to_2, 32);

            int NUM_verts = 0;

            if (if_trucncated == false)
                NUM_verts = 4;
            else
                NUM_verts = verts[i].NumVertsTruncated;

            cuDFNsys::Vector3<T> *Frac_verts = (cuDFNsys::Vector3<T> *)malloc(NUM_verts * sizeof(cuDFNsys::Vector3<T>));
            // = verts[i].Verts3D;

            if (if_trucncated == false)
                memcpy(Frac_verts, verts[i].Verts3D, sizeof(cuDFNsys::Vector3<T>) * NUM_verts);
            else
                memcpy(Frac_verts, verts[i].Verts3DTruncated, sizeof(cuDFNsys::Vector3<T>) * NUM_verts);

            cuDFNsys::Vector3<T> Center_ = verts[i].Center;
            cuDFNsys::Vector2<T> *Frac_verts_2D = (cuDFNsys::Vector2<T> *)malloc(NUM_verts * sizeof(cuDFNsys::Vector2<T>));

            for (int k = 0; k < NUM_verts; ++k)
            {
                Frac_verts[k] = cuDFNsys::MakeVector3<T>(Frac_verts[k].x - Center_.x,
                                                         Frac_verts[k].y - Center_.y,
                                                         Frac_verts[k].z - Center_.z);
                cuDFNsys::Vector3<T> HJ = cuDFNsys::ProductSquare3Vector3<T>(Mat_3_to_2, Frac_verts[k]);
                Frac_verts_2D[k] = cuDFNsys::MakeVector2<T>(HJ.x, HJ.y);
            }

            int NUM_verts_j = 0;
            if (if_trucncated == false)
                NUM_verts_j = 4;
            else
                NUM_verts_j = verts[j].NumVertsTruncated;

            cuDFNsys::Vector3<T> *Frac_verts_j = (cuDFNsys::Vector3<T> *)malloc(sizeof(cuDFNsys::Vector3<T>) * NUM_verts_j); // = verts[j].verts_3D;

            if (if_trucncated == false)
                memcpy(Frac_verts_j, verts[j].Verts3D, sizeof(cuDFNsys::Vector3<T>) * NUM_verts_j);
            else
                memcpy(Frac_verts_j, verts[j].Verts3DTruncated, sizeof(cuDFNsys::Vector3<T>) * NUM_verts_j);

            for (int k = 0; k < NUM_verts_j; ++k)
            {
                Frac_verts_j[k] = cuDFNsys::MakeVector3<T>(Frac_verts_j[k].x - Center_.x,
                                                           Frac_verts_j[k].y - Center_.y,
                                                           Frac_verts_j[k].z - Center_.z);
                Frac_verts_j[k] = cuDFNsys::ProductSquare3Vector3<T>(Mat_3_to_2,
                                                                     Frac_verts_j[k]);
            }

            //---------------start to identify
            cuDFNsys::Vector3<T> Intersection_xyPlane_Poly2[2];
            bool ik = cuDFNsys::Intersection3DPolyXYPlane<T>(Frac_verts_j,
                                                             NUM_verts_j,
                                                             Intersection_xyPlane_Poly2,
                                                             _TOL_Intersection3DPolyXYPlane);
            if (ik == false)
            {
                free(Frac_verts);
                free(Frac_verts_j);
                free(Frac_verts_2D);
                continue;
            }
            cuDFNsys::Vector2<T> LineInf[2];
            LineInf[0] = cuDFNsys::MakeVector2<T>(Intersection_xyPlane_Poly2[0].x, Intersection_xyPlane_Poly2[0].y);
            LineInf[1] = cuDFNsys::MakeVector2<T>(Intersection_xyPlane_Poly2[1].x, Intersection_xyPlane_Poly2[1].y);

            cuDFNsys::Vector2<T> Intersection_Poly1_InfLine[2];
            bool jt = cuDFNsys::Intersection2DLine2DPoly<T>(Frac_verts_2D,
                                                            NUM_verts,
                                                            LineInf,
                                                            Intersection_Poly1_InfLine,
                                                            _TOL_Intersection2DLine2DPoly);
            if (jt == false)
            {
                free(Frac_verts);
                free(Frac_verts_j);
                free(Frac_verts_2D);
                continue;
            }
            int sign_ = -1;
            cuDFNsys::Vector2<T> Intersection_two_segs[2];
            bool kh = cuDFNsys::IntersectionTwoCollinear2DSegs<T>(LineInf,
                                                                  Intersection_Poly1_InfLine,
                                                                  Intersection_two_segs,
                                                                  &sign_,
                                                                  _TOL_IntersectionTwoCollinearSegs);

            if (kh == true)
            {
                if (sign_ == 2)
                {
                    cuDFNsys::Vector3<T> Intersection_f[2];
                    Intersection_f[0] = cuDFNsys::MakeVector3<T>(Intersection_two_segs[0].x,
                                                                 Intersection_two_segs[0].y, 0);
                    Intersection_f[1] = cuDFNsys::MakeVector3<T>(Intersection_two_segs[1].x,
                                                                 Intersection_two_segs[1].y, 0);

                    for (int k = 0; k < 2; ++k)
                    {
                        Intersection_f[k] = cuDFNsys::ProductSquare3Vector3<T>(Mat_2_to_3,
                                                                               Intersection_f[k]);
                        Intersection_f[k] = cuDFNsys::MakeVector3<T>(Intersection_f[k].x + Center_.x, Intersection_f[k].y + Center_.y, Intersection_f[k].z + Center_.z);
                    }

                    pair<size_t, size_t> key_ = std::make_pair((size_t)i, (size_t)j);
                    pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>> p;
                    p.first = Intersection_f[0];
                    p.second = Intersection_f[1];
                    pair<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> element_ =
                        std::make_pair(key_, p);
                    Intersection_map.insert(element_);

                    free(Frac_verts);
                    free(Frac_verts_j);
                    free(Frac_verts_2D);
                    continue;
                }
                else
                {

                    free(Frac_verts);
                    free(Frac_verts_j);
                    free(Frac_verts_2D);
                    continue;
                }
            }
            else
            {

                free(Frac_verts);
                free(Frac_verts_j);
                free(Frac_verts_2D);
                continue;
            }
        }
    }
}; // IdentifyIntersection
template cuDFNsys::IdentifyIntersection<double>::IdentifyIntersection(thrust::host_vector<cuDFNsys::Fracture<double>> verts,
                                                                      const bool &if_trucncated,
                                                                      std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<double>, cuDFNsys::Vector3<double>>> &Intersection_map);
template cuDFNsys::IdentifyIntersection<float>::IdentifyIntersection(thrust::host_vector<cuDFNsys::Fracture<float>> verts,
                                                                     const bool &if_trucncated,
                                                                     std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<float>, cuDFNsys::Vector3<float>>> &Intersection_map);

// ====================================================
// NAME:        IdentifyIntersection
// DESCRIPTION: Identify Intersections of fractures
//              in a  DFN with GPU
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
template <typename T>
cuDFNsys::IdentifyIntersection<T>::IdentifyIntersection(const size_t &Fracsize,
                                                        cuDFNsys::Fracture<T> *Frac_verts_device_ptr,
                                                        const bool &if_trucncated,
                                                        std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> &Intersection_map)
{
    Intersection_map.clear();

    int DSIZE = (int)Fracsize;
    int NUM_frac_pairs = DSIZE * floor((DSIZE - 1) / 2) + (DSIZE - 1) % 2 * DSIZE * 0.5;

    thrust::host_vector<int3> FracPairHost(NUM_frac_pairs);
    int HowManyPairesOneTime = 256 * 6;

    for (int i = 0; i <= NUM_frac_pairs / HowManyPairesOneTime; ++i)
    {
        int First_index = i * HowManyPairesOneTime;
        int Last_index = First_index + ((First_index + HowManyPairesOneTime) < NUM_frac_pairs ? HowManyPairesOneTime : (NUM_frac_pairs - First_index));

        if (Last_index == First_index)
            break;

        int Span_ = Last_index - First_index;
        //cout << "First_index: " << First_index << ", Last_index: " << Last_index << endl;
        thrust::device_vector<int3> FracPairDevice(Span_);
        int3 *ptrFracPairDevice = thrust::raw_pointer_cast(FracPairDevice.data());
        //cout << "NUM_frac_pairs: " << NUM_frac_pairs << endl;
        //cout << "\tSpherical test started\n";
        cuDFNsys::IdentifyFracPairSphericalDetection<T><<<Span_ / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                                  ptrFracPairDevice,
                                                                                  First_index,
                                                                                  Span_);
        cudaDeviceSynchronize();
        //cout << "\tSpherical test finished\n";
        //FracPairHost = FracPairDevice;
        //FracPairDevice.resize(0);
        //FracPairDevice.shrink_to_fit();
        thrust::copy(FracPairDevice.begin(), FracPairDevice.end(), FracPairHost.begin() + First_index);
    };
    thrust::host_vector<cuDFNsys::Intersection<T>> IntersectionHost;
    IntersectionHost.reserve(NUM_frac_pairs);

    for (size_t i = 0; i < NUM_frac_pairs; ++i)
    {
        if (FracPairHost[i].z != 0)
        {
            cuDFNsys::Intersection<T> TMP;
            TMP.FracIDPair.x = FracPairHost[i].x;
            TMP.FracIDPair.y = FracPairHost[i].y;
            IntersectionHost.push_back(TMP);
            //cout << TMP.FracIDPair.x << ", " << TMP.FracIDPair.y << endl;
        }
    }
    //cout << "IntersectionHost.size = " << IntersectionHost.size() << endl;
    FracPairHost.clear();
    FracPairHost.shrink_to_fit();
    IntersectionHost.shrink_to_fit();
    NUM_frac_pairs = IntersectionHost.size();

    HowManyPairesOneTime = 256 * 3;

    //cout << "NUM_pairs = " << NUM_frac_pairs << endl;
    for (int i = 0; i <= NUM_frac_pairs / HowManyPairesOneTime; ++i)
    {

        int First_index = i * HowManyPairesOneTime;
        int Last_index = First_index + ((First_index + HowManyPairesOneTime) < NUM_frac_pairs ? HowManyPairesOneTime : (NUM_frac_pairs - First_index));

        if (Last_index == First_index)
            break;
        //cout << "First_index: " << First_index << ", Last_index: " << Last_index << endl;
        thrust::device_vector<cuDFNsys::Intersection<T>> IntersectionDevice{IntersectionHost.begin() + First_index, IntersectionHost.begin() + Last_index};
        cuDFNsys::Intersection<T> *ptrIntersection = thrust::raw_pointer_cast(IntersectionDevice.data());

        //cout << "\tIdentifyIntersectionKernel started\n";
        //cout << "Number of cores: " << IntersectionDevice.size() << endl;
        cuDFNsys::IdentifyIntersectionKernel<T><<<IntersectionDevice.size() / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                                              IntersectionDevice.size(),
                                                                                              ptrIntersection,
                                                                                              if_trucncated);
        cudaDeviceSynchronize();

        //cout << "\tIdentifyIntersectionKernel finished\n";
        thrust::copy(IntersectionDevice.begin(), IntersectionDevice.end(), IntersectionHost.begin() + First_index);
    };

    for (size_t i = 0; i < IntersectionHost.size(); ++i)
    {
        if (IntersectionHost[i].FracIDPair.x != -1)
        {
            pair<size_t, size_t> key_ = std::make_pair((size_t)IntersectionHost[i].FracIDPair.x,
                                                       (size_t)IntersectionHost[i].FracIDPair.y);
            pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>> p;
            p.first = IntersectionHost[i].Coord[0];
            p.second = IntersectionHost[i].Coord[1];
            pair<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> element_ =
                std::make_pair(key_, p);
            Intersection_map.insert(element_);
        }
    };
    //cout << "IdentifyIntersection GPU finished\n";
}; // IdentifyIntersection
template cuDFNsys::IdentifyIntersection<double>::IdentifyIntersection<double>(const size_t &Fracsize,
                                                                              cuDFNsys::Fracture<double> *Frac_verts_device_ptr,
                                                                              const bool &if_trucncated,
                                                                              std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<double>, cuDFNsys::Vector3<double>>> &Intersection_map);
template cuDFNsys::IdentifyIntersection<float>::IdentifyIntersection<float>(const size_t &Fracsize,
                                                                            cuDFNsys::Fracture<float> *Frac_verts_device_ptr,
                                                                            const bool &if_trucncated,
                                                                            std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<float>, cuDFNsys::Vector3<float>>> &Intersection_map);