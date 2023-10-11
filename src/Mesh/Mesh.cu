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

#include "Mesh/Mesh.cuh"

// ====================================================
// NAME:        Mesh
// DESCRIPTION: Constructor of Mesh
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
cuDFNsys::Mesh<T>::Mesh(const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs,
                        const std::vector<pair<int, int>> &IntersectionPair_percol,
                        std::vector<size_t> *Fracs_percol,
                        const T &min_ele_edge,
                        const T &max_ele_edge,
                        const int &dir_,
                        const T &L,
                        double3 DomainDimensionRatio)
{
    this->Dir = dir_;
    this->FracID = Fracs_percol;
    size_t NUM_frac = (*Fracs_percol).size();
    if (NUM_frac < 1)
    {
        MeshSuccess = false;
        return;
    }

    try
    {
        gmsh::initialize();
        cout << "\tmesh init" << endl;
        //gmsh::option::setNumber("General.NumThreads", Nproc);
        gmsh::option::setNumber("General.Verbosity", 2); // default level is 5
        gmsh::model::add("t2");
        gmsh::model::mesh::clear();
        map<size_t, size_t> Corre_Tag;

        std::vector<size_t>::iterator it_fracID = Fracs_percol->begin();

        thrust::host_vector<cuDFNsys::Fracture<T>> Fracs_ss(NUM_frac);

        double istart = cuDFNsys::CPUSecond();

        for (size_t i = 0; i < NUM_frac; ++i)
        {
            size_t FracID = (*it_fracID);

            it_fracID++;

            int NUM_Points = Fracs[FracID].NumVertsTruncated;

            std::vector<int> Pointloop(NUM_Points);

            for (int k = 0; k < NUM_Points; ++k)
                Pointloop[k] = gmsh::model::occ::addPoint(Fracs[FracID].Verts3DTruncated[k].x,
                                                          Fracs[FracID].Verts3DTruncated[k].y,
                                                          Fracs[FracID].Verts3DTruncated[k].z,
                                                          0);
            std::vector<int> curveloop(NUM_Points);

            for (int k = 0; k < NUM_Points; ++k)
                curveloop[k] = gmsh::model::occ::addLine(Pointloop[k],
                                                         Pointloop[(k + 1) % NUM_Points]);

            int CurveLoopID = gmsh::model::occ::addCurveLoop(curveloop);

            //int SurfaceID = gmsh::model::occ::addPlaneSurface({CurveLoopID});
            int SurfaceID = gmsh::model::occ::addSurfaceFilling(CurveLoopID);
            // cout << "SurfaceID: " << SurfaceID << endl;
            Corre_Tag.insert(std::make_pair(FracID, SurfaceID));
            Fracs_ss[i] = Fracs[FracID];
        }
        gmsh::model::occ::synchronize();
        cout << "\t\tadded entities, running time: " << cuDFNsys::CPUSecond() - istart << " sec\n";

        //cout << "mesh fragment" << endl;
        size_t pair_size = IntersectionPair_percol.size();
        std::vector<std::pair<int, int>> object_entity(pair_size), tool_entity(pair_size);

        istart = cuDFNsys::CPUSecond();

        // gmsh::fltk::run();

        std::vector<std::vector<std::pair<int, int>>> outmap;

        if (NUM_frac > 1)
        {
            double istart1 = cuDFNsys::CPUSecond();
            // cout << "pair_size: " << pair_size << endl;
            for (size_t i = 0; i < pair_size; ++i)
            {
                object_entity[i].first = 2;
                tool_entity[i].first = 2;

                object_entity[i].second = Corre_Tag[IntersectionPair_percol[i].first];
                tool_entity[i].second = Corre_Tag[IntersectionPair_percol[i].second];
                // cout << IntersectionPair_percol[i].first << ", " << IntersectionPair_percol[i].second << endl;
            }
            cout << "\t\tadded object-tool pairs; running time: " << cuDFNsys::CPUSecond() - istart1 << " sec\n";

            //bool IfFrag = IntersectionPair_percol.size() == 0 ? false : true;

            if (IntersectionPair_percol.size() != 0)
            {
                double istart2 = cuDFNsys::CPUSecond();
                std::vector<std::pair<int, int>> out;

                gmsh::model::occ::fragment(object_entity, tool_entity, out, outmap);
                gmsh::model::occ::synchronize();
                // gmsh::write("exod.brep");
                cout << "\t\tfragmented entities, running time: " << cuDFNsys::CPUSecond() - istart2 << " sec\n";
                std::set<gmsh::vectorpair> set_map;
                pair<std::set<gmsh::vectorpair>::iterator, bool> kkit;

                for (size_t i = 0; i < outmap.size(); ++i)
                    kkit = set_map.insert(outmap[i]);
                outmap.clear();
                outmap.assign(set_map.begin(), set_map.end());
            }
            else
            {
                gmsh::vectorpair IO;
                gmsh::model::getEntities(IO,
                                         2);
                outmap.resize(IO.size());
                for (int u = 0; u < IO.size(); ++u)
                {
                    outmap[u].resize(1);
                    outmap[u][0] = IO[u];
                }
            }
            // cout << "outDimTags\n";
            // for (size_t i = 0; i < out.size(); ++i)
            //     cout << out[i].first << ",  " << out[i].second << endl;
            // cout << "\n NUM_frac, size = " << NUM_frac << "\n";
            // cout << "\noutDimTagsMap, size = " << outmap.size() << "\n";
            //for (size_t i = 0; i < outmap.size(); ++i)
            //{
            //    for (size_t j = 0; j < outmap[i].size(); ++j)
            //        cout << outmap[i][j].first << ",  " << outmap[i][j].second << "; ";
            //    cout << endl;
            //}

            // cout << "\nafter erase, outDimTagsMap, size = " << outmap.size() << "\n";
            // for (size_t i = 0; i < outmap.size(); ++i)
            // {
            //     for (size_t j = 0; j < outmap[i].size(); ++j)
            //         cout << outmap[i][j].first << ",  " << outmap[i][j].second << "; ";
            //     cout << endl
            //          << endl;
            // }
            // exit(0);
            // gmsh::fltk::run();
        }
        else
        {
            outmap.resize(1);
            outmap[0].resize(1);
            outmap[0][0] = std::make_pair(2, 1);
        };

        istart = cuDFNsys::CPUSecond();
        gmsh::option::setNumber("Mesh.MeshSizeMin", min_ele_edge);
        gmsh::option::setNumber("Mesh.MeshSizeMax", max_ele_edge);

        gmsh::option::setNumber("Mesh.Algorithm", 1);
        cout << "\tmesh ing" << endl;
        gmsh::model::mesh::generate(2);
        /// gmsh::fltk::run();
        cout << "\t\tmeshing finished, running time: " << cuDFNsys::CPUSecond() - istart << " sec\n";

        //gmsh::write("exod.msh");

        cout << "\t\tGeting coordinates" << endl;
        this->GetCoordinates();
        cout << "\t\tGeting elements" << endl;
        this->GetElements(Fracs_ss, outmap);
        //cudaDeviceReset();
        cout << "\t\tNumbering edges\n";
        istart = cuDFNsys::CPUSecond();
        this->NumberingEdges(L, DomainDimensionRatio /*, Fracs*/);
        cout << "\t\tFinished! Running time: " << cuDFNsys::CPUSecond() - istart << " sec\n";

        gmsh::model::mesh::clear();
        gmsh::clear();
        gmsh::finalize();
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        cout << "\t" << e.what() << endl;
        MeshSuccess = false;
        gmsh::clear();
        gmsh::finalize();
        throw;
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << "\tPause now!\n"
             << e.what() << endl;
        MeshSuccess = false;
        gmsh::clear();
        gmsh::finalize();
        exit(0);
    }
    catch (std::exception &e)
    {
        cout << "\033[31m\t" << e.what() << "\033[0m\n";
        MeshSuccess = false;
        gmsh::clear();
        gmsh::finalize();
        throw;
    };
}; // Mesh
template cuDFNsys::Mesh<double>::Mesh(const thrust::host_vector<cuDFNsys::Fracture<double>> &Fracs,
                                      const std::vector<pair<int, int>> &IntersectionPair_percol,
                                      std::vector<size_t> *Fracs_percol,
                                      const double &min_ele_edge,
                                      const double &max_ele_edge,
                                      const int &dir_,
                                      const double &L, double3 DomainDimensionRatio);
template cuDFNsys::Mesh<float>::Mesh(const thrust::host_vector<cuDFNsys::Fracture<float>> &Fracs,
                                     const std::vector<pair<int, int>> &IntersectionPair_percol,
                                     std::vector<size_t> *Fracs_percol,
                                     const float &min_ele_edge,
                                     const float &max_ele_edge,
                                     const int &dir_,
                                     const float &L, double3 DomainDimensionRatio);

// ====================================================
// NAME:        MatlabPlot
// DESCRIPTION: Plot mesh
//              (generate command file and data file)
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
double cuDFNsys::Mesh<T>::MatlabPlot(const string &mat_key,
                                     const string &command_key,
                                     thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                                     const T &L,
                                     const bool &if_check_2D_coordinates,
                                     const bool &if_check_edge_Numbering,
                                     bool if_python_visualization,
                                     string PythonName_Without_suffix, double3 DomainDimensionRatio)
{
    //cuDFNsys::MatlabAPI M1;

    cuDFNsys::HDF5API h5gg;

    h5gg.NewFile(mat_key);

    size_t node_num = this->Coordinate3D.size();
    T *ptr_coordinates_3D;
    ptr_coordinates_3D = new T[node_num * 3];

    if (ptr_coordinates_3D == NULL)
    {
        string AS = "Alloc error in Mesh::MatlabPlot\n";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    for (size_t i = 0; i < node_num; ++i)
    {
        ptr_coordinates_3D[i] = this->Coordinate3D[i].x;
        ptr_coordinates_3D[i + node_num] = this->Coordinate3D[i].y;
        ptr_coordinates_3D[i + node_num * 2] = this->Coordinate3D[i].z;
    }

    uint2 dim_f = make_uint2(3, node_num);
    h5gg.AddDataset(mat_key, "N", "coordinate_3D", ptr_coordinates_3D, dim_f);

    dim_f = make_uint2(1, 1);
    T yu[1] = {L};
    h5gg.AddDataset(mat_key, "N", "L_m", yu, dim_f);

    // M1.WriteMat(mat_key, "w", node_num * 3,
    //             node_num, 3, ptr_coordinates_3D, "coordinate_3D");

    delete[] ptr_coordinates_3D;
    ptr_coordinates_3D = NULL;

    uint2 dim_ds = make_uint2(3, 1);
    double DomainDimensionRatio_D[3] = {DomainDimensionRatio.x, DomainDimensionRatio.y, DomainDimensionRatio.z};

    h5gg.AddDataset(mat_key, "N", "DomainDimensionRatio", DomainDimensionRatio_D, dim_ds);

    size_t ele_num = this->Element3D.size();
    T *ptr_element_3D = new T[ele_num * 3];
    if (ptr_element_3D == NULL)
    {
        string AS = "Alloc error in Mesh::MatlabPlot\n";
        throw cuDFNsys::ExceptionsPause(AS);
    }

    double Area_characteristic = 0;

    for (size_t i = 0; i < ele_num; ++i)
    {
        ptr_element_3D[i] = this->Element3D[i].x;
        ptr_element_3D[i + ele_num] = this->Element3D[i].y;
        ptr_element_3D[i + ele_num * 2] = this->Element3D[i].z;

        Area_characteristic += double(cuDFNsys::Triangle3DArea<T>(this->Coordinate3D[this->Element3D[i].x - 1],
                                                                  this->Coordinate3D[this->Element3D[i].y - 1],
                                                                  this->Coordinate3D[this->Element3D[i].z - 1]));
    };
    Area_characteristic = Area_characteristic / ele_num;

    //M1.WriteMat(mat_key, "u", ele_num * 3,
    //            ele_num, 3, ptr_element_3D, "element_3D");
    dim_f = make_uint2(3, ele_num);
    h5gg.AddDataset(mat_key, "N", "element_3D", ptr_element_3D, dim_f);

    delete[] ptr_element_3D;
    ptr_element_3D = NULL;

    if (if_check_2D_coordinates == true)
    {
        size_t frac_tag_num = this->ElementFracTag.size();

        T *ptr_element_Frac_Tag = new T[frac_tag_num];
        if (ptr_element_Frac_Tag == NULL)
        {
            string AS = "Alloc error in Mesh::MatlabPlot\n";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        for (size_t i = 0; i < frac_tag_num; ++i)
            ptr_element_Frac_Tag[i] = this->ElementFracTag[i] + 1;

        //M1.WriteMat(mat_key, "u", frac_tag_num * 1,
        //            frac_tag_num, 1, ptr_element_Frac_Tag, "element_Frac_Tag");
        dim_f = make_uint2(1, frac_tag_num);
        h5gg.AddDataset(mat_key, "N", "element_Frac_Tag", ptr_element_Frac_Tag, dim_f);

        delete[] ptr_element_Frac_Tag;
        ptr_element_Frac_Tag = NULL;

        //thrust::host_vector<ELE_COOR> coordinate_2D;

        T *verts2D = new T[frac_tag_num * 6];
        if (verts2D == NULL)
        {
            string AS = "Alloc error in Mesh::MatlabPlot\n";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        for (size_t i = 0; i < frac_tag_num; ++i)
        {
            verts2D[i] = this->Coordinate2D[i].x[0];
            verts2D[i + frac_tag_num] = this->Coordinate2D[i].y[0];

            verts2D[i + 2 * frac_tag_num] = this->Coordinate2D[i].x[1];
            verts2D[i + 3 * frac_tag_num] = this->Coordinate2D[i].y[1];

            verts2D[i + 4 * frac_tag_num] = this->Coordinate2D[i].x[2];
            verts2D[i + 5 * frac_tag_num] = this->Coordinate2D[i].y[2];
        }

        // M1.WriteMat(mat_key, "u", frac_tag_num * 6,
        //             frac_tag_num, 6, verts2D, "coordinate_2D");
        dim_f = make_uint2(6, frac_tag_num);
        h5gg.AddDataset(mat_key, "N", "coordinate_2D", verts2D, dim_f);

        delete[] verts2D;
        verts2D = NULL;

        // std::vector<size_t> *Frac_ID
        T *fracs = new T[(*this->FracID).size() * 4 * 2];
        int vb = 0;
        for (std::vector<size_t>::iterator it = this->FracID->begin(); it != this->FracID->end(); it++)
        {
            size_t tag = *(it);

            cuDFNsys::Vector2<T> verts2D_[4];
            Fracs[tag].Generate2DVerts(verts2D_, 4, false);

            fracs[vb] = verts2D_[0].x;
            fracs[vb + (*this->FracID).size() * 4] = verts2D_[0].y;
            vb++;
            fracs[vb] = verts2D_[1].x;
            fracs[vb + (*this->FracID).size() * 4] = verts2D_[1].y;
            vb++;
            fracs[vb] = verts2D_[2].x;
            fracs[vb + (*this->FracID).size() * 4] = verts2D_[2].y;
            vb++;
            fracs[vb] = verts2D_[3].x;
            fracs[vb + (*this->FracID).size() * 4] = verts2D_[3].y;
            vb++;
        }
        if (fracs == NULL)
        {
            string AS = "Alloc error in Mesh::MatlabPlot\n";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        //M1.WriteMat(mat_key, "u", (*this->FracID).size() * 4 * 2, (*this->FracID).size() * 4, 2, fracs, "fracs_2D");
        dim_f = make_uint2(2, (*this->FracID).size() * 4);
        h5gg.AddDataset(mat_key, "N", "fracs_2D", fracs, dim_f);
        delete[] fracs;
        fracs = NULL;
    };

    if (if_check_edge_Numbering == true)
    {
        size_t NUM_ele = this->Element3D.size();

        T *edge_attri = new T[NUM_ele * 6];
        if (edge_attri == NULL)
        {
            string AS = "Alloc error in Mesh::MatlabPlot\n";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        for (size_t i = 0; i < NUM_ele; ++i)
        {
            edge_attri[i] = this->Element3D[i].x;
            edge_attri[i + NUM_ele] = this->Element3D[i].y;
            edge_attri[i + 2 * NUM_ele] = this->Element3D[i].z;

            edge_attri[i + 3 * NUM_ele] = this->EdgeAttri[i].e[0];
            edge_attri[i + 4 * NUM_ele] = this->EdgeAttri[i].e[1];
            edge_attri[i + 5 * NUM_ele] = this->EdgeAttri[i].e[2];
        }
        //M1.WriteMat(mat_key, "u", NUM_ele * 6,
        //            NUM_ele, 6, edge_attri, "edge_attri");

        dim_f = make_uint2(6, NUM_ele);
        h5gg.AddDataset(mat_key, "N", "edge_attri", edge_attri, dim_f);

        delete[] edge_attri;
        edge_attri = NULL;
    }

    if (if_python_visualization)
    {
        std::ofstream oss(PythonName_Without_suffix + ".py", ios::out);
        oss << "import h5py\n";
        oss << "import numpy as np\n";
        oss << "from mayavi import mlab as ML\n";
        oss << "f = h5py.File('" << mat_key << "')\n";
        oss << "coordinate_3D = np.array(f['coordinate_3D'][:])\n";
        oss << "element_3D = np.array(f['element_3D'][:])\n";
        if (if_check_edge_Numbering == true)
            oss << "edge_attri = np.transpose(np.array(f['edge_attri'][:]))\n";
        oss << "L_m = f['L_m'][:][0]\n";
        oss << "DomainDimensionRatio = f['DomainDimensionRatio'][:]\n";

        oss << "f.close()\n";
        oss << "ML.triangular_mesh(coordinate_3D[0, :], coordinate_3D[1, :], coordinate_3D[2, :], np.transpose(element_3D - 1), scalars=coordinate_3D[2, :], opacity=0.8)\n";
        oss << "ML.triangular_mesh(coordinate_3D[0, :], coordinate_3D[1, :], coordinate_3D[2, :], np.transpose(element_3D-1), representation='wireframe', color=(0, 0, 0), line_width=1.0)\n";
        oss << "ML.outline(extent=[-0.5 * DomainDimensionRatio[0] * L_m, 0.5 * DomainDimensionRatio[0] * L_m, -0.5 * DomainDimensionRatio[1] * L_m, 0.5 * DomainDimensionRatio[1] * L_m, -0.5 * DomainDimensionRatio[2] * L_m, 0.5 * DomainDimensionRatio[2] * L_m])\n";
        oss << "ML.axes()\n";
        oss << "ML.colorbar(orientation='vertical')\n";
        oss << "ML.xlabel('x (m)')\n";
        oss << "ML.ylabel('y (m)')\n";
        oss << "ML.zlabel('z (m)')\n";
        oss << "ML.show()\n";
        if (if_check_edge_Numbering == true)
        {
            oss << "kk = np.zeros([edge_attri.shape[0] * 3, 3], dtype=int)\n";
            oss << "kk[[list(range(0, kk.shape[0], 3))], :] = np.concatenate((edge_attri[:, [0, 1]], edge_attri[:, [3]]), axis=1)\n";
            oss << "kk[[list(range(1, kk.shape[0], 3))], :] = np.concatenate((edge_attri[:, [1, 2]], edge_attri[:, [4]]), axis=1)\n";
            oss << "kk[[list(range(2, kk.shape[0], 3))], :] = np.concatenate((edge_attri[:, [2, 0]], edge_attri[:, [5]]), axis=1)\n";
            oss << "inlet_loc = np.where(kk[:, 2] == 0)\n";
            oss << "outlet_loc = np.where(kk[:, 2] == 1)\n";
            oss << "neumann_loc = np.where(kk[:, 2] == 2)\n";
            oss << "interior_loc = np.where(kk[:, 2] == 3)\n";
            oss << "ML.triangular_mesh(coordinate_3D[0, :], coordinate_3D[1, :], coordinate_3D[2, :], np.transpose(element_3D - 1), representation='wireframe', color=(0, 0, 0), line_width=1.0)\n";
            oss << "ML.outline(extent=[-0.5 * DomainDimensionRatio[0] * L_m, 0.5 * DomainDimensionRatio[0] * L_m, -0.5 * DomainDimensionRatio[1] * L_m, 0.5 * DomainDimensionRatio[1] * L_m, -0.5 * DomainDimensionRatio[2] * L_m, 0.5 * DomainDimensionRatio[2] * L_m])\n";
            oss << "ML.axes()\n";
            oss << "ML.colorbar(orientation='vertical')\n";
            oss << "ML.xlabel('x (m)')\n";
            oss << "ML.ylabel('y (m)')\n";
            oss << "ML.zlabel('z (m)')\n";
            oss << "src1 = ML.pipeline.scalar_scatter(coordinate_3D[0, :], coordinate_3D[1, :], coordinate_3D[2, :])\n";
            oss << "src1.mlab_source.dataset.lines = (kk[inlet_loc, 0:2][0] - 1).tolist()\n";
            oss << "src1.update()\n";
            oss << "lines = ML.pipeline.stripper(src1)\n";
            oss << "ML.pipeline.surface(lines, color=(1, 0, 0), line_width=1.2, opacity=1)\n";
            oss << "src2 = ML.pipeline.scalar_scatter(coordinate_3D[0, :], coordinate_3D[1, :], coordinate_3D[2, :])\n";
            oss << "src2.mlab_source.dataset.lines = (kk[outlet_loc, 0:2][0] - 1).tolist()\n";
            oss << "src2.update()\n";
            oss << "lines = ML.pipeline.stripper(src2)\n";
            oss << "ML.pipeline.surface(lines, color=(1, 0, 0), line_width=1.2, opacity=1)\n";
            oss << "src3 = ML.pipeline.scalar_scatter(coordinate_3D[0, :], coordinate_3D[1, :], coordinate_3D[2, :])\n";
            oss << "src3.mlab_source.dataset.lines = (kk[neumann_loc, 0:2][0] - 1).tolist()\n";
            oss << "src3.update()\n";
            oss << "lines = ML.pipeline.stripper(src3)\n";
            oss << "ML.pipeline.surface(lines, color=(0, 0, 1), line_width=1.2, opacity=1)\n";
            oss << "src4 = ML.pipeline.scalar_scatter(coordinate_3D[0, :], coordinate_3D[1, :], coordinate_3D[2, :])\n";
            oss << "src4.mlab_source.dataset.lines = (kk[interior_loc, 0:2][0] - 1).tolist()\n";
            oss << "src4.update()\n";
            oss << "lines = ML.pipeline.stripper(src4)\n";
            oss << "ML.pipeline.surface(lines, color=(0, 1, 0), line_width=1.2, opacity=1)\n";
            oss << "ML.show()\n";
        }

        oss.close();

        if (if_check_2D_coordinates)
        {
            std::ofstream oss("CHECK_2D_fac_" + PythonName_Without_suffix + ".py", ios::out);
            oss << "import numpy as np\n";
            oss << "import matplotlib.pyplot as plt\n";
            oss << "import h5py\n";
            oss << "f = h5py.File('DFN_mesh_1.h5')\n";
            oss << "element_Frac_Tag  = np.array(f['element_Frac_Tag'][:])\n";
            oss << "coordinate_2D = np.array(f['coordinate_2D'][:])\n";
            oss << "fracs_2D = np.array(f['fracs_2D'][:])\n";
            oss << "NumFracs = int((fracs_2D.shape[1]) / 4)\n";
            oss << "print('Number of fracures is', NumFracs)\n";
            oss << "f.close()\n";
            oss << "\n";
            oss << "NumFrac_show = input('How many fractures and the associating meshes you want to visualize?\n')\n";
            oss << "NumFrac_show = int(NumFrac_show)\n";
            oss << "if (NumFrac_show > NumFracs):\n";
            oss << "    NumFrac_show = NumFracs;\n";
            oss << "    print('Sorry, the maximum number is', NumFracs)\n";
            oss << "\n";
            oss << "for i in range(NumFrac_show):\n";
            oss << "    a = np.where(element_Frac_Tag == i + 1)\n";
            oss << "    vertsLocal = coordinate_2D[:, a[1]]\n";
            oss << "    vertsLocal = np.concatenate((vertsLocal[[0, 1], :], vertsLocal[[2, 3], :], vertsLocal[[4, 5], :]), axis=1)\n";
            oss << "    NUM_ele_local = int(vertsLocal.shape[1] / 3)\n";
            oss << "    element_local = np.transpose(np.array([range(0, NUM_ele_local, 1)]))\n";
            oss << "    element_local = np.concatenate((element_local, element_local + NUM_ele_local, element_local + 2 * NUM_ele_local), axis=1)\n";
            oss << "    plt.triplot(vertsLocal[0, :], vertsLocal[1, :], element_local, 'g-')\n";
            oss << "    tmy = [i * 4 + j for j in range(4)]\n";
            oss << "    tmy = tmy + [tmy[0]]\n";
            oss << "    plt.plot(fracs_2D[0, tmy], fracs_2D[1, tmy], 'k-')\n";
            oss << "    plt.title('Frac NO. ' + str(i + 1))\n";
            oss << "    plt.show()\n";

            oss.close();
        }
    }
    //---------------------------------------------------------
    if (command_key != "N")
    {
        std::ofstream oss(command_key, ios::out);
        oss << "clc;\nclose all;\nclear all;\n";
        //oss << "load('" << mat_key << "');\n";
        oss << "currentPath = fileparts(mfilename('fullpath'));\n";
        oss << "coordinate_3D = h5read([currentPath, '/" << mat_key << "'], '/coordinate_3D');\n";
        oss << "element_3D = h5read([currentPath, '/" << mat_key << "'], '/element_3D');\n";
        oss << "L = h5read([currentPath, '/" << mat_key << "'], '/L_m');\n";
        oss << "DomainDimensionRatio = h5read([currentPath, '/" << mat_key << "'], '/DomainDimensionRatio');\n";
        oss << "cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];\n";
        oss << "cube_frame(:, 1) = 0.5 .* cube_frame(:, 1) .* DomainDimensionRatio(1); ";
        oss << "cube_frame(:, 2) = 0.5 .* cube_frame(:, 2) .* DomainDimensionRatio(2); ";
        oss << "cube_frame(:, 3) = 0.5 .* cube_frame(:, 3) .* DomainDimensionRatio(3);\n";

        oss << "figure(1); view(3); title('DFN mesh'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on\n";
        oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on\n";
        oss << "patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', coordinate_3D(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1); hold on\n";

        oss << "axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * DomainDimensionRatio(1) * L, -1.1 / 2 * DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) * L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * DomainDimensionRatio(3) * L]);\n";
        oss << "pbaspect([DomainDimensionRatio]); hold on\n";

        if (if_check_edge_Numbering == true)
        {
            oss << "edge_attri = h5read([currentPath, '/" << mat_key << "'], '/edge_attri');\n";
            oss << "figure(2); view(3); title('check edge numbering'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on\n";
            oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on\n";

            oss << "axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * DomainDimensionRatio(1) * L, -1.1 / 2 * DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) * L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * DomainDimensionRatio(3) * L]);\n";
            oss << "pbaspect([DomainDimensionRatio]); hold on\n";

            oss << "kk = zeros(size(edge_attri, 1) * 3, 3);kk([1:3:end], :) = [edge_attri(:, [1, 2]), edge_attri(:, 4)];kk([2:3:end], :) = [edge_attri(:, [2, 3]), edge_attri(:, 5)];kk([3:3:end], :) = [edge_attri(:, [3, 1]), edge_attri(:, 6)];\n\n";

            oss << "\ninlet_loc = find(kk(:, 3)==0);\n";
            oss << "patch('Vertices', coordinate_3D, 'Faces', kk(inlet_loc, [1, 2]), 'FaceVertexCData', zeros(size(coordinate_3D, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'b'); hold on\n";
            oss << "outlet_loc = find(kk(:, 3)==1);\n";
            oss << "patch('Vertices', coordinate_3D, 'Faces', kk(outlet_loc, [1, 2]), 'FaceVertexCData', zeros(size(coordinate_3D, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'b'); hold on\n";
            oss << "neumann_loc = find(kk(:, 3)==2);\n";
            oss << "patch('Vertices', coordinate_3D, 'Faces', kk(neumann_loc, [1, 2]), 'FaceVertexCData', zeros(size(coordinate_3D, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'r'); hold on\n";
            oss << "interior_loc = find(kk(:, 3)==3);\n";
            oss << "patch('Vertices', coordinate_3D, 'Faces', kk(interior_loc, [1, 2]), 'FaceVertexCData', zeros(size(coordinate_3D, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'k'); hold on\n";
        }
        oss.close();

        if (if_check_2D_coordinates == true)
        {
            std::ofstream oss_e("CHECK_2D_fac_" + command_key, ios::out);
            oss_e << "clc;\nclose all;\nclear all;\n";
            //oss_e << "load('" << mat_key << "');\n";
            oss_e << "currentPath = fileparts(mfilename('fullpath'));\n";
            oss_e << "element_Frac_Tag = h5read([currentPath, '/" << mat_key << "'], '/element_Frac_Tag');\n";
            oss_e << "fracs_2D = h5read([currentPath, '/" << mat_key << "'], '/fracs_2D');\n";
            oss_e << "coordinate_2D = h5read([currentPath, '/" << mat_key << "'], '/coordinate_2D');\n";
            oss_e << "num_frac = " << this->Element2D.size() << ";\n\n";

            oss_e << "tmp_cc = 1;\n";
            oss_e << "disp(['num_frac = ', num2str(num_frac), ';']);\n";
            oss_e << "pause(1.5);\n";

            oss_e << "for i = 1:num_frac\n";
            oss_e << "\ta = find(element_Frac_Tag == i);\n";
            oss_e << "\tfigure(1); view(2);\n";
            oss_e << "\tvertsLocal = coordinate_2D([a], :);\n";
            oss_e << "\tvertsLocal = [vertsLocal(:, [1 2]); vertsLocal(:, [3 4]); vertsLocal(:, [5 6])];\n";
            oss_e << "\tNUM_ele_local = size(vertsLocal, 1) / 3;\n";
            oss_e << "\telement_local = [[1:NUM_ele_local]', [1:NUM_ele_local]' + NUM_ele_local, [1:NUM_ele_local]' + 2 * NUM_ele_local];\n";
            oss_e << "\tpatch('Vertices', vertsLocal, 'Faces', element_local, 'FaceVertexCData', zeros(size(vertsLocal, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'k'); hold on\n";
            oss_e << "\tpatch('Vertices', fracs_2D, 'Faces', [1 2 3 4] + (i - 1) * 4, 'FaceVertexCData', zeros(size(fracs_2D, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'r'); hold on\n";

            oss_e << "\tdisp(i);\n";
            oss_e << "\tpause();\n";
            oss_e << "\tclose 1\n";
            oss_e << "end\n";
            oss_e.close();
        }
    }
    return Area_characteristic;
}; // MatlabPlot
template double cuDFNsys::Mesh<double>::MatlabPlot(const string &mat_key,
                                                   const string &command_key,
                                                   thrust::host_vector<cuDFNsys::Fracture<double>> Fracs,
                                                   const double &L,
                                                   const bool &if_check_2D_coordinates,
                                                   const bool &if_check_edge_Numbering,
                                                   bool if_python_visualization,
                                                   string PythonName_Without_suffix, double3 DomainDimensionRatio);
template double cuDFNsys::Mesh<float>::MatlabPlot(const string &mat_key,
                                                  const string &command_key,
                                                  thrust::host_vector<cuDFNsys::Fracture<float>> Fracs,
                                                  const float &L,
                                                  const bool &if_check_2D_coordinates,
                                                  const bool &if_check_edge_Numbering,
                                                  bool if_python_visualization,
                                                  string PythonName_Without_suffix, double3 DomainDimensionRatio);

// ====================================================
// NAME:        GetCoordinates
// DESCRIPTION: Get Coordinates of 3D element
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
void cuDFNsys::Mesh<T>::GetCoordinates()
{
    std::vector<size_t> nodes;
    std::vector<double> coord, coordParam;
    gmsh::model::mesh::getNodes(nodes, coord, coordParam);
    size_t NUM_nodes = coord.size() / 3;

    this->Coordinate3D.resize(NUM_nodes);

    for (size_t i = 0; i < coord.size(); i += 3)
        Coordinate3D[i / 3] = cuDFNsys::MakeVector3((T)coord[i],
                                                    (T)coord[i + 1],
                                                    (T)coord[i + 2]);
}; // GetCoordinates
template void cuDFNsys::Mesh<double>::GetCoordinates();
template void cuDFNsys::Mesh<float>::GetCoordinates();

// ====================================================
// NAME:        GetElements
// DESCRIPTION: Get elements of 2D / 3D element
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
void cuDFNsys::Mesh<T>::GetElements(const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs_s, const std::vector<std::vector<std::pair<int, int>>> &outmap)
{
    thrust::host_vector<thrust::host_vector<uint3>> elementEntities_2D;
    thrust::host_vector<uint> Largest_ele;

    this->GetEntitiesElements(elementEntities_2D, Largest_ele, outmap);
    cout << "\t\t\tGot entity elements" << endl;
    //cout << 1 << endl;

    // for (size_t i = 0; i < elementEntities_2D.size(); ++i)
    // cout << "this entity: " << elementEntities_2D[i].size() << " elements, the largest ele: " << Largest_ele[i] << endl;

    thrust::host_vector<uint3> One_entity_one_ele(elementEntities_2D.size());
    //cout << "node:\n";
    for (size_t i = 0; i < One_entity_one_ele.size(); ++i)
    {
        One_entity_one_ele[i] = elementEntities_2D[i][Largest_ele[i] - 1];
        //cout << One_entity_one_ele[i].x << ", " << One_entity_one_ele[i].y << ", " << One_entity_one_ele[i].z << endl;
    }

    thrust::device_vector<uint3> One_entity_one_ele_dev(elementEntities_2D.size());
    thrust::device_vector<cuDFNsys::Vector3<T>> coordinate_3D_dev;
    thrust::device_vector<int> Elements_Frac_dev(elementEntities_2D.size());
    thrust::host_vector<int> Elements_Frac_host;

    One_entity_one_ele_dev = One_entity_one_ele;
    coordinate_3D_dev = this->Coordinate3D;

    int *Elements_Frac_dev_ptr = thrust::raw_pointer_cast(Elements_Frac_dev.data());

    uint3 *One_entity_one_ele_dev_ptr = thrust::raw_pointer_cast(One_entity_one_ele_dev.data());
    cuDFNsys::Vector3<T> *coordinate_3D_dev_ptr = thrust::raw_pointer_cast(coordinate_3D_dev.data());

    thrust::device_vector<cuDFNsys::Fracture<T>> Fracss;
    Fracss = Fracs_s;
    cuDFNsys::Fracture<T> *Fracturesss_vert_dev_ptr = thrust::raw_pointer_cast(Fracss.data());

    cuDFNsys::IdentifyEleFrac<T><<<elementEntities_2D.size() / 256 + 1, 256>>>(One_entity_one_ele_dev_ptr,
                                                                               coordinate_3D_dev_ptr,
                                                                               Fracturesss_vert_dev_ptr,
                                                                               Elements_Frac_dev_ptr,
                                                                               elementEntities_2D.size(),
                                                                               Fracs_s.size(),
                                                                               _TOL_IdentifyEleFrac);
    cudaDeviceSynchronize();
    cout << "\t\t\tIdentified element's fracID" << endl;
    //------------------------------------------------------------------------------------------
    // cout << "elementEntities_2D.size: " << elementEntities_2D.size() << endl;
    // for (size_t i = 0; i < elementEntities_2D.size(); ++i)
    //     cout << Elements_Frac_dev[i] << endl;

    Elements_Frac_host = Elements_Frac_dev;

    for (int i = 0; i < elementEntities_2D.size(); ++i)
    {
        if (Elements_Frac_host[i] == -1)
        {
            string AS = "entity: " + std::to_string(i) + ", ";
            AS += "cannot find which fracture that the enetity is belonging to!\n";
            //cout << AS;
            throw cuDFNsys::ExceptionsIgnore(AS);
        }
    }

    this->Element2D.resize((*(this->FracID)).size());

    for (int i = 0; i < Elements_Frac_host.size(); ++i)
    {
        int FracIDLocal = Elements_Frac_host[i];
        // cout << "local ID: " << FracIDLocal << ", elementEntities_2D.size = " << elementEntities_2D[i].size() << endl;
        this->Element2D[FracIDLocal].insert(this->Element2D[FracIDLocal].end(),
                                            elementEntities_2D[i].begin(),
                                            elementEntities_2D[i].end());
    }

    for (int i = 0; i < this->Element2D.size(); ++i)
    {
        if (this->Element2D[i].size() == 0)
        {
            string AS = "one frac has no grids!\n";
            //cout << AS;
            //exit(0);
            throw cuDFNsys::ExceptionsIgnore(AS);
        }
        this->Element3D.insert(this->Element3D.end(),
                               this->Element2D[i].begin(),
                               this->Element2D[i].end());

        thrust::host_vector<uint> KK(this->Element2D[i].size(),
                                     i);

        this->ElementFracTag.insert(this->ElementFracTag.end(),
                                    KK.begin(),
                                    KK.end());
    };
    //------------------
    /*
        for (int i = 0; i < this->Element2D.size(); ++i)
        {
            cout << "Frac " << i + 1 << endl;
            for (int j = 0; j < this->Element2D[i].size(); ++j)
            {
                cout << this->Element2D[i][j].x << ", " << this->Element2D[i][j].y << ", " << this->Element2D[i][j].z << "\n";
            }
            cout << endl;
        }*/

    thrust::device_vector<uint3> element_3D_dev;
    element_3D_dev = this->Element3D;
    uint3 *element_3D_dev_ptr = thrust::raw_pointer_cast(element_3D_dev.data());

    thrust::device_vector<cuDFNsys::EleCoor<T>> coordinate_2D_dev(this->Element3D.size());
    cuDFNsys::EleCoor<T> *coordinate_2D_dev_ptr = thrust::raw_pointer_cast(coordinate_2D_dev.data());

    thrust::device_vector<uint> element_Frac_Tag_dev;
    element_Frac_Tag_dev = ElementFracTag;
    uint *element_Frac_Tag_dev_ptr = thrust::raw_pointer_cast(element_Frac_Tag_dev.data());

    cuDFNsys::GetLocalCoordiates<T><<<this->Element3D.size() / 256 + 1, 256>>>(element_3D_dev_ptr,
                                                                               Fracturesss_vert_dev_ptr,
                                                                               element_Frac_Tag_dev_ptr,
                                                                               coordinate_2D_dev_ptr,
                                                                               coordinate_3D_dev_ptr,
                                                                               this->Element3D.size());
    cudaDeviceSynchronize();
    this->Element3D = element_3D_dev;
    this->Coordinate2D = coordinate_2D_dev;
}; // GetElements
template void cuDFNsys::Mesh<double>::GetElements(const thrust::host_vector<cuDFNsys::Fracture<double>> &Fracs_s, const std::vector<std::vector<std::pair<int, int>>> &outmap);
template void cuDFNsys::Mesh<float>::GetElements(const thrust::host_vector<cuDFNsys::Fracture<float>> &Fracs_s, const std::vector<std::vector<std::pair<int, int>>> &outmap);

// ====================================================
// NAME:        NumberingEdges
// DESCRIPTION: Numbering edges of mesh
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
void cuDFNsys::Mesh<T>::NumberingEdges(const T L, double3 DomainDimensionRatio)
{
    uint NUM_ele = this->Element3D.size();
    uint NUM_node = this->Coordinate3D.size();

    this->EdgeAttri.resize(NUM_ele);

    //SparseMatrix<size_t> Shared_edge_global_NO(NUM_node, NUM_node);
    //Shared_edge_global_NO.reserve(VectorXi::Constant(NUM_node, 4));
    UMapEdge Shared_edge_global_NO;
    Shared_edge_global_NO.reserve(3 * NUM_node);

    size_t edge_shared = 1;

    size_t Sep_edge_NO_in = 1;
    size_t Sep_edge_NO_out = 1;
    size_t Sep_edge_NO_neumann = 1;

    this->InletEdgeNOLen.reserve(NUM_ele * 3);
    this->OutletEdgeNOLen.reserve(NUM_ele * 3);

    //----
    //SparseMatrix<float> shared_edge_len(NUM_node, NUM_node);
    //shared_edge_len.reserve(VectorXi::Constant(NUM_node, 4));
    //double istart = 0;
    //int element_count = -1;
    for (size_t i = 0; i < this->Element2D.size(); ++i)
    {
        bool if_change_ori = false;
        size_t eleNO_tmp_1 = this->GetElementID(i, 0);

        if (this->Element2D[i][0].y != this->Element3D[eleNO_tmp_1].y &&
            this->Element2D[i][0].z != this->Element3D[eleNO_tmp_1].z)
            if_change_ori = true;

        //istart = cuDFNsys::CPUSecond();

        UMapEdge node2element_local = this->SparseMatEdgeAttri(i, if_change_ori);

        //cout << "i: " << i << " gen Matrix " << cuDFNsys::CPUSecond() - istart << " sec\n";

        //istart = cuDFNsys::CPUSecond();

        /// get 2D coordinates of the fracture at present
        /// get 2D coordinates of the fracture at present
        /// get 2D coordinates of the fracture at present
        //size_t FracID__ = this->ElementFracTag[element_count + 1];
        // cuDFNsys::Vector2<T> verts2DDD[Fracs[FracID__].NumVertsTruncated];
        // cuDFNsys::Fracture<T> FL = Fracs[FracID__];
        // FL.Generate2DVerts(verts2DDD, FL.NumVertsTruncated, true);
        // cuDFNsys::Vector1<T> tmp_R_1[3][3];
        // FL.RoationMatrix(tmp_R_1, 32);

        for (size_t j = 0; j < this->Element2D[i].size(); ++j)
        {
            // element_count++;

            size_t tmp_ele_NO = this->GetElementID(i, j);

            uint element_list_e[3] = {this->Element2D[i][j].x,
                                      this->Element2D[i][j].y,
                                      this->Element2D[i][j].z};

            for (size_t k = 0; k < 3; ++k)
            {
                size_t Sep_NO = tmp_ele_NO * 3 + k + 1;

                this->EdgeAttri[tmp_ele_NO].e[k] = -1;
                this->EdgeAttri[tmp_ele_NO].no[k] = -1;

                size_t node1 = element_list_e[k],
                       node2 = element_list_e[(k + 1) % 3];

                pair<size_t, size_t> key_ = make_pair(node1 < node2 ? node1 : node2, node1 > node2 ? node1 : node2);

                if (node2element_local[key_] > 0)
                {
                    this->EdgeAttri[tmp_ele_NO].e[k] = 3; // interior

                    if (Shared_edge_global_NO.find(key_) == Shared_edge_global_NO.end())
                    {
                        this->EdgeAttri[tmp_ele_NO].no[k] = edge_shared;
                        Shared_edge_global_NO[key_] = edge_shared;
                        //Shared_edge_global_NO.coeffRef(node2 - 1, node1 - 1) = edge_shared;
                        edge_shared++;

                        //shared_edge_len.coeffRef(node1 - 1, node2 - 1) = len;
                        //shared_edge_len.coeffRef(node2 - 1, node1 - 1) = len;
                    }
                    else
                    {
                        this->EdgeAttri[tmp_ele_NO].no[k] = Shared_edge_global_NO[key_]; // Sep_NO;

                        //cout << "compare !!! ";
                        //printf("len: %lf, pre: %lf, error: %lf%\n", len, shared_edge_len.coeffRef(node1 - 1, node2 - 1), abs(len - shared_edge_len.coeffRef(node1 - 1, node2 - 1)) * 100);
                    }
                }
                else
                {
                    pair<bool, string> if_d = this->IfTwoEndsDirchlet(node1,
                                                                      node2,
                                                                      L, DomainDimensionRatio);

                    if (if_d.first == true)
                    {
                        cuDFNsys::Vector3<T> vert1 = this->Coordinate3D[node1 - 1];
                        cuDFNsys::Vector3<T> vert2 = this->Coordinate3D[node2 - 1];
                        cuDFNsys::Vector3<T> vect = cuDFNsys::MakeVector3(vert1.x - vert2.x,
                                                                          vert1.y - vert2.y,
                                                                          vert1.z - vert2.z);
                        T len = sqrt(vect.x * vect.x + vect.y * vect.y + vect.z * vect.z);

                        if (if_d.second == "in")
                        {
                            this->EdgeAttri[tmp_ele_NO].e[k] = 0;
                            this->EdgeAttri[tmp_ele_NO].no[k] = Sep_NO;
                            this->InletEdgeNOLen.push_back(cuDFNsys::MakeVector2((T)Sep_NO, len));

                            Sep_edge_NO_in++;
                        }
                        else
                        {
                            this->EdgeAttri[tmp_ele_NO].e[k] = 1;
                            this->EdgeAttri[tmp_ele_NO].no[k] = Sep_NO;
                            this->OutletEdgeNOLen.push_back(cuDFNsys::MakeVector2((T)Sep_NO, len));

                            Sep_edge_NO_out++;
                        }
                    }
                    else
                    {
                        this->EdgeAttri[tmp_ele_NO].e[k] = 2;
                        this->EdgeAttri[tmp_ele_NO].no[k] = Sep_NO;

                        // determine if two nodes are inside the fracture or not,
                        // if inside, then the edge must be problematic: it is not a neumman edge!!!

                        //if (element_count + 1 == 254555)
                        //{
                        //    // if (element_count != 0)
                        //    // exit(0);
                        //
                        //    cuDFNsys::Vector3<T> coord_1_s = cuDFNsys::MakeVector3(this->Coordinate3D[node1 - 1].x /*- FL.Center.x*/,
                        //                                                           this->Coordinate3D[node1 - 1].y /*- FL.Center.y*/,
                        //                                                           this->Coordinate3D[node1 - 1].z /*- FL.Center.z*/);
                        //    cuDFNsys::Vector3<T> coord_2_s = cuDFNsys::MakeVector3(this->Coordinate3D[node2 - 1].x /*- FL.Center.x*/,
                        //                                                           this->Coordinate3D[node2 - 1].y /*- FL.Center.y*/,
                        //                                                           this->Coordinate3D[node2 - 1].z /*- FL.Center.z*/);
                        //
                        //    // coord_1_s = cuDFNsys::ProductSquare3Vector3<T>(tmp_R_1, coord_1_s);
                        //    // coord_2_s = cuDFNsys::ProductSquare3Vector3<T>(tmp_R_1, coord_2_s);
                        //    // cuDFNsys::Vector2<T> coord_1_k = cuDFNsys::MakeVector2(coord_1_s.x, coord_1_s.y);
                        //    // cuDFNsys::Vector2<T> coord_2_k = cuDFNsys::MakeVector2(coord_2_s.x, coord_2_s.y);
                        //    // bool pls_1 = cuDFNsys::IfPntInside2DConvexPoly<T>(coord_1_k, verts2DDD, FL.NumVertsTruncated);
                        //    // bool pls_2 = cuDFNsys::IfPntInside2DConvexPoly<T>(coord_2_k, verts2DDD, FL.NumVertsTruncated);
                        //
                        //    for (int uuo = 0; uuo < Fracs[FracID__].NumVertsTruncated; ++uuo)
                        //    {
                        //        bool pls_1 = false, pls_2 = false;
                        //
                        //        cuDFNsys::Vector3<T> end_1 = Fracs[FracID__].Verts3DTruncated[uuo];
                        //        cuDFNsys::Vector3<T> end_2 = Fracs[FracID__].Verts3DTruncated[(uuo + 1) % Fracs[FracID__].NumVertsTruncated];
                        //
                        //        cuDFNsys::Vector3<T> A1, B1, A2, B2;
                        //
                        //        A1.x = end_1.x - coord_1_s.x,
                        //        A1.y = end_1.y - coord_1_s.y,
                        //        A1.z = end_1.z - coord_1_s.z;
                        //        B1.x = end_2.x - coord_1_s.x,
                        //        B1.y = end_2.y - coord_1_s.y,
                        //        B1.z = end_2.z - coord_1_s.z;
                        //
                        //        A2.x = end_1.x - coord_2_s.x,
                        //        A2.y = end_1.y - coord_2_s.y,
                        //        A2.z = end_1.z - coord_2_s.z;
                        //        B2.x = end_2.x - coord_2_s.x,
                        //        B2.y = end_2.y - coord_2_s.y,
                        //        B2.z = end_2.z - coord_2_s.z;
                        //
                        //        T norm_A1 = sqrt(A1.x * A1.x + A1.y * A1.y + A1.z * A1.z);
                        //        T norm_A2 = sqrt(A2.x * A2.x + A2.y * A2.y + A2.z * A2.z);
                        //        T norm_B1 = sqrt(B1.x * B1.x + B1.y * B1.y + B1.z * B1.z);
                        //        T norm_B2 = sqrt(B2.x * B2.x + B2.y * B2.y + B2.z * B2.z);
                        //
                        //        T normError = 1e-5;
                        //        T degreeError = 0.5;
                        //
                        //        if (norm_A1 < normError || norm_B1 < normError)
                        //            pls_1 = true;
                        //
                        //        if (norm_A2 < normError || norm_B2 < normError)
                        //            pls_2 = true;
                        //
                        //        T theta1, theta2;
                        //        if (norm_A1 > normError && norm_B1 > normError)
                        //        {
                        //            theta1 = acos((A1.x * B1.x + A1.y * B1.y + A1.z * B1.z) / (norm_A1 * norm_B1)) * 180.0 / M_PI;
                        //            if (abs(theta1 - 180.0) < degreeError)
                        //                pls_1 = true;
                        //        }
                        //
                        //        if (norm_A2 > normError && norm_B2 > normError)
                        //        {
                        //            theta2 = acos((A2.x * B2.x + A2.y * B2.y + A2.z * B2.z) / (norm_A2 * norm_B2)) * 180.0 / M_PI;
                        //            if (abs(theta2 - 180.0) < degreeError)
                        //                pls_2 = true;
                        //        }
                        //
                        //        cout << "Edge node: [" << end_1.x << ", " << end_1.y << ", " << end_1.z << "], [";
                        //        cout << end_2.x << ", " << end_2.y << ", " << end_2.z << "\n";
                        //        cout << "Iden: " << pls_1 << ", " << pls_2 << "; ";
                        //        cout << "norm_A1: " << norm_A1 << ", ";
                        //        cout << "norm_A2: " << norm_A2 << ", ";
                        //        cout << "norm_B1: " << norm_B1 << ", ";
                        //        cout << "norm_B2: " << norm_B2 << ", ";
                        //        cout << "theta1: " << theta1 << ", ";
                        //        cout << "theta2: " << theta2 << "\n";
                        //        if (pls_1 && pls_2)
                        //            break;
                        //
                        //        if ((!pls_1 || !pls_2) && uuo == Fracs[FracID__].NumVertsTruncated - 1)
                        //        {
                        //
                        //            cout << "\t\tFracid: " << FracID__ << ", boundID: " << uuo << ", Found a problematic neumman edge in element ID: " << element_count + 1 << "; ";
                        //            cout << "node1: " << this->Coordinate3D[node1 - 1].x << ", ";
                        //            cout << this->Coordinate3D[node1 - 1].y << ", " << this->Coordinate3D[node1 - 1].z << "; node2: ";
                        //            cout << this->Coordinate3D[node2 - 1].x << ", " << this->Coordinate3D[node2 - 1].y << ", " << this->Coordinate3D[node2 - 1].z << endl;
                        //        }
                        //    }
                        //}

                        Sep_edge_NO_neumann++;
                    }
                }
            }
        }
        //cout << "i: " << i << " identify edge " << cuDFNsys::CPUSecond() - istart << " sec\n";
    };

    this->InletEdgeNOLen.shrink_to_fit();
    this->OutletEdgeNOLen.shrink_to_fit();

    NumInteriorEdges = edge_shared - 1;
    NumInletEdges = Sep_edge_NO_in - 1;
    NumOutletEdges = Sep_edge_NO_out - 1;
    NumNeumannEdges = Sep_edge_NO_neumann - 1;
}; // NumberingEdges
template void cuDFNsys::Mesh<double>::NumberingEdges(const double L, double3 DomainDimensionRatio /*, const thrust::host_vector<cuDFNsys::Fracture<double>> &Fracs*/);
template void cuDFNsys::Mesh<float>::NumberingEdges(const float L, double3 DomainDimensionRatio /*, const thrust::host_vector<cuDFNsys::Fracture<float>> &Fracs*/);

// ====================================================
// NAME:        GetEntitiesElements
// DESCRIPTION: Get elements in each entity (GMSH)
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
void cuDFNsys::Mesh<T>::GetEntitiesElements(thrust::host_vector<thrust::host_vector<uint3>> &elementEntities_2D,
                                            thrust::host_vector<uint> &Largest_ele, const std::vector<std::vector<std::pair<int, int>>> &outmap)
{
    elementEntities_2D.resize(outmap.size());
    Largest_ele.resize(outmap.size());
    vector<size_t> Zero_element_entityNO;

    for (size_t i = 0; i < outmap.size(); ++i)
    {
        Largest_ele[i] = 0;
        T area_ll = 0;

        std::vector<std::size_t> elemNodeTags;
        for (size_t j = 0; j < outmap[i].size(); ++j)
        {
            std::vector<int> elemTypes;
            std::vector<std::vector<std::size_t>> elemTags, elemNodeTags_hh;
            gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags_hh, 2, outmap[i][j].second);

            elemNodeTags.insert(elemNodeTags.end(), elemNodeTags_hh[0].begin(), elemNodeTags_hh[0].end());
        }
        size_t NUM_ele = elemNodeTags.size() / 3;

        elementEntities_2D[i].reserve(NUM_ele);
        for (size_t j = 0; j < NUM_ele * 3; j += 3)
        {
            int node1 = (size_t)elemNodeTags[j];
            int node2 = (size_t)elemNodeTags[j + 1];
            int node3 = (size_t)elemNodeTags[j + 2];

            bool skinny_if = cuDFNsys::If3DTriangleSkinny<T>(this->Coordinate3D[node1 - 1],
                                                             this->Coordinate3D[node2 - 1],
                                                             this->Coordinate3D[node3 - 1],
                                                             _TOL_If3DTriangleSkinny);

            if (skinny_if == false)
            {
                elementEntities_2D[i].push_back(make_uint3(node1, node2, node3));
                //cout << "YY " << area << "; node " << RowVector3d(node1, node2, node3) << endl;

                T area = cuDFNsys::Triangle3DArea<T>(this->Coordinate3D[node1 - 1],
                                                     this->Coordinate3D[node2 - 1],
                                                     this->Coordinate3D[node3 - 1]);
                if (area > area_ll)
                {
                    area_ll = area;
                    Largest_ele[i] = elementEntities_2D[i].size();
                }
            }
        }
        elementEntities_2D[i].shrink_to_fit();

        if (elementEntities_2D[i].size() < 1)
        {
            Zero_element_entityNO.push_back(i);
            //string AS = "Error! One entity has zero element!\n";
            //throw Error_throw_ignore(AS);
        }
    };

    if (Zero_element_entityNO.size() > 0)
    {
        cout << "\t\t\033[33mFound zero element entities! I am removing these entities\033[0m\n";
        thrust::host_vector<thrust::host_vector<uint3>> tmp_e(elementEntities_2D.size() - Zero_element_entityNO.size());
        thrust::host_vector<uint> tmp_l(elementEntities_2D.size() - Zero_element_entityNO.size());

        for (int i = 0, j = 0; i < elementEntities_2D.size(); ++i)
        {
            if (find(Zero_element_entityNO.begin(), Zero_element_entityNO.end(), i) == Zero_element_entityNO.end())
            {
                tmp_e[j] = elementEntities_2D[i];
                tmp_l[j] = Largest_ele[i];
                j++;
            }
        }
        cout << "\t\t\033[33mremoved!\033[0m\n";
        elementEntities_2D.clear();
        elementEntities_2D.resize(tmp_e.size());
        elementEntities_2D = tmp_e;

        Largest_ele.clear();
        Largest_ele.resize(tmp_l.size());
        Largest_ele = tmp_l;
    }
}; // GetEntitiesElements
template void cuDFNsys::Mesh<double>::GetEntitiesElements(thrust::host_vector<thrust::host_vector<uint3>> &elementEntities_2D,
                                                          thrust::host_vector<uint> &Largest_ele, const std::vector<std::vector<std::pair<int, int>>> &outmap);
template void cuDFNsys::Mesh<float>::GetEntitiesElements(thrust::host_vector<thrust::host_vector<uint3>> &elementEntities_2D,
                                                         thrust::host_vector<uint> &Largest_ele, const std::vector<std::vector<std::pair<int, int>>> &outmap);

// ====================================================
// NAME:        SparseMatEdgeAttri
// DESCRIPTION: generate a sparse matrix (Eigen)
//              to check edge attribute
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
UMapEdge cuDFNsys::Mesh<T>::SparseMatEdgeAttri(uint i, bool if_change_ori)
{
    size_t Dim = this->Coordinate3D.size();

    UMapEdge umap_s;
    umap_s.reserve(this->Element2D[i].size() * 3);

    for (int j = 0; j < this->Element2D[i].size(); ++j)
    {
        uint node_list_e[3];

        if (if_change_ori == false)
        {
            node_list_e[0] = this->Element2D[i][j].x;
            node_list_e[1] = this->Element2D[i][j].y;
            node_list_e[2] = this->Element2D[i][j].z;
        }
        else
        {
            uint tmp_k = this->Element2D[i][j].y;
            this->Element2D[i][j].y = this->Element2D[i][j].z;
            this->Element2D[i][j].z = tmp_k;

            node_list_e[0] = this->Element2D[i][j].x;
            node_list_e[1] = this->Element2D[i][j].y;
            node_list_e[2] = this->Element2D[i][j].z;
        }

        for (size_t k = 0; k < 3; ++k)
        {
            size_t node1 = node_list_e[k];
            size_t node2 = node_list_e[(k + 1) % 3];

            pair<size_t, size_t> key_ = make_pair(node1 < node2 ? node1 : node2, node1 > node2 ? node1 : node2);

            if (umap_s.find(key_) == umap_s.end())
                umap_s[key_] = 0;
            else
                umap_s[key_]++;
        }
    }

    return umap_s;
}; // SparseMatEdgeAttri
template UMapEdge cuDFNsys::Mesh<double>::SparseMatEdgeAttri(uint i, bool if_change_ori);
template UMapEdge cuDFNsys::Mesh<float>::SparseMatEdgeAttri(uint i, bool if_change_ori);

// ====================================================
// NAME:        GetElementID
// DESCRIPTION: get elementID
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
size_t cuDFNsys::Mesh<T>::GetElementID(size_t i, size_t j)
{
    size_t ID = 0;

    if (i == 0)
        return ID + j;

    for (size_t y = 0; y < i; ++y)
        ID += this->Element2D[y].size();

    return ID + j;
}; // GetElementID
template size_t cuDFNsys::Mesh<double>::GetElementID(size_t i, size_t j);
template size_t cuDFNsys::Mesh<float>::GetElementID(size_t i, size_t j);

// ====================================================
// NAME:        IfTwoEndsDirchlet
// DESCRIPTION: Identify dirchlet boundary
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
pair<bool, string> cuDFNsys::Mesh<T>::IfTwoEndsDirchlet(const size_t node1,
                                                        const size_t node2,
                                                        const T L, double3 DomainDimensionRatio)
{
    pair<bool, string> kk = std::make_pair(false, "N");

    T coord_1[3] = {this->Coordinate3D[node1 - 1].x,
                    this->Coordinate3D[node1 - 1].y,
                    this->Coordinate3D[node1 - 1].z};

    T coord_2[3] = {this->Coordinate3D[node2 - 1].x,
                    this->Coordinate3D[node2 - 1].y,
                    this->Coordinate3D[node2 - 1].z};

    double *DomainDimensionRatio_rr = &DomainDimensionRatio.x;

    if (abs(coord_1[this->Dir] - DomainDimensionRatio_rr[this->Dir] * L * 0.5) < _TOL_IfTwoEndsDirchlet &&
        abs(coord_2[this->Dir] - DomainDimensionRatio_rr[this->Dir] * L * 0.5) < _TOL_IfTwoEndsDirchlet)
    {
        kk.first = true;
        kk.second = "in";
        return kk;
    }

    if (abs(coord_1[this->Dir] - DomainDimensionRatio_rr[this->Dir] * L * (-0.5)) < _TOL_IfTwoEndsDirchlet &&
        abs(coord_2[this->Dir] - DomainDimensionRatio_rr[this->Dir] * L * (-0.5)) < _TOL_IfTwoEndsDirchlet)
    {
        kk.first = true;
        kk.second = "out";
        return kk;
    }

    return kk;
}; // IfTwoEndsDirchlet
template pair<bool, string> cuDFNsys::Mesh<double>::IfTwoEndsDirchlet(const size_t node1,
                                                                      const size_t node2,
                                                                      const double L, double3 DomainDimensionRatio);
template pair<bool, string> cuDFNsys::Mesh<float>::IfTwoEndsDirchlet(const size_t node1,
                                                                     const size_t node2,
                                                                     const float L, double3 DomainDimensionRatio);