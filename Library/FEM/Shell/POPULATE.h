#pragma once
#include <FEM/Shell/ORIGAMI.h>
// for populate a certain smocking pattern on a garment stage

namespace py = pybind11;
namespace JGSL {

template<class T, int dim = 3>
void populate_pattern_hinge(MESH_NODE<T, dim>& X_0, std::vector<VECTOR<int, 3>>& rodHinge, 
    std::vector<VECTOR<T, 3>>& rodHingeInfo, MESH_ELEM<dim - 1>& Smock_pattern,
    int start_idx, int end_idx, 
    int stage_start, int stage_size, 
    T scale, VECTOR<T, 3>& x_axis, VECTOR<T, 3>& y_axis, T stiffness, T thickness)
{   
    x_axis = x_axis.Normalized();
    y_axis = y_axis.Normalized();
    // determine the boundary
    const VECTOR<T, dim>& root = std::get<0>(X_0.Get_Unchecked(start_idx)); // pop root
    std::cout << "The root pos is: " << root[0] << " " << root[1] << std::endl;
    const VECTOR<T, dim>& terminal = std::get<0>(X_0.Get_Unchecked(end_idx)); // pop terminal
    std::cout << "The terminal pos is: " << terminal[0] << " " << terminal[1] << std::endl;
    int x_times = std::abs((terminal - root).dot(x_axis) / (3.0 * scale));
    int y_times = std::abs((terminal - root).dot(y_axis) / (2.0 * scale));
    std::cout << "The repeat info is: " << x_times << " " << y_times << std::endl;
    
    int smock_hinge_size = 0;
    // braid pattern
    for(int i  = 0; i < x_times; i++){
        for(int j = 0; j < y_times; j++){
            // four smocking nodes for an unit braid pattern
            VECTOR<T, dim> target1 = root + (scale + i * 3 * scale) * x_axis - (j * 2 * scale) * y_axis; // VECTOR<T, dim>(scale + i * 3 * scale, - j * 2 * scale, 0.0);
            VECTOR<T, dim> target2 = root + (i * 3 * scale) * x_axis - (scale + j * 2 * scale) * y_axis; // VECTOR<T, dim>(i * 3 * scale, -(scale + j * 2 * scale), 0.0);

            VECTOR<T, dim> target3 = root + (scale + i * 3 * scale) * x_axis - (scale + j * 2 * scale) * y_axis; // VECTOR<T, dim>(scale + i * 3 * scale, -(scale + j * 2 * scale), 0.0);
            
            VECTOR<T, dim> target4 = root + (scale * 2 + i * 3 * scale) * x_axis - (scale * 2 + j * 2 * scale) * y_axis; //VECTOR<T, dim>(scale * 2 + i * 3 * scale, -(scale * 2 + j * 2 * scale), 0.0);
            VECTOR<T, dim> target5 = root + (scale + i * 3 * scale) * x_axis - (scale * 2 + j * 2 * scale) * y_axis; //VECTOR<T, dim>(scale * 2 + i * 3 * scale, -(scale * 2 + j * 2 * scale), 0.0);

            double dist1 = 100.0;
            double dist2 = 100.0;
            double dist3 = 100.0;
            double dist4 = 100.0;
            double dist5 = 100.0;

            int found_node1, found_node2, found_node3, found_node4, found_node5;
            VECTOR<T, 3> hingeInfo;
            hingeInfo[0] = stiffness;
            hingeInfo[1] = 2 * scale;
            hingeInfo[2] = thickness;
            
            VECTOR<int, 3> hinge1, hinge2;
            VECTOR<int, 3> smock1, smock2;

            // find target within the stage
            for(int k = 0; k < stage_size; k++){
                int cur_idx = stage_start + k;
                const VECTOR<T, dim>& cur_node = std::get<0>(X_0.Get_Unchecked(cur_idx)); 
                double cur_dist1 = (target1 - cur_node).length();
                double cur_dist2 = (target2 - cur_node).length();
                double cur_dist3 = (target3 - cur_node).length();
                double cur_dist4 = (target4 - cur_node).length();
                double cur_dist5 = (target5 - cur_node).length();

                if(cur_dist1 < dist1){
                    dist1 = cur_dist1;
                    found_node1 = cur_idx;
                }  

                if(cur_dist2 < dist2){
                    dist2 = cur_dist2;
                    found_node2 = cur_idx;
                }
                
                if(cur_dist3 < dist3){
                    dist3 = cur_dist3;
                    found_node3 = cur_idx;
                } 
                
                if(cur_dist4 < dist4){
                    dist4 = cur_dist4;
                    found_node4 = cur_idx;
                }

                if(cur_dist5 < dist5){
                    dist5 = cur_dist5;
                    found_node5 = cur_idx;
                }
            }

            hinge1[0] = found_node1;
            hinge1[1] = found_node3;
            hinge1[2] = found_node2;
            
            hinge2[0] = found_node3;
            hinge2[1] = found_node5;
            hinge2[2] = found_node4;

            printf("Cur hinge1 idx: h1_0 = %d, h1_1 = %d, h1_2 = %d\n", hinge1[0], hinge1[1],hinge1[2]);
            printf("Cur hinge2 idx: h2_0 = %d, h2_1 = %d, h2_2 = %d\n", hinge2[0], hinge2[1],hinge2[2]);

            smock1(0) = found_node1;
            smock1(1) = found_node2;
            smock1(2) = -1;

            smock2(0) = found_node3;
            smock2(1) = found_node4;
            smock2(2) = -1;

            rodHinge.push_back(hinge1);
            rodHinge.push_back(hinge2);

            rodHingeInfo.push_back(hingeInfo);
            rodHingeInfo.push_back(hingeInfo);

            Smock_pattern.Append(smock1);
            Smock_pattern.Append(smock2);

            smock_hinge_size += 2;
        }
    }
    std::cout << "The size of the smock rod hinge is: " << smock_hinge_size << std::endl;
}

template<class T, int dim = 3>
void populate_pattern(MESH_NODE<T, dim>& X_0, std::vector<VECTOR<int, 3>>& smock_stitch, 
    std::vector<T>& stitchRatio, MESH_ELEM<dim - 1>& Smock_pattern,
    int start_idx, int end_idx, 
    int stage_start, int stage_size, 
    T scale, VECTOR<T, 3>& x_axis, VECTOR<T, 3>& y_axis)
{   
    x_axis = x_axis.Normalized();
    y_axis = y_axis.Normalized();
    // determine the boundary
    const VECTOR<T, dim>& root = std::get<0>(X_0.Get_Unchecked(start_idx)); // pop root
    std::cout << "The root pos is: " << root[0] << " " << root[1] << std::endl;
    const VECTOR<T, dim>& terminal = std::get<0>(X_0.Get_Unchecked(end_idx)); // pop terminal
    std::cout << "The terminal pos is: " << terminal[0] << " " << terminal[1] << std::endl;
    int x_times = std::abs((terminal - root).dot(x_axis) / (3.0 * scale));
    int y_times = std::abs((terminal - root).dot(y_axis) / (2.0 * scale));
    std::cout << "The repeat info is: " << x_times << " " << y_times << std::endl;
    
    int smock_stitch_size = 0;
    // braid pattern
    for(int i  = 0; i < x_times; i++){
        for(int j = 0; j < y_times; j++){
            // four smocking nodes for an unit braid pattern
            VECTOR<T, dim> target1 = root + (scale + i * 3 * scale) * x_axis - (j * 2 * scale) * y_axis; // VECTOR<T, dim>(scale + i * 3 * scale, - j * 2 * scale, 0.0);
            VECTOR<T, dim> target2 = root + (i * 3 * scale) * x_axis - (scale + j * 2 * scale) * y_axis; // VECTOR<T, dim>(i * 3 * scale, -(scale + j * 2 * scale), 0.0);

            VECTOR<T, dim> target3 = root + (scale + i * 3 * scale) * x_axis - (scale + j * 2 * scale) * y_axis; // VECTOR<T, dim>(scale + i * 3 * scale, -(scale + j * 2 * scale), 0.0);
            VECTOR<T, dim> target4 = root + (scale * 2 + i * 3 * scale) * x_axis - (scale * 2 + j * 2 * scale) * y_axis; //VECTOR<T, dim>(scale * 2 + i * 3 * scale, -(scale * 2 + j * 2 * scale), 0.0);

            double dist1 = 100.0;
            double dist2 = 100.0;
            double dist3 = 100.0;
            double dist4 = 100.0;

            int found_node1, found_node2, found_node3, found_node4;
            VECTOR<int, 3> stitch1, stitch2;
            VECTOR<int, 3> smock1, smock2;

            // find target within the stage
            for(int k = 0; k < stage_size; k++){
                int cur_idx = stage_start + k;
                const VECTOR<T, dim>& cur_node = std::get<0>(X_0.Get_Unchecked(cur_idx)); 
                double cur_dist1 = (target1 - cur_node).length();
                double cur_dist2 = (target2 - cur_node).length();
                double cur_dist3 = (target3 - cur_node).length();
                double cur_dist4 = (target4 - cur_node).length();

                if(cur_dist1 < dist1){
                    dist1 = cur_dist1;
                    found_node1 = cur_idx;
                }  

                if(cur_dist2 < dist2){
                    dist2 = cur_dist2;
                    found_node2 = cur_idx;
                }
                
                if(cur_dist3 < dist3){
                    dist3 = cur_dist3;
                    found_node3 = cur_idx;
                } 
                
                if(cur_dist4 < dist4){
                    dist4 = cur_dist4;
                    found_node4 = cur_idx;
                }
            }

            // also need to add smock pattern

            stitch1[0] = found_node1;
            stitch1[1] = found_node2;
            stitch1[2] = found_node2;
            
            stitch2[0] = found_node3;
            stitch2[1] = found_node4;
            stitch2[2] = found_node4;

            smock1(0) = found_node1;
            smock1(1) = found_node2;
            smock1(2) = -1;

            smock2(0) = found_node3;
            smock2(1) = found_node4;
            smock2(2) = -1;

            smock_stitch.push_back(stitch1);
            smock_stitch.push_back(stitch2);

            stitchRatio.push_back(1.0);
            stitchRatio.push_back(1.0);

            Smock_pattern.Append(smock1);
            Smock_pattern.Append(smock2);

            smock_stitch_size += 2;
        }
    }
    std::cout << "The size of the smock stitch is: " << smock_stitch_size << std::endl;
}

template <class T, int dim = 3>
void offset_stitching(T offset, MESH_NODE<T, dim>& X, MESH_ELEM<dim - 1>& graph_Elem){
    // std::cout << "The size of all stitch is: " << stitchNodes.size() << std::endl;
    graph_Elem.Each([&](int id, auto data) {
        auto &[elemVInd] = data;
        VECTOR<T, dim>& X1 = std::get<0>(X.Get_Unchecked(elemVInd[0]));
        VECTOR<T, dim>& X2 = std::get<0>(X.Get_Unchecked(elemVInd[1]));
        VECTOR<T, dim>& X3 = std::get<0>(X.Get_Unchecked(elemVInd[2]));

        if (X1[2] > -0.0001)
            X1[2] -= offset;

        if (X2[2] > -0.0001)
            X2[2] -= offset;

        if (X3[2] > -0.0001)
            X3[2] -= offset;

    });
}

template <class T, int dim = 3>
void populate_coarse_graph(MESH_NODE<T, dim>& X_0, MESH_NODE<T, dim>& X_load, MESH_ELEM<dim - 1>& graph_Elem,
    int start_idx, int end_idx, 
    int stage_start, int stage_size, 
    T scale, VECTOR<T, 3>& x_axis, VECTOR<T, 3>& y_axis, std::vector<VECTOR<int, 4>>& edgeStencil,
    std::vector<VECTOR<T, 3>>& edgeInfo, T alignmult)
{
    x_axis = x_axis.Normalized();
    y_axis = y_axis.Normalized();
    const VECTOR<T, dim>& root = std::get<0>(X_0.Get_Unchecked(start_idx)); // pop root
    const VECTOR<T, dim>& terminal = std::get<0>(X_0.Get_Unchecked(end_idx)); // pop terminal
    int x_times = std::abs((terminal - root).dot(x_axis) / scale) + 1;
    int y_times = std::abs((terminal - root).dot(y_axis) / scale) + 1;
    int x_units = std::abs((terminal - root).dot(x_axis) / (3.0 * scale));
    int y_units = std::abs((terminal - root).dot(y_axis) / (2.0 * scale));

    std::vector<int> graph_nodes(x_times * y_times);

    for(int i  = 0; i < x_times; i++){
        for(int j = 0; j < y_times; j++){

            VECTOR<T, dim> target = root + (i * scale) * x_axis - (j * scale) * y_axis; //VECTOR<T, dim>(i*scale, -j*scale, 0.0);
            double dist = 100.0;
            int found_node;
            
            for(int k = 0; k < stage_size; k++){
                int cur_idx = stage_start + k;
                const VECTOR<T, dim>& cur_node = std::get<0>(X_0.Get_Unchecked(cur_idx)); 
                double cur_dist = (target - cur_node).length();

                if(cur_dist < dist){
                    dist = cur_dist;
                    found_node = cur_idx;
                }  
            }
            graph_nodes[j*x_times + i] = found_node; // row-major storage
        }
    }

    std::vector<int> underlay_idx;

    for(int j = 0; j < y_times-1; j++){
        for(int i = 0; i < x_times-1; i++){
            // TODO: Delete pleat elems!!!
            VECTOR<int, 3> tri;
            
            if(i%3 == 1 && j%2 == 0 && i/3 < x_units && j/2 < y_units){
                tri(0) = graph_nodes[j*x_times + i];
                tri(1) = graph_nodes[j*x_times + i + 1];
                tri(2) = graph_nodes[j*x_times + i + x_times];
                graph_Elem.Append(tri);
                underlay_idx.push_back(graph_Elem.size - 1);
            }
            
            else if(i%3 == 0 && j%2 == 1 && i/3 < x_units && j/2 < y_units){
                tri(0) = graph_nodes[j*x_times + i];
                tri(1) = graph_nodes[j*x_times + i + 1];
                tri(2) = graph_nodes[j*x_times + i + x_times + 1];
                graph_Elem.Append(tri);
                underlay_idx.push_back(graph_Elem.size - 1);
            }
        }
    }
    std::cout << "The size of the graph elements are: " << graph_Elem.size << std::endl;
    add_underlay_edges(X_load, graph_Elem, underlay_idx, edgeStencil, edgeInfo, x_units, alignmult);
    
}


void Export_Populate(py::module& m) {
    py::module smock_m = m.def_submodule("smock", "A submodule of JGSL for FEM smock simulation");

    smock_m.def("Populate_pattern", &populate_pattern<double>);
    smock_m.def("Populate_pattern_hinge", &populate_pattern_hinge<double>);
    smock_m.def("Offset_stitching", &offset_stitching<double>);
    smock_m.def("Populate_coarse_graph", &populate_coarse_graph<double>);
}



}