#pragma once

// for populate a certain smocking pattern on a garment stage

namespace py = pybind11;
namespace JGSL {

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
void offset_stitching(T offset, MESH_NODE<T, dim>& X, std::vector<VECTOR<int, 3>>& stitchNodes, const int& smock_stitch_size){
    // std::cout << "The size of all stitch is: " << stitchNodes.size() << std::endl;
    int start = stitchNodes.size() - smock_stitch_size;
    for(int i = start; i < stitchNodes.size(); i++){
        int node_0_idx = stitchNodes[i][0];
        int node_1_idx = stitchNodes[i][1];
        // std::cout << "The offsetting node is: " << node_0_idx << "" << node_1_idx << std::endl;

        VECTOR<T, dim>& X1_0 = std::get<0>(X.Get_Unchecked(node_0_idx));
        VECTOR<T, dim>& X2_0 = std::get<0>(X.Get_Unchecked(node_1_idx));

        if(X1_0[2] > 0.0){
            X1_0[2] += offset;
            X2_0[2] += offset;
        }

        else{
            X1_0[2] -= offset;
            X2_0[2] -= offset;
        }
        
    }
}

template <class T, int dim = 3>
void populate_coarse_graph(MESH_NODE<T, dim>& X_0, MESH_ELEM<dim - 1>& graph_Elem,
    int start_idx, int end_idx, 
    int stage_start, int stage_size, 
    T scale, VECTOR<T, 3>& x_axis, VECTOR<T, 3>& y_axis)
{
    x_axis = x_axis.Normalized();
    y_axis = y_axis.Normalized();
    const VECTOR<T, dim>& root = std::get<0>(X_0.Get_Unchecked(start_idx)); // pop root
    const VECTOR<T, dim>& terminal = std::get<0>(X_0.Get_Unchecked(end_idx)); // pop terminal
    int x_times = std::abs((terminal - root).dot(x_axis) / scale) + 1;
    int y_times = std::abs((terminal - root).dot(y_axis) / scale) + 1;

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
            graph_nodes[j*x_times + i] = found_node;
        }
    }

    for(int i = 0; i < x_times-1; i++){
        for(int j = 0; j < y_times-1; j++){
            VECTOR<int, 3> tri0,tri1,tri2,tri3;
            tri0(0) = graph_nodes[j*x_times + i];
            tri0(1) = graph_nodes[j*x_times + i + 1];
            tri0(2) = graph_nodes[j*x_times + i + x_times];
            
            tri1(0) = graph_nodes[j*x_times + i + 1];
            tri1(1) = graph_nodes[j*x_times + i + x_times];
            tri1(2) = graph_nodes[j*x_times + i + x_times + 1];

            tri2(0) = graph_nodes[j*x_times + i];
            tri2(1) = graph_nodes[j*x_times + i + 1];
            tri2(2) = graph_nodes[j*x_times + i + x_times + 1];

            tri3(0) = graph_nodes[j*x_times + i] ;
            tri3(1) = graph_nodes[j*x_times + i + x_times] ;
            tri3(2) = graph_nodes[j*x_times + i + x_times + 1] ;

            graph_Elem.Append(tri0);
            graph_Elem.Append(tri1);
            graph_Elem.Append(tri2);
            graph_Elem.Append(tri3);
        }
    }

    std::cout << "The size of the graph elements are: " << graph_Elem.size << std::endl;
}


void Export_Populate(py::module& m) {
    py::module smock_m = m.def_submodule("smock", "A submodule of JGSL for FEM smock simulation");

    smock_m.def("Populate_pattern", &populate_pattern<double>);
    smock_m.def("Offset_stitching", &offset_stitching<double>);
    smock_m.def("Populate_coarse_graph", &populate_coarse_graph<double>);
}



}