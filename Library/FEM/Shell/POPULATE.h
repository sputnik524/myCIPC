#pragma once

// for populate a certain smocking pattern on a garment stage

namespace py = pybind11;
namespace JGSL {

template<class T, int dim = 3>
void populate_pattern(MESH_NODE<T, dim>& X_0, std::vector<VECTOR<int, 3>>& smock_stitch, std::vector<T>& stitchRatio,
    int start_idx, int end_idx, 
    int stage_start, int stage_size, 
    T scale, int& smock_stitch_size)
{
    // determine the boundary
    const VECTOR<T, dim>& root = std::get<0>(X_0.Get_Unchecked(start_idx)); // pop root
    std::cout << "The root pos is: " << root[0] << " " << root[1] << std::endl;
    const VECTOR<T, dim>& terminal = std::get<0>(X_0.Get_Unchecked(end_idx)); // pop terminal
    std::cout << "The terminal pos is: " << terminal[0] << " " << terminal[1] << std::endl;
    int x_times = std::abs((terminal - root)[0] / (3.0 * scale));
    int y_times = std::abs((terminal - root)[1] / (2.0 * scale));
    std::cout << "The repeat info is: " << x_times << " " << y_times << std::endl;
    
    smock_stitch_size = 0;
    // braid pattern
    for(int i  = 0; i < x_times; i++){
        for(int j = 0; j < y_times; j++){
            // four smocking nodes for an unit braid pattern
            VECTOR<T, dim> target1 = root + VECTOR<T, dim>(scale + i * 3 * scale, - j * 2 * scale, 0.0);
            VECTOR<T, dim> target2 = root + VECTOR<T, dim>(i * 3 * scale, -(scale + j * 2 * scale), 0.0);

            VECTOR<T, dim> target3 = root + VECTOR<T, dim>(scale + i * 3 * scale, -(scale + j * 2 * scale), 0.0);
            VECTOR<T, dim> target4 = root + VECTOR<T, dim>(scale * 2 + i * 3 * scale, -(scale * 2 + j * 2 * scale), 0.0);

            double dist1 = 100.0;
            double dist2 = 100.0;
            double dist3 = 100.0;
            double dist4 = 100.0;

            int found_node1, found_node2, found_node3, found_node4;
            VECTOR<int, 3> stitch1, stitch2;

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

            // std::cout << "Distance1: " << dist1 << "with found node:" << found_node1 << std::endl;

            stitch1[0] = found_node1;
            stitch1[1] = found_node2;
            stitch1[2] = found_node2;
            
            stitch2[0] = found_node3;
            stitch2[1] = found_node4;
            stitch2[2] = found_node4;

            smock_stitch.push_back(stitch1);
            smock_stitch.push_back(stitch2);

            stitchRatio.push_back(1.0);
            stitchRatio.push_back(1.0);

            smock_stitch_size += 2;
        }
    }
    std::cout << "The size of the smock stitch is: " << smock_stitch_size << std::endl;
}

template <class T, int dim = 3>
void offset_stitching(T offset, MESH_NODE<T, dim>& X, std::vector<VECTOR<int, 3>>& stitchNodes, const int& smock_stitch_size){
    std::cout << "The size of all stitch is: " << stitchNodes.size() << std::endl;
    int start = stitchNodes.size() - smock_stitch_size;
    for(int i = start; i < stitchNodes.size(); i++){
        int node_0_idx = stitchNodes[i][0];
        int node_1_idx = stitchNodes[i][1];
        // std::cout << "The offsetting node is: " << node_0_idx << "" << node_1_idx << std::endl;

        VECTOR<T, dim>& X1_0 = std::get<0>(X.Get_Unchecked(node_0_idx));
        VECTOR<T, dim>& X2_0 = std::get<0>(X.Get_Unchecked(node_1_idx));

        X1_0[2] += offset;
        X2_0[2] += offset;
    }
}


void Export_Populate(py::module& m) {
    py::module smock_m = m.def_submodule("smock", "A submodule of JGSL for FEM smock simulation");

    smock_m.def("Populate_pattern", &populate_pattern<double>);
    smock_m.def("Offset_stitching", &offset_stitching<double>);
}



}