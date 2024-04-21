#pragma once

// #include <Physics/FIXED_COROTATED.h>

namespace JGSL {

template <int dim = 3>
void graph_construct(MESH_ELEM<dim - 1>& graph_Elem, int mesh_size = 13){
    for(int i = 0; i < mesh_size-1; i++){
        for(int j = 0; j < mesh_size-1; j++){
            VECTOR<int, 3> tri0,tri1,tri2,tri3;
            tri0(0) = i*mesh_size + j;
            tri0(1) = i*mesh_size + j + 1;
            tri0(2) = i*mesh_size + j + mesh_size;
            
            tri1(0) = i*mesh_size + j + 1;
            tri1(1) = i*mesh_size + j + mesh_size;
            tri1(2) = i*mesh_size + j + mesh_size + 1;

            tri2(0) = i*mesh_size + j;
            tri2(1) = i*mesh_size + j + 1;
            tri2(2) = i*mesh_size + j + mesh_size + 1;

            tri3(0) = i*mesh_size + j;
            tri3(1) = i*mesh_size + j + mesh_size;
            tri3(2) = i*mesh_size + j + mesh_size + 1;

            graph_Elem.Append(tri0);
            graph_Elem.Append(tri1);
            graph_Elem.Append(tri2);
            graph_Elem.Append(tri3);
        }
    }

    std::cout << "The size of the graph elements are: " << graph_Elem.size << std::endl; 
}

template<class T, int dim = 3>
double compute_edge_length(int v1, int v2, int pleat_1, int pleat_2, MESH_NODE<T, dim>& X, MESH_ELEM<dim - 1>& Smock_pattern){
    if(pleat_1 < 0 && pleat_2 < 0){
        const VECTOR<T, dim>& X1 = std::get<0>(X.Get_Unchecked(v1));
        const VECTOR<T, dim>& X2 = std::get<0>(X.Get_Unchecked(v2));
        return (X2-X1).length();
    }
    
    else if(pleat_1 >= 0 && pleat_2 >= 0){
        std::cout << "Degenerated edge:" << pleat_1 << " " << pleat_2 << std::endl;  
        if (pleat_1 == pleat_2)
            return 0.0;
        else{
            double shortest = 100.0;
            //return the closest length btw pair of smocking lines
            // TODO: consider smock pattern with length of 2, which means -1 for the third elem, now the for loop only iter over 2
            const VECTOR<int, dim>& L1 = std::get<0>(Smock_pattern.Get_Unchecked(pleat_1));
            const VECTOR<int, dim>& L2 = std::get<0>(Smock_pattern.Get_Unchecked(pleat_2));
            for(int i = 0; i < 2; i++){     
                const VECTOR<T, dim>& X1 = std::get<0>(X.Get_Unchecked(L1[i]));
                for(int j = 0; j < 2; j++){
                    const VECTOR<T, dim>& X2 = std::get<0>(X.Get_Unchecked(L2[j]));
                    double length_ = (X2-X1).length(); 
                    if(length_ < shortest) 
                        shortest = length_; 
                }
            }
            return shortest;
        }
    }

    else{
        double shortest = 100.0;   
        if(pleat_1 >= 0){
            const VECTOR<int, dim>& L1 = std::get<0>(Smock_pattern.Get_Unchecked(pleat_1));
            const VECTOR<T, dim>& X2 = std::get<0>(X.Get_Unchecked(v2));
            for(int i = 0; i < 2; i++){
                // compute the closest length btw node and edge
                const VECTOR<T, dim>& X1 = std::get<0>(X.Get_Unchecked(L1[i]));
                double length_ = (X2-X1).length(); 
                    if(length_ < shortest) 
                        shortest = length_; 
            }
        }

        else{
            const VECTOR<int, dim>& L2 = std::get<0>(Smock_pattern.Get_Unchecked(pleat_2));
            const VECTOR<T, dim>& X1 = std::get<0>(X.Get_Unchecked(v1));
            for(int i = 0; i < 2; i++){
                // compute the closest length btw node and edge
                const VECTOR<T, dim>& X2 = std::get<0>(X.Get_Unchecked(L2[i]));
                double length_ = (X2-X1).length(); 
                    if(length_ < shortest) 
                        shortest = length_; 
            }
        }
        return shortest;
    }

    return 1.0;
}

template <class T, int dim = 3>
MATRIX<T, dim - 1> compute_ref(const VECTOR<int, 3>& tri_elem, MESH_NODE<T, dim>& X, MESH_ELEM<dim - 1>& Smock_pattern, double scale = 0.02){
    
    std::cout << "Current graph elem: " << tri_elem[0] << " " << tri_elem[1] << " " << tri_elem[2] << std::endl;
    
    MATRIX<T, dim - 1> IB;
    VECTOR<int, 3> tri_elem_ref;
    VECTOR<T, 3> tri_elem_ref_length;
    tri_elem_ref[0] = -1; // negtive for pleat nodes
    tri_elem_ref[1] = -1;
    tri_elem_ref[2] = -1;

    Smock_pattern.Each([&](int id, auto data) {
        auto &[elemVInd] = data;
        for(int i = 0; i < 3; i++){     
            for(int j = 0; j < 3; j++){
                if(tri_elem[i] == elemVInd[j] && elemVInd[j] >= 0)
                    tri_elem_ref[i] = id; // record the idx of the smocking line the current node belongs to
            }
        }
    });

    // three edges
    double length01,length12,length20;
    length01 = compute_edge_length(tri_elem[0], tri_elem[1], tri_elem_ref[0], tri_elem_ref[1],X,Smock_pattern);
    length12 = compute_edge_length(tri_elem[1], tri_elem[2], tri_elem_ref[1], tri_elem_ref[2],X,Smock_pattern);
    length20 = compute_edge_length(tri_elem[2], tri_elem[0], tri_elem_ref[2], tri_elem_ref[0],X,Smock_pattern);

    std::cout << "Current edge length: " << length01 << " " << length12 << " " << length20 << std::endl;

    if(length01 < INT_MIN || length12 < INT_MIN || length20 < INT_MIN)
    {
        //degenerated edge
        IB(0,0) = IB(1,1) = 1.0;
        IB(1,0) = IB(0,1) = 1.0;
    }

    else
    {
        IB(0,0) = length01 * length01;
        IB(1,1) = length20 * length20;
        IB(0,1) = IB(1,0) = 0.5 * (length01 * length01 - length12 * length12 + length20 * length20);// given a,b,c, want ab*cos(C),law of cosine
    }

    return IB;
}




}