#pragma once

namespace JGSL {

template <class T, int dim = 3>
void add_underlay_edges(MESH_NODE<T, dim>& X, MESH_ELEM<dim - 1>& graph_Elem, const std::vector<int>& underlay_idx,
    std::vector<VECTOR<int ,4>>& edgeStencil, std::vector<VECTOR<T, 3>>& edgeInfo, int x_units, T alignmult)
{   
    std::cout << "Underlay size: " << underlay_idx.size() << std::endl;
    for(int i = 0; i < underlay_idx.size() - x_units; i++){
        // int graph_offset = graph_Elem.size - underlay_idx.size();
        // std::cout << "Current underlay idx: " << underlay_idx[i] << std::endl;
        const VECTOR<int, 3>& underlay_tri_1 = std::get<0>(graph_Elem.Get_Unchecked(underlay_idx[i]));
        const VECTOR<int, 3>& underlay_tri_2 = std::get<0>(graph_Elem.Get_Unchecked(underlay_idx[i + x_units]));
        
        int v3_0, v3_1;

        if((i/x_units)%2 == 0){ 
            v3_0 = underlay_tri_2[2];
            v3_1 = underlay_tri_1[1];
            edgeStencil.emplace_back(VECTOR<int, 4>(underlay_tri_1[1], underlay_tri_1[2], underlay_tri_1[0], v3_0));
            edgeStencil.emplace_back(VECTOR<int, 4>(underlay_tri_2[2], underlay_tri_2[0], underlay_tri_2[1], v3_1));
        }

        else if((i/x_units)%2 == 1){
            v3_0 = underlay_tri_2[2];
            v3_1 = underlay_tri_1[0];
            edgeStencil.emplace_back(VECTOR<int, 4>(underlay_tri_1[0], underlay_tri_1[1], underlay_tri_1[2], v3_0));
            edgeStencil.emplace_back(VECTOR<int, 4>(underlay_tri_2[2], underlay_tri_2[0], underlay_tri_2[1], v3_1));
        }

        // std::cout << "EdgeStencil size: " << edgeStencil.size() << std::endl;

        for(int i = 0; i < 2; i++){
            VECTOR<int, 4> cur_edgeStencil = edgeStencil[edgeStencil.size() - 2 + i];
            const VECTOR<T, dim>& X0 = std::get<0>(X.Get_Unchecked(cur_edgeStencil[0]));
            const VECTOR<T, dim>& X1 = std::get<0>(X.Get_Unchecked(cur_edgeStencil[1]));
            const VECTOR<T, dim>& X2 = std::get<0>(X.Get_Unchecked(cur_edgeStencil[2]));
            const VECTOR<T, dim>& X3 = std::get<0>(X.Get_Unchecked(cur_edgeStencil[3]));

            Eigen::Matrix<T, dim, 1> X0e(X0.data), X1e(X1.data), X2e(X2.data), X3e(X3.data);

            edgeInfo.resize(edgeInfo.size() + 1);
            Compute_Dihedral_Angle(X0e, X1e, X2e, X3e, edgeInfo.back()[0]);
            std::cout << "Cur dihedral angle: " << edgeInfo.back()[0] << std::endl;
            edgeInfo.back()[0] = 0.0;
            edgeInfo.back()[1] = (X1 - X2).length();
            VECTOR<T, 3> n1 = cross(X1 - X0, X2 - X0);
            VECTOR<T, 3> n2 = cross(X2 - X3, X1 - X3);
            edgeInfo.back()[2] = (n1.length() + n2.length()) / (edgeInfo.back()[1] * 6); 
            edgeInfo.back()[1] *= alignmult;
        }
    }
    std::cout << edgeStencil.size() << "Added smocking hinges" << std::endl;
}

}