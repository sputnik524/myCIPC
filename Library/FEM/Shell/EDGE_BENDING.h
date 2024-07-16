#pragma once

#include <Math/DIHEDRAL_ANGLE.h>
#include <FEM/Shell/Rod/ROD_BENDING_DERIVATIVES.h>
#include <FEM/Shell/UTILS.h>

namespace JGSL {

template<class T, int dim = 3>
void Compute_Edge_Bending_Energy(
    T h,
    const std::vector<VECTOR<int, 3>>& rodHinge,
    const std::vector<VECTOR<T, 3>>& rodHingeInfo,
    MESH_NODE<T, dim>& X, // mid-surface node coordinates
    T& E)
{
    TIMER_FLAG("Compute_Edge_Bending_Energy");

    int i = 0;
    for (const auto& hingeI : rodHinge) {
        const VECTOR<T, dim>& x0 = std::get<0>(X.Get_Unchecked(hingeI[0]));
        const VECTOR<T, dim>& x1 = std::get<0>(X.Get_Unchecked(hingeI[1]));
        const VECTOR<T, dim>& x2 = std::get<0>(X.Get_Unchecked(hingeI[2]));

        const VECTOR<T, dim> e0 = x1 - x0;
        const VECTOR<T, dim> x2_ref = 2.0 * x1 - x2;
        const VECTOR<T, dim> e1 = x2_ref - x1;
        // const VECTOR<T, dim> e1 = x2 - x1;

        const VECTOR<T, dim> kappa = cross(e0, e1) * 2 / (sqrt(e0.length2() * e1.length2()) + e0.dot(e1));
        const T alpha = rodHingeInfo[i][0] * std::pow(rodHingeInfo[i][2], 4) * M_PI / 64;
        E += h * h * alpha * kappa.length2() / rodHingeInfo[i][1];

        ++i;
    }

}

template<class T, int dim = 3>
void Compute_Edge_Bending_Gradient(
    T h,
    const std::vector<VECTOR<int, 3>>& rodHinge,
    const std::vector<VECTOR<T, 3>>& rodHingeInfo,
    MESH_NODE<T, dim>& X, // mid-surface node coordinates
    MESH_NODE_ATTR<T, dim>& nodeAttr)
{
    TIMER_FLAG("Compute_Edge_Bending_Gradient");
    int i = 0;
    for (const auto& hingeI : rodHinge) {
        const VECTOR<T, dim>& x0 = std::get<0>(X.Get_Unchecked(hingeI[0]));
        const VECTOR<T, dim>& x1 = std::get<0>(X.Get_Unchecked(hingeI[1]));
        const VECTOR<T, dim>& x2 = std::get<0>(X.Get_Unchecked(hingeI[2]));

        const VECTOR<T, dim> x2_ref = 2.0 * x1 - x2;

        T g_kappa2[9];
        g_RB(x0[0], x0[1], x0[2], x1[0], x1[1], x1[2], x2_ref[0], x2_ref[1], x2_ref[2], g_kappa2);
        // g_RB(x0[0], x0[1], x0[2], x1[0], x1[1], x1[2], x2[0], x2[1], x2[2], g_kappa2);


        for (int j = 0; j < dim; ++j)
        {
            g_kappa2[3+j] += g_kappa2[6+j] * 2.0;
            g_kappa2[6+j] *= -1.0;
        }

        if (hingeI[0] == 1956){
            std::cout << "First closing pair pos: " << x0[0]<< " " << x0[1]<< " " << x0[2] << std::endl;
            std::cout  << x1[0]<< " " << x1[1]<< " " << x1[2] << std::endl;
            std::cout  << x2[0]<< " " << x2[1]<< " " << x2[2] << std::endl;
            std::cout << "g_kappa2 after transform: " << g_kappa2[0]<< " " << g_kappa2[1]<< " " << g_kappa2[2] << std::endl;
            std::cout  << g_kappa2[3]<< " " << g_kappa2[4]<< " " << g_kappa2[5] << std::endl;
            std::cout  << g_kappa2[6]<< " " << g_kappa2[7]<< " " << g_kappa2[8] << std::endl;
        }

        const T w = h * h / rodHingeInfo[i][1] *
            rodHingeInfo[i][0] * std::pow(rodHingeInfo[i][2], 4) * M_PI / 64;
        for (int endI = 0; endI < 3; ++endI) {
            VECTOR<T, dim>& grad = std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(nodeAttr.Get_Unchecked(hingeI[endI]));
            for (int dimI = 0; dimI < dim; ++dimI) {
                grad[dimI] += w * g_kappa2[endI * dim + dimI];
            }
        }

        ++i;
    }
}

template<class T, int dim = 3>
void Compute_Edge_Bending_Hessian(
    T h, bool projectSPD,
    const std::vector<VECTOR<int, 3>>& rodHinge,
    const std::vector<VECTOR<T, 3>>& rodHingeInfo,
    MESH_NODE<T, dim>& X,
    std::vector<Eigen::Triplet<T>>& triplets)
{
    TIMER_FLAG("Compute_Edge_Bending_Hessian");

    BASE_STORAGE<int> threads(rodHinge.size());
    for (int i = 0; i < rodHinge.size(); ++i) {
        threads.Append(triplets.size() + i * 81);
    }

    triplets.resize(triplets.size() + rodHinge.size() * 81);
    threads.Par_Each([&](int i, auto data) {
        const auto& [tripletStartInd] = data;
        const VECTOR<int, 3>& hingeI = rodHinge[i];
        const VECTOR<T, dim>& x0 = std::get<0>(X.Get_Unchecked(hingeI[0]));
        const VECTOR<T, dim>& x1 = std::get<0>(X.Get_Unchecked(hingeI[1]));
        const VECTOR<T, dim>& x2 = std::get<0>(X.Get_Unchecked(hingeI[2]));

        const VECTOR<T, dim> x2_ref = 2.0 * x1 - x2;

        Eigen::Matrix<T, 9, 9> hessian;
        H_RB(x0[0], x0[1], x0[2], x1[0], x1[1], x1[2], x2_ref[0], x2_ref[1], x2_ref[2], hessian.data());
        // H_RB(x0[0], x0[1], x0[2], x1[0], x1[1], x1[2], x2[0], x2[1], x2[2], hessian.data());

        hessian.template block<3,3>(3,0) += 2.0 * (hessian.template block<3,3>(6,0));
        hessian.template block<3,3>(0,3) += 2.0 * (hessian.template block<3,3>(0,6));

        hessian.template block<3,3>(3,3) += 2.0 * (hessian.template block<3,3>(3,6)) + 2.0 * (hessian.template block<3,3>(6,3)) + 4.0 * (hessian.template block<3,3>(6,6));

        hessian.template block<3,3>(6,3) *= -1.0;
        hessian.template block<3,3>(6,3) -= 2.0 * hessian.template block<3,3>(6,6);
        hessian.template block<3,3>(3,6) *= -1.0;
        hessian.template block<3,3>(3,6) -= 2.0 * hessian.template block<3,3>(6,6);

        hessian.template block<3,3>(6,0) *= -1.0;
        hessian.template block<3,3>(0,6) *= -1.0;

        if (projectSPD) {
            makePD(hessian);
        }

        int globalInd[9] = { 
            hingeI[0] * dim,
            hingeI[0] * dim + 1,
            hingeI[0] * dim + 2,
            hingeI[1] * dim,
            hingeI[1] * dim + 1,
            hingeI[1] * dim + 2,
            hingeI[2] * dim,
            hingeI[2] * dim + 1,
            hingeI[2] * dim + 2,
        };
        const T w = h * h / rodHingeInfo[i][1] *
            rodHingeInfo[i][0] * std::pow(rodHingeInfo[i][2], 4) * M_PI / 64;
        for (int rowI = 0; rowI < 9; ++rowI) {
            for (int colI = 0; colI < 9; ++colI) {
                triplets[tripletStartInd + rowI * 9 + colI] = Eigen::Triplet<T>(
                    globalInd[rowI], globalInd[colI], w * hessian(rowI, colI)
                );
            }
        }
    });
}

template <class T, int dim>
void Check_Edge_Bending_Gradient(
    T h,
    const std::vector<VECTOR<int, 3>>& rodHinge,
    const std::vector<VECTOR<T, 3>>& rodHingeInfo,
    MESH_NODE<T, dim>& X, // mid-surface node coordinates
    MESH_NODE_ATTR<T, dim>& nodeAttr)
{
    T eps = 1.0e-6;

    T E0 = 0;
    Compute_Edge_Bending_Energy(h, rodHinge, rodHingeInfo, X, E0);
    nodeAttr.template Fill<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(VECTOR<T, dim>(0));
    Compute_Edge_Bending_Gradient(h, rodHinge, rodHingeInfo, X, nodeAttr);

    std::vector<T> grad_FD(X.size * dim);
    for (int i = 0; i < X.size * dim; ++i) {
        MESH_NODE<T, dim> Xperturb;
        Append_Attribute(X, Xperturb);
        std::get<0>(Xperturb.Get_Unchecked(i / dim))[i % dim] += eps;
        
        T E = 0;
        Compute_Edge_Bending_Energy(h, rodHinge, rodHingeInfo, Xperturb, E);
        grad_FD[i] = (E - E0) / eps;
    }

    T err = 0.0, norm = 0.0;
    int error_num = 0;
    nodeAttr.Each([&](int id, auto data) {
        auto &[x0, v, g, m] = data;

        err += std::pow(grad_FD[id * dim] - g[0], 2);
        err += std::pow(grad_FD[id * dim + 1] - g[1], 2);

        norm += std::pow(grad_FD[id * dim], 2);
        norm += std::pow(grad_FD[id * dim + 1], 2);

        if constexpr (dim == 3) {
            err += std::pow(grad_FD[id * dim + 2] - g[2], 2);
            norm += std::pow(grad_FD[id * dim + 2], 2);
        }
        
        if(std::pow(grad_FD[id * dim + 1] - g[1], 2) > 0.01){
            // T v0_err = grad_FD[id * dim] - g[0];
            // T v1_err = grad_FD[id * dim + 1] - g[1];
            // T v2_err = grad_FD[id * dim + 2] - g[2];
            
            // const VECTOR<T, dim>& x0 = std::get<0>(X.Get_Unchecked(id));
            // printf("Cur error idx: %10d \n", id);
            // printf("Cur node pos: x0 = %le, x1 = %le, x2 = %le\n", x0[0], x0[1],x0[2]);
            // printf("v0_err = %le, v1_err = %le, v2_err = %le\n", v0_err, v1_err, v2_err);
            // printf("v0_FD = %le, v1_FD = %le, v2_FD = %le\n", grad_FD[id * dim], grad_FD[id * dim + 1], grad_FD[id * dim + 2]);
            // printf("g0 = %le, g1 = %le, g2 = %le\n", g[0], g[1], g[2]);
            error_num ++;
        }
        
    });
    std::cout << error_num << std::endl;
    printf("err_abs = %le, sqnorm_FD = %le, err_rel = %le\n", err, norm, err / norm);
}

template <class T, int dim>
void Check_Edge_Bending_Hessian(
    T h,
    const std::vector<VECTOR<int, 3>>& rodHinge,
    const std::vector<VECTOR<T, 3>>& rodHingeInfo,
    MESH_NODE<T, dim>& X, // mid-surface node coordinates
    MESH_NODE_ATTR<T, dim>& nodeAttr)
{
    T eps = 1.0e-6;

    MESH_NODE_ATTR<T, dim> nodeAttr0;
    nodeAttr.deep_copy_to(nodeAttr0);
    nodeAttr0.template Fill<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(VECTOR<T, dim>(0));
    Compute_Edge_Bending_Gradient(h, rodHinge, rodHingeInfo, X, nodeAttr0);
    std::vector<Eigen::Triplet<T>> HStriplets;
    Compute_Edge_Bending_Hessian(h, false, rodHinge, rodHingeInfo, X, HStriplets);
    CSR_MATRIX<T> HS;
    HS.Construct_From_Triplet(X.size * dim, X.size * dim, HStriplets);

    std::vector<Eigen::Triplet<T>> HFDtriplets;
    HFDtriplets.reserve(HStriplets.size());
    for (int i = 0; i < X.size * dim; ++i) {
        MESH_NODE<T, dim> Xperturb;
        Append_Attribute(X, Xperturb);
        std::get<0>(Xperturb.Get_Unchecked(i / dim))[i % dim] += eps;
        
        nodeAttr.template Fill<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(VECTOR<T, dim>(0));
        Compute_Edge_Bending_Gradient(h, rodHinge, rodHingeInfo, Xperturb, nodeAttr);
        for (int vI = 0; vI < X.size; ++vI) {
            const VECTOR<T, dim>& g = std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(nodeAttr.Get_Unchecked(vI));
            const VECTOR<T, dim>& g0 = std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(nodeAttr0.Get_Unchecked(vI));
            const VECTOR<T, dim> hFD = (g - g0) / eps;
            if (hFD.length2() != 0) {
                HFDtriplets.emplace_back(i, vI * dim, hFD[0]);
                HFDtriplets.emplace_back(i, vI * dim + 1, hFD[1]);
                if constexpr (dim == 3) {
                    HFDtriplets.emplace_back(i, vI * dim + 2, hFD[2]);
                }
            }
        }
    }
    CSR_MATRIX<T> HFD;
    HFD.Construct_From_Triplet(X.size * dim, X.size * dim, HFDtriplets);

    T err = (HS.Get_Matrix() - HFD.Get_Matrix()).squaredNorm(), norm = HFD.Get_Matrix().squaredNorm();
    printf("err_abs = %le, sqnorm_FD = %le, err_rel = %le\n", err, norm, err / norm);
}

}