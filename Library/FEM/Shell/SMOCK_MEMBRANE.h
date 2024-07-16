#pragma once

#include <Physics/FIXED_COROTATED.h>

namespace JGSL {

template<class T>
void compute_pFpX(const Eigen::Matrix2d& Dm_iInv, Eigen::Matrix<T, 6, 9>& pF_ipX){
    T b00,b01,b10,b11,b00_10,b01_11;
    b00 = Dm_iInv(0,0);
    b01 = Dm_iInv(0,1);
    b10 = Dm_iInv(1,0);
    b11 = Dm_iInv(1,1);

    b00_10 = -(b00 + b10);
    b01_11 = -(b01 + b11);

    pF_ipX(0,0) = b00_10;
    pF_ipX(3,0) = b01_11;
    // vec(pF_ipX_0y)
    pF_ipX(1,1) = b00_10;
    pF_ipX(4,1) = b01_11;
    // vec(pF_ipX_0z)
    pF_ipX(2,2) = b00_10;
    pF_ipX(5,2) = b01_11;

    // vec(pF_ipX_1x)
    pF_ipX(0,3) = b00;
    pF_ipX(3,3) = b01;
    // vec(pF_ipX_1y)
    pF_ipX(1,4) = b00;
    pF_ipX(4,4) = b01;
    // vec(pF_ipX_1z)
    pF_ipX(2,5) = b00;
    pF_ipX(5,5) = b01;

    // vec(pF_ipX_2x)
    pF_ipX(0,6) = b10;
    pF_ipX(3,6) = b11;
    // vec(pF_ipX_2y)
    pF_ipX(1,7) = b10;
    pF_ipX(4,7) = b11;
    // vec(pF_ipX_2z)
    pF_ipX(2,8) = b10;
    pF_ipX(5,8) = b11;
}

template<class T, int dim, bool useNH = true, bool useARAP = false>
void Compute_Smock_Membrane_Energy(
    MESH_ELEM<dim - 1>& Elem, MESH_ELEM<dim - 1>& Elem_smock, T h,
    const std::vector<bool>& DBCb,
    MESH_NODE<T, dim>& X, // mid-surface node coordinates
    MESH_NODE_ATTR<T, dim>& nodeAttr,
    MESH_ELEM_ATTR<T, dim - 1>& elemAttr,
    MESH_ELEM_ATTR<T, dim - 1>& elemAttr_smock,
    FIXED_COROTATED<T, dim - 1>& elasticityAttr,
    FIXED_COROTATED<T, dim - 1>& elasticityAttr_smock,
    T& E)
{
    TIMER_FLAG("Compute_Membrane_Energy");
    if constexpr (dim == 2) {
        //TODO
    }
    else {
        //TODO: parallelize
        Elem.Join(elasticityAttr).Each([&](int id, auto data) {
            auto &[elemVInd, F, vol, lambda, mu] = data;
            if (!(DBCb[elemVInd[0]] && DBCb[elemVInd[1]] && DBCb[elemVInd[2]])) {
                const VECTOR<T, dim>& x1 = std::get<0>(X.Get_Unchecked(elemVInd[0]));
                const VECTOR<T, dim>& x2 = std::get<0>(X.Get_Unchecked(elemVInd[1]));
                const VECTOR<T, dim>& x3 = std::get<0>(X.Get_Unchecked(elemVInd[2]));
                MATRIX<T, dim - 1> IB = std::get<FIELDS<MESH_ELEM_ATTR<T, dim>>::IB>(elemAttr.Get_Unchecked(id));
                IB.invert();

                const VECTOR<T, dim> e01 = x2 - x1;
                const VECTOR<T, dim> e02 = x3 - x1;
                MATRIX<T, dim - 1> A;
                A(0, 0) = e01.length2();
                A(1, 0) = A(0, 1) = e01.dot(e02);
                A(1, 1) = e02.length2();

                if constexpr (useNH) {
                    const T lnJ = std::log(A.determinant() * IB.determinant()) / 2;
                    T e = h * h * vol * (mu / 2 * ((IB * A).trace() - 2 - 2 * lnJ) + lambda / 2 * lnJ * lnJ);
                    if (!std::isnan(e)){
                        // if(e < 0.0){
                        //     std::cout << "E negative with elem id: " << id << std::endl;
                        //     std::cout << "A det: " << A.determinant() << std::endl;
                        //     std::cout << "x1: " << x1[0] << " " << x1[1] << " " << x1[2] << std::endl;
                        //     std::cout << "x2: " << x2[0] << " " << x2[1] << " " << x2[2] << std::endl;
                        //     std::cout << "x3: " << x3[0] << " " << x3[1] << " " << x3[2] << std::endl;
                        //     std::cout << "A00: " << A(0, 0) << std::endl;
                        //     std::cout << "A10, A01: " << A(1, 0) << std::endl;
                        //     std::cout << "A11: " << A(1, 1) << std::endl;
                        //     std::cout << "IB det: " << IB.determinant() << std::endl;
                        // }
                        // else
                        E += e;
                        
                    }
                        
                    else{
                        std::cout << "Nan detected with elem id: " << id << std::endl;
                        std::cout << "A det: " << A.determinant() << std::endl;
                        std::cout << "x1: " << x1[0] << " " << x1[1] << " " << x1[2] << std::endl;
                        std::cout << "x2: " << x2[0] << " " << x2[1] << " " << x2[2] << std::endl;
                        std::cout << "x3: " << x3[0] << " " << x3[1] << " " << x3[2] << std::endl;
                        std::cout << "A00: " << A(0, 0) << std::endl;
                        std::cout << "A10, A01: " << A(1, 0) << std::endl;
                        std::cout << "A11: " << A(1, 1) << std::endl;
                        std::cout << "IB det: " << IB.determinant() << std::endl;
                    }
                }
                else {
                    MATRIX<T, dim - 1> M = IB * A; M(0, 0) -= 1; M(1, 1) -= 1;
                    E += h * h * vol / 4 * (0.5 * lambda * pow(M.trace(), 2) + mu * (M * M).trace());
                }
            }
        });

        T E_0 = E;

        Elem_smock.Join(elasticityAttr_smock).Each([&](int id, auto data) {
            auto &[elemVInd, F, vol, lambda, mu] = data;
            if (!(DBCb[elemVInd[0]] && DBCb[elemVInd[1]] && DBCb[elemVInd[2]])) {
                const VECTOR<T, dim>& x1 = std::get<0>(X.Get_Unchecked(elemVInd[0]));
                const VECTOR<T, dim>& x2 = std::get<0>(X.Get_Unchecked(elemVInd[1]));
                const VECTOR<T, dim>& x3 = std::get<0>(X.Get_Unchecked(elemVInd[2]));
                Eigen::Matrix<T, dim, 1> x1e(x1.data), x2e(x2.data), x3e(x3.data);
                MATRIX<T, dim - 1> IB = std::get<FIELDS<MESH_ELEM_ATTR<T, dim>>::IB>(elemAttr_smock.Get_Unchecked(id));
                IB.invert();

                if constexpr (useARAP) {
                    // right Cauchy-Green tensor C == F^T @ F == IB * A, 2 by 2 for membrane
                    // TODO: IB is now Dm_inv!
                    Eigen::Matrix2d Dm_iInv;
                    Dm_iInv(0,0) = IB(0,0);
                    Dm_iInv(0,1) = IB(0,1);
                    Dm_iInv(1,0) = IB(1,0);
                    Dm_iInv(1,1) = IB(1,1);

                    Eigen::Matrix<T, dim, 2> Ds, tri_F;
                    Ds << x2e - x1e, x3e - x1e;
                    tri_F = Ds * Dm_iInv;


                    T i1 = tri_F.squaredNorm();
                    T i3 = (tri_F.transpose() * tri_F).determinant();

                    T t1, t2;
                    t1 = -std::sqrt(i1 + 2 * std::sqrt(i3));
                    t2 = -t1;
                    
                    T e = h * h * vol * mu / 2 * (i1 - 2 * t2 + 2);
                    E += e;
                }

                else {
                    const VECTOR<T, dim> e01 = x2 - x1;
                    const VECTOR<T, dim> e02 = x3 - x1;
                    MATRIX<T, dim - 1> A;
                    A(0, 0) = e01.length2();
                    A(1, 0) = A(0, 1) = e01.dot(e02);
                    A(1, 1) = e02.length2();

                    if constexpr (useNH) {
                        const T lnJ = std::log(A.determinant() * IB.determinant()) / 2;
                        T e = h * h * vol * (mu / 2 * ((IB * A).trace() - 2 - 2 * lnJ) + lambda / 2 * lnJ * lnJ); 
                        if (!std::isnan(e))
                            E += e;
                            
                    }
                    else {
                        MATRIX<T, dim - 1> M = IB * A; M(0, 0) -= 1; M(1, 1) -= 1;
                        E += h * h * vol / 4 * (0.5 * lambda * pow(M.trace(), 2) + mu * (M * M).trace());
                    }
                }
            }
        });

        std::cout << "Underlay distortion energy change: " << E - E_0 << std::endl;

    }
}

template<class T, int dim, bool useNH = true, bool useARAP = false>
void Compute_Smock_Membrane_Gradient(
    MESH_ELEM<dim - 1>& Elem, MESH_ELEM<dim - 1>& Elem_smock, T h,
    const std::vector<bool>& DBCb,
    MESH_NODE<T, dim>& X, // mid-surface node coordinates
    MESH_NODE_ATTR<T, dim>& nodeAttr,
    MESH_ELEM_ATTR<T, dim - 1>& elemAttr,MESH_ELEM_ATTR<T, dim - 1>& elemAttr_smock,
    FIXED_COROTATED<T, dim - 1>& elasticityAttr, FIXED_COROTATED<T, dim - 1>& elasticityAttr_smock) 
{
    TIMER_FLAG("Compute_Membrane_Gradient");
    if constexpr (dim == 2) {
        //TODO
    }
    else {
        //TODO: parallelize
        Elem.Join(elasticityAttr).Each([&](int id, auto data) {
            auto &[elemVInd, F, vol, lambda, mu] = data;
            if (!(DBCb[elemVInd[0]] && DBCb[elemVInd[1]] && DBCb[elemVInd[2]])) {
                const VECTOR<T, dim>& x1 = std::get<0>(X.Get_Unchecked(elemVInd[0]));
                const VECTOR<T, dim>& x2 = std::get<0>(X.Get_Unchecked(elemVInd[1]));
                const VECTOR<T, dim>& x3 = std::get<0>(X.Get_Unchecked(elemVInd[2]));
                Eigen::Matrix<T, dim, 1> x1e(x1.data), x2e(x2.data), x3e(x3.data);
                MATRIX<T, dim - 1> IB = std::get<FIELDS<MESH_ELEM_ATTR<T, dim>>::IB>(elemAttr.Get_Unchecked(id));
                IB.invert();

                const VECTOR<T, dim> e01 = x2 - x1;
                const VECTOR<T, dim> e02 = x3 - x1;
                MATRIX<T, dim - 1> A;
                A(0, 0) = e01.length2();
                A(1, 0) = A(0, 1) = e01.dot(e02);
                A(1, 1) = e02.length2();

                MATRIX<T, dim - 1> temp;
                if constexpr (useNH) {
                    MATRIX<T, dim - 1> IA = A; IA.invert();
                    const T lnJ = std::log(A.determinant() * IB.determinant()) / 2;
                    temp = h * h * vol * (mu / 2 * IB + (-mu + lambda * lnJ) / 2 * IA);
                }
                else {
                    MATRIX<T, dim - 1> M = IB * A; M(0, 0) -= 1; M(1, 1) -= 1;
                    temp = h * h * vol / 4 * (lambda * M.trace() * IB + 2 * mu * M * IB);
                }

                Eigen::Matrix<T, 4, 9> dA_div_dx;
                dA_div_dx.setZero();
                dA_div_dx.template block<1, 3>(0, 3) += 2.0 * (x2e - x1e).transpose();
                dA_div_dx.template block<1, 3>(0, 0) -= 2.0 * (x2e - x1e).transpose();
                dA_div_dx.template block<1, 3>(1, 6) += (x2e - x1e).transpose();
                dA_div_dx.template block<1, 3>(1, 3) += (x3e - x1e).transpose();
                dA_div_dx.template block<1, 3>(1, 0) += -(x2e - x1e).transpose() - (x3e - x1e).transpose();
                dA_div_dx.template block<1, 3>(2, 6) += (x2e - x1e).transpose();
                dA_div_dx.template block<1, 3>(2, 3) += (x3e - x1e).transpose();
                dA_div_dx.template block<1, 3>(2, 0) += -(x2e - x1e).transpose() - (x3e - x1e).transpose();
                dA_div_dx.template block<1, 3>(3, 6) += 2.0 * (x3e - x1e).transpose();
                dA_div_dx.template block<1, 3>(3, 0) -= 2.0 * (x3e - x1e).transpose();

                for (int endI = 0; endI < dim; ++endI) {
                    VECTOR<T, dim>& g = std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(nodeAttr.Get_Unchecked(elemVInd[endI]));
                    for (int dimI = 0; dimI < dim; ++dimI) {
                        int i = endI * dim + dimI;
                        g[dimI] += dA_div_dx(0, i) * temp(0, 0) + dA_div_dx(1, i) * temp(1, 0) + 
                            dA_div_dx(2, i) * temp(0, 1) + dA_div_dx(3, i) * temp(1, 1);
                    }
                }
            }
        });

        int degenerated_smock_elems = 0;
        Elem_smock.Join(elasticityAttr_smock).Each([&](int id, auto data) {
            auto &[elemVInd, F, vol, lambda, mu] = data;
            if (!(DBCb[elemVInd[0]] && DBCb[elemVInd[1]] && DBCb[elemVInd[2]])) {
                const VECTOR<T, dim>& x1 = std::get<0>(X.Get_Unchecked(elemVInd[0]));
                const VECTOR<T, dim>& x2 = std::get<0>(X.Get_Unchecked(elemVInd[1]));
                const VECTOR<T, dim>& x3 = std::get<0>(X.Get_Unchecked(elemVInd[2]));
                Eigen::Matrix<T, dim, 1> x1e(x1.data), x2e(x2.data), x3e(x3.data);
                MATRIX<T, dim - 1> IB = std::get<FIELDS<MESH_ELEM_ATTR<T, dim>>::IB>(elemAttr_smock.Get_Unchecked(id));
                IB.invert();
                // TODO: Adding a Dm.invert() at init!

                if constexpr (useARAP) {
                    Eigen::Matrix<T, 6, 9> pF_ipX;
                    pF_ipX.setZero();
                    Eigen::Matrix2d Dm_iInv;
                    Dm_iInv(0,0) = IB(0,0);
                    Dm_iInv(0,1) = IB(0,1);
                    Dm_iInv(1,0) = IB(1,0);
                    Dm_iInv(1,1) = IB(1,1);

                    compute_pFpX(Dm_iInv, pF_ipX);

                    Eigen::Matrix<T, dim, 2> Ds, tri_F;
                    Ds << x2e - x1e, x3e - x1e;
                    tri_F = Ds * Dm_iInv;

                    Eigen::Matrix<T, 6, 1> g1, g2, g3;
                    g1.template block<3, 1>(0, 0) = 2 * tri_F.col(0);
                    g1.template block<3, 1>(3, 0) = 2 * tri_F.col(1);

                    Eigen::Matrix<T, dim, 2> mat_g2;
                    g2.template block<3, 1>(0, 0) = mat_g2.col(0);
                    g2.template block<3, 1>(3, 0) = mat_g2.col(1);

                    Eigen::Matrix2d C_F = tri_F.transpose() * tri_F;
                    T det_C = C_F.determinant(); 

                    T a,b,c,d;
                    a = C_F(0,0);
                    b = C_F(0,1);
                    c = C_F(1,0);
                    d = C_F(1,1);
                    g3(0,0) = 2.0 * d * tri_F(0,0) - (b+c) * tri_F(0,1);
                    g3(1,0) = 2.0 * d * tri_F(1,0) - (b+c) * tri_F(1,1);
                    g3(2,0) = 2.0 * d * tri_F(2,0) - (b+c) * tri_F(2,1);

                    g3(3,0) = 2.0 * a * tri_F(0,1) - (b+c) * tri_F(0,0);
                    g3(4,0) = 2.0 * a * tri_F(1,1) - (b+c) * tri_F(1,0);
                    g3(5,0) = 2.0 * a * tri_F(2,1) - (b+c) * tri_F(2,0);

                    T i1 = tri_F.squaredNorm();
                    T i2 = (tri_F.transpose() * tri_F).squaredNorm();
                    T i3 = det_C;
                    T t1, t2;
                    t1 = -sqrt(i1 + 2 * sqrt(i3));
                    t2 = -t1; // pos root

                    //according to eq.22 in [iARAP] supp.
                    T f1 = 0.5 / t2;
                    T f2 = 0.;
                    T f3 = 0.5 / (t2 * sqrt(i3));
                    Eigen::Matrix<T, 6, 1> pEpf = h * h * vol * 0.5 * ((1.0 - 2 * f1) * g1 - 2.0 * f3 * g3); // TODO: check 0.5 * 
                    Eigen::Matrix<T, 9, 1> pEpX = pF_ipX.transpose() * pEpf;

                    for (int endI = 0; endI < dim; ++endI) {
                        VECTOR<T, dim>& g = std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(nodeAttr.Get_Unchecked(elemVInd[endI]));
                        for (int dimI = 0; dimI < dim; ++dimI) {
                            int i = endI * dim + dimI;
                            g[dimI] += pEpX(0, i);// TODO: check g & idx mapping
                        }
                    }
                }

                else {
                    const VECTOR<T, dim> e01 = x2 - x1;
                    const VECTOR<T, dim> e02 = x3 - x1;
                    MATRIX<T, dim - 1> A;
                    A(0, 0) = e01.length2();
                    A(1, 0) = A(0, 1) = e01.dot(e02);
                    A(1, 1) = e02.length2();

                    MATRIX<T, dim - 1> temp;

                    if constexpr (useNH) {
                        MATRIX<T, dim - 1> IA = A; 
                        IA.invert();
                        const T lnJ = std::log(A.determinant() * IB.determinant()) / 2;
                        temp = h * h * vol * (mu / 2 * IB + (-mu + lambda * lnJ) / 2 * IA);
                    }
                    else {
                        MATRIX<T, dim - 1> M = IB * A; M(0, 0) -= 1; M(1, 1) -= 1;
                        temp = h * h * vol / 4 * (lambda * M.trace() * IB + 2 * mu * M * IB);
                    }

                    Eigen::Matrix<T, 4, 9> dA_div_dx;
                    dA_div_dx.setZero();
                    dA_div_dx.template block<1, 3>(0, 3) += 2.0 * (x2e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(0, 0) -= 2.0 * (x2e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(1, 6) += (x2e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(1, 3) += (x3e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(1, 0) += -(x2e - x1e).transpose() - (x3e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(2, 6) += (x2e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(2, 3) += (x3e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(2, 0) += -(x2e - x1e).transpose() - (x3e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(3, 6) += 2.0 * (x3e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(3, 0) -= 2.0 * (x3e - x1e).transpose();

                    for (int endI = 0; endI < dim; ++endI) {
                        VECTOR<T, dim>& g = std::get<FIELDS<MESH_NODE_ATTR<T, dim>>::g>(nodeAttr.Get_Unchecked(elemVInd[endI]));
                        for (int dimI = 0; dimI < dim; ++dimI) {
                            int i = endI * dim + dimI;
                            g[dimI] += dA_div_dx(0, i) * temp(0, 0) + dA_div_dx(1, i) * temp(1, 0) + 
                                dA_div_dx(2, i) * temp(0, 1) + dA_div_dx(3, i) * temp(1, 1);
                        }
                    }
                }
            }
        });
        std::cout << "Number of degenerated_smock_elems: " << degenerated_smock_elems << std::endl;
    }
}

template<class T, int dim, bool useNH = true, bool useARAP = false>
void Compute_Smock_Membrane_Hessian(
    MESH_ELEM<dim - 1>& Elem, 
    MESH_ELEM<dim - 1>& Elem_smock, 
    T h, bool projectSPD,
    const std::vector<bool>& DBCb,
    MESH_NODE<T, dim>& X, // mid-surface node coordinates
    MESH_NODE_ATTR<T, dim>& nodeAttr,
    MESH_ELEM_ATTR<T, dim - 1>& elemAttr,
    MESH_ELEM_ATTR<T, dim - 1>& elemAttr_smock,
    FIXED_COROTATED<T, dim - 1>& elasticityAttr,
    FIXED_COROTATED<T, dim - 1>& elasticityAttr_smock,
    std::vector<Eigen::Triplet<T>>& triplets)
{
    TIMER_FLAG("Compute_Membrane_Hessian");
    if constexpr (dim == 2) {
        //TODO
    }
    else {
        Eigen::Matrix<T, 9, 9> ahess[4];
        for (int i = 0; i < 4; i++) {
            ahess[i].setZero();
        }
        Eigen::Matrix<T, 3, 3> I = Eigen::Matrix<T, 3, 3>::Identity();
        ahess[0].template block<3, 3>(0, 0) += 2.0*I;
        ahess[0].template block<3, 3>(3, 3) += 2.0*I;
        ahess[0].template block<3, 3>(0, 3) -= 2.0*I;
        ahess[0].template block<3, 3>(3, 0) -= 2.0*I;

        ahess[1].template block<3, 3>(3, 6) += I;
        ahess[1].template block<3, 3>(6, 3) += I;
        ahess[1].template block<3, 3>(0, 3) -= I;
        ahess[1].template block<3, 3>(0, 6) -= I;
        ahess[1].template block<3, 3>(3, 0) -= I;
        ahess[1].template block<3, 3>(6, 0) -= I;
        ahess[1].template block<3, 3>(0, 0) += 2.0*I;

        ahess[2].template block<3, 3>(3, 6) += I;
        ahess[2].template block<3, 3>(6, 3) += I;
        ahess[2].template block<3, 3>(0, 3) -= I;
        ahess[2].template block<3, 3>(0, 6) -= I;
        ahess[2].template block<3, 3>(3, 0) -= I;
        ahess[2].template block<3, 3>(6, 0) -= I;
        ahess[2].template block<3, 3>(0, 0) += 2.0*I;

        ahess[3].template block<3, 3>(0, 0) += 2.0*I;
        ahess[3].template block<3, 3>(6, 6) += 2.0*I;
        ahess[3].template block<3, 3>(0, 6) -= 2.0*I;
        ahess[3].template block<3, 3>(6, 0) -= 2.0*I;

        std::vector<int> tripletStartInd(Elem.size + Elem_smock.size);
        int nonDBCElemCount = 0;
        Elem.Each([&](int id, auto data) {
            auto &[elemVInd] = data;
            if (!(DBCb[elemVInd[0]] && DBCb[elemVInd[1]] && DBCb[elemVInd[2]])) {
                tripletStartInd[id] = triplets.size() + nonDBCElemCount * 81;
                ++nonDBCElemCount;
            }
            else {
                tripletStartInd[id] = -1;
            }
        });

        Elem_smock.Each([&](int id, auto data) {
            auto &[elemVInd] = data;
            id += Elem.size;
            if (!(DBCb[elemVInd[0]] && DBCb[elemVInd[1]] && DBCb[elemVInd[2]])) {
                tripletStartInd[id] = triplets.size() + nonDBCElemCount * 81;
                ++nonDBCElemCount;
            }
            else {
                tripletStartInd[id] = -1;
            }
        });

        std::cout << "Starting computing hessian elems with elemCount: " << nonDBCElemCount << std::endl;

        triplets.resize(triplets.size() + nonDBCElemCount * 81);
        Elem.Join(elasticityAttr).Par_Each([&](int id, auto data) {
            if (tripletStartInd[id] >= 0) {
                auto &[elemVInd, F, vol, lambda, mu] = data;
                const VECTOR<T, dim>& x1 = std::get<0>(X.Get_Unchecked(elemVInd[0]));
                const VECTOR<T, dim>& x2 = std::get<0>(X.Get_Unchecked(elemVInd[1]));
                const VECTOR<T, dim>& x3 = std::get<0>(X.Get_Unchecked(elemVInd[2]));
                Eigen::Matrix<T, dim, 1> x1e(x1.data), x2e(x2.data), x3e(x3.data);
                MATRIX<T, dim - 1> IB = std::get<FIELDS<MESH_ELEM_ATTR<T, dim>>::IB>(elemAttr.Get_Unchecked(id));
                IB.invert();

                const VECTOR<T, dim> e01 = x2 - x1;
                const VECTOR<T, dim> e02 = x3 - x1;
                MATRIX<T, dim - 1> A;
                A(0, 0) = e01.length2();
                A(1, 0) = A(0, 1) = e01.dot(e02);
                A(1, 1) = e02.length2();

                Eigen::Matrix<T, 4, 9> dA_div_dx;
                dA_div_dx.setZero();
                dA_div_dx.template block<1, 3>(0, 3) += 2.0 * (x2e - x1e).transpose();
                dA_div_dx.template block<1, 3>(0, 0) -= 2.0 * (x2e - x1e).transpose();
                dA_div_dx.template block<1, 3>(1, 6) += (x2e - x1e).transpose();
                dA_div_dx.template block<1, 3>(1, 3) += (x3e - x1e).transpose();
                dA_div_dx.template block<1, 3>(1, 0) += -(x2e - x1e).transpose() - (x3e - x1e).transpose();
                dA_div_dx.template block<1, 3>(2, 6) += (x2e - x1e).transpose();
                dA_div_dx.template block<1, 3>(2, 3) += (x3e - x1e).transpose();
                dA_div_dx.template block<1, 3>(2, 0) += -(x2e - x1e).transpose() - (x3e - x1e).transpose();
                dA_div_dx.template block<1, 3>(3, 6) += 2.0 * (x3e - x1e).transpose();
                dA_div_dx.template block<1, 3>(3, 0) -= 2.0 * (x3e - x1e).transpose();

                Eigen::Matrix<T, 9, 9> hessian;
                if constexpr (useNH) {
                    MATRIX<T, dim - 1> IA = A; IA.invert();
                    Eigen::Matrix<T, 1, 9> ainvda;
                    for (int endI = 0; endI < dim; ++endI) {
                        for (int dimI = 0; dimI < dim; ++dimI) {
                            int i = endI * dim + dimI;
                            ainvda[i] = dA_div_dx(0, i) * IA(0, 0) + dA_div_dx(1, i) * IA(1, 0) + 
                                dA_div_dx(2, i) * IA(0, 1) + dA_div_dx(3, i) * IA(1, 1);
                        }
                    }
                    const T deta = A.determinant();
                    const T lnJ = std::log(deta * IB.determinant()) / 2;
                    const T term1 = (-mu + lambda * lnJ) / 2;
                    hessian = (-term1 + lambda / 4) * ainvda.transpose() * ainvda;

                    Eigen::Matrix<T, 4, 9> aderivadj;
                    aderivadj << dA_div_dx.row(3), -dA_div_dx.row(1), -dA_div_dx.row(2), dA_div_dx.row(0);
                    hessian += term1 / deta * aderivadj.transpose() * dA_div_dx;

                    for(int i = 0; i < 2; ++i) {
                        for(int j = 0; j < 2; ++j) {
                            hessian += (term1 * IA(i, j) + mu / 2 * IB(i, j)) * ahess[i + j * 2];
                        }
                    }
                    hessian *= h * h * vol;
                }
                else {
                    Eigen::Matrix<T, 1, 9> inner;
                    for (int i = 0; i < 9; ++i) {
                        inner[i] = dA_div_dx(0, i) * IB(0, 0) + dA_div_dx(1, i) * IB(1, 0) + 
                            dA_div_dx(2, i) * IB(0, 1) + dA_div_dx(3, i) * IB(1, 1);
                    }
                    hessian = lambda * inner.transpose() * inner;

                    MATRIX<T, dim - 1> M = IB * A; M(0, 0) -= 1; M(1, 1) -= 1;
                    MATRIX<T, dim - 1> Mainv = M * IB;
                    for(int i = 0; i < 2; ++i) {
                        for(int j = 0; j < 2; ++j) {
                            hessian += (lambda * M.trace() * IB(i, j) + 2 * mu * Mainv(i, j)) * ahess[i + j * 2];
                        }
                    }

                    Eigen::Matrix<T, 1, 9> inner00 = IB(0, 0) * dA_div_dx.row(0) + IB(0, 1) * dA_div_dx.row(2);
                    Eigen::Matrix<T, 1, 9> inner01 = IB(0, 0) * dA_div_dx.row(1) + IB(0, 1) * dA_div_dx.row(3);
                    Eigen::Matrix<T, 1, 9> inner10 = IB(1, 0) * dA_div_dx.row(0) + IB(1, 1) * dA_div_dx.row(2);
                    Eigen::Matrix<T, 1, 9> inner11 = IB(1, 0) * dA_div_dx.row(1) + IB(1, 1) * dA_div_dx.row(3);
                    hessian += 2 * mu * inner00.transpose() * inner00;
                    hessian += 2 * mu * (inner01.transpose() * inner10  + inner10.transpose() * inner01);
                    hessian += 2 * mu * inner11.transpose() * inner11;

                    hessian *= h * h * vol / 4;
                }

                if (projectSPD) {
                    makePD(hessian);
                }

                int indMap[9] = {
                    elemVInd[0] * dim,
                    elemVInd[0] * dim + 1,
                    elemVInd[0] * dim + 2,
                    elemVInd[1] * dim,
                    elemVInd[1] * dim + 1,
                    elemVInd[1] * dim + 2,
                    elemVInd[2] * dim,
                    elemVInd[2] * dim + 1,
                    elemVInd[2] * dim + 2
                };
                for (int i = 0; i < 9; ++i) {
                    int tripletIStart = tripletStartInd[id] + i * 9;
                    for (int j = 0; j < 9; ++j) {
                        triplets[tripletIStart + j] = std::move(Eigen::Triplet<T>(indMap[i], indMap[j], hessian(i, j)));
                    }
                }
            }
        });

        std::cout << "Starting computing hessian smocking elems" << std::endl;
        Elem_smock.Join(elasticityAttr_smock).Par_Each([&](int id, auto data) {
            
            if (tripletStartInd[id + Elem.size] >= 0) { // TODO: reset id 
                auto &[elemVInd, F, vol, lambda, mu] = data;
                const VECTOR<T, dim>& x1 = std::get<0>(X.Get_Unchecked(elemVInd[0]));
                const VECTOR<T, dim>& x2 = std::get<0>(X.Get_Unchecked(elemVInd[1]));
                const VECTOR<T, dim>& x3 = std::get<0>(X.Get_Unchecked(elemVInd[2]));
                Eigen::Matrix<T, dim, 1> x1e(x1.data), x2e(x2.data), x3e(x3.data);
                MATRIX<T, dim - 1> IB = std::get<FIELDS<MESH_ELEM_ATTR<T, dim>>::IB>(elemAttr_smock.Get_Unchecked(id));
                IB.invert();

                Eigen::Matrix<T, 9, 9> hessian;

                if constexpr (useARAP){
                    Eigen::Matrix2d Dm_iInv;
                    Dm_iInv(0,0) = IB(0,0);
                    Dm_iInv(0,1) = IB(0,1);
                    Dm_iInv(1,0) = IB(1,0);
                    Dm_iInv(1,1) = IB(1,1);

                    Eigen::Matrix<T, dim, 2> Ds, tri_F;
                    Ds << x2e - x1e, x3e - x1e;
                    tri_F = Ds * Dm_iInv;

                    T i1 = tri_F.squaredNorm();
                    T i2 = (tri_F.transpose() * tri_F).squaredNorm();
                    T i3 = (tri_F.transpose() * tri_F).determinant();
                    T t1, t2;
                    t1 = -sqrt(i1 + 2 * sqrt(i3));
                    t2 = -t1; // pos root

                    T f1 = 0.5 / t2;
                    T f2 = 0.;
                    T f3 = 0.5 / (t2 * sqrt(i3));
                    Eigen::Matrix<T, 3, 2> g1 = 2 * tri_F;
                    Eigen::Matrix<T, 3, 2> g2 = 4 * tri_F * tri_F.transpose() * tri_F;
                    Eigen::Matrix<T, 3, 2> g3;
                    Eigen::Matrix<T, 2, 2> C_F = tri_F.transpose() * tri_F;
                    T det_C = C_F.determinant(); 

                    T a,b,c,d;
                    a = C_F(0,0);
                    b = C_F(0,1);
                    c = C_F(1,0);
                    d = C_F(1,1);
                    g3(0,0) = 2.0 * d * tri_F(0,0) - (b+c) * tri_F(0,1);
                    g3(1,0) = 2.0 * d * tri_F(1,0) - (b+c) * tri_F(1,1);
                    g3(2,0) = 2.0 * d * tri_F(2,0) - (b+c) * tri_F(2,1);

                    g3(3,0) = 2.0 * a * tri_F(0,1) - (b+c) * tri_F(0,0);
                    g3(4,0) = 2.0 * a * tri_F(1,1) - (b+c) * tri_F(1,0);
                    g3(5,0) = 2.0 * a * tri_F(2,1) - (b+c) * tri_F(2,0);

                    Eigen::Matrix<T, dim, 2> R = f1 * g1 + f2 * g2 + f3 * g3;
                    Eigen::Matrix2d S = R.transpose() * tri_F;

                    Eigen::Matrix2d V,sigmas;
                    T sigma1, sigma2;

                    Eigen::EigenSolver<Eigen::Matrix2d> es(S);
                    sigmas = es.pseudoEigenvalueMatrix();
                    V = es.pseudoEigenvectors();
                    sigma1 = sigmas(0,0);
                    sigma2 = sigmas(1,1);

                    Eigen::Matrix<T,3,2> U_thin;
                    Eigen::Matrix3d U; // thrid col of U is the deformed triMesh surface normal n;
                    U_thin = R * V;
                    U.col(0) = U_thin.col(0);
                    U.col(1) = U_thin.col(1);
                    VECTOR<T, 3> deformed_n = cross(x2 - x1, x3 - x1);
                    U.col(2)[0] = deformed_n[0];
                    U.col(2)[1] = deformed_n[1];
                    U.col(2)[2] = deformed_n[2];

                    Eigen::MatrixXd T0(3,2);
                    T0 << 0,-1,
                        1,0,
                        0,0;
                    T0 = (1 / sqrt(2)) * U * T0 * V.transpose();

                    Eigen::MatrixXd T1(3,2);
                    T1 << 0,0,
                        0,0,
                        1,0;
                    T1 = U * T1 * V.transpose();

                    Eigen::MatrixXd T2(3,2);
                    T2 << 0,0,
                        0,0,
                        0,1;
                    T2 = U * T2 * V.transpose();

                    Eigen::Matrix<T, 6, 1> vec_T0, vec_T1, vec_T2;
                    vec_T0.template block<3, 1>(0, 0) = T0.col(0);
                    vec_T0.template block<3, 1>(3, 0) = T0.col(1);
                    vec_T1.template block<3, 1>(0, 0) = T1.col(0);
                    vec_T1.template block<3, 1>(3, 0) = T1.col(1);
                    vec_T2.template block<3, 1>(0, 0) = T2.col(0);
                    vec_T2.template block<3, 1>(3, 0) = T2.col(1);

                    T lambda0,lambda1,lambda2;
                    lambda0 = 2.0 / (sigma1 + sigma2);
                    lambda1 = 1.0 / sigma1;
                    lambda2 = 1.0 / sigma2;

                    if (sigma1 + sigma2 < 2)lambda0 = 1;
                    if (sigma1 < 1)lambda1 = 1;
                    if (sigma2 < 1)lambda2 = 1;

                    Eigen::Matrix<T, 6, 6> SH = Eigen::Matrix<T, 6, 6>::Identity();
                    SH -= lambda0 * (vec_T0 * vec_T0.transpose());
                    SH -= lambda1 * (vec_T1 * vec_T1.transpose());
                    SH -= lambda2 * (vec_T2 * vec_T2.transpose());

                    Eigen::Matrix<T, 6, 9> pF_ipX;
                    pF_ipX.setZero();
                    compute_pFpX(Dm_iInv,pF_ipX);

                    hessian = pF_ipX.transpose() * SH * pF_ipX;

                    hessian *= 0.5 * h * h * vol;
                }


                else{
                    const VECTOR<T, dim> e01 = x2 - x1;
                    const VECTOR<T, dim> e02 = x3 - x1;
                    MATRIX<T, dim - 1> A;
                    A(0, 0) = e01.length2();
                    A(1, 0) = A(0, 1) = e01.dot(e02);
                    A(1, 1) = e02.length2();

                    Eigen::Matrix<T, 4, 9> dA_div_dx;
                    dA_div_dx.setZero();
                    dA_div_dx.template block<1, 3>(0, 3) += 2.0 * (x2e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(0, 0) -= 2.0 * (x2e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(1, 6) += (x2e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(1, 3) += (x3e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(1, 0) += -(x2e - x1e).transpose() - (x3e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(2, 6) += (x2e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(2, 3) += (x3e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(2, 0) += -(x2e - x1e).transpose() - (x3e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(3, 6) += 2.0 * (x3e - x1e).transpose();
                    dA_div_dx.template block<1, 3>(3, 0) -= 2.0 * (x3e - x1e).transpose();

                    
                    if constexpr (useNH) {
                        MATRIX<T, dim - 1> IA = A; IA.invert();
                        Eigen::Matrix<T, 1, 9> ainvda;
                        for (int endI = 0; endI < dim; ++endI) {
                            for (int dimI = 0; dimI < dim; ++dimI) {
                                int i = endI * dim + dimI;
                                ainvda[i] = dA_div_dx(0, i) * IA(0, 0) + dA_div_dx(1, i) * IA(1, 0) + 
                                    dA_div_dx(2, i) * IA(0, 1) + dA_div_dx(3, i) * IA(1, 1);
                            }
                        }
                        const T deta = A.determinant();
                        const T lnJ = std::log(deta * IB.determinant()) / 2;
                        const T term1 = (-mu + lambda * lnJ) / 2;
                        hessian = (-term1 + lambda / 4) * ainvda.transpose() * ainvda;

                        Eigen::Matrix<T, 4, 9> aderivadj;
                        aderivadj << dA_div_dx.row(3), -dA_div_dx.row(1), -dA_div_dx.row(2), dA_div_dx.row(0);
                        hessian += term1 / deta * aderivadj.transpose() * dA_div_dx;

                        for(int i = 0; i < 2; ++i) {
                            for(int j = 0; j < 2; ++j) {
                                hessian += (term1 * IA(i, j) + mu / 2 * IB(i, j)) * ahess[i + j * 2];
                            }
                        }
                        hessian *= h * h * vol;
                    }
                    else {
                        Eigen::Matrix<T, 1, 9> inner;
                        for (int i = 0; i < 9; ++i) {
                            inner[i] = dA_div_dx(0, i) * IB(0, 0) + dA_div_dx(1, i) * IB(1, 0) + 
                                dA_div_dx(2, i) * IB(0, 1) + dA_div_dx(3, i) * IB(1, 1);
                        }
                        hessian = lambda * inner.transpose() * inner;

                        MATRIX<T, dim - 1> M = IB * A; M(0, 0) -= 1; M(1, 1) -= 1;
                        MATRIX<T, dim - 1> Mainv = M * IB;
                        for(int i = 0; i < 2; ++i) {
                            for(int j = 0; j < 2; ++j) {
                                hessian += (lambda * M.trace() * IB(i, j) + 2 * mu * Mainv(i, j)) * ahess[i + j * 2];
                            }
                        }

                        Eigen::Matrix<T, 1, 9> inner00 = IB(0, 0) * dA_div_dx.row(0) + IB(0, 1) * dA_div_dx.row(2);
                        Eigen::Matrix<T, 1, 9> inner01 = IB(0, 0) * dA_div_dx.row(1) + IB(0, 1) * dA_div_dx.row(3);
                        Eigen::Matrix<T, 1, 9> inner10 = IB(1, 0) * dA_div_dx.row(0) + IB(1, 1) * dA_div_dx.row(2);
                        Eigen::Matrix<T, 1, 9> inner11 = IB(1, 0) * dA_div_dx.row(1) + IB(1, 1) * dA_div_dx.row(3);
                        hessian += 2 * mu * inner00.transpose() * inner00;
                        hessian += 2 * mu * (inner01.transpose() * inner10  + inner10.transpose() * inner01);
                        hessian += 2 * mu * inner11.transpose() * inner11;

                        hessian *= h * h * vol / 4;
                    }
                }

                if (projectSPD) {
                    makePD(hessian);
                }

                int indMap[9] = {
                    elemVInd[0] * dim,
                    elemVInd[0] * dim + 1,
                    elemVInd[0] * dim + 2,
                    elemVInd[1] * dim,
                    elemVInd[1] * dim + 1,
                    elemVInd[1] * dim + 2,
                    elemVInd[2] * dim,
                    elemVInd[2] * dim + 1,
                    elemVInd[2] * dim + 2
                };
                id += Elem.size;
                for (int i = 0; i < 9; ++i) {
                    int tripletIStart = tripletStartInd[id] + i * 9;
                    for (int j = 0; j < 9; ++j) {
                        triplets[tripletIStart + j] = std::move(Eigen::Triplet<T>(indMap[i], indMap[j], hessian(i, j)));
                    }
                }
            }
        });
    }
}



}