/*****************************************************************************
* Copyright (C) 2011 Adrien Maglo
*
* This file is part of PPMC.
*
* PPMC is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* PPMC is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with PPMC.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include "../PPMC/frenetRotation.h"


using namespace CGAL;


// Determine the Frenet coordinate frame.
void determineFrenetFrame(const MyMesh::Halfedge_handle &heh_gate,
                          const Vector &normal,
                          Vector &t1, Vector &t2)
{
    Point P = heh_gate->vertex()->point();
    Point Q = heh_gate->opposite()->vertex()->point();

    Vector gateVector = P - Q;

    // Determine t1.
    t1 = gateVector - (gateVector * normal) * normal;
    float f_t1Length = sqrt(t1 * t1);
    if (f_t1Length != 0)
        t1 = t1 / f_t1Length;

    // Determine t2;
    t2 = CGAL::cross_product(normal,t1);
}


inline int signe(const float x)
{
    return (x < 0) ? -1 : 1;
}

// Matrix product.
Matrix matProd(Matrix const &m, Matrix const &n)
{
    float mT[3][3];
    mT[0][0] = m.r0().x(); mT[0][1] = m.r0().y(); mT[0][2] = m.r0().z();
    mT[1][0] = m.r1().x(); mT[1][1] = m.r1().y(); mT[1][2] = m.r1().z();
    mT[2][0] = m.r2().x(); mT[2][1] = m.r2().y(); mT[2][2] = m.r2().z();

    float nT[3][3];
    nT[0][0] = n.r0().x(); nT[0][1] = n.r0().y(); nT[0][2] = n.r0().z();
    nT[1][0] = n.r1().x(); nT[1][1] = n.r1().y(); nT[1][2] = n.r1().z();
    nT[2][0] = n.r2().x(); nT[2][1] = n.r2().y(); nT[2][2] = n.r2().z();

    float C[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

    for(int i=0; i< 3 ; i++)
    {
        for(int j=0; j<3; j++)
        {
            for(int k=0;k<3;k++)
                C[i][j] += mT[i][k]*nT[k][j];
        }
    }

    return Matrix(C[0][0], C[0][1], C[0][2],
                  C[1][0], C[1][1], C[1][2],
                  C[2][0], C[2][1], C[2][2]);
}


// Description : To find a bijection through a rotation transformation in 3D with only integer coordinates.
std::vector<Matrix> compRotationMat(const Vector &T1, const Vector &T2,
                                    const Vector &normal)
{
    Matrix rotMat(T1.x(), T2.x(), normal.x(),
                  T1.y(), T2.y(), normal.y(),
                  T1.z(), T2.z(), normal.z());
    Matrix M = rotMat;

    Matrix D1(1, 0, 0,
              0, 1, 0,
              0, 0, 1);
    Matrix D2(0, 1, 0,
              1, 0, 0,
              0, 0, 1);
    Matrix D3(-1,  0, 0,
               0, -1, 0,
               0,  0, 1);
    Matrix D4(1, 0, 0,
              0, 0, 1,
              0, 1, 0);
    Matrix D5(1,  0,  0,
              0, -1,  0,
              0,  0, -1);

    std::vector<Matrix> S(15, NULL_MATRIX);

    // Verify in order to find the smallest rotation angle.
    if (abs(M.r0().z()) > abs(M.r1().z()))
        S[0] = D2;
    else
        S[0] = D1;

    M = matProd(S[0], M);

    if (M.r1().z() < 0)
        S[1] = D3;
    else
        S[1] = D1;

    M = matProd(S[1], M);


    // Determine first rotation angle phi.
    float len = square(M.r0().z()) + square(M.r1().z());
    float phi = len == 0 ? 0 : signe(-1 * M.r0().z()) * acos(M.r1().z() / sqrt(len));

    S[2] = Matrix(1, -tan(phi / 2), 0,
                  0,             1, 0,
                  0,             0, 1);
    S[3] = Matrix(1,        0, 0,
                  sin(phi), 1, 0,
                  0,        0, 1);
    S[4] = S[2];

    Matrix R1inv( cos(phi), sin(phi), 0,
                 -sin(phi), cos(phi), 0,
                         0,        0, 1);

    M = matProd(R1inv, M);

    if (abs(M.r1().z()) > abs(M.r2().z()))
        S[5] = D4;
    else
        S[5] = D1;

    M = matProd(S[5], M);

    if(M.r2().z() < 0)
        S[6] = D5;
    else
        S[6] = D1;

    M = matProd(S[6], M);


    // Determine second rotation angle psi.
    len = square(M.r1().z()) + square(M.r2().z());
    float psi = len == 0 ? 0 : signe(-1 * M.r1().z()) * acos(M.r2().z() / sqrt(len));

    S[7] = Matrix(1, 0,           0,
                  0, 1, -tan(psi/2),
                  0, 0,           1);
    S[8] = Matrix(1,        0, 0,
                  0,        1, 0,
                  0, sin(psi), 1);
    S[9] = S[7];

    Matrix R2inv(1,         0,        0,
                 0,  cos(psi), sin(psi),
                 0, -sin(psi), cos(psi));

    M = matProd(R2inv, M);

    if(abs(M.r0().y()) > abs(M.r1().y()))
        S[10] = D2;
    else
        S[10] = D1;

    M = matProd(S[10], M);

    if(M.r1().y() < 0)
        S[11] = D3;
    else
        S[11] = D1;

    M = matProd(S[11], M);

    // Determine last rotation angle theta.
    len = square(M.r0().y()) + square(M.r1().y());
    float theta = len == 0 ? 0 : signe(-1 * M.r0().y()) * std::acos(M.r1().y() / sqrt(len));

    S[12] = Matrix(1, -tan(theta/2), 0,
                   0,             1, 0,
                   0,             0, 1);
    S[13] = Matrix(         1, 0, 0,
                   sin(theta), 1, 0,
                            0, 0, 1);
    S[14] = S[12];

    return S;
}


VectorInt frenetRotation(const VectorInt &Dist, const Vector &T1,
                         const Vector &T2, const Vector &normal)
{
    std::vector<Matrix> S = compRotationMat(T1,T2, normal);

    VectorDouble u(Dist.x(),Dist.y(),Dist.z());
    Matrix m_inter(NULL_MATRIX);

    // Procedure of the bijection.
    for (int i = 0; i < 15; i++)
    {
        if (i == 0 || i == 1 || i == 5
            || i == 6 || i == 10 || i == 11)
            m_inter = S[i];
        else
            m_inter = Matrix(S[i].r0().x(), -S[i].r0().y(), -S[i].r0().z(),
                             -S[i].r1().x(), S[i].r1().y(), -S[i].r1().z(),
                             -S[i].r2().x(), -S[i].r2().y(), S[i].r2().z());
        u = m_inter * u;

        int x = round(u.x());
        int y = round(u.y());
        int z = round(u.z());

        u = VectorDouble(x, y, z);
    }

    return VectorInt(u.x(), u.y(), u.z());
}


VectorInt invFrenetRotation(const VectorInt &Frenet, const Vector &T1,
                            const Vector &T2, const Vector &normal)
{
    std::vector<Matrix> S = compRotationMat(T1,T2, normal);

    VectorDouble u(Frenet.x(),Frenet.y(),Frenet.z());
    Matrix m_inter(NULL_MATRIX);

    for (int i = 14; i > -1; i--)
    {
        m_inter = S[i];
        u =  m_inter * u;

        int x = round(u.x());
        int y = round(u.y());
        int z = round(u.z());

        u = VectorDouble(x, y, z);
    }

    return VectorInt(u.x(), u.y(), u.z());
}
