#ifndef SFC_GEN_HPP
#define SFC_GEN_HPP

#define _USE_MATH_DEFINES

#include "quickhull.hpp"
#include "geoutils.hpp"

#include <random>
#include <deque>
#include <Eigen/Eigen>

#include <cfloat>
#include <iostream>
#include <set>

inline bool extractV(const Eigen::MatrixXd &hPoly,
                     Eigen::Matrix3Xd &vPoly)
{
    auto hPolyON = hPoly;
    hPolyON.topRows<3>().array() *= -1.0;

    return geoutils::enumerateVs(hPolyON, vPoly);
}

inline void simplifySFC(std::vector<Eigen::MatrixXd> &hPolys)
{
    std::vector<Eigen::MatrixXd> origPolys = hPolys;
    hPolys.clear();
    int M = origPolys.size();
    Eigen::MatrixXd hPoly;
    Eigen::Matrix3Xd vPoly;
    bool overlap;
    std::deque<int> idices;
    idices.push_front(M - 1);
    for (int i = M - 1; i >= 0; i--)
    {
        for (int j = 0; j < i; j++)
        {
            hPoly.resize(6, origPolys[i].cols() + origPolys[j].cols());
            hPoly << origPolys[i], origPolys[j];
            extractV(hPoly, vPoly);
            overlap = extractV(hPoly, vPoly);
            if (overlap)
            {
                idices.push_front(j);
                i = j + 1;
                break;
            }
        }
    }
    for (const auto &ele : idices)
    {
        hPolys.push_back(origPolys[ele]);
    }
}
inline Eigen::MatrixXd genPolyH(const int &nf,
                                const double &rmin,
                                const double &rmax)
{
    static std::mt19937_64 gen;
    static std::uniform_real_distribution<double> uniformReal(0.0, 1.0);

    const int dmf = nf + 4;
    const double smax = 16.0;
    const double cmin = sqrt(6.0) / 12.0 / smax;
    const double crange = 1.0 - cmin;
    Eigen::MatrixXd dualPolyIR(3, dmf);

    dualPolyIR.col(0) << sqrt(3.0) / 3.0, 0.0, -sqrt(6.0) / 12.0;
    dualPolyIR.col(1) << -sqrt(3.0) / 6.0, 0.5, -sqrt(6.0) / 12.0;
    dualPolyIR.col(2) << -sqrt(3.0) / 6.0, -0.5, -sqrt(6.0) / 12.0;
    dualPolyIR.col(3) << 0.0, 0.0, sqrt(6.0) / 4.0;
    dualPolyIR.leftCols(4) /= smax;

    double cur;
    for (int i = 4; i < dmf; i++)
    {
        cur = (2.0 * uniformReal(gen) - 1.0) * crange;
        cur += cur > 0 ? cmin : -cmin;
        dualPolyIR.col(i)(0) = cur;
        cur = (2.0 * uniformReal(gen) - 1.0) * crange;
        cur += cur > 0 ? cmin : -cmin;
        dualPolyIR.col(i)(1) = cur;
        cur = (2.0 * uniformReal(gen) - 1.0) * crange;
        cur += cur > 0 ? cmin : -cmin;
        dualPolyIR.col(i)(2) = cur;
    }

    Eigen::MatrixXd dualPolyHR(6, dmf);
    for (int i = 0; i < dmf; i++)
    {
        dualPolyHR.col(i) << -dualPolyIR.col(i),
            dualPolyIR.col(i) / dualPolyIR.col(i).squaredNorm();
    }
    Eigen::Matrix3Xd polyI;
    extractV(dualPolyHR, polyI);

    int mf = polyI.cols();
    Eigen::MatrixXd polyH(6, mf);
    for (int i = 0; i < mf; i++)
    {
        polyH.col(i) << -polyI.col(i).normalized(),
            polyI.col(i) / polyI.col(i).squaredNorm();
    }
    Eigen::Matrix3Xd polyV;
    extractV(polyH, polyV);

    double curd = (rmax - rmin) * uniformReal(gen) + rmin;
    curd /= polyV.colwise().norm().maxCoeff();
    polyH.bottomRows<3>() *= curd;

    return polyH;
}

inline double intersectDist(const Eigen::Vector3d &o,
                            const Eigen::Vector3d &ray,
                            const Eigen::MatrixXd &hPoly)
{
    const Eigen::Vector3d r = ray.normalized();
    double md = INFINITY;
    double a, d;
    for (int i = 0; i < hPoly.cols(); i++)
    {
        a = hPoly.col(i).head<3>().dot(r);
        if (a != 0.0)
        {
            d = hPoly.col(i).head<3>().dot(hPoly.col(i).tail<3>() - o) / a;
            md = (d >= 0 && md > d) ? d : md;
        }
    }
    return md;
}

inline std::vector<Eigen::MatrixXd> genPolySFC(const int &num,
                                               const Eigen::Vector3d &offset,
                                               const double &rmin,
                                               const double &rmax,
                                               const double &drange,
                                               const int &nf,
                                               const double &vturn,
                                               const double &hturn,
                                               Eigen::MatrixXd &inifin)
{
    static std::mt19937_64 gen;
    static std::uniform_real_distribution<double> uniformReal(0.0, 1.0);

    std::vector<Eigen::MatrixXd> hPolys;
    Eigen::MatrixXd hPoly;
    hPolys.reserve(num);
    inifin.resize(3, 2);

    double theta, psi, ldm, cdm, dmax, dmin, d;
    Eigen::Vector3d lo, co, fd;
    Eigen::Vector3d e1(1.0, 0.0, 0.0), e2(0.0, 1.0, 0.0), e3(0.0, 0.0, 1.0);

    co = offset;
    inifin.col(0) = co;
    hPoly = genPolyH(nf, rmin, rmax);
    hPoly.bottomRows<3>().colwise() += co;
    hPolys.push_back(hPoly);
    for (int i = 1; i < num; i++)
    {
        lo = co;
        theta = (uniformReal(gen) - 0.5) * M_PI * vturn;
        psi = (uniformReal(gen) - 0.5) * M_PI * hturn;
        fd = e1 + e2 * tan(psi) + e3 * tan(theta);
        fd.normalize();

        hPoly = genPolyH(nf, rmin, rmax);

        ldm = intersectDist(lo, fd, hPolys.back());
        cdm = intersectDist(Eigen::Vector3d::Zero(), -fd, hPoly);
        dmax = ldm + cdm;
        dmin = std::max(ldm, cdm);
        d = uniformReal(gen) * (dmax - dmin) * drange + dmin;
        co = lo + fd * d;

        hPoly.bottomRows<3>().colwise() += co;
        hPolys.push_back(hPoly);
    }
    inifin.col(1) = co;

    return hPolys;
}

inline Eigen::MatrixXd genGatePoly(const Eigen::Vector3d &gatePos,
                                   const Eigen::Vector3d &gateDir,
                                   const double &gateWidth,
                                   const double &gateLength) {
    double hlfLength = gateLength / 2;
    double hlfWidth = gateWidth / 2;
    Eigen::Vector3d unitZ(0, 0, 1);
    Eigen::Vector3d unitL;
    unitL = unitZ.cross(gateDir);
    unitL.normalize();

    Eigen::MatrixXd polyH(6, 6);

    /* front and back */
    polyH.col(0) << gatePos + gateDir * hlfLength, gateDir;
    polyH.col(1) << gatePos - gateDir * hlfLength, -gateDir;

    /* up and down */
    polyH.col(2) << gatePos + unitZ * hlfWidth, unitZ;
    polyH.col(3) << gatePos - unitZ * hlfWidth, -unitZ;

    /* left and right */
    polyH.col(4) << gatePos + unitL * hlfWidth, unitL;
    polyH.col(5) << gatePos - unitL * hlfWidth, -unitL;

    return polyH;
}

inline Eigen::MatrixXd genBetweenGatePoly(const Eigen::Vector3d &lpos,
                                          const Eigen::Vector3d &ldir,
                                          const Eigen::Vector3d &npos,
                                          const Eigen::Vector3d &ndir,
                                          const double &gateWidth,
                                          const int &alpha) {
    Eigen::MatrixXd polyH(6, 6);
    double zLength = npos[2] - lpos[2];

    Eigen::Vector3d unitZ(0, 0, 1);
    double zSigned = copysign(1, zLength);

    /* up and down */
    double zShift = zLength + zSigned * gateWidth;
    polyH.col(0) << npos + zShift * unitZ, zSigned * unitZ;
    polyH.col(1) << lpos - zShift * unitZ, zSigned * unitZ;

    /* front and back */
    polyH.col(2) << npos, ndir;
    polyH.col(3) << lpos, -ldir;

    /* left and right */
    Eigen::Vector3d nUnitLeft, lUnitLeft;
    nUnitLeft = unitZ.cross(ndir);
    lUnitLeft = unitZ.cross(ldir);
    double hShift = alpha * gateWidth;  // horizontal shift

    Eigen::Vector3d lpl, lpr, npl, npr;  // last left, new right, etc...
    lpl = lpos + hShift * lUnitLeft;
    lpr = lpos - hShift * lUnitLeft;
    npl = npos + hShift * nUnitLeft;
    npr = npos - hShift * nUnitLeft;

    Eigen::Vector3d leftEdge = npl - lpl;
    leftEdge = unitZ.cross(leftEdge);
    leftEdge.normalize();

    polyH.col(4) << npl, leftEdge;

    Eigen::Vector3d rightEdge = npr - lpr;
    rightEdge = rightEdge.cross(unitZ);
    rightEdge.normalize();

    polyH.col(5) << npr, rightEdge;

    return polyH;
}

inline std::vector<Eigen::MatrixXd> genGateSFC(const int &num,
                                               const Eigen::Vector3d &offset,
                                               const std::vector<Eigen::Matrix<double, 3, 2>> &gates,
                                               const double &gateWidth,
                                               const double &gateLength,
                                               Eigen::MatrixXd &inifin) {

    static std::mt19937_64 gen;
    static std::uniform_real_distribution<double> uniformReal(0.0, 1.0);
    
    std::vector<Eigen::MatrixXd> hPolys;
    Eigen::MatrixXd hPoly;
    hPolys.reserve(2 * num + 1);
    inifin.resize(3, 2);

    inifin.col(0) << 0, 0, 0;
    inifin.col(1) << 40, 0, 0;

    Eigen::Vector3d gatePos, gateAngle;
    Eigen::Vector3d lgatePos, lgateAngle;  // last value

    lgatePos = offset;
    lgateAngle << 1, 0, 0;

    
    for (int i = 0; i < num; i++) {
        gatePos   = gates[i].col(0);
        gateAngle = gates[i].col(1);

        hPoly = genBetweenGatePoly(lgatePos, lgateAngle,
                                   gatePos, gateAngle,
                                   gateWidth, 2);
        hPolys.push_back(hPoly);

        hPoly = genGatePoly(gatePos, gateAngle, gateWidth, gateLength);
        hPolys.push_back(hPoly);

        lgatePos = gatePos;
        lgateAngle = gateAngle;        
    }

    gatePos = inifin.col(1);
    gateAngle << -1, 0, 0;
    hPoly = genBetweenGatePoly(lgatePos, lgateAngle,
                               gatePos, gateAngle,
                               gateWidth, 2);
    hPolys.push_back(hPoly);

    return hPolys;
}

#endif
