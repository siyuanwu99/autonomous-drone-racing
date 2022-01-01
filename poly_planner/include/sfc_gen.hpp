#ifndef SFC_GEN_HPP
#define SFC_GEN_HPP

#define _USE_MATH_DEFINES

#include "gcopter/quickhull.hpp"
#include "gcopter/geoutils.hpp"

#include <random>
#include <deque>
#include <Eigen/Eigen>

#include <cfloat>
#include <iostream>
#include <string>
#include <fstream>
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
        psi   = (uniformReal(gen) - 0.5) * M_PI * hturn;
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
    polyH.col(0) << gateDir, gatePos + gateDir * hlfLength;
    polyH.col(1) << -gateDir, gatePos - gateDir * hlfLength;

    /* up and down */
    polyH.col(2) << unitZ, gatePos + unitZ * hlfWidth;
    polyH.col(3) << -unitZ, gatePos - unitZ * hlfWidth;

    /* left and right */
    polyH.col(4) << unitL, gatePos + unitL * hlfWidth;
    polyH.col(5) << -unitL, gatePos - unitL * hlfWidth;

    std::cout << polyH << std::endl;

    return polyH;
}

/**
 * @brief generate polygon between two gates
 * 
 * Modified Dec. 22, 2021: remove gates, overlap between gate polygons
 * 
 * @param lpos 
 * @param ldir 
 * @param npos 
 * @param ndir 
 * @param gateWidth 
 * @param gateLength 
 * @param alpha 
 * @return Eigen::MatrixXd 
 */
inline Eigen::MatrixXd genBetweenGatePoly(const Eigen::Vector3d &lpos,
                                          const Eigen::Vector3d &ldir,
                                          const Eigen::Vector3d &npos,
                                          const Eigen::Vector3d &ndir,
                                          const double &gateWidth,
                                          const double &gateLength,
                                          const int &alpha) {

        Eigen::MatrixXd polyH(6, 6);
        double zLength = npos[2] - lpos[2];

        Eigen::Vector3d unitZ(0, 0, 1);
        double zSigned = copysign(1, zLength);

        /* up and down */
        double zShift = 0.25 * zLength + 0.5 * zSigned * gateWidth;
        polyH.col(0) << zSigned * unitZ, npos + zShift * unitZ;
        polyH.col(1) << -zSigned * unitZ, lpos - zShift * unitZ;

        /* front and back */
        double delta = gateLength / 2 - 0.1;
        polyH.col(2) << ndir, npos - delta * ndir;
        polyH.col(3) << -ldir,lpos + delta * ldir;

        /* left and right */
        Eigen::Vector3d nUnitLeft, lUnitLeft;
        nUnitLeft = unitZ.cross(ndir);
        lUnitLeft = unitZ.cross(ldir);
        double hShift = alpha * gateWidth;  // horizontal shift

    if (ndir.dot(ldir) > 0) {  /* if alpha < 180 degree */

        Eigen::Vector3d lpl, lpr, npl, npr;  // last left, new right, etc...
        lpl = lpos + hShift * lUnitLeft;
        lpr = lpos - hShift * lUnitLeft;
        npl = npos + hShift * nUnitLeft;
        npr = npos - hShift * nUnitLeft;

        Eigen::Vector3d leftEdge = npl - lpl;
        leftEdge = unitZ.cross(leftEdge);
        leftEdge.normalize();

        polyH.col(4) << leftEdge, npl;

        Eigen::Vector3d rightEdge = npr - lpr;
        rightEdge = rightEdge.cross(unitZ);
        rightEdge.normalize();

        polyH.col(5) << rightEdge, npr;
    } else { /* if alpha > 180 degree */
        Eigen::Vector3d d = npos - lpos;

        double lSigned = copysign(1, d.dot(nUnitLeft));

        Eigen::Vector3d lpl, npl;
        npl = npos + lSigned * hShift * nUnitLeft;
        lpl = lpos + lSigned * hShift * lUnitLeft;

        polyH.col(4) << -ndir, lpl;
        polyH.col(5) << ldir, npl;
    }

    std::cout << "----- out ------" << std::endl;
    std::cout << polyH << std::endl;
    std::cout << "----- in  ------" << std::endl;

    return polyH;
}

/**
 * @brief generate safe corridors
 * 
 * @param num 
 * @param offset 
 * @param gates 
 * @param gateWidth 
 * @param gateLength 
 * @param inifin 
 * @return std::vector<Eigen::Matrix<double, 6, -1>> 
 */
inline std::vector<Eigen::Matrix<double, 6, -1>> genGateSFC(const int &num,
                                               const Eigen::Vector3d &offset,
                                               const std::vector<Eigen::Matrix<double, 3, 2>> &gates,
                                               const double &gateWidth,
                                               const double &gateLength,
                                               const Eigen::MatrixXd &inifin) {


    
    std::vector<Eigen::Matrix<double, 6, -1>> hPolys;
    Eigen::Matrix<double, 6, -1> hPoly;
    hPolys.reserve(2 * num + 1);


    Eigen::Vector3d gatePos, gateAngle;
    Eigen::Vector3d lgatePos, lgateAngle;  // last value

    lgatePos = inifin.col(0);
    lgateAngle << 1, 0, 0;

    
    for (int i = 0; i < num; i++) {
        gatePos   = gates[i].col(0);
        gateAngle = gates[i].col(1);

        hPoly = genBetweenGatePoly(lgatePos, lgateAngle,
                                   gatePos, gateAngle,
                                   gateWidth, gateLength, 2);
        hPolys.push_back(hPoly);

        // hPoly = genGatePoly(gatePos, gateAngle, gateWidth, gateLength);
        // hPolys.push_back(hPoly);

        lgatePos = gatePos;
        lgateAngle = gateAngle;        
    }

    gatePos = inifin.col(1);
    gateAngle << 1, 0, 0;
    hPoly = genBetweenGatePoly(lgatePos, lgateAngle,
                               gatePos, gateAngle,
                               gateWidth, gateLength, 2);
    hPolys.push_back(hPoly);

    return hPolys;
}


/**
 * @brief generate gates
 * 
 * @param num 
 * @param zMax 
 * @param yMax 
 * @param xMax 
 * @param tMax 
 * @param inifin 
 * @return std::vector<Eigen::Matrix<double, 3, 2>> 
 */
inline std::vector<Eigen::Matrix<double, 3, 2>> genGate(const int &num,
                                                        const double zMax,
                                                        const double yMax,
                                                        const double xMax,
                                                        const double tMax,
                                                        Eigen::MatrixXd &inifin
                                                        ) {
    static std::mt19937_64 gen;
    static std::uniform_real_distribution<double> uniformReal(0.0, 1.0);
    
    std::vector<Eigen::Matrix<double, 3, 2>> gates;
    Eigen::Matrix<double, 3, 2> gate;
    gates.reserve(num);

    inifin.resize(3, 2);
    inifin.col(0) << 0, 0, 0;
    inifin.col(1) << xMax, 0, 0;

    double xGap = xMax / (num + 1);

    for (int i = 0; i < num; i++) {
        double theta = tMax * M_PI / 180 * (uniformReal(gen) - 0.5);

        gate.col(0) << xGap * (i + 1),
                       yMax * (uniformReal(gen) - 0.5),
                       zMax * uniformReal(gen);
        gate.col(1) << cos(theta), sin(theta), 0;

        gates.push_back(gate);

        std::cout << "===============" << std::endl;
        std::cout << "Gate" << i + 1 << '\n' << gate << std::endl;
    }

    return gates;

}

/**
 * @brief read gates data (location and direction) from filepath
 * 
 *      first line: start and end position
 *      folloing lines: gates position and direction angle
 * @param filepath 
 * @param inifin 
 * @return std::vector<Eigen::Matrix<double, 3, 2>> 
 */
inline std::vector<Eigen::Matrix<double, 3, 2>> readGate(const std::string filepath,
                                                         Eigen::MatrixXd &inifin) {
    inifin.resize(3, 2);

    std::vector<Eigen::Matrix<double, 3, 2>> gates;
    Eigen::Matrix<double, 3, 2> gate;
    std::vector<double> buffer;
    buffer.reserve(4);

    std::ifstream fin(filepath);
    std::string line_info, input_rst;
    if (fin) {
        for (int i = 0; getline(fin, line_info); i++) {
            std::stringstream input(line_info);
            for (int j = 0; input >> input_rst; ++j) {
                std::string::size_type size;
                if (i == 0) {
                    inifin(j) = std::stof(input_rst, &size);
                } else {
                    buffer[j] = std::stof(input_rst, &size);
                }
            }

            if (i == 0 ) {
                std::cout << inifin << std::endl;
                continue;
            }
            gate.col(0) << buffer[0], buffer[1], buffer[2];
            double theta = buffer[3] * M_PI / 180;
            buffer.clear();
            gate.col(1) << cos(theta), sin(theta), 0;
            gates.push_back(gate);
            std::cout << "===== Gate " << i << " =====" << std::endl; 
            std::cout << gate << std::endl;
        }

    } else {
        ROS_ERROR("[Reader]: Failed to read gates from files, file doesnot exist.");
    }
    return gates;
}

#endif
