#ifndef GCOPTER_HPP
#define GCOPTER_HPP

#include "root_finder.hpp"
#include "lbfgs.hpp"
#include "geoutils.hpp"

#include <Eigen/Eigen>

#include <iostream>
#include <cmath>
#include <cfloat>
#include <vector>

typedef Eigen::Matrix<double, 3, 6> CoefficientMat;
typedef Eigen::Matrix<double, 3, 5> VelCoefficientMat;
typedef Eigen::Matrix<double, 3, 4> AccCoefficientMat;

class Piece
{
private:
    double duration;
    CoefficientMat coeffMat;

public:
    Piece() = default;

    Piece(double dur, const CoefficientMat &cMat)
        : duration(dur), coeffMat(cMat) {}

    inline int getDim() const
    {
        return 3;
    }

    inline int getOrder() const
    {
        return 5;
    }

    inline double getDuration() const
    {
        return duration;
    }

    inline const CoefficientMat &getCoeffMat() const
    {
        return coeffMat;
    }

    inline Eigen::Vector3d getPos(const double &t) const
    {
        Eigen::Vector3d pos(0.0, 0.0, 0.0);
        double tn = 1.0;
        for (int i = 5; i >= 0; i--)
        {
            pos += tn * coeffMat.col(i);
            tn *= t;
        }
        return pos;
    }

    inline Eigen::Vector3d getVel(const double &t) const
    {
        Eigen::Vector3d vel(0.0, 0.0, 0.0);
        double tn = 1.0;
        int n = 1;
        for (int i = 4; i >= 0; i--)
        {
            vel += n * tn * coeffMat.col(i);
            tn *= t;
            n++;
        }
        return vel;
    }

    inline Eigen::Vector3d getAcc(const double &t) const
    {
        Eigen::Vector3d acc(0.0, 0.0, 0.0);
        double tn = 1.0;
        int m = 1;
        int n = 2;
        for (int i = 3; i >= 0; i--)
        {
            acc += m * n * tn * coeffMat.col(i);
            tn *= t;
            m++;
            n++;
        }
        return acc;
    }

    inline Eigen::Vector3d getJer(const double &t) const
    {
        Eigen::Vector3d jer(0.0, 0.0, 0.0);
        double tn = 1.0;
        int l = 1;
        int m = 2;
        int n = 3;
        for (int i = 2; i >= 0; i--)
        {
            jer += l * m * n * tn * coeffMat.col(i);
            tn *= t;
            l++;
            m++;
            n++;
        }
        return jer;
    }

    inline CoefficientMat normalizePosCoeffMat() const
    {
        CoefficientMat nPosCoeffsMat;
        double t = 1.0;
        for (int i = 5; i >= 0; i--)
        {
            nPosCoeffsMat.col(i) = coeffMat.col(i) * t;
            t *= duration;
        }
        return nPosCoeffsMat;
    }

    inline VelCoefficientMat normalizeVelCoeffMat() const
    {
        VelCoefficientMat nVelCoeffMat;
        int n = 1;
        double t = duration;
        for (int i = 4; i >= 0; i--)
        {
            nVelCoeffMat.col(i) = n * coeffMat.col(i) * t;
            t *= duration;
            n++;
        }
        return nVelCoeffMat;
    }

    inline AccCoefficientMat normalizeAccCoeffMat() const
    {
        AccCoefficientMat nAccCoeffMat;
        int n = 2;
        int m = 1;
        double t = duration * duration;
        for (int i = 3; i >= 0; i--)
        {
            nAccCoeffMat.col(i) = n * m * coeffMat.col(i) * t;
            n++;
            m++;
            t *= duration;
        }
        return nAccCoeffMat;
    }

    inline double getMaxVelRate() const
    {
        VelCoefficientMat nVelCoeffMat = normalizeVelCoeffMat();
        Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                                RootFinder::polySqr(nVelCoeffMat.row(1)) +
                                RootFinder::polySqr(nVelCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++)
        {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON)
        {
            return getVel(0.0).norm();
        }
        else
        {
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON)
            {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON)
            {
                r = 0.5 * (r + 1.0);
            }
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);
            candidates.insert(0.0);
            candidates.insert(1.0);
            double maxVelRateSqr = -INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
                 it != candidates.end();
                 it++)
            {
                if (0.0 <= *it && 1.0 >= *it)
                {
                    tempNormSqr = getVel((*it) * duration).squaredNorm();
                    maxVelRateSqr = maxVelRateSqr < tempNormSqr ? tempNormSqr : maxVelRateSqr;
                }
            }
            return sqrt(maxVelRateSqr);
        }
    }

    inline double getMaxAccRate() const
    {
        AccCoefficientMat nAccCoeffMat = normalizeAccCoeffMat();
        Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                                RootFinder::polySqr(nAccCoeffMat.row(1)) +
                                RootFinder::polySqr(nAccCoeffMat.row(2));
        int N = coeff.size();
        int n = N - 1;
        for (int i = 0; i < N; i++)
        {
            coeff(i) *= n;
            n--;
        }
        if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON)
        {
            return getAcc(0.0).norm();
        }
        else
        {
            double l = -0.0625;
            double r = 1.0625;
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON)
            {
                l = 0.5 * l;
            }
            while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON)
            {
                r = 0.5 * (r + 1.0);
            }
            std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                      FLT_EPSILON / duration);
            candidates.insert(0.0);
            candidates.insert(1.0);
            double maxAccRateSqr = -INFINITY;
            double tempNormSqr;
            for (std::set<double>::const_iterator it = candidates.begin();
                 it != candidates.end();
                 it++)
            {
                if (0.0 <= *it && 1.0 >= *it)
                {
                    tempNormSqr = getAcc((*it) * duration).squaredNorm();
                    maxAccRateSqr = maxAccRateSqr < tempNormSqr ? tempNormSqr : maxAccRateSqr;
                }
            }
            return sqrt(maxAccRateSqr);
        }
    }

    inline bool checkMaxVelRate(const double &maxVelRate) const
    {
        double sqrMaxVelRate = maxVelRate * maxVelRate;
        if (getVel(0.0).squaredNorm() >= sqrMaxVelRate ||
            getVel(duration).squaredNorm() >= sqrMaxVelRate)
        {
            return false;
        }
        else
        {
            VelCoefficientMat nVelCoeffMat = normalizeVelCoeffMat();
            Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                                    RootFinder::polySqr(nVelCoeffMat.row(1)) +
                                    RootFinder::polySqr(nVelCoeffMat.row(2));
            double t2 = duration * duration;
            coeff.tail<1>()(0) -= sqrMaxVelRate * t2;
            return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
        }
    }

    inline bool checkMaxAccRate(const double &maxAccRate) const
    {
        double sqrMaxAccRate = maxAccRate * maxAccRate;
        if (getAcc(0.0).squaredNorm() >= sqrMaxAccRate ||
            getAcc(duration).squaredNorm() >= sqrMaxAccRate)
        {
            return false;
        }
        else
        {
            AccCoefficientMat nAccCoeffMat = normalizeAccCoeffMat();
            Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                                    RootFinder::polySqr(nAccCoeffMat.row(1)) +
                                    RootFinder::polySqr(nAccCoeffMat.row(2));
            double t2 = duration * duration;
            double t4 = t2 * t2;
            coeff.tail<1>()(0) -= sqrMaxAccRate * t4;
            return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
        }
    }
};

class Trajectory
{
private:
    typedef std::vector<Piece> Pieces;
    Pieces pieces;

public:
    Trajectory() = default;

    Trajectory(const std::vector<double> &durs,
               const std::vector<CoefficientMat> &cMats)
    {
        int N = std::min(durs.size(), cMats.size());
        pieces.reserve(N);
        for (int i = 0; i < N; i++)
        {
            pieces.emplace_back(durs[i], cMats[i]);
        }
    }

    inline int getPieceNum() const
    {
        return pieces.size();
    }

    inline Eigen::VectorXd getDurations() const
    {
        int N = getPieceNum();
        Eigen::VectorXd durations(N);
        for (int i = 0; i < N; i++)
        {
            durations(i) = pieces[i].getDuration();
        }
        return durations;
    }

    inline double getTotalDuration() const
    {
        int N = getPieceNum();
        double totalDuration = 0.0;
        for (int i = 0; i < N; i++)
        {
            totalDuration += pieces[i].getDuration();
        }
        return totalDuration;
    }

    inline Eigen::Matrix3Xd getPositions() const
    {
        int N = getPieceNum();
        Eigen::Matrix3Xd positions(3, N + 1);
        for (int i = 0; i < N; i++)
        {
            positions.col(i) = pieces[i].getCoeffMat().col(5);
        }
        positions.col(N) = pieces[N - 1].getPos(pieces[N - 1].getDuration());
        return positions;
    }

    inline const Piece &operator[](int i) const
    {
        return pieces[i];
    }

    inline Piece &operator[](int i)
    {
        return pieces[i];
    }

    inline void clear(void)
    {
        pieces.clear();
        return;
    }

    inline Pieces::const_iterator begin() const
    {
        return pieces.begin();
    }

    inline Pieces::const_iterator end() const
    {
        return pieces.end();
    }

    inline Pieces::iterator begin()
    {
        return pieces.begin();
    }

    inline Pieces::iterator end()
    {
        return pieces.end();
    }

    inline void reserve(const int &n)
    {
        pieces.reserve(n);
        return;
    }

    inline void emplace_back(const Piece &piece)
    {
        pieces.emplace_back(piece);
        return;
    }

    inline void emplace_back(const double &dur,
                             const CoefficientMat &cMat)
    {
        pieces.emplace_back(dur, cMat);
        return;
    }

    inline void append(const Trajectory &traj)
    {
        pieces.insert(pieces.end(), traj.begin(), traj.end());
        return;
    }

    inline int locatePieceIdx(double &t) const
    {
        int N = getPieceNum();
        int idx;
        double dur;
        for (idx = 0;
             idx < N &&
             t > (dur = pieces[idx].getDuration());
             idx++)
        {
            t -= dur;
        }
        if (idx == N)
        {
            idx--;
            t += pieces[idx].getDuration();
        }
        return idx;
    }

    inline Eigen::Vector3d getPos(double t) const
    {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getPos(t);
    }

    inline Eigen::Vector3d getVel(double t) const
    {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getVel(t);
    }

    inline Eigen::Vector3d getAcc(double t) const
    {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getAcc(t);
    }

    inline Eigen::Vector3d getJer(double t) const
    {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getJer(t);
    }

    inline Eigen::Vector3d getJuncPos(int juncIdx) const
    {
        if (juncIdx != getPieceNum())
        {
            return pieces[juncIdx].getCoeffMat().col(5);
        }
        else
        {
            return pieces[juncIdx - 1].getPos(pieces[juncIdx - 1].getDuration());
        }
    }

    inline Eigen::Vector3d getJuncVel(int juncIdx) const
    {
        if (juncIdx != getPieceNum())
        {
            return pieces[juncIdx].getCoeffMat().col(4);
        }
        else
        {
            return pieces[juncIdx - 1].getVel(pieces[juncIdx - 1].getDuration());
        }
    }

    inline Eigen::Vector3d getJuncAcc(int juncIdx) const
    {
        if (juncIdx != getPieceNum())
        {
            return pieces[juncIdx].getCoeffMat().col(3) * 2.0;
        }
        else
        {
            return pieces[juncIdx - 1].getAcc(pieces[juncIdx - 1].getDuration());
        }
    }

    inline double getMaxVelRate() const
    {
        int N = getPieceNum();
        double maxVelRate = -INFINITY;
        double tempNorm;
        for (int i = 0; i < N; i++)
        {
            tempNorm = pieces[i].getMaxVelRate();
            maxVelRate = maxVelRate < tempNorm ? tempNorm : maxVelRate;
        }
        return maxVelRate;
    }

    inline double getMaxAccRate() const
    {
        int N = getPieceNum();
        double maxAccRate = -INFINITY;
        double tempNorm;
        for (int i = 0; i < N; i++)
        {
            tempNorm = pieces[i].getMaxAccRate();
            maxAccRate = maxAccRate < tempNorm ? tempNorm : maxAccRate;
        }
        return maxAccRate;
    }

    inline bool checkMaxVelRate(const double &maxVelRate) const
    {
        int N = getPieceNum();
        bool feasible = true;
        for (int i = 0; i < N && feasible; i++)
        {
            feasible = feasible && pieces[i].checkMaxVelRate(maxVelRate);
        }
        return feasible;
    }

    inline bool checkMaxAccRate(const double &maxAccRate) const
    {
        int N = getPieceNum();
        bool feasible = true;
        for (int i = 0; i < N && feasible; i++)
        {
            feasible = feasible && pieces[i].checkMaxAccRate(maxAccRate);
        }
        return feasible;
    }
};

// The banded system class is used for solving
// banded linear system Ax=b efficiently.
// A is an N*N band matrix with lower band width lowerBw
// and upper band width upperBw.
// Banded LU factorization has O(N) time complexity.
class BandedSystem
{
public:
    // The size of A, as well as the lower/upper
    // banded width p/q are needed
    inline void create(const int &n, const int &p, const int &q)
    {
        // In case of re-creating before destroying
        destroy();
        N = n;
        lowerBw = p;
        upperBw = q;
        int actualSize = N * (lowerBw + upperBw + 1);
        ptrData = new double[actualSize];
        std::fill_n(ptrData, actualSize, 0.0);
        return;
    }

    inline void destroy()
    {
        if (ptrData != nullptr)
        {
            delete[] ptrData;
            ptrData = nullptr;
        }
        return;
    }

private:
    int N;
    int lowerBw;
    int upperBw;
    // Compulsory nullptr initialization here
    double *ptrData = nullptr;

public:
    // Reset the matrix to zero
    inline void reset(void)
    {
        std::fill_n(ptrData, N * (lowerBw + upperBw + 1), 0.0);
        return;
    }

    // The band matrix is stored as suggested in "Matrix Computation"
    inline const double &operator()(const int &i, const int &j) const
    {
        return ptrData[(i - j + upperBw) * N + j];
    }

    inline double &operator()(const int &i, const int &j)
    {
        return ptrData[(i - j + upperBw) * N + j];
    }

    // This function conducts banded LU factorization in place
    // Note that NO PIVOT is applied on the matrix "A" for efficiency!!!
    inline void factorizeLU()
    {
        int iM, jM;
        double cVl;
        for (int k = 0; k <= N - 2; k++)
        {
            iM = std::min(k + lowerBw, N - 1);
            cVl = operator()(k, k);
            for (int i = k + 1; i <= iM; i++)
            {
                if (operator()(i, k) != 0.0)
                {
                    operator()(i, k) /= cVl;
                }
            }
            jM = std::min(k + upperBw, N - 1);
            for (int j = k + 1; j <= jM; j++)
            {
                cVl = operator()(k, j);
                if (cVl != 0.0)
                {
                    for (int i = k + 1; i <= iM; i++)
                    {
                        if (operator()(i, k) != 0.0)
                        {
                            operator()(i, j) -= operator()(i, k) * cVl;
                        }
                    }
                }
            }
        }
        return;
    }

    // This function solves Ax=b, then stores x in b
    // The input b is required to be N*m, i.e.,
    // m vectors to be solved.
    template <typename EIGENMAT>
    inline void solve(EIGENMAT &b) const
    {
        int iM;
        for (int j = 0; j <= N - 1; j++)
        {
            iM = std::min(j + lowerBw, N - 1);
            for (int i = j + 1; i <= iM; i++)
            {
                if (operator()(i, j) != 0.0)
                {
                    b.row(i) -= operator()(i, j) * b.row(j);
                }
            }
        }
        for (int j = N - 1; j >= 0; j--)
        {
            b.row(j) /= operator()(j, j);
            iM = std::max(0, j - upperBw);
            for (int i = iM; i <= j - 1; i++)
            {
                if (operator()(i, j) != 0.0)
                {
                    b.row(i) -= operator()(i, j) * b.row(j);
                }
            }
        }
        return;
    }

    // This function solves ATx=b, then stores x in b
    // The input b is required to be N*m, i.e.,
    // m vectors to be solved.
    template <typename EIGENMAT>
    inline void solveAdj(EIGENMAT &b) const
    {
        int iM;
        for (int j = 0; j <= N - 1; j++)
        {
            b.row(j) /= operator()(j, j);
            iM = std::min(j + upperBw, N - 1);
            for (int i = j + 1; i <= iM; i++)
            {
                if (operator()(j, i) != 0.0)
                {
                    b.row(i) -= operator()(j, i) * b.row(j);
                }
            }
        }
        for (int j = N - 1; j >= 0; j--)
        {
            iM = std::max(0, j - lowerBw);
            for (int i = iM; i <= j - 1; i++)
            {
                if (operator()(j, i) != 0.0)
                {
                    b.row(i) -= operator()(j, i) * b.row(j);
                }
            }
        }
        return;
    }
};

class MinJerkOpt
{
public:
    MinJerkOpt() = default;
    ~MinJerkOpt() { A.destroy(); }

private:
    int N;
    Eigen::Matrix3d headPVA;
    Eigen::Matrix3d tailPVA;
    Eigen::VectorXd T1;
    BandedSystem A;
    Eigen::MatrixX3d b;

    // Temp variables
    Eigen::VectorXd T2;
    Eigen::VectorXd T3;
    Eigen::VectorXd T4;
    Eigen::VectorXd T5;
    Eigen::MatrixX3d gdC;
    double smoothEps;

private:
    inline void addGradJbyT(Eigen::VectorXd &gdT) const
    {
        for (int i = 0; i < N; i++)
        {
            gdT(i) += 36.0 * b.row(6 * i + 3).squaredNorm() +
                      288.0 * b.row(6 * i + 4).dot(b.row(6 * i + 3)) * T1(i) +
                      576.0 * b.row(6 * i + 4).squaredNorm() * T2(i) +
                      720.0 * b.row(6 * i + 5).dot(b.row(6 * i + 3)) * T2(i) +
                      2880.0 * b.row(6 * i + 5).dot(b.row(6 * i + 4)) * T3(i) +
                      3600.0 * b.row(6 * i + 5).squaredNorm() * T4(i);
        }
        return;
    }

    inline void addGradJbyC(Eigen::MatrixX3d &gdC) const
    {
        for (int i = 0; i < N; i++)
        {
            gdC.row(6 * i + 5) += 240.0 * b.row(6 * i + 3) * T3(i) +
                                  720.0 * b.row(6 * i + 4) * T4(i) +
                                  1440.0 * b.row(6 * i + 5) * T5(i);
            gdC.row(6 * i + 4) += 144.0 * b.row(6 * i + 3) * T2(i) +
                                  384.0 * b.row(6 * i + 4) * T3(i) +
                                  720.0 * b.row(6 * i + 5) * T4(i);
            gdC.row(6 * i + 3) += 72.0 * b.row(6 * i + 3) * T1(i) +
                                  144.0 * b.row(6 * i + 4) * T2(i) +
                                  240.0 * b.row(6 * i + 5) * T3(i);
        }
        return;
    }

    inline void solveAdjGradC(Eigen::MatrixX3d &gdC) const
    {
        A.solveAdj(gdC);
        return;
    }

    inline void addPropCtoT(const Eigen::MatrixX3d &adjGdC, Eigen::VectorXd &gdT) const
    {

        Eigen::Matrix<double, 6, 3> B1;
        Eigen::Matrix3d B2;

        for (int i = 0; i < N - 1; i++)
        {
            // negative velocity
            B1.row(2) = -(b.row(i * 6 + 1) +
                          2.0 * T1(i) * b.row(i * 6 + 2) +
                          3.0 * T2(i) * b.row(i * 6 + 3) +
                          4.0 * T3(i) * b.row(i * 6 + 4) +
                          5.0 * T4(i) * b.row(i * 6 + 5));
            B1.row(3) = B1.row(2);

            // negative acceleration
            B1.row(4) = -(2.0 * b.row(i * 6 + 2) +
                          6.0 * T1(i) * b.row(i * 6 + 3) +
                          12.0 * T2(i) * b.row(i * 6 + 4) +
                          20.0 * T3(i) * b.row(i * 6 + 5));

            // negative jerk
            B1.row(5) = -(6.0 * b.row(i * 6 + 3) +
                          24.0 * T1(i) * b.row(i * 6 + 4) +
                          60.0 * T2(i) * b.row(i * 6 + 5));

            // negative snap
            B1.row(0) = -(24.0 * b.row(i * 6 + 4) +
                          120.0 * T1(i) * b.row(i * 6 + 5));

            // negative crackle
            B1.row(1) = -120.0 * b.row(i * 6 + 5);

            gdT(i) += B1.cwiseProduct(adjGdC.block<6, 3>(6 * i + 3, 0)).sum();
        }

        // negative velocity
        B2.row(0) = -(b.row(6 * N - 5) +
                      2.0 * T1(N - 1) * b.row(6 * N - 4) +
                      3.0 * T2(N - 1) * b.row(6 * N - 3) +
                      4.0 * T3(N - 1) * b.row(6 * N - 2) +
                      5.0 * T4(N - 1) * b.row(6 * N - 1));

        // negative acceleration
        B2.row(1) = -(2.0 * b.row(6 * N - 4) +
                      6.0 * T1(N - 1) * b.row(6 * N - 3) +
                      12.0 * T2(N - 1) * b.row(6 * N - 2) +
                      20.0 * T3(N - 1) * b.row(6 * N - 1));

        // negative jerk
        B2.row(2) = -(6.0 * b.row(6 * N - 3) +
                      24.0 * T1(N - 1) * b.row(6 * N - 2) +
                      60.0 * T2(N - 1) * b.row(6 * N - 1));

        gdT(N - 1) += B2.cwiseProduct(adjGdC.block<3, 3>(6 * N - 3, 0)).sum();

        return;
    }

    inline void addPropCtoP(const Eigen::MatrixX3d &adjGdC, Eigen::Matrix3Xd &gdInP) const
    {
        for (int i = 0; i < N - 1; i++)
        {
            gdInP.col(i) += adjGdC.row(6 * i + 5).transpose();
        }
        return;
    }

    inline void positiveSmoothedL1(const double &x, double &f, double &df) const
    {
        const double pe = smoothEps;
        const double half = 0.5 * pe;
        const double f3c = 1.0 / (pe * pe);
        const double f4c = -0.5 * f3c / pe;
        const double d2c = 3.0 * f3c;
        const double d3c = 4.0 * f4c;

        if (x < pe)
        {
            f = (f4c * x + f3c) * x * x * x;
            df = (d3c * x + d2c) * x * x;
        }
        else
        {
            f = x - half;
            df = 1.0;
        }

        return;
    }

    inline void addTimeIntPenalty(const Eigen::VectorXi cons,
                                  const Eigen::VectorXi &idxHs,
                                  const std::vector<Eigen::Matrix<double, 6, -1>> &cfgHs,
                                  const double vmax,
                                  const double amax,
                                  const Eigen::Vector3d ci,
                                  double &cost,
                                  Eigen::VectorXd &gdT,
                                  Eigen::MatrixX3d &gdC) const
    {
        double pena = 0.0;
        const double vmaxSqr = vmax * vmax;
        const double amaxSqr = amax * amax;

        Eigen::Vector3d pos, vel, acc, jer;
        double step, alpha;
        double s1, s2, s3, s4, s5;
        Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3;
        Eigen::Vector3d outerNormal;
        int K;
        double violaPos, violaVel, violaAcc;
        double violaPosPenaD, violaVelPenaD, violaAccPenaD;
        double violaPosPena, violaVelPena, violaAccPena;
        Eigen::Matrix<double, 6, 3> gradViolaVc, gradViolaAc;
        double gradViolaVt, gradViolaAt;
        double omg;

        int innerLoop, idx;
        for (int i = 0; i < N; i++)
        {
            const auto &c = b.block<6, 3>(i * 6, 0);
            step = T1(i) / cons(i);
            s1 = 0.0;
            innerLoop = cons(i) + 1;
            for (int j = 0; j < innerLoop; j++)
            {
                s2 = s1 * s1;
                s3 = s2 * s1;
                s4 = s2 * s2;
                s5 = s4 * s1;
                beta0(0) = 1.0, beta0(1) = s1, beta0(2) = s2, beta0(3) = s3, beta0(4) = s4, beta0(5) = s5;
                beta1(0) = 0.0, beta1(1) = 1.0, beta1(2) = 2.0 * s1, beta1(3) = 3.0 * s2, beta1(4) = 4.0 * s3, beta1(5) = 5.0 * s4;
                beta2(0) = 0.0, beta2(1) = 0.0, beta2(2) = 2.0, beta2(3) = 6.0 * s1, beta2(4) = 12.0 * s2, beta2(5) = 20.0 * s3;
                beta3(0) = 0.0, beta3(1) = 0.0, beta3(2) = 0.0, beta3(3) = 6.0, beta3(4) = 24.0 * s1, beta3(5) = 60.0 * s2;
                alpha = 1.0 / cons(i) * j;
                pos = c.transpose() * beta0;
                vel = c.transpose() * beta1;
                acc = c.transpose() * beta2;
                jer = c.transpose() * beta3;
                violaVel = vel.squaredNorm() - vmaxSqr;
                violaAcc = acc.squaredNorm() - amaxSqr;

                omg = (j == 0 || j == innerLoop - 1) ? 0.5 : 1.0;

                idx = idxHs(i);
                K = cfgHs[idx].cols();
                for (int k = 0; k < K; k++)
                {
                    outerNormal = cfgHs[idx].col(k).head<3>();
                    violaPos = outerNormal.dot(pos - cfgHs[idx].col(k).tail<3>());
                    if (violaPos > 0.0)
                    {
                        positiveSmoothedL1(violaPos, violaPosPena, violaPosPenaD);
                        gdC.block<6, 3>(i * 6, 0) += omg * step * ci(0) * violaPosPenaD * beta0 * outerNormal.transpose();
                        gdT(i) += omg * (ci(0) * violaPosPenaD * alpha * outerNormal.dot(vel) * step +
                                         ci(0) * violaPosPena / cons(i));
                        pena += omg * step * ci(0) * violaPosPena;
                    }
                }

                if (violaVel > 0.0)
                {
                    positiveSmoothedL1(violaVel, violaVelPena, violaVelPenaD);
                    gradViolaVc = 2.0 * beta1 * vel.transpose();
                    gradViolaVt = 2.0 * alpha * vel.transpose() * acc;
                    gdC.block<6, 3>(i * 6, 0) += omg * step * ci(1) * violaVelPenaD * gradViolaVc;
                    gdT(i) += omg * (ci(1) * violaVelPenaD * gradViolaVt * step +
                                     ci(1) * violaVelPena / cons(i));
                    pena += omg * step * ci(1) * violaVelPena;
                }

                if (violaAcc > 0.0)
                {
                    positiveSmoothedL1(violaAcc, violaAccPena, violaAccPenaD);
                    gradViolaAc = 2.0 * beta2 * acc.transpose();
                    gradViolaAt = 2.0 * alpha * acc.transpose() * jer;
                    gdC.block<6, 3>(i * 6, 0) += omg * step * ci(2) * violaAccPenaD * gradViolaAc;
                    gdT(i) += omg * (ci(2) * violaAccPenaD * gradViolaAt * step +
                                     ci(2) * violaAccPena / cons(i));
                    pena += omg * step * ci(2) * violaAccPena;
                }

                s1 += step;
            }
        }

        cost += pena;
        return;
    }

public:
    inline void reset(const Eigen::Matrix3d &headState,
                      const Eigen::Matrix3d &tailState,
                      const double &smoEps,
                      const int &pieceNum)
    {
        N = pieceNum;
        headPVA = headState;
        tailPVA = tailState;
        T1.resize(N);
        A.create(6 * N, 6, 6);
        b.resize(6 * N, 3);
        gdC.resize(6 * N, 3);
        smoothEps = smoEps;
        return;
    }

    inline void generate(const Eigen::Matrix3Xd &inPs,
                         const Eigen::VectorXd &ts)
    {
        T1 = ts;
        T2 = T1.cwiseProduct(T1);
        T3 = T2.cwiseProduct(T1);
        T4 = T2.cwiseProduct(T2);
        T5 = T4.cwiseProduct(T1);

        A.reset();
        b.setZero();

        A(0, 0) = 1.0;
        A(1, 1) = 1.0;
        A(2, 2) = 2.0;
        b.row(0) = headPVA.col(0).transpose();
        b.row(1) = headPVA.col(1).transpose();
        b.row(2) = headPVA.col(2).transpose();

        for (int i = 0; i < N - 1; i++)
        {
            A(6 * i + 3, 6 * i + 3) = 6.0;
            A(6 * i + 3, 6 * i + 4) = 24.0 * T1(i);
            A(6 * i + 3, 6 * i + 5) = 60.0 * T2(i);
            A(6 * i + 3, 6 * i + 9) = -6.0;
            A(6 * i + 4, 6 * i + 4) = 24.0;
            A(6 * i + 4, 6 * i + 5) = 120.0 * T1(i);
            A(6 * i + 4, 6 * i + 10) = -24.0;
            A(6 * i + 5, 6 * i) = 1.0;
            A(6 * i + 5, 6 * i + 1) = T1(i);
            A(6 * i + 5, 6 * i + 2) = T2(i);
            A(6 * i + 5, 6 * i + 3) = T3(i);
            A(6 * i + 5, 6 * i + 4) = T4(i);
            A(6 * i + 5, 6 * i + 5) = T5(i);
            A(6 * i + 6, 6 * i) = 1.0;
            A(6 * i + 6, 6 * i + 1) = T1(i);
            A(6 * i + 6, 6 * i + 2) = T2(i);
            A(6 * i + 6, 6 * i + 3) = T3(i);
            A(6 * i + 6, 6 * i + 4) = T4(i);
            A(6 * i + 6, 6 * i + 5) = T5(i);
            A(6 * i + 6, 6 * i + 6) = -1.0;
            A(6 * i + 7, 6 * i + 1) = 1.0;
            A(6 * i + 7, 6 * i + 2) = 2 * T1(i);
            A(6 * i + 7, 6 * i + 3) = 3 * T2(i);
            A(6 * i + 7, 6 * i + 4) = 4 * T3(i);
            A(6 * i + 7, 6 * i + 5) = 5 * T4(i);
            A(6 * i + 7, 6 * i + 7) = -1.0;
            A(6 * i + 8, 6 * i + 2) = 2.0;
            A(6 * i + 8, 6 * i + 3) = 6 * T1(i);
            A(6 * i + 8, 6 * i + 4) = 12 * T2(i);
            A(6 * i + 8, 6 * i + 5) = 20 * T3(i);
            A(6 * i + 8, 6 * i + 8) = -2.0;

            b.row(6 * i + 5) = inPs.col(i).transpose();
        }

        A(6 * N - 3, 6 * N - 6) = 1.0;
        A(6 * N - 3, 6 * N - 5) = T1(N - 1);
        A(6 * N - 3, 6 * N - 4) = T2(N - 1);
        A(6 * N - 3, 6 * N - 3) = T3(N - 1);
        A(6 * N - 3, 6 * N - 2) = T4(N - 1);
        A(6 * N - 3, 6 * N - 1) = T5(N - 1);
        A(6 * N - 2, 6 * N - 5) = 1.0;
        A(6 * N - 2, 6 * N - 4) = 2 * T1(N - 1);
        A(6 * N - 2, 6 * N - 3) = 3 * T2(N - 1);
        A(6 * N - 2, 6 * N - 2) = 4 * T3(N - 1);
        A(6 * N - 2, 6 * N - 1) = 5 * T4(N - 1);
        A(6 * N - 1, 6 * N - 4) = 2;
        A(6 * N - 1, 6 * N - 3) = 6 * T1(N - 1);
        A(6 * N - 1, 6 * N - 2) = 12 * T2(N - 1);
        A(6 * N - 1, 6 * N - 1) = 20 * T3(N - 1);

        b.row(6 * N - 3) = tailPVA.col(0).transpose();
        b.row(6 * N - 2) = tailPVA.col(1).transpose();
        b.row(6 * N - 1) = tailPVA.col(2).transpose();

        A.factorizeLU();
        A.solve(b);

        return;
    }

    inline double getTrajJerkCost() const
    {
        double objective = 0.0;
        for (int i = 0; i < N; i++)
        {
            objective += 36.0 * b.row(6 * i + 3).squaredNorm() * T1(i) +
                         144.0 * b.row(6 * i + 4).dot(b.row(6 * i + 3)) * T2(i) +
                         192.0 * b.row(6 * i + 4).squaredNorm() * T3(i) +
                         240.0 * b.row(6 * i + 5).dot(b.row(6 * i + 3)) * T3(i) +
                         720.0 * b.row(6 * i + 5).dot(b.row(6 * i + 4)) * T4(i) +
                         720.0 * b.row(6 * i + 5).squaredNorm() * T5(i);
        }
        return objective;
    }

    inline void evalTrajCostGrad(const Eigen::VectorXi &cons,
                                 const Eigen::VectorXi &idxHs,
                                 const std::vector<Eigen::Matrix<double, 6, -1>> &cfgHs,
                                 const double &vmax,
                                 const double &amax,
                                 const Eigen::Vector3d &ci,
                                 double &cost,
                                 Eigen::VectorXd &gdT,
                                 Eigen::Matrix3Xd &gdInPs)
    {
        gdT.setZero();
        gdInPs.setZero();
        gdC.setZero();

        cost = getTrajJerkCost();
        addGradJbyT(gdT);
        addGradJbyC(gdC);

        addTimeIntPenalty(cons, idxHs, cfgHs, vmax, amax, ci, cost, gdT, gdC);

        solveAdjGradC(gdC);
        addPropCtoT(gdC, gdT);
        addPropCtoP(gdC, gdInPs);
    }

    inline Trajectory getTraj(void) const
    {
        Trajectory traj;
        traj.reserve(N);
        for (int i = 0; i < N; i++)
        {
            traj.emplace_back(T1(i), b.block<6, 3>(6 * i, 0).transpose().rowwise().reverse());
        }
        return traj;
    }
};

class GCOPTER
{
private:
    // Use C2 or Cinf diffeo
    bool c2dfm;

    // Use retraction for hyperball
    bool retraction;

    // Use soft time or not
    bool softT;

    // Weight for time regularization term
    double rho;

    // Fixed total time
    double sumT;

    //Minimum Jerk Optimizer
    MinJerkOpt jerkOpt;

    // Temp variables for problem solving
    Eigen::Matrix3d iState;
    Eigen::Matrix3d fState;

    // Each col of cfgHs denotes a facet (outter_normal^T,point^T)^T
    std::vector<Eigen::Matrix3Xd> cfgVs;
    std::vector<Eigen::Matrix<double, 6, -1>> cfgHs;
    Eigen::Matrix3Xd gdInPs;

    // Piece num for each polytope
    Eigen::VectorXi intervals;
    // Assignment vector for point in V-polytope
    Eigen::VectorXi idxVs;
    // Assignment vector for piece in H-polytope
    Eigen::VectorXi idxHs;

    int coarseN;
    int fineN;
    int dimFreeT;
    int dimFreeP;
    Eigen::VectorXd coarseT;
    Eigen::VectorXd fineT;
    Eigen::Matrix3Xd innerP;

    // Params for constraints
    Eigen::VectorXi cons;
    double smoothEps;
    Eigen::Vector3d chi;
    double vmax;
    double amax;

    // L-BFGS Solver Parameters
    lbfgs::lbfgs_parameter_t lbfgs_params;

private:
    template <typename EIGENVEC>
    static inline void forwardT(const EIGENVEC &t,
                                Eigen::VectorXd &vecT,
                                bool soft,
                                const double &sT,
                                bool c2)
    {
        if (soft)
        {
            if (c2)
            {
                int M = vecT.size();
                for (int i = 0; i < M; i++)
                {
                    vecT(i) = t(i) > 0.0
                                  ? ((0.5 * t(i) + 1.0) * t(i) + 1.0)
                                  : 1.0 / ((0.5 * t(i) - 1.0) * t(i) + 1.0);
                }
            }
            else
            {
                vecT = t.array().exp();
            }
        }
        else
        {
            if (c2)
            {
                int Ms1 = t.size();
                for (int i = 0; i < Ms1; i++)
                {
                    vecT(i) = t(i) > 0.0
                                  ? ((0.5 * t(i) + 1.0) * t(i) + 1.0)
                                  : 1.0 / ((0.5 * t(i) - 1.0) * t(i) + 1.0);
                }
                vecT(Ms1) = 0.0;
                vecT /= 1.0 + vecT.sum();
                vecT(Ms1) = 1.0 - vecT.sum();
                vecT *= sT;
            }
            else
            {
                int Ms1 = t.size();
                vecT.head(Ms1) = t.array().exp();
                vecT(Ms1) = 0.0;
                vecT /= 1.0 + vecT.sum();
                vecT(Ms1) = 1.0 - vecT.sum();
                vecT *= sT;
            }
        }
        return;
    }

    template <typename EIGENVEC>
    static inline void backwardT(const Eigen::VectorXd &vecT,
                                 EIGENVEC &t,
                                 bool soft,
                                 bool c2)
    {
        if (soft)
        {
            if (c2)
            {
                int M = vecT.size();
                for (int i = 0; i < M; i++)
                {
                    t(i) = vecT(i) > 1.0
                               ? (sqrt(2.0 * vecT(i) - 1.0) - 1.0)
                               : (1.0 - sqrt(2.0 / vecT(i) - 1.0));
                }
            }
            else
            {
                t = vecT.array().log();
            }
        }
        else
        {
            if (c2)
            {
                int Ms1 = t.size();
                t = vecT.head(Ms1) / vecT(Ms1);
                for (int i = 0; i < Ms1; i++)
                {
                    t(i) = t(i) > 1.0
                               ? (sqrt(2.0 * t(i) - 1.0) - 1.0)
                               : (1.0 - sqrt(2.0 / t(i) - 1.0));
                }
            }
            else
            {
                int Ms1 = t.size();
                t = (vecT.head(Ms1) / vecT(Ms1)).array().log();
            }
        }
        return;
    }

    template <typename EIGENVEC>
    static inline void forwardP(const EIGENVEC &p,
                                const Eigen::VectorXi &idVs,
                                const std::vector<Eigen::Matrix3Xd> &cfgPolyVs,
                                Eigen::Matrix3Xd &inP,
                                bool retract)
    {
        if (!retract)
        {
            int M = inP.cols();
            Eigen::VectorXd q;
            int j = 0, k, idx;
            for (int i = 0; i < M; i++)
            {
                idx = idVs(i);
                k = cfgPolyVs[idx].cols() - 1;
                q = 2.0 / (1.0 + p.segment(j, k).squaredNorm()) * p.segment(j, k);
                inP.col(i) = cfgPolyVs[idx].rightCols(k) * q.cwiseProduct(q) +
                             cfgPolyVs[idx].col(0);
                j += k;
            }
        }
        else
        {
            int M = inP.cols();
            Eigen::VectorXd q;
            double normInv;
            int j = 0, k, idx;
            for (int i = 0; i < M; i++)
            {
                idx = idVs(i);
                k = cfgPolyVs[idx].cols() - 1;
                normInv = 1.0 / p.segment(j, k + 1).norm();
                q = p.segment(j, k + 1).head(k) * normInv;
                inP.col(i) = cfgPolyVs[idx].rightCols(k) * q.cwiseProduct(q) +
                             cfgPolyVs[idx].col(0);
                j += k + 1;
            }
        }
        return;
    }

    static inline double objectiveNLS(void *ptrPOBs,
                                      const double *x,
                                      double *grad,
                                      const int n)
    {
        const Eigen::Matrix3Xd &pobs = *(Eigen::Matrix3Xd *)ptrPOBs;
        Eigen::Map<const Eigen::VectorXd> p(x, n);
        Eigen::Map<Eigen::VectorXd> gradp(grad, n);

        double qnsqr = p.squaredNorm();
        double qnsqrp1 = qnsqr + 1.0;
        double qnsqrp1sqr = qnsqrp1 * qnsqrp1;
        Eigen::VectorXd r = 2.0 / qnsqrp1 * p;

        Eigen::Vector3d delta = pobs.rightCols(n) * r.cwiseProduct(r) +
                                pobs.col(1) - pobs.col(0);
        double cost = delta.squaredNorm();
        Eigen::Vector3d gradR3 = 2 * delta;

        Eigen::VectorXd gdr = pobs.rightCols(n).transpose() * gradR3;
        gdr = gdr.array() * r.array() * 2.0;
        gradp = gdr * 2.0 / qnsqrp1 -
                p * 4.0 * gdr.dot(p) / qnsqrp1sqr;

        return cost;
    }

    static inline double objectiveRetractNLS(void *ptrPOBs,
                                             const double *x,
                                             double *grad,
                                             const int n)
    {
        const Eigen::Matrix3Xd &pobs = *(Eigen::Matrix3Xd *)ptrPOBs;
        Eigen::Map<const Eigen::VectorXd> p(x, n);
        Eigen::Map<Eigen::VectorXd> gradp(grad, n);

        double sqrnormp = p.squaredNorm();
        double normInv = 1.0 / sqrt(sqrnormp);
        Eigen::VectorXd unitp = p * normInv;
        Eigen::VectorXd r = unitp.head(n - 1);

        Eigen::Vector3d delta = pobs.rightCols(n - 1) * r.cwiseProduct(r) +
                                pobs.col(1) - pobs.col(0);
        double cost = delta.squaredNorm();
        Eigen::Vector3d gradR3 = 2 * delta;

        Eigen::VectorXd gdr = pobs.rightCols(n - 1).transpose() * gradR3;
        gdr = gdr.array() * r.array() * 2.0;

        gradp.head(n - 1) = gdr;
        gradp(n - 1) = 0.0;
        gradp = (gradp - unitp * unitp.dot(gradp)).eval() * normInv;

        const double normViola = sqrnormp - 1.0;
        if (normViola > 0.0)
        {
            double c = normViola * normViola;
            const double dc = 3.0 * c;
            c *= normViola;
            cost += c;
            gradp += dc * 2.0 * p;
        }

        return cost;
    }

    template <typename EIGENVEC>
    static inline void backwardP(const Eigen::Matrix3Xd &inP,
                                 const Eigen::VectorXi &idVs,
                                 const std::vector<Eigen::Matrix3Xd> &cfgPolyVs,
                                 EIGENVEC &p,
                                 bool retract)
    {
        int M = inP.cols();
        int j = 0, k, idx;

        // Parameters for tiny nonlinear least squares
        double minSqrD;
        lbfgs::lbfgs_parameter_t nls_params;
        lbfgs::lbfgs_load_default_parameters(&nls_params);
        nls_params.g_epsilon = FLT_EPSILON;
        nls_params.max_iterations = 128;

        Eigen::Matrix3Xd pobs;
        if (!retract)
        {
            for (int i = 0; i < M; i++)
            {
                idx = idVs(i);
                k = cfgPolyVs[idx].cols() - 1;
                p.segment(j, k).setConstant(1.0 / (sqrt(k + 1.0) + 1.0));
                pobs.resize(3, k + 2);
                pobs.col(0) = inP.col(i);
                pobs.rightCols(k + 1) = cfgPolyVs[idx];
                lbfgs::lbfgs_optimize(k,
                                      p.data() + j,
                                      &minSqrD,
                                      &GCOPTER::objectiveNLS,
                                      nullptr,
                                      nullptr,
                                      &pobs,
                                      &nls_params);

                j += k;
            }
        }
        else
        {
            for (int i = 0; i < M; i++)
            {
                idx = idVs(i);
                k = cfgPolyVs[idx].cols() - 1;
                p.segment(j, k + 1).setConstant(1.0 / sqrt(k + 1.0));
                pobs.resize(3, k + 2);
                pobs.col(0) = inP.col(i);
                pobs.rightCols(k + 1) = cfgPolyVs[idx];
                lbfgs::lbfgs_optimize(k + 1,
                                      p.data() + j,
                                      &minSqrD,
                                      &GCOPTER::objectiveRetractNLS,
                                      nullptr,
                                      nullptr,
                                      &pobs,
                                      &nls_params);

                j += k + 1;
            }
        }

        return;
    }

    static inline void addLayerTGrad(const Eigen::VectorXd &t,
                                     Eigen::VectorXd &gradT,
                                     bool soft,
                                     const double &sT,
                                     bool c2)
    {
        if (soft)
        {
            if (c2)
            {
                int M = t.size();
                double denSqrt;
                for (int i = 0; i < M; i++)
                {
                    if (t(i) > 0)
                    {
                        gradT(i) *= t(i) + 1.0;
                    }
                    else
                    {
                        denSqrt = (0.5 * t(i) - 1.0) * t(i) + 1.0;
                        gradT(i) *= (1.0 - t(i)) / (denSqrt * denSqrt);
                    }
                }
            }
            else
            {
                int M = t.size();
                gradT.head(M).array() *= t.array().exp();
            }
        }
        else
        {
            if (c2)
            {
                int Ms1 = t.size();
                Eigen::VectorXd gFree = sT * gradT.head(Ms1);
                double gTail = sT * gradT(Ms1);
                Eigen::VectorXd dExpTau(Ms1);
                double expTauSum = 0.0, gFreeDotExpTau = 0.0;
                double denSqrt, expTau;
                for (int i = 0; i < Ms1; i++)
                {
                    if (t(i) > 0)
                    {
                        expTau = (0.5 * t(i) + 1.0) * t(i) + 1.0;
                        dExpTau(i) = t(i) + 1.0;
                        expTauSum += expTau;
                        gFreeDotExpTau += expTau * gFree(i);
                    }
                    else
                    {
                        denSqrt = (0.5 * t(i) - 1.0) * t(i) + 1.0;
                        expTau = 1.0 / denSqrt;
                        dExpTau(i) = (1.0 - t(i)) / (denSqrt * denSqrt);
                        expTauSum += expTau;
                        gFreeDotExpTau += expTau * gFree(i);
                    }
                }
                denSqrt = expTauSum + 1.0;
                gradT.head(Ms1) = (gFree.array() - gTail) * dExpTau.array() / denSqrt -
                                  (gFreeDotExpTau - gTail * expTauSum) * dExpTau.array() / (denSqrt * denSqrt);
                gradT(Ms1) = 0.0;
            }
            else
            {
                int Ms1 = t.size();
                Eigen::VectorXd gFree = sT * gradT.head(Ms1);
                double gTail = sT * gradT(Ms1);
                Eigen::VectorXd expTau = t.array().exp();
                double expTauSum = expTau.sum();
                double denom = expTauSum + 1.0;
                gradT.head(Ms1) = (gFree.array() - gTail) * expTau.array() / denom -
                                  (gFree.dot(expTau) - gTail * expTauSum) * expTau.array() / (denom * denom);
                gradT(Ms1) = 0.0;
            }
        }
        return;
    }

    template <typename EIGENVEC>
    static inline void addLayerPGrad(const Eigen::VectorXd &p,
                                     const Eigen::VectorXi &idVs,
                                     const std::vector<Eigen::Matrix3Xd> &cfgPolyVs,
                                     const Eigen::Matrix3Xd &gradInPs,
                                     EIGENVEC &grad,
                                     bool retract)
    {
        int M = gradInPs.cols();

        if (!retract)
        {
            int j = 0, k, idx;
            double qnsqr, qnsqrp1, qnsqrp1sqr;
            Eigen::VectorXd q, r, gdr;
            for (int i = 0; i < M; i++)
            {
                idx = idVs(i);
                k = cfgPolyVs[idx].cols() - 1;

                q = p.segment(j, k);
                qnsqr = q.squaredNorm();
                qnsqrp1 = qnsqr + 1.0;
                qnsqrp1sqr = qnsqrp1 * qnsqrp1;
                r = 2.0 / qnsqrp1 * q;
                gdr = cfgPolyVs[idx].rightCols(k).transpose() * gradInPs.col(i);
                gdr = gdr.array() * r.array() * 2.0;

                grad.segment(j, k) = gdr * 2.0 / qnsqrp1 -
                                     q * 4.0 * gdr.dot(q) / qnsqrp1sqr;

                j += k;
            }
        }
        else
        {
            int j = 0, k, idx;
            double normInv;
            Eigen::VectorXd q, r, gdr, gradq, unitq;
            for (int i = 0; i < M; i++)
            {
                idx = idVs(i);
                k = cfgPolyVs[idx].cols() - 1;

                q = p.segment(j, k + 1);
                normInv = 1.0 / q.norm();
                unitq = q * normInv;
                r = unitq.head(k);
                gdr = cfgPolyVs[idx].rightCols(k).transpose() * gradInPs.col(i);
                gdr = gdr.array() * r.array() * 2.0;

                gradq.resize(k + 1);
                gradq.head(k) = gdr;
                gradq(k) = 0.0;
                grad.segment(j, k + 1) = (gradq - unitq * unitq.dot(gradq)) * normInv;

                j += k + 1;
            }
        }

        return;
    }

    template <typename EIGENVEC>
    static inline void restrictNorm(const Eigen::VectorXd &p,
                                    const Eigen::VectorXi &idVs,
                                    const std::vector<Eigen::Matrix3Xd> &cfgPolyVs,
                                    double &cost,
                                    EIGENVEC &grad)
    {
        int M = idVs.size();

        int j = 0, k, idx;
        double sqrnormq, normViola, c, dc;
        Eigen::VectorXd q;
        for (int i = 0; i < M; i++)
        {
            idx = idVs(i);
            k = cfgPolyVs[idx].cols() - 1;

            q = p.segment(j, k + 1);
            sqrnormq = q.squaredNorm();
            normViola = sqrnormq - 1.0;
            if (normViola > 0.0)
            {
                c = normViola * normViola;
                dc = 3.0 * c;
                c *= normViola;
                cost += c;
                grad.segment(j, k + 1) += dc * 2.0 * q;
            }

            j += k + 1;
        }

        return;
    }

    static inline void splitToFineT(const Eigen::VectorXd &cT,
                                    const Eigen::VectorXi &intervs,
                                    Eigen::VectorXd &fT)
    {
        int M = intervs.size();
        int offset = 0;
        int inverv;
        for (int i = 0; i < M; i++)
        {
            inverv = intervs(i);
            fT.segment(offset, inverv).setConstant(cT(i) / inverv);
            offset += inverv;
        }
        return;
    }

    static inline void mergeToCoarseGradT(const Eigen::VectorXi &intervs,
                                          Eigen::VectorXd &fineGdT)
    {
        int M = intervs.size();
        int offset = 0;
        int inverv;
        for (int i = 0; i < M; i++)
        {
            inverv = intervs(i);
            fineGdT(i) = fineGdT.segment(offset, inverv).mean();
            offset += inverv;
        }
        return;
    }

    static inline double objectiveFunc(void *ptrObj,
                                       const double *x,
                                       double *grad,
                                       const int n)
    {
        GCOPTER &obj = *(GCOPTER *)ptrObj;
        const int dimT = obj.dimFreeT;
        const int dimP = obj.dimFreeP;
        const double rh = obj.rho;
        Eigen::Map<const Eigen::VectorXd> t(x, dimT);
        Eigen::Map<const Eigen::VectorXd> p(x + dimT, dimP);
        Eigen::Map<Eigen::VectorXd> gradt(grad, dimT);
        Eigen::VectorXd proxyGradT(obj.fineN);
        Eigen::Map<Eigen::VectorXd> gradp(grad + dimT, dimP);

        forwardT(t, obj.coarseT, obj.softT, obj.sumT, obj.c2dfm);
        splitToFineT(obj.coarseT, obj.intervals, obj.fineT);
        forwardP(p, obj.idxVs, obj.cfgVs, obj.innerP, obj.retraction);

        double cost;

        obj.jerkOpt.generate(obj.innerP, obj.fineT);
        obj.jerkOpt.evalTrajCostGrad(obj.cons,
                                     obj.idxHs, obj.cfgHs,
                                     obj.vmax, obj.amax,
                                     obj.chi, cost,
                                     proxyGradT, obj.gdInPs);

        cost += rh * obj.coarseT.sum();
        proxyGradT.array() += rh;

        mergeToCoarseGradT(obj.intervals, proxyGradT);
        addLayerTGrad(t, proxyGradT, obj.softT, obj.sumT, obj.c2dfm);
        addLayerPGrad(p, obj.idxVs, obj.cfgVs, obj.gdInPs, gradp, obj.retraction);
        if (obj.retraction)
        {
            restrictNorm(p, obj.idxVs, obj.cfgVs, cost, gradp);
        }

        gradt = proxyGradT.head(dimT);

        return cost;
    }

public:
    inline void gridMesh(const Eigen::Matrix3d &iState,
                         const Eigen::Matrix3d &fState,
                         const std::vector<Eigen::Matrix3Xd> &cfgPolyVs,
                         const double &gridResolution,
                         Eigen::VectorXi &intervalsVec) const
    {
        int M = intervalsVec.size();

        int curInterval, k;
        Eigen::Vector3d lastP, curP;
        curP = iState.col(0);
        for (int i = 0; i < M - 1; i++)
        {
            lastP = curP;
            k = cfgPolyVs[2 * i + 1].cols() - 1;
            curP = cfgPolyVs[2 * i + 1].rightCols(k).rowwise().sum() / (1.0 + k) +
                   cfgPolyVs[2 * i + 1].col(0);
            curInterval = ceil((curP - lastP).norm() / gridResolution);
            intervalsVec(i) = curInterval > 0 ? curInterval : 1;
        }
        lastP = curP;
        curP = fState.col(0);
        curInterval = ceil((curP - lastP).norm() / gridResolution);
        intervalsVec(M - 1) = curInterval > 0 ? curInterval : 1;

        return;
    }

    inline bool extractVs(const std::vector<Eigen::Matrix<double, 6, -1>> &hPs,
                          std::vector<Eigen::Matrix3Xd> &vPs) const
    {
        const int M = hPs.size() - 1;

        vPs.clear();
        vPs.reserve(2 * M + 1);

        int nv;
        Eigen::Matrix<double, 6, -1> curIH;
        Eigen::Matrix3Xd curIV, curIOB;
        for (int i = 0; i < M; i++)
        {
            if (!geoutils::enumerateVs(hPs[i], curIV))
            {
                return false;
            }
            nv = curIV.cols();
            curIOB.resize(3, nv);
            curIOB.col(0) = curIV.col(0);
            curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
            vPs.push_back(curIOB);

            curIH.resize(6, hPs[i].cols() + hPs[i + 1].cols());
            curIH.leftCols(hPs[i].cols()) = hPs[i];
            curIH.rightCols(hPs[i + 1].cols()) = hPs[i + 1];
            if (!geoutils::enumerateVs(curIH, curIV))
            {
                return false;
            }
            nv = curIV.cols();
            curIOB.resize(3, nv);
            curIOB.col(0) = curIV.col(0);
            curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
            vPs.push_back(curIOB);
        }

        if (!geoutils::enumerateVs(hPs.back(), curIV))
        {
            return false;
        }
        nv = curIV.cols();
        curIOB.resize(3, nv);
        curIOB.col(0) = curIV.col(0);
        curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
        vPs.push_back(curIOB);

        return true;
    }

    inline bool setup(const double &rh,
                      const double &st,
                      const Eigen::Matrix3d &iniState,
                      const Eigen::Matrix3d &finState,
                      const std::vector<Eigen::Matrix<double, 6, -1>> &cfgPolyHs,
                      const double &gridRes,
                      const int &itgSpaces,
                      const double &vm,
                      const double &am,
                      const double &smoEps,
                      const Eigen::Vector3d &w,
                      bool c2diffeo,
                      bool retract)
    {
        // Setup for optimization parameters
        c2dfm = c2diffeo;

        retraction = retract;

        softT = rh > 0;
        if (softT)
        {
            rho = rh;
            sumT = 1.0;
        }
        else
        {
            rho = 0.0;
            sumT = st;
        }

        iState = iniState;
        fState = finState;

        cfgHs = cfgPolyHs;
        coarseN = cfgHs.size();
        for (int i = 0; i < coarseN; i++)
        {
            cfgHs[i].topRows<3>().colwise().normalize();
        }
        if (!extractVs(cfgHs, cfgVs))
        {
            return false;
        }

        intervals.resize(coarseN);
        gridMesh(iState, fState, cfgVs, gridRes, intervals);
        fineN = intervals.sum();
        cons.resize(fineN);
        cons.setConstant(itgSpaces);

        idxVs.resize(fineN - 1);
        idxHs.resize(fineN);
        dimFreeT = softT ? coarseN : coarseN - 1;
        dimFreeP = 0;
        int offset = 0, interval;
        for (int i = 0; i < coarseN; i++)
        {
            interval = intervals(i);
            for (int j = 0; j < interval; j++)
            {
                if (j < interval - 1)
                {
                    idxVs(offset) = 2 * i;
                    dimFreeP += cfgVs[2 * i].cols() - 1;
                    if (retraction)
                    {
                        dimFreeP++;
                    }
                }
                else if (i < coarseN - 1)
                {
                    idxVs(offset) = 2 * i + 1;
                    dimFreeP += cfgVs[2 * i + 1].cols() - 1;
                    if (retraction)
                    {
                        dimFreeP++;
                    }
                }
                idxHs(offset) = i;
                offset++;
            }
        }

        smoothEps = smoEps;
        chi = w;
        vmax = vm;
        amax = am;

        // Make all conditions legal
        double tempNorm;
        tempNorm = iState.col(1).norm();
        iState.col(1) *= tempNorm > vmax ? (vmax / tempNorm) : 1.0;
        tempNorm = fState.col(1).norm();
        fState.col(1) *= tempNorm > vmax ? (vmax / tempNorm) : 1.0;
        tempNorm = iState.col(2).norm();
        iState.col(2) *= tempNorm > amax ? (amax / tempNorm) : 1.0;
        tempNorm = fState.col(2).norm();
        fState.col(2) *= tempNorm > amax ? (amax / tempNorm) : 1.0;

        // Setup for L-BFGS solver
        lbfgs::lbfgs_load_default_parameters(&lbfgs_params);

        // Allocate temp variables
        coarseT.resize(coarseN);
        fineT.resize(fineN);
        innerP.resize(3, fineN - 1);
        gdInPs.resize(3, fineN - 1);
        jerkOpt.reset(iniState, finState, smoothEps, fineN);

        return true;
    }

    inline void setInitial(const std::vector<Eigen::Matrix3Xd> &cfgPolyVs,
                           const Eigen::VectorXi &intervs,
                           Eigen::VectorXd &vecT,
                           Eigen::Matrix3Xd &vecInP) const
    {
        constexpr double maxSpeedForAllocation = 10.0;
        const double allocationSpeed = std::min(vmax, maxSpeedForAllocation);

        int M = vecT.size();
        Eigen::Vector3d lastP, curP, delta;
        int offset, interv, k;

        offset = 0;
        curP = iState.col(0);
        for (int i = 0; i < M - 1; i++)
        {
            lastP = curP;
            interv = intervs(i);
            k = cfgPolyVs[2 * i + 1].cols() - 1;
            curP = cfgPolyVs[2 * i + 1].rightCols(k).rowwise().sum() / (1.0 + k) +
                   cfgPolyVs[2 * i + 1].col(0);
            delta = curP - lastP;
            vecT(i) = delta.norm() / allocationSpeed;
            delta /= interv;
            for (int j = 0; j < interv; j++)
            {
                vecInP.col(offset++) = (j + 1) * delta + lastP;
            }
        }
        interv = intervs(M - 1);
        lastP = curP;
        curP = fState.col(0);
        delta = curP - lastP;
        vecT(M - 1) = delta.norm() / allocationSpeed;
        delta /= interv;
        for (int j = 0; j < interv - 1; j++)
        {
            vecInP.col(offset++) = (j + 1) * delta + lastP;
        }

        return;
    }

    inline double optimize(Trajectory &traj,
                           const double &relCostTol)
    {
        double *x = new double[dimFreeT + dimFreeP];
        Eigen::Map<Eigen::VectorXd> t(x, dimFreeT);
        Eigen::Map<Eigen::VectorXd> p(x + dimFreeT, dimFreeP);

        setInitial(cfgVs, intervals, coarseT, innerP);

        backwardT(coarseT, t, softT, c2dfm);
        backwardP(innerP, idxVs, cfgVs, p, retraction);

        double minObjectivePenalty;
        lbfgs_params.mem_size = 64;
        lbfgs_params.past = 3;
        lbfgs_params.g_epsilon = 1.0e-16;
        lbfgs_params.min_step = 1.0e-32;
        lbfgs_params.abs_curv_cond = 0;
        lbfgs_params.delta = relCostTol;

        int retCode = lbfgs::lbfgs_optimize(dimFreeT + dimFreeP,
                                            x,
                                            &minObjectivePenalty,
                                            &GCOPTER::objectiveFunc,
                                            nullptr,
                                            nullptr,
                                            this,
                                            &lbfgs_params);

        std::cout << "L-BFGS: " << lbfgs::lbfgs_strerror(retCode) << std::endl;

        forwardT(t, coarseT, softT, sumT, c2dfm);
        splitToFineT(coarseT, intervals, fineT);
        forwardP(p, idxVs, cfgVs, innerP, retraction);

        jerkOpt.generate(innerP, fineT);
        traj = jerkOpt.getTraj();

        delete[] x;
        return jerkOpt.getTrajJerkCost();
    }
};

#endif
