/**
 * @file corridor_minisnap.h
 * @author Siyuan Wu (siyuanwu99@gmail.com)
 * @brief
 * @version 1.0
 * @date 2021-12-22
 *
 * @copyright Copyright (c) 2021
 *
 */

#if !defined(CORRIDOR_MINI_SNAP_H_) x
#define CORRIDOR_MINI_SNAP_H_
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <iosqp.hpp>
#include <iostream>
#include <vector>

#define N_ORDER 7  // order of polynomial trajectory
#define DIM 3      // number of dimensions in Cspace

namespace trajopt {

/**
 * @brief A piece of polynomial trajectory
 * use relative time t_ = t / T
 */
class PolyPiece {
 private:
  double _duration;
  Eigen::Matrix<double, 3, 8> _coeffs;

 public:
  PolyPiece() {}
  ~PolyPiece() {}
  inline Eigen::Vector3d getPos(double t);
  inline Eigen::Vector3d getVel(double t);
  inline Eigen::Vector3d getAcc(double t);
  inline double getDuration();
};

class Trajectory {
 private:
  typedef std::vector<PolyPiece> Pieces;
  int N;
  Pieces _pieces;

 public:
  Trajectory() {}
  ~Trajectory() {}
  void setDuration(const std::vector<double> &t);
  void setTrajectory(const Eigen::VectorXd &x);
  inline void locatePiece(const double &t0, double &t, int &idx);
  inline Eigen::Vector3d getPos(double t);
  inline Eigen::Vector3d getVel(double t);
  inline Eigen::Vector3d getAcc(double t);
};

class MiniSnap {
 private:
  int N;                     // number of pieces
  Eigen::Matrix3d _headPVA;  // head's pos, vel, acc
  Eigen::Matrix3d _tailPVA;  // tail's pos, vel, acc
  Eigen::MatrixXd _Q;        //  matrix, J = x^T Q x
  Eigen::VectorXd _x;        // solutions
  Eigen::MatrixXd _A;
  Eigen::VectorXd _b;
  std::vector<Eigen::Vector3d> _waypoints;
  std::vector<double> _timeAlloc;

 public:
  MiniSnap() {}
  ~MiniSnap() {}
  void reset(const Eigen::Matrix3d &head, const Eigen::Matrix3d &tail,
             const std::vector<Eigen::Vector3d> &waypoints,
             const std::vector<double> &timeAlloc);
  bool solve();
  bool solveQP();
  void getCostFunc();
  void getContinuityConstraint();
  void getWaypointsConstraint();
  void getTrajectory(Trajectory *traj);
};

class CorridorMiniSnap {
 private:
  int N;
  Eigen::Vector3d _headPVA;
  Eigen::Vector3d _tailPVA;
};

};  // namespace trajopt

#endif  // CORRIDOR_MINI_SNAP_H_
