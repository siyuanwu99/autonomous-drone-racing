/**
 * @file corridor_minisnap.cpp
 * @author Siyuan Wu (siyuanwu99@gmail.com)
 * @brief
 * @version 1.0
 * @date 2021-12-22
 *
 * @copyright Copyright (c) 2021
 *
 */

#include <CorridorMiniSnap/corridor_minisnap.h>

using namespace trajopt;

/**
 * @brief get position
 * @param t relative time stamp
 * @return Eigen::Vector3d
 */
inline Eigen::Vector3d PolyPiece::getPos(double t) {
  Eigen::Vector3d pos = Eigen::Vector3d::Zero();
  Eigen::VectorXd T(N_ORDER);
  T(0) = 1;
  for (int i = 1; i <= N_ORDER; i++) {
    T(i) = pow(t, i);
  }
  pos = _coeffs * T;
  return pos;
}
inline Eigen::Vector3d PolyPiece::getVel(double t) {
  Eigen::VectorXd T(6);
  T(0) = 1;
  for (int i = 1; i <= 6; i++) {
    T(i) = pow(t, i) * (i + 1);
  }
  Eigen::Vector3d vel = _coeffs.block<3, 6>(1, 2) * T;
  return vel;
}

inline Eigen::Vector3d PolyPiece::getAcc(double t) {
  Eigen::VectorXd T(5);
  T(0) = 1;
  for (int i = 1; i <= 5; i++) {
    T(i) = pow(t, i) * (i + 1) * (i + 2);
  }
  Eigen::Vector3d acc = _coeffs.block<3, 5>(1, 3) * T;
  return acc;
}

inline double PolyPiece::getDuration() { return _duration; }

void Trajectory::setDuration(const std::vector<double> &t) {
  N = t.size();
  for (int i = 0; i < N; i++) {
  }
}

void Trajectory::setTrajectory(const Eigen::VectorXd &x) {}

/**
 * @brief
 *
 * @param t0 absolute time stamp
 * @param t  return t \in [0, 1)
 * @param idx
 */
inline void Trajectory::locatePiece(const double &t0, double &t, int &idx) {
  idx = N;
  double tmp = t0;
  for (int i = 0; i < N; i++) {
    double Ti = _pieces[i].getDuration();
    if (tmp > Ti) {
      tmp -= Ti;
    } else {
      idx = i;
      t = tmp / Ti;
    }
  }
  /* if t0 is longer than all durations */
  if (idx == N) {
    idx = N - 1;
    t = 1;
  }
}

/**
 * @brief
 *
 * @param t absolute time stamp
 * @return Eigen::Vector3d
 */
inline Eigen::Vector3d Trajectory::getPos(double t) {
  double relative_time;
  int index;
  locatePiece(t, relative_time, index);
  return _pieces[index].getPos(relative_time);
}

inline Eigen::Vector3d Trajectory::getVel(double t) {
  double relative_time;
  int index;
  locatePiece(t, relative_time, index);
  return _pieces[index].getVel(relative_time);
}

inline Eigen::Vector3d Trajectory::getAcc(double t) {
  double relative_time;
  int index;
  locatePiece(t, relative_time, index);
  return _pieces[index].getAcc(relative_time);
}

void MiniSnap::reset(const Eigen::Matrix3d &head, const Eigen::Matrix3d &tail,
                     const std::vector<Eigen::Vector3d> &waypoints,
                     const std::vector<double> &timeAlloc) {
  _headPVA = head;
  _tailPVA = tail;
  _waypoints = waypoints;
  _timeAlloc = timeAlloc;
  N = timeAlloc.size();
  int S = N * (N_ORDER + 1) * DIM;
  _x.resize(S);
  _Q.resize(S, S);
  int M = 15 * N + 9;
  _A.resize(M, S);
  _b = Eigen::VectorXd::Zero(M);  // number of all constraints
}

/* T = 1 */
void MiniSnap::getCostFunc() {
  /* for single piece, single dimension */
  int D = N_ORDER + 1;  // size of matrix Q
  Eigen::Matrix<double, N_ORDER + 1, N_ORDER + 1> Q;
  for (int i = 0; i <= N_ORDER; i++) {
    for (int j = 0; j = N_ORDER; j++) {
      if (i < 4 || j < 4) {
        Q(i, j) = 0;
      }
      if (i + j > N_ORDER) {
        Q(i, j) = i * (i - 1) * (i - 2) * (i - 3) * j * (j - 1) * (j - 2) *
                  (j - 3) / (i + j - N_ORDER);
      }
    }
  }
  /* iterate all dimensions and all pieces */
  for (int i = 0; i < N * DIM; i++) {
    _Q.block(i * D, i * D, D, D) = Q;
  }
}

void MiniSnap::getWaypointsConstraint() {
  /* constraints for starting and ending states*/
  for (int i = 0; i < DIM; i++) {
    _A(0 + 4 * i, 0) = 1;
    _b(0 + 4 * i) = _headPVA(i, 0);
    _A(1 + 4 * i, 1) = 1;
    _b(1 + 4 * i) = _headPVA(i, 1);
    _A(2 + 4 * i, 2) = 2;
    _b(2 + 4 * i) = _headPVA(i, 2);
    _A(3 + 4 * i, 3) = 6;
    _b(3 + 4 * i) = 0;
  }

  int M = N * N_ORDER * DIM;
  Eigen::Matrix<double, 1, N_ORDER + 1> pos_1d;
  Eigen::Matrix<double, 1, N_ORDER + 1> vel_1d;
  Eigen::Matrix<double, 1, N_ORDER + 1> acc_1d;
  Eigen::Matrix<double, 1, N_ORDER + 1> jer_1d;
  pos_1d << 1, 1, 1, 1, 1, 1, 1, 1;
  vel_1d << 0, 1, 2, 3, 4, 5, 6, 7;
  acc_1d << 0, 0, 2, 6, 12, 20, 30, 42;
  jer_1d << 0, 0, 0, 6, 24, 60, 120, 210;

  for (int i = 0; i < DIM; i++) {
    _A.block(12 + 0 + 4 * i, M, 1, N_ORDER + 1) = pos_1d;
    _A.block(12 + 1 + 4 * i, M, 1, N_ORDER + 1) = vel_1d;
    _A.block(12 + 2 + 4 * i, M, 1, N_ORDER + 1) = acc_1d;
    _A.block(12 + 3 + 4 * i, M, 1, N_ORDER + 1) = jer_1d;
    _b(12 + 0 + 4 * i) = _tailPVA(i, 0);
    _b(12 + 1 + 4 * i) = _tailPVA(i, 1);
    _b(12 + 2 + 4 * i) = _tailPVA(i, 2);
    _b(12 + 3 + 4 * i) = 0;
  }

  /* constraints for medium states*/
  for (int j = 0; j < N - 1; j++) {
    for (int i = 0; i < DIM; i++) {
      _A.block(24 + i + 3 * j, M, 1, N_ORDER + 1) = pos_1d;
      _b(24 + i + 3 * j) = _waypoints[j](0);
    }
  }
}

void MiniSnap::getContinuityConstraint() {
  int M = N_ORDER + 1;
  int K = DIM * M;  // 3*8
  int P = 24 + 3 * (N - 1);
  Eigen::Matrix<double, 1, N_ORDER + 1> pos_1d;
  Eigen::Matrix<double, 1, N_ORDER + 1> vel_1d;
  Eigen::Matrix<double, 1, N_ORDER + 1> acc_1d;
  Eigen::Matrix<double, 1, N_ORDER + 1> jer_1d;
  pos_1d << 1, 1, 1, 1, 1, 1, 1, 1;
  vel_1d << 0, 1, 2, 3, 4, 5, 6, 7;
  acc_1d << 0, 0, 2, 6, 12, 20, 30, 42;
  jer_1d << 0, 0, 0, 6, 24, 60, 120, 210;

  for (int j = 0; j < N; j++) {
    for (int i = 0; i < DIM; i++) {
      _A.block(P + 0 + 4 * i + 12 * j, i * M + j * K, 1, N_ORDER + 1) = pos_1d;
      _A.block(P + 1 + 4 * i + 12 * j, i * M + j * K, 1, N_ORDER + 1) = vel_1d;
      _A.block(P + 2 + 4 * i + 12 * j, i * M + j * K, 1, N_ORDER + 1) = acc_1d;
      _A.block(P + 3 + 4 * i + 12 * j, i * M + j * K, 1, N_ORDER + 1) = jer_1d;
      _A(P + 0 + 4 * i + 12 * j, 0 + K * (j + 1) + i * M) = 1;
      _A(P + 1 + 4 * i + 12 * j, 1 + K * (j + 1) + i * M) = 1;
      _A(P + 2 + 4 * i + 12 * j, 2 + K * (j + 1) + i * M) = 2;
      _A(P + 3 + 4 * i + 12 * j, 3 + K * (j + 1) + i * M) = 6;
    }
  }
}

bool MiniSnap::solveQP() {
  IOSQP solver;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(N * (N_ORDER + 1) * DIM);
  Eigen::VectorXd l = Eigen::VectorXd::Zero(0);
  Eigen::SparseMatrix<double> Q = _Q.sparseView();
  Eigen::SparseMatrix<double> A = _A.sparseView();
  solver.setMats(Q, q, A, l, _b, 1e-3, 1e-3);
  solver.solve();
  _x = solver.getPrimalSol();
}

bool MiniSnap::solve() {
  getCostFunc();
  getWaypointsConstraint();
  getContinuityConstraint();
  bool isSuccess = solveQP();
  return isSuccess;
}


void MiniSnap::getTrajectory(Trajectory *traj) {
  traj->setDuration(_timeAlloc);
  traj->setTrajectory(_x);
}
