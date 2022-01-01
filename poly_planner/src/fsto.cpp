#include <fsto/fsto.h>

using namespace std;
using namespace ros;
using namespace Eigen;

MavGlobalPlanner::MavGlobalPlanner(Config &conf, NodeHandle &nh_)
    : config(conf), nh(nh_), visualization(config, nh) {
  targetSub =
      nh.subscribe(config.targetTopic, 1, &MavGlobalPlanner::targetCallBack,
                   this, TransportHints().tcpNoDelay());
  trajPub = nh.advertise<quadrotor_msgs::PolynomialTrajectory>(
      config.trajectoryTopic, 1);
}

inline std::vector<Eigen::Matrix<double, 6, -1>> flipNormal(
    const std::vector<Eigen::MatrixXd> &hPolys) {
  std::vector<Eigen::Matrix<double, 6, -1>> hPolysON;
  hPolysON.reserve(hPolys.size());
  for (const Eigen::MatrixXd &ele : hPolys) {
    hPolysON.emplace_back(6, ele.cols());
    hPolysON.back() << ele;
    hPolysON.back().topRows<3>().array() *= -1.0;
  }
  return hPolysON;
}

inline Eigen::Vector3d jetColor(double a) {
  double s = a * 4;
  Eigen::Vector3d c;  // [r, g, b]
  switch ((int)floor(s)) {
    case 0:
      c << 0, 0, s;
      break;
    case 1:
      c << 0, s - 1, 1;
      break;
    case 2:
      c << s - 2, 1, 3 - s;
      break;
    case 3:
      c << 1, 4 - s, 0;
      break;
    default:
      c << 0.5, 0.5, 0.5;
      break;
  }
  return c;
}

void MavGlobalPlanner::targetCallBack(
    const geometry_msgs::PoseStamped::ConstPtr &msg) {
  int M = config.num;
  Eigen::Vector3d zeroVec(0.0, 0.0, 0.0);

  vector<Vector3d> route;
  Eigen::MatrixXd inifin;
  // std::vector<Eigen::Matrix<double, 3, 2>> gates = genGate(M, config.zBounds,
  // config.yBounds,
  //                                                          config.xBounds,
  //                                                          config.tBounds,
  //                                                          inifin);
  std::vector<Eigen::Matrix<double, 3, 2>> gates =
      readGate(config.filepath, inifin);
  M = gates.size();
  std::vector<Eigen::Matrix<double, 6, -1>> hPolys =
      genGateSFC(M, zeroVec, gates, 0.2, 0.05, inifin);

  Eigen::Matrix3d iniState;
  Eigen::Matrix3d finState;
  iniState << inifin.col(0), zeroVec, zeroVec;
  finState << inifin.col(1), zeroVec, zeroVec;
  double vmax = config.maxVelRate, amax = config.maxAccRate;
  Eigen::Vector3d chi(config.chiVec[0], config.chiVec[1],
                      config.chiVec[2]); /* weights of p, v, a */
  double smoothEps = config.smoothEps;
  double res = config.res;
  int itg = config.itg;

  std::chrono::high_resolution_clock::time_point tic =
      std::chrono::high_resolution_clock::now();
  GCOPTER nonlinOpt;
  Trajectory traj;

  if (!nonlinOpt.setup(config.weightT, 1.0, iniState, finState, hPolys, res,
                       itg, vmax, amax, smoothEps, chi, config.c2Diffeo,
                       config.retraction)) {
    vec_E<Polyhedron3D> polyhedra;
    polyhedra.reserve(hPolys.size());
    for (const auto &ele : hPolys) {
      Polyhedron3D hPoly;
      for (int i = 0; i < ele.cols(); i++) {
        hPoly.add(Hyperplane3D(ele.col(i).tail<3>(), ele.col(i).head<3>()));
      }
      polyhedra.push_back(hPoly);
    }
    visualization.visualizePolyH(polyhedra, ros::Time::now());
    ROS_WARN("Planner cannot find a feasible solution in the current problem.");
    return;
  }

  nonlinOpt.optimize(traj, config.relCostTol);
  std::chrono::high_resolution_clock::time_point toc =
      std::chrono::high_resolution_clock::now();
  double compTime =
      std::chrono::duration_cast<std::chrono::microseconds>(toc - tic).count() *
      1.0e-3;
  std::cout << std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                       .count() *
                   1.0e-3
            << "ms" << std::endl;
  std::cout << "Max Vel Rate: " << traj.getMaxVelRate() << std::endl;
  std::cout << "Total Time: " << traj.getTotalDuration() << std::endl;

  if (traj.getPieceNum() > 0) {
    quadrotor_msgs::PolynomialTrajectory trajMsg;
    Time stamp = ros::Time::now();
    polynomialTrajConverter(traj, trajMsg, Eigen::Isometry3d::Identity(),
                            stamp);
    trajPub.publish(trajMsg);
    visualization.visualize(traj, route, ros::Time::now(), compTime,
                            traj.getMaxVelRate(), traj.getTotalDuration());

    vec_E<Polyhedron3D> polyhedra;
    polyhedra.reserve(hPolys.size());
    for (const auto &ele : hPolys) {
      Polyhedron3D hPoly;
      for (int i = 0; i < ele.cols(); i++) {
        hPoly.add(Hyperplane3D(ele.col(i).tail<3>(), ele.col(i).head<3>()));
      }
      polyhedra.push_back(hPoly);
    }
    visualization.visualizePolyH(polyhedra, ros::Time::now());
    ROS_WARN("PolyH has been published");
  }
}

void MavGlobalPlanner::polynomialTrajConverter(
    const Trajectory &traj, quadrotor_msgs::PolynomialTrajectory &trajMsg,
    Eigen::Isometry3d tfR2L, Time &iniStamp) {
  trajMsg.header.stamp = iniStamp;
  static uint32_t traj_id = 0;
  traj_id++;
  trajMsg.trajectory_id = traj_id;
  trajMsg.action = quadrotor_msgs::PolynomialTrajectory::ACTION_ADD;
  trajMsg.num_order = traj[0].getOrder();
  trajMsg.num_segment = traj.getPieceNum();
  Eigen::Vector3d initialVel, finalVel;
  initialVel = tfR2L * traj.getVel(0.0);
  finalVel = tfR2L * traj.getVel(traj.getTotalDuration());
  trajMsg.start_yaw = atan2(initialVel(1), initialVel(0));
  trajMsg.final_yaw = atan2(finalVel(1), finalVel(0));

  for (size_t p = 0; p < (size_t)traj.getPieceNum(); p++) {
    trajMsg.time.push_back(traj[p].getDuration());
    trajMsg.order.push_back(traj[p].getCoeffMat().cols() - 1);

    Eigen::VectorXd linearTr(2);
    linearTr << 0.0, trajMsg.time[p];
    std::vector<Eigen::VectorXd> linearTrCoeffs;
    linearTrCoeffs.emplace_back(1);
    linearTrCoeffs[0] << 1;
    for (size_t k = 0; k < trajMsg.order[p]; k++) {
      linearTrCoeffs.push_back(
          RootFinder::polyConv(linearTrCoeffs[k], linearTr));
    }

    Eigen::MatrixXd coefMat(3, traj[p].getCoeffMat().cols());
    for (int i = 0; i < coefMat.cols(); i++) {
      coefMat.col(i) =
          tfR2L.rotation() *
          traj[p].getCoeffMat().col(coefMat.cols() - i - 1).head<3>();
    }
    coefMat.col(0) = (coefMat.col(0) + tfR2L.translation()).eval();

    for (int i = 0; i < coefMat.cols(); i++) {
      double coefx(0.0), coefy(0.0), coefz(0.0);
      for (int j = i; j < coefMat.cols(); j++) {
        coefx += coefMat(0, j) * linearTrCoeffs[j](i);
        coefy += coefMat(1, j) * linearTrCoeffs[j](i);
        coefz += coefMat(2, j) * linearTrCoeffs[j](i);
      }
      trajMsg.coef_x.push_back(coefx);
      trajMsg.coef_y.push_back(coefy);
      trajMsg.coef_z.push_back(coefz);
    }
  }

  trajMsg.mag_coeff = 1.0;
  trajMsg.debug_info = "";
}

Visualization::Visualization(Config &conf, NodeHandle &nh_)
    : config(conf), nh(nh_) {
  routePub =
      nh.advertise<visualization_msgs::Marker>("/visualization/route", 1);
  wayPointsPub =
      nh.advertise<visualization_msgs::Marker>("/visualization/waypoints", 1);
  appliedTrajectoryPub = nh.advertise<visualization_msgs::Marker>(
      "/visualization/applied_trajectory", 1);
  hPolyPub = nh.advertise<decomp_ros_msgs::PolyhedronArray>(
      "/visualization/polyhedra", 1);
  textPub = nh.advertise<visualization_msgs::Marker>("/visualization/text", 1);
}

void Visualization::visualize(const Trajectory &appliedTraj,
                              const vector<Vector3d> &route, Time timeStamp,
                              double compT, double maxV, double totalT) {
  visualization_msgs::Marker routeMarker, wayPointsMarker, appliedTrajMarker;

  routeMarker.id = 0;
  routeMarker.type = visualization_msgs::Marker::LINE_LIST;
  routeMarker.header.stamp = timeStamp;
  routeMarker.header.frame_id = config.odomFrame;
  routeMarker.pose.orientation.w = 1.00;
  routeMarker.action = visualization_msgs::Marker::ADD;
  routeMarker.ns = "route";
  routeMarker.color.r = 1.00;
  routeMarker.color.g = 0.00;
  routeMarker.color.b = 0.00;
  routeMarker.color.a = 1.00;
  routeMarker.scale.x = 0.05;

  wayPointsMarker = routeMarker;
  wayPointsMarker.type = visualization_msgs::Marker::SPHERE_LIST;
  wayPointsMarker.ns = "waypoints";
  wayPointsMarker.color.r = 1.00;
  wayPointsMarker.color.g = 0.00;
  wayPointsMarker.color.b = 0.00;
  wayPointsMarker.scale.x = 0.20;
  wayPointsMarker.scale.y = 0.20;
  wayPointsMarker.scale.z = 0.20;

  appliedTrajMarker = routeMarker;
  appliedTrajMarker.header.frame_id = config.odomFrame;
  appliedTrajMarker.id = 0;
  appliedTrajMarker.ns = "applied_trajectory";
  appliedTrajMarker.color.r = 0.00;
  appliedTrajMarker.color.g = 0.50;
  appliedTrajMarker.color.b = 1.00;
  appliedTrajMarker.scale.x = 0.10;

  if (route.size() > 0) {
    bool first = true;
    Vector3d last;
    for (auto it : route) {
      if (first) {
        first = false;
        last = it;
        continue;
      }
      geometry_msgs::Point point;

      point.x = last(0);
      point.y = last(1);
      point.z = last(2);
      routeMarker.points.push_back(point);
      point.x = it(0);
      point.y = it(1);
      point.z = it(2);
      routeMarker.points.push_back(point);
      last = it;
    }

    routePub.publish(routeMarker);
  }

  if (route.size() > 0) {
    Eigen::Matrix3Xd wps = appliedTraj.getPositions();
    for (int i = 0; i < wps.cols(); i++) {
      geometry_msgs::Point point;
      point.x = wps.col(i)(0);
      point.y = wps.col(i)(1);
      point.z = wps.col(i)(2);
      wayPointsMarker.points.push_back(point);
    }

    wayPointsPub.publish(wayPointsMarker);
  }

  if (appliedTraj.getPieceNum() > 0) {
    double T = 0.01;
    Vector3d lastX = appliedTraj.getPos(0.0);
    for (double t = T; t < appliedTraj.getTotalDuration(); t += T) {
      std_msgs::ColorRGBA c;
      Eigen::Vector3d jets =
          jetColor(appliedTraj.getVel(t).norm() / config.maxVelRate);
      c.r = jets[0];
      c.g = jets[1];
      c.b = jets[2];

      geometry_msgs::Point point;
      Vector3d X = appliedTraj.getPos(t);
      point.x = lastX(0);
      point.y = lastX(1);
      point.z = lastX(2);
      appliedTrajMarker.points.push_back(point);
      appliedTrajMarker.colors.push_back(c);
      point.x = X(0);
      point.y = X(1);
      point.z = X(2);
      appliedTrajMarker.points.push_back(point);
      appliedTrajMarker.colors.push_back(c);
      lastX = X;
    }
    appliedTrajectoryPub.publish(appliedTrajMarker);
  }

  if (appliedTraj.getPieceNum() > 0) {
    visualization_msgs::Marker textMarker;
    textMarker.header.frame_id = config.odomFrame;
    textMarker.header.stamp = timeStamp;
    textMarker.ns = "text";
    textMarker.id = 1;
    textMarker.type = visualization_msgs::Marker::TEXT_VIEW_FACING;
    textMarker.action = visualization_msgs::Marker::ADD;

    textMarker.pose.position.x = -9;
    textMarker.pose.position.y = 0.0;
    textMarker.pose.position.z = 6.0;
    textMarker.pose.orientation.x = 0.0;
    textMarker.pose.orientation.y = 0.0;
    textMarker.pose.orientation.z = 0.0;
    textMarker.pose.orientation.w = 1.0;
    textMarker.scale.x = 1.0;
    textMarker.scale.y = 1.0;
    textMarker.scale.z = 1.0;
    textMarker.color.r = 1.0;
    textMarker.color.g = 0.0;
    textMarker.color.b = 0.0;
    textMarker.color.a = 1.0;
    textMarker.text = "Comp: ";
    textMarker.text += to_string((int)(compT));
    textMarker.text += ".";
    textMarker.text += to_string((int)(compT * 10) % 10);
    textMarker.text += "ms\n";
    textMarker.text += "Max speed: ";
    textMarker.text += to_string((int)(maxV));
    textMarker.text += ".";
    textMarker.text += to_string((int)(maxV * 100) % 100);
    textMarker.text += "m/s\n";
    textMarker.text += "Total time: ";
    textMarker.text += to_string((int)(totalT));
    textMarker.text += ".";
    textMarker.text += to_string((int)(totalT * 100) % 100);
    textMarker.text += "s\n";
    textPub.publish(textMarker);
  }
}

void Visualization::visualizePolyH(const vec_E<Polyhedron3D> &polyhedra,
                                   ros::Time timeStamp) {
  decomp_ros_msgs::PolyhedronArray poly_msg =
      DecompROS::polyhedron_array_to_ros(polyhedra);
  poly_msg.header.frame_id = config.odomFrame;
  poly_msg.header.stamp = timeStamp;
  hPolyPub.publish(poly_msg);
}
