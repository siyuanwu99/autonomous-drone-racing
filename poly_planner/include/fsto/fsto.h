#ifndef FSTO_H
#define FSTO_H

#include "fsto/config.h"
#include "quadrotor_msgs/PolynomialTrajectory.h"
#include "decomp_util/ellipsoid_decomp.h"
#include "decomp_ros_utils/data_ros_utils.h"
#include "gcopter/gcopter.hpp"
#include "gcopter/sfc_gen.hpp"

#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>

#include <ros/ros.h>
#include <ros/console.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

class Visualization
{
public:
    Visualization(Config &conf, ros::NodeHandle &nh_);

    Config config;
    ros::NodeHandle nh;

    ros::Publisher routePub;
    ros::Publisher wayPointsPub;
    ros::Publisher appliedTrajectoryPub;
    ros::Publisher hPolyPub;
    ros::Publisher textPub;

    void visualize(const Trajectory &appliedTraj, const std::vector<Eigen::Vector3d> &route, ros::Time timeStamp, double compT);
    void visualizePolyH(const vec_E<Polyhedron3D> &polyhedra, ros::Time timeStamp);
};

class MavGlobalPlanner
{
public:
    MavGlobalPlanner(Config &conf, ros::NodeHandle &nh_);

    Config config;
    ros::NodeHandle nh;

    ros::Subscriber targetSub;
    void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg);
    ros::Publisher trajPub;

    static void polynomialTrajConverter(const Trajectory &traj,
                                        quadrotor_msgs::PolynomialTrajectory &trajMsg,
                                        Eigen::Isometry3d tfR2L, ros::Time &iniStamp);

    EllipsoidDecomp3D cvxDecomp;
    Visualization visualization;
};
#endif