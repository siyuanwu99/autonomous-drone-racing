#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>

#include <ros/ros.h>

struct Config
{
    std::string targetTopic;
    std::string trajectoryTopic;
    std::string odomFrame;
    double maxAccRate;
    double maxVelRate;
    double weightT;
    double smoothEps;
    std::vector<double> chiVec;
    double relCostTol;
    bool c2Diffeo;
    bool retraction;

    int num;
    double xBounds;
    double yBounds;
    double zBounds;
    double tBounds; 

    static void loadParameters(Config &conf, const ros::NodeHandle &nh_priv)
    {
        nh_priv.getParam("TargetTopic", conf.targetTopic);
        nh_priv.getParam("TrajectoryTopic", conf.trajectoryTopic);
        nh_priv.getParam("OdomFrame", conf.odomFrame);
        nh_priv.getParam("MaxAccRate", conf.maxAccRate);
        nh_priv.getParam("MaxVelRate", conf.maxVelRate);
        nh_priv.getParam("WeightT", conf.weightT);
        nh_priv.getParam("SmoothEps", conf.smoothEps);
        nh_priv.getParam("ChiVec", conf.chiVec);
        nh_priv.getParam("RelCostTol", conf.relCostTol);
        nh_priv.getParam("C2Diffeo", conf.c2Diffeo);
        nh_priv.getParam("Retraction", conf.retraction);
        nh_priv.getParam("NumGates", conf.num);
        nh_priv.getParam("MaxX", conf.xBounds);
        nh_priv.getParam("MaxY", conf.yBounds);
        nh_priv.getParam("MaxZ", conf.zBounds);
        nh_priv.getParam("MaxTheta", conf.tBounds);
    }
};

#endif