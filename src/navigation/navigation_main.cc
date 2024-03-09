//========================================================================
//  This software is free: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License Version 3,
//  as published by the Free Software Foundation.
//
//  This software is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  Version 3 in the file COPYING that came with this distribution.
//  If not, see <http://www.gnu.org/licenses/>.
//========================================================================
/*!
\file    navigation_main.cc
\brief   Main entry point for reference Navigation implementation
\author  Joydeep Biswas, (C) 2019
*/
//========================================================================

#include <inttypes.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "amrl_msgs/Localization2DMsg.h"
#include "config_reader/config_reader.h"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include "geometry_msgs/Pose2D.h"
#include "geometry_msgs/PoseArray.h"
#include "geometry_msgs/PoseStamped.h"
#include "geometry_msgs/PoseWithCovarianceStamped.h"
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "nav_msgs/Odometry.h"
#include "navigation.h"
#include "ros/ros.h"
#include "sensor_msgs/LaserScan.h"
#include "shared/math/math_util.h"
#include "shared/ros/ros_helpers.h"
#include "shared/util/timer.h"
#include "std_msgs/Bool.h"
#include "std_msgs/String.h"
#include "visualization_msgs/Marker.h"
#include "visualization_msgs/MarkerArray.h"

using amrl_msgs::Localization2DMsg;
using Eigen::Vector2f;
using math_util::DegToRad;
using math_util::RadToDeg;
using navigation::Navigation;
using ros::Time;
using ros_helpers::Eigen2DToRosPoint;
using ros_helpers::Eigen3DToRosPoint;
using ros_helpers::RosPoint;
using ros_helpers::SetRosVector;
using std::string;
using std::vector;

// Create command line arguments
DEFINE_string(nav_config, "config/navigation.lua", "Navigation config file");
DEFINE_string(laser_topic, "scan", "Name of ROS topic for LIDAR data");
DEFINE_string(odom_topic, "odom", "Name of ROS topic for odometry data");
DEFINE_string(loc_topic, "localization", "Name of ROS topic for localization");
DEFINE_string(init_topic, "initialpose", "Name of ROS topic for initialization");
DEFINE_string(map, "GDC1", "Name of vector map file");

bool run_ = true;
sensor_msgs::LaserScan last_laser_msg_;
Navigation* navigation_ = nullptr;

void LaserCallback(const sensor_msgs::LaserScan& msg) {
  if (FLAGS_v > 0) {
    cout << "=============== [Navigation Main] LaserCallback ==============" << endl;
    printf("Laser t=%f, dt=%f\n",
           msg.header.stamp.toSec(),
           GetWallTime() - msg.header.stamp.toSec());
    cout << "==============================================================\n" << endl;
  }
  // Location of the laser on the robot. Assumes the laser is forward-facing.
  const Vector2f kLaserLoc(0.2, 0);

  static vector<Vector2f> point_cloud_;
  point_cloud_.clear();

  // Convert the LaserScan to a point cloud
  // The LaserScan parameters are accessible as follows:
  // msg.angle_increment // Angular increment between subsequent rays
  // msg.angle_max // Angle of the first ray
  // msg.angle_min // Angle of the last ray
  // msg.range_max // Maximum observable range
  // msg.range_min // Minimum observable range
  // msg.ranges[i] // The range of the i'th ray

  int angle_range = msg.angle_max - msg.angle_min;
  int num_rays = angle_range / msg.angle_increment;

  // Convert the LaserScan to a point cloud in the base_link frame
  for (int i = 0; i < num_rays; i++) {
    double r = msg.ranges[i];
    if (r < msg.range_max && r > msg.range_min) {
      double theta = msg.angle_min + i * msg.angle_increment;
      Vector2f point = {r * cos(theta), r * sin(theta)};
      point += kLaserLoc;
      point_cloud_.push_back(point);
    }
  }

  navigation_->ObservePointCloud(point_cloud_, msg.header.stamp.toSec());
  last_laser_msg_ = msg;
}

void OdometryCallback(const nav_msgs::Odometry& msg) {
  if (FLAGS_v > 0) {
    cout << "============= [Navigation Main] OdometryCallback =============" << endl;
    printf("Odometry t=%f\n", msg.header.stamp.toSec());
    printf("Position: (%f,%f)\n", msg.pose.pose.position.x, msg.pose.pose.position.y);
    cout << "==============================================================\n" << endl;
  }
  navigation_->UpdateOdometry(
      Vector2f(msg.pose.pose.position.x, msg.pose.pose.position.y),
      2.0 * atan2(msg.pose.pose.orientation.z, msg.pose.pose.orientation.w),
      Vector2f(msg.twist.twist.linear.x, msg.twist.twist.linear.y),
      msg.twist.twist.angular.z,
      msg.header.stamp.toSec());
}

void GoToCallback(const geometry_msgs::PoseStamped& msg) {
  const Vector2f loc(msg.pose.position.x, msg.pose.position.y);
  const float angle = 2.0 * atan2(msg.pose.orientation.z, msg.pose.orientation.w);
  printf("Goal: (%f,%f) %f\u00b0\n", loc.x(), loc.y(), angle);
  navigation_->SetNavGoal(loc, angle);
}

void AutonomyCallback(const std_msgs::Bool& msg) { navigation_->SetAutonomy(msg.data); }

void SignalHandler(int) {
  if (!run_) {
    printf("Force Exit.\n");
    exit(0);
  }
  printf("Exiting.\n");
  run_ = false;
}

void LocalizationCallback(const amrl_msgs::Localization2DMsg msg) {
  if (FLAGS_v > 0) {
    cout << "=========== [Navigation Main] LocalizationCallback ============" << endl;
    printf("Localization t=%f\n", GetWallTime());
    printf("Position: (%f,%f)\n", msg.pose.x, msg.pose.y);
    cout << "==============================================================\n" << endl;
  }
  navigation_->UpdateLocation(Vector2f(msg.pose.x, msg.pose.y), msg.pose.theta);
}

void StringCallback(const std_msgs::String& msg) { std::cout << msg.data << "\n"; }

void LoadConfig(navigation::NavigationParams& params) {
#define REAL_PARAM(x) CONFIG_DOUBLE(x, "NavigationParameters." #x);
  REAL_PARAM(dt);
  REAL_PARAM(system_latency);
  REAL_PARAM(max_linear_accel);
  REAL_PARAM(max_linear_deccel);
  REAL_PARAM(max_linear_speed);
  REAL_PARAM(max_angular_accel);
  REAL_PARAM(max_angular_deccel);
  REAL_PARAM(max_angular_speed);
  REAL_PARAM(max_curvature);
  REAL_PARAM(max_path_length);
  REAL_PARAM(max_clearance);
  REAL_PARAM(clearance_weight);
  REAL_PARAM(arc_length_weight);
  REAL_PARAM(distance_weight);
  REAL_PARAM(goal_tolerance);
  REAL_PARAM(robot_length);
  REAL_PARAM(robot_width);
  REAL_PARAM(robot_wheelbase);
  REAL_PARAM(obstacle_margin);
  REAL_PARAM(lidar_offset);

  config_reader::ConfigReader reader({FLAGS_nav_config});
  params.dt = CONFIG_dt;
  params.system_latency = CONFIG_system_latency;
  params.linear_limits = navigation::MotionLimits(
      CONFIG_max_linear_accel, CONFIG_max_linear_deccel, CONFIG_max_linear_speed);
  params.angular_limits = navigation::MotionLimits(
      CONFIG_max_angular_accel, CONFIG_max_angular_deccel, CONFIG_max_angular_speed);
  params.max_curvature = CONFIG_max_curvature;
  params.max_path_length = CONFIG_max_path_length;
  params.max_clearance = CONFIG_max_clearance;
  params.clearance_weight = CONFIG_clearance_weight;
  params.arc_length_weight = CONFIG_arc_length_weight;
  params.distance_weight = CONFIG_distance_weight;
  params.goal_tolerance = CONFIG_goal_tolerance;
  params.robot_length = CONFIG_robot_length;
  params.robot_width = CONFIG_robot_width;
  params.robot_wheelbase = CONFIG_robot_wheelbase;
  params.obstacle_margin = CONFIG_obstacle_margin;
  params.lidar_offset = CONFIG_lidar_offset;
}

int main(int argc, char** argv) {
  google::ParseCommandLineFlags(&argc, &argv, false);
  signal(SIGINT, SignalHandler);
  // Initialize ROS.
  ros::init(argc, argv, "navigation", ros::init_options::NoSigintHandler);
  ros::NodeHandle n;

  // Load Configurations
  navigation::NavigationParams params;
  LoadConfig(params);

  // Construct navigation
  navigation_ = new Navigation(FLAGS_map, params, &n);

  ros::Subscriber string_sub = n.subscribe("string_topic", 1, &StringCallback);
  ros::Subscriber velocity_sub = n.subscribe(FLAGS_odom_topic, 1, &OdometryCallback);
  ros::Subscriber localization_sub =
      n.subscribe(FLAGS_loc_topic, 1, &LocalizationCallback);
  ros::Subscriber laser_sub = n.subscribe(FLAGS_laser_topic, 1, &LaserCallback);
  ros::Subscriber goto_sub = n.subscribe("/move_base_simple/goal", 1, &GoToCallback);
  ros::Subscriber autonomy_sub = n.subscribe("/autonomy_enabler", 1, &AutonomyCallback);

  RateLoop loop(20.0);
  while (run_ && ros::ok()) {
    ros::spinOnce();

    navigation_->Run();
    loop.Sleep();
  }
  delete navigation_;
  return 0;
}
