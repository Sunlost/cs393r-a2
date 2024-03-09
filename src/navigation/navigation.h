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
\file    navigation.h
\brief   Interface for reference Navigation class.
\author  Joydeep Biswas, (C) 2019
*/
//========================================================================
#ifndef NAVIGATION_H
#define NAVIGATION_H

#include <deque>
#include <vector>

#include "ackermann_motion_primitives.h"
#include "amrl_msgs/AckermannCurvatureDriveMsg.h"
#include "amrl_msgs/VisualizationMsg.h"
#include "constant_curvature_arc.h"
#include "eigen3/Eigen/Dense"
#include "navigation_params.h"
#include "ros/ros.h"
#include "ros/ros_helpers.h"
#include "vector_map/vector_map.h"

using amrl_msgs::AckermannCurvatureDriveMsg;
using amrl_msgs::VisualizationMsg;
using Eigen::Rotation2Df;
using Eigen::Vector2f;
using std::string;
using std::vector;

using namespace math_util;
using namespace ros_helpers;

namespace ros {
class NodeHandle;
}  // namespace ros

namespace navigation {

struct PathOption {
  float curvature;
  float clearance;
  float free_path_length;
  Eigen::Vector2f obstruction;
  Eigen::Vector2f closest_point;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

struct Command {
  double time;
  // Eigen::Vector2f velocity;
  // float omega;
  AckermannCurvatureDriveMsg drive_msg;
};

class Navigation {
 public:
  // Constructor
  explicit Navigation(const string& map_name,
                      NavigationParams& params,
                      ros::NodeHandle* n);

  // Used in callback from localization to update position.
  void UpdateLocation(const Eigen::Vector2f& loc, float angle);

  void PruneCommandQueue();

  // Used in callback for odometry messages to update based on odometry.
  void UpdateOdometry(const Eigen::Vector2f& loc,
                      float angle,
                      const Eigen::Vector2f& vel,
                      float ang_vel,
                      double time);

  // Updates based on an observed laser scan
  void ObservePointCloud(const std::vector<Eigen::Vector2f>& cloud, double time);

  // Used to predict the robot's state forward in time.
  void UpdateCommandHistory(const AckermannCurvatureDriveMsg& drive_msg);

  void ForwardPredict(double time);

  // Main function called continously from main
  void Run();
  
  // Used to set the next target pose.
  void SetNavGoal(const Eigen::Vector2f& loc, float angle);
  // Used to set autonomy status
  void SetAutonomy(bool autonomy_enabled);

 private:
  // Map of the environment.
  vector_map::VectorMap map_;
  // Navigation parameters
  navigation::NavigationParams params_;

  // Autonomy
  bool autonomy_enabled_;
  
  bool odom_initialized_;
  bool localization_initialized_;

  // odom states
  Eigen::Vector2f odom_start_loc_;
  float odom_start_angle_;
  Eigen::Vector2f odom_loc_;
  float odom_angle_;
  double t_odom_;

  // Current robot velocity.
  Eigen::Vector2f robot_vel_;
  float robot_omega_;

  // localization states
  Eigen::Vector2f robot_start_loc_;
  float robot_start_angle_;
  Eigen::Vector2f robot_loc_;
  float robot_angle_;
  double t_localization_;

  // command history
  std::deque<Command> command_history_;

  // Latest observed point cloud.
  std::vector<Eigen::Vector2f> point_cloud_;  // base_link frame
  std::vector<Eigen::Vector2f> fp_point_cloud_; // forward predicted base_link frame
  double t_point_cloud_;

  // Whether navigation is complete.
  bool nav_complete_;
  // Navigation goal location.
  Eigen::Vector2f nav_goal_loc_;
  // Navigation goal angle.
  float nav_goal_angle_;

  // Ackermann Path Sampler object
  motion_primitives::AckermannSampler ackermann_sampler_;

  

  void test1DTOC();

  void testSamplePaths(AckermannCurvatureDriveMsg& drive_msg);
};

}  // namespace navigation

#endif  // NAVIGATION_H
