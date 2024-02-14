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
\file    navigation.cc
\brief   Starter code for navigation.
\author  Joydeep Biswas, (C) 2019
*/
//========================================================================

#include "gflags/gflags.h"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include "amrl_msgs/AckermannCurvatureDriveMsg.h"
#include "amrl_msgs/Pose2Df.h"
#include "amrl_msgs/VisualizationMsg.h"
#include "glog/logging.h"
#include "ros/ros.h"
#include "ros/package.h"
#include "shared/math/math_util.h"
#include "shared/util/timer.h"
#include "shared/ros/ros_helpers.h"
#include "navigation.h"
#include "visualization/visualization.h"
#include<cmath>

using Eigen::Vector2f;
using amrl_msgs::AckermannCurvatureDriveMsg;
using amrl_msgs::VisualizationMsg;
using std::string;
using std::vector;
using std::min;

using namespace math_util;
using namespace ros_helpers;

namespace {
ros::Publisher drive_pub_;
ros::Publisher viz_pub_;
VisualizationMsg local_viz_msg_;
VisualizationMsg global_viz_msg_;
AckermannCurvatureDriveMsg drive_msg_;
// Epsilon value for handling limited numerical precision.
const float kEpsilon = 1e-5;
} //namespace

namespace navigation {

string GetMapFileFromName(const string& map) {
  string maps_dir_ = ros::package::getPath("amrl_maps");
  return maps_dir_ + "/" + map + "/" + map + ".vectormap.txt";
}

Navigation::Navigation(const string& map_name, ros::NodeHandle* n) :
    odom_initialized_(false),
    localization_initialized_(false),
    robot_loc_(0, 0),
    robot_angle_(0),
    robot_vel_(0, 0),
    robot_omega_(0),
    nav_complete_(true),
    nav_goal_loc_(0, 0),
    nav_goal_angle_(0) {
  map_.Load(GetMapFileFromName(map_name));
  drive_pub_ = n->advertise<AckermannCurvatureDriveMsg>(
      "ackermann_curvature_drive", 1);
  viz_pub_ = n->advertise<VisualizationMsg>("visualization", 1);
  local_viz_msg_ = visualization::NewVisualizationMessage(
      "base_link", "navigation_local");
  global_viz_msg_ = visualization::NewVisualizationMessage(
      "map", "navigation_global");
  InitRosHeader("base_link", &drive_msg_.header);
}

void Navigation::SetNavGoal(const Vector2f& loc, float angle) {
  nav_goal_loc_ = loc;
  nav_goal_angle_ = angle;
}

void Navigation::UpdateLocation(const Eigen::Vector2f& loc, float angle) {
  localization_initialized_ = true;
  robot_loc_ = loc;
  robot_angle_ = angle;
}

void Navigation::UpdateOdometry(const Vector2f& loc,
                                float angle,
                                const Vector2f& vel,
                                float ang_vel) {
  robot_omega_ = ang_vel;
  robot_vel_ = vel;
  if (!odom_initialized_) {
    odom_start_angle_ = angle;
    odom_start_loc_ = loc;
    odom_initialized_ = true;
    odom_loc_ = loc;
    odom_angle_ = angle;
    return;
  }
  odom_loc_ = loc;
  odom_angle_ = angle;
}

void Navigation::ObservePointCloud(const vector<Vector2f>& cloud,
                                   double time) {
  point_cloud_ = cloud;                                     
}

void Navigation::Run() {
  // This function gets called 20 times a second to form the control loop.
  
  // Clear previous visualizations.
  visualization::ClearVisualizationMsg(local_viz_msg_);
  visualization::ClearVisualizationMsg(global_viz_msg_);

  // If odometry has not been initialized, we can't do anything.
  if (!odom_initialized_) return;

  // here is where we should be able to see the obstacle and the goal. what does that mean?
  // The control iteration goes here. 
  // Feel free to make helper functions to structure the control appropriately.
  
  // The latest observed point cloud is accessible via "point_cloud_"
  // drive_msg_.velocity = 1;
  // drive_msg_.curvature = 1;


  // Eventually, you will have to set the control values to issue drive commands:
  // drive_msg_.curvature = ...;
  // drive_msg_.velocity = 1;

  // For a fixed arc, implement obstacle detection and integrate it with 1-D TOC to drive up to the observed obstacle.
  pick_arc();

  // Add timestamps to all messages.
  local_viz_msg_.header.stamp = ros::Time::now();
  global_viz_msg_.header.stamp = ros::Time::now();
  drive_msg_.header.stamp = ros::Time::now();
  // Publish messages.
  viz_pub_.publish(local_viz_msg_);
  viz_pub_.publish(global_viz_msg_);
  drive_pub_.publish(drive_msg_);
}

float magnitude(double x, double y) {
  return sqrt(pow(x, 2) + pow(y, 2));
}

PathOption Navigation::pick_arc() {
  // for loop of arcs
  // for each arc, calc score
  // return the best arc
  float arc_score = 0.0; 
  float best_arc_score = -1;
  // float clearance = 0.0;
  double temp_fpl = 100;
  double h = 0.535 + .1; // add .1 for safety margin
  double w = 0.281 / 2 + .1;
  vector<PathOption> path_options;
  PathOption best_path_option;

  int index = 0;

  // uncomment for debugging - I need to figure out how to set a nav target
  Eigen::Vector2f goal(5, 3);
  // visualization::DrawCross(goal, .3, 0x239847, local_viz_msg_);

  // curvature options from right to left
  // max curvature is 1
  for(double i = -1; i <= 1; i += 0.1) {

    path_options.push_back(PathOption());
    double radius = 1 / (i + 1e-6); // adding small value to account for 0 curvature
    path_options[index].free_path_length = 100; // init to some high value
    
    Eigen::Vector2f center(0, radius); // right = negative value
    double goal_mag = magnitude(goal.x() - center.x(), goal.y() - center.y());
    Eigen::Vector2f closest_point(
          center.x() + (goal.x() - center.x()) / goal_mag * abs(radius),
          center.y() + (goal.y() - center.y()) / goal_mag * abs(radius)
        );

    // uncomment for debugging
    // visualization::DrawCross(P, .3, 0xab4865, local_viz_msg_);
    // visualization::DrawLine(goal, P, 0, local_viz_msg_);

    // check for potential collisions with all points in the point cloud
    for (Vector2f point : point_cloud_) {
      // visualize point cloud, uncomment for debugging
      // visualization::DrawPoint(point, 0xB92348, local_viz_msg_);

      double mag = magnitude(point.x() - center.x(), point.y() - center.y());
      
      // calc arcs for rear left and front right car corners
      double r_1 = radius - w;
      double r_2 = magnitude(radius + w, h);
      // whether point is in front or behind vehicle
      double theta = atan2(point.x(), radius - point.y());
      // if radius is negative things get weird
      if (radius < 0) {
        r_1 = w - radius;
        r_2 = magnitude(abs(radius) + w, h);
        theta = atan2(point.x(), point.y() - radius);
      }

      if (mag >= r_1 && mag <= r_2 && theta > 0) {
        temp_fpl = min(
          radius * (theta - atan2(h, radius - w)),
          2 * abs(radius) * asin(magnitude(closest_point.x(), closest_point.y()) / abs(2 * radius))
        );
        if(radius < 0) temp_fpl = min(
          abs(radius * (theta - atan2(h, abs(radius) + w))),
          2 * abs(radius) * asin(magnitude(closest_point.x(), closest_point.y()) / (2 * radius))
        );

        // only save the smallest free path length for each curvature
        if (temp_fpl < path_options[index].free_path_length) {
          path_options[index].free_path_length = temp_fpl;
          path_options[index].curvature = i;
          path_options[index].obstruction = point;
          path_options[index].closest_point = closest_point;

        }

        // uncomment for debugging
        // visualization::DrawPathOption(
        //               i,
        //               path_options[index].free_path_length,
        //               0,
        //               0,
        //               true,
        //               local_viz_msg_
        //               );
      }
    }
    
    // calculate clearance around obstacle
    double dtgoal = magnitude(goal.x(), goal.y());
    arc_score = path_options[index].free_path_length + dtgoal*1;
    if (arc_score > best_arc_score) {
      best_path_option = path_options[index];
      best_arc_score = arc_score;
    }

  }    
    index++;
    return best_path_option;
  }


}  // namespace navigation

