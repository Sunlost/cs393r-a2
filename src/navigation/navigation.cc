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
#include <algorithm>
#include <cstdint>
#include <eigen3/Eigen/src/Core/GenericPacketMath.h>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <sys/types.h>
#include <tuple>
#include "simple_queue.h"

using Eigen::Vector2f;
using amrl_msgs::AckermannCurvatureDriveMsg;
using amrl_msgs::VisualizationMsg;
using std::string;
using std::vector;

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

float d_curr;
float d_max;

float v_max;
float a_max;
float decel_max;

float x_pred;
float y_pred;
float theta_pred;
float v_i_pred;
float d_curr_pred;

uint8_t cycles_per_second;
float cycle_time;
uint64_t cycle_num;

Vector2f prev_loc;
SimpleQueue<uint64_t, float, float> toc_queue;
uint8_t toc_queue_size;

bool debug_print;

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

  // phase of 1dTOC we are currently in
  phase = PHASE_ACCEL;

  // distance we have travelled so far
  d_curr = 0;
  // distance we want to go
  d_max = 3.65;

  // max velocity: 1.0 m/s
  v_max = 1.0;
  // max acceleration: 4.0 m/s^2
  a_max = 4.0;
  // max deceleration: 4.0 m/s^2
  decel_max = -4.0;

  x_pred = 0;
  y_pred = 0;
  theta_pred = 0;

  cycles_per_second = 20;
  cycle_time = (float) 1 / cycles_per_second;
  cycle_num = 0;

  prev_loc = Vector2f(0, 0);
  toc_queue_size = 1; // assume 0.15s latency @ 0.05s/cycle = 3 cycles

  debug_print = true;
}

void Navigation::SetNavGoal(const Vector2f& loc, float angle) {
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
  if(debug_print) printf("new odometry! (x,y,v_x,v_y): %f, %f, %f, %f\n", loc.x(), loc.y(), vel.x(), vel.y());
  robot_omega_ = ang_vel;
  robot_vel_ = vel;
  if (!odom_initialized_) {
    odom_start_angle_ = angle;
    odom_start_loc_ = loc;
    odom_initialized_ = true;
    odom_loc_ = loc;
    prev_loc = loc;
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

  // The control iteration goes here. 
  // Feel free to make helper functions to structure the control appropriately.
  
  // The latest observed point cloud is accessible via "point_cloud_"

  // Eventually, you will have to set the control values to issue drive commands:
  // drive_msg_.curvature = ...;
  // drive_msg_.velocity = ...;

  cycle_num++;

  // predict current position, odometry
  position_prediction();

  // invoke curve-generator

  // pass curves into cost function

  // handle 1d toc
  toc1dstraightline();

  // Add timestamps to all messages.
  local_viz_msg_.header.stamp = ros::Time::now();
  global_viz_msg_.header.stamp = ros::Time::now();
  drive_msg_.header.stamp = ros::Time::now();
  // Publish messages.
  viz_pub_.publish(local_viz_msg_);
  viz_pub_.publish(global_viz_msg_);
  drive_pub_.publish(drive_msg_);
}

void Navigation::position_prediction() {
  // TODO: add additional past command logging to reconcile with LIDAR delay

  // 1. get/set actual car movement/velocity
  float v_i = hypot(robot_vel_.x(), robot_vel_.y());
  float d_travelled = sqrt(pow((odom_loc_.x() - prev_loc.x()), 2) + pow((odom_loc_.y() - prev_loc.y()), 2));
  d_curr = d_curr + d_travelled;
  x_pred = odom_loc_.x();
  y_pred = odom_loc_.y();
  theta_pred = odom_angle_;

  if(debug_print) printf("\n");
  if(debug_print) printf("prev_loc(x,y): %f, %f\n", prev_loc.x(), prev_loc.y());
  if(debug_print) printf("odom_loc_(x,y): %f, %f\n", odom_loc_.x(), odom_loc_.y());
  if(debug_print) printf("d_travelled: %f, d_curr %f\n", d_travelled, d_curr);
  if(debug_print) printf("v_i is now %f\n", v_i);

  // 2. predict our future position, building off of actual movement/velocity
  d_curr_pred = d_curr;
  v_i_pred = v_i;
  if(debug_print) printf("cycle_num: %ld, toc_queue_size + 0x1UL: %ld, actual queue size: %d\n", cycle_num, toc_queue_size + 0x1UL, toc_queue.Size());
  if(cycle_num > toc_queue_size + 0x1UL) {
    if(debug_print) printf("POPPED! cycle %ld\n", cycle_num);
    toc_queue.Pop();
  }
  for(unsigned i = 0; i < toc_queue.Size(); i++) {
    // get values out of queue
    float v_delta = std::get<1>(toc_queue.values_[i]);
    // float c = std::get<2>_(toc_queue.values[i]);
    // predict new velocity
    float new_v_f = v_delta + v_i_pred;
    if(new_v_f < 0) new_v_f = 0;
    if(new_v_f > 1) new_v_f = 1;
    // predict new distance
    float d_delta;
    if(v_delta > 0) d_delta = (pow(new_v_f, 2) - pow(v_i_pred, 2)) / (2 * a_max);
    else if(v_delta == 0) d_delta = new_v_f * cycle_time;
    else d_delta = (pow(new_v_f, 2) - pow(v_i_pred, 2)) / (2 * decel_max);
    // // predict new x, y, angle
    // double radius = (1 / (path.curvature + 1e-6));
    // double radians = radius * d_delta;
    // double a_x = 
    // update v, d predictions
    v_i_pred = new_v_f;
    d_curr_pred += d_delta;
    if(debug_print) printf("pred v_delta = %f, new_v_f now = %f\n", v_delta, new_v_f);
    if(debug_print) printf("pred d_delta = %f, d_curr_pred now = %f\n", d_delta, d_curr_pred);
  }
}

// calculate what phase of ToC we are in, update state
// TODO: add curvature as a parameter
void Navigation::toc1dstraightline() {
  // formulas used:
    // v_f = v_i + at
    // d = (v_f^2 - v_i^2) / (2a)
    // t = (V_f - V_i) / a
    // d = vt

  // initial velocity
  float v_i = v_i_pred;
  // final velocity
  float v_f = 0;
  // distance we will travel in this cycle
  float d_this_cycle = 0;
  // total distance after this cycle ends
  float d_total_after_this_cycle = 0;
  // total distance after we fully decel to zero
  float d_total_after_decel_to_zero = 0;

  // 1. calculate which phase we're in
  if(phase != PHASE_DECEL) phase = (v_i == v_max) ? PHASE_CRUISE : PHASE_ACCEL;
  
  // 2. predict future state
  tocPhases new_phase = phase;
  float v_i2 = 0;
  float v_f2 = 0;
  switch(phase) {
    case PHASE_ACCEL:
      if(debug_print) printf("ACCEL PHASE\n");
      v_f = v_i + (a_max * cycle_time);
      float d_accel;
      float d_at_max_vel;
      if(v_f <= v_max) {
        d_accel = (pow(v_f, 2) - pow(v_i, 2)) / (2 * a_max);
        d_at_max_vel = 0;
      } else {
        v_f = v_max;
        d_accel = (pow(v_f, 2) - pow(v_i, 2)) / (2 * a_max);
        float t_accel = (v_f - v_i) / a_max;
        float t_at_max_vel = cycle_time - t_accel; // max vel for rest of cycle
        d_at_max_vel = v_f * t_at_max_vel;
      }
      d_this_cycle = d_accel + d_at_max_vel;
      d_total_after_this_cycle = d_curr_pred + d_this_cycle;

      v_i2 = v_f;
      v_f2 = 0;
      d_total_after_decel_to_zero = (v_f2 - pow(v_i2, 2)) / (2 * decel_max)
                                    + d_total_after_this_cycle;
      if(d_total_after_decel_to_zero > d_max) new_phase = PHASE_DECEL;
    break;

    case PHASE_CRUISE:
      if(debug_print) printf("CRUISE PHASE\n");
      v_f = v_max;
      d_this_cycle = (v_max / cycles_per_second);
      d_total_after_this_cycle = d_curr_pred + d_this_cycle;
      d_total_after_decel_to_zero = (0 - pow(v_i, 2)) / (2 * decel_max) 
                + d_total_after_this_cycle;
      if(d_total_after_decel_to_zero > d_max) new_phase = PHASE_DECEL;
    break;

    case PHASE_DECEL:
      if(debug_print) printf("ORG DECEL PHASE\n");
      v_f = v_i + (decel_max * cycle_time);
      if(v_f < 0) v_f = 0;
      d_this_cycle = (pow(v_f, 2) - pow(v_i, 2)) / (2 * decel_max);
      d_total_after_this_cycle = d_this_cycle + d_curr_pred;
      d_total_after_decel_to_zero = (0 - pow(v_i, 2)) / (2 * decel_max)
                                + d_total_after_this_cycle;
    break;

    default:
      assert(0); // should never occur
    break;
  } 

  // 3. check if our prediction changed to decel
  if(phase != new_phase) {
    if(debug_print) printf("SWAPPED TO DECEL PHASE\n");
    v_f = v_i + (decel_max * cycle_time);
    if(v_f < 0) v_f = 0;
    d_this_cycle = (pow(v_f, 2) - pow(v_i, 2)) / (2 * decel_max);
    d_total_after_this_cycle = d_this_cycle + d_curr_pred;
    d_total_after_decel_to_zero = (0 - pow(v_i, 2)) / (2 * decel_max)
                              + d_total_after_this_cycle;
    phase = new_phase;
  }

  // 4. act on predictions, update internal state
  switch(phase) {
    case PHASE_ACCEL:
      drive_msg_.velocity = v_f;
    break;

    case PHASE_CRUISE:
      drive_msg_.velocity = v_f;
    break;

    case PHASE_DECEL:
      drive_msg_.velocity = v_f;
    break;

    default:
      assert(0); // should never occur
    break;
  }

  // 5. save past state
  prev_loc = odom_loc_;
  // TODO: replace 0 with set curvature value
  toc_queue.Push(cycle_num, v_f - v_i, 0);
  if(debug_print) printf("pushed %f to queue with value %ld. now size %d.\n", v_f - v_i, cycle_num, toc_queue.Size());

  return;
}


}  // namespace navigation
