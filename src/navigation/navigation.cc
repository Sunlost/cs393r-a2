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
#include <algorithm>
#include <cstdint>
#include <eigen3/Eigen/src/Core/GenericPacketMath.h>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Geometry/Rotation2D.h>
#include <eigen3/Eigen/src/Geometry/Transform.h>
#include <eigen3/Eigen/src/Geometry/Translation.h>
#include <math.h>
#include <sys/types.h>
#include <tuple>
#include "simple_queue.h"

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

bool debug_print; //maybe should make separate flags for convenience.

} //namespace

namespace navigation {

PathOption prev_path;
double prev_score;
unsigned prev_path_id;

PathOption curr_path;
double curr_score;
unsigned curr_path_id;

const double safety_margin = .05;
const double h = 0.4295 + safety_margin;
const double w = (0.281 / 2) + safety_margin;
const Eigen::Vector2f map_goal(20, 0);
const double clearance_cap = .005;
//                  curvature clearance   fpl     obstruction                 closest
const PathOption empty = {  0, -INFINITY, INFINITY, Eigen::Vector2f(0,0), Eigen::Vector2f(0,0) };

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

  debug_print = false;
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

  curr_path = pick_arc();
  visualization::DrawPathOption(curr_path.curvature,
                                curr_path.free_path_length,
                                curr_path.clearance,
                                0x3EB489,
                                true,
                                local_viz_msg_);

  if (debug_print) printf("\nchosen path's fpl %f\n", curr_path.free_path_length);
  if (debug_print) printf("chosen path's clearance %f\n", curr_path.clearance);
  if (debug_print) printf("chosen path's curvature %f\n", curr_path.curvature);
  if (debug_print) printf("chosen path's optimal cutoff %f, %f\n", curr_path.closest_point.x(), curr_path.closest_point.y());
  if (debug_print) printf("chosen path's obstruction %f, %f\n\n", curr_path.obstruction.x(), curr_path.obstruction.y());

  // predict current position, odometry
  position_prediction();

  // going to make a change, save prev path vars.
  if (curr_path.free_path_length != INFINITY) {
    d_max = curr_path.free_path_length;
    d_curr = 0;
    drive_msg_.curvature = curr_path.curvature;
    if (debug_print) printf("curvature set to %f\n", drive_msg_.curvature);
    d_curr_pred = 0;
    phase = PHASE_ACCEL;
    // prepare for next cycle
    prev_path = curr_path;
    prev_score = curr_score;
    prev_path_id = curr_path_id;
    curr_path = empty;
    curr_score = -1;
    curr_path_id = 17 ; //def out of bounds for our arcs
  }

  cycle_num++;

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

/* calculate magnitude */
float magnitude(float x, float y) {
  return sqrt(pow(x, 2) + pow(y, 2));
}

/* for loop of arcs, for each arc, calc score, return the best arc */
PathOption Navigation::pick_arc() {
  // for loop of arcs
  // for each arc, calc score
  // return the best arc
  float path_score = 0.0; 
  float best_path_score = -1;
  unsigned best_path_id = 0;
  double temp_fpl = INFINITY;
  Vector2f temp_endpt(INFINITY, INFINITY);
  vector<double> curves;
  vector<PathOption> path_options;
  PathOption best_path_option = PathOption();
  best_path_option.free_path_length = INFINITY;

  Eigen::Vector2f robot_rel_goal = map_goal;

  // curvature options from right to left
  // -ve y is the right direction, +ve y is the left direction
  // max curvature is 1

  for(double i = -1; i <= 1; i += 0.1) {
    bool mirrored = false;

    PathOption path_i = PathOption();
    double radius = 1 / (i + 1e-6); // adding small value to account for 0 curvature

    // calculate all arcs based on positive curvature
    if (radius < 0.0) {
      mirrored = true;
      radius = -1 * radius; 
    } 

    Eigen::Vector2f center(0, radius);
    double goal_mag = magnitude(robot_rel_goal.x() - center.x(), robot_rel_goal.y() - center.y());

    // fpl = f(c, p) if c > 0
    // fpl = f(-c, Mirrored(p)) if c < 0, we flip our curve back into the +ve, and flip points in cloud y

    // Straight path init vals
    Eigen::Vector2f optimal_endpt(robot_rel_goal.x(), 0);
    path_i.free_path_length = robot_rel_goal.x();
    // path is curved.
    optimal_endpt.x() = center.x() + (robot_rel_goal.x() - center.x()) / goal_mag * radius;
    optimal_endpt.y() = center.y() + (robot_rel_goal.y() - center.y()) / goal_mag * radius;
    path_i.free_path_length = (2 * radius) * asin(magnitude(optimal_endpt.x(), optimal_endpt.y()) / (2 * radius)); // init to some high value
    path_i.clearance = INFINITY;
    path_i.curvature = i; 

    // uncomment for debugging goal
    // visualization::DrawCross(optimal_endpt, .3, 0xab4865, local_viz_msg_);
    // visualization::DrawLine(robot_rel_goal, optimal_endpt, 0, local_viz_msg_);

    // check for potential collisions with all points in the point cloud
    for (Vector2f point : point_cloud_) {
      // visualize point cloud, uncomment for debugging
      // visualization::DrawPoint(point, 0xB92348, local_viz_msg_);
      Eigen::Vector2f eval_point(point.x(), point.y());

      // flip the point cloud across the x axis
      if (mirrored) eval_point.y() = -1 * point.y();

      // now the math should work as we know it should, for curves
      double r_1 = radius - w;
      double r_2 = magnitude(radius + w, h);
      double omega = atan2(h, radius - w);
      double mag = magnitude(eval_point.x() - center.x(), eval_point.y() - center.y());
      double theta = atan2(eval_point.x(), radius - eval_point.y());
      double phi = (theta - omega);

      // rotate r from (0,r) by phi radians to find the endpoint
      Eigen::Affine2f rotate_phi = Eigen::Translation2f(0, radius) * Eigen::Rotation2Df(phi);
      Eigen::Vector2f circle_center(0, -radius);
      Eigen::Vector2f obstructed_endpt = rotate_phi * circle_center; // this wrong, need to use affines.

      // this point is an obstruction for this path
      if ((mag >= r_1 && mag <= r_2) && theta > 0) {
        // check the optimal fpl math.
        double obstructed_fpl = radius * phi;
        Eigen::Affine2f translate_center = Eigen::Translation2f(0, -radius) * Eigen::Rotation2Df(0);
        Eigen::Vector2f center_optimal_endpt = translate_center * optimal_endpt;
        double optimal_central_angle = abs(atan(center_optimal_endpt.x() / center_optimal_endpt.y()));
        double optimal_fpl = optimal_central_angle * radius;

        // truncate free path length if necessary
        if (obstructed_fpl < optimal_fpl) {
          temp_fpl = obstructed_fpl;
          temp_endpt = obstructed_endpt;
        } else {
          temp_fpl = optimal_fpl;
          temp_endpt = optimal_endpt;
        }

        if (temp_fpl < path_i.free_path_length) {
          path_i.free_path_length = min(path_i.free_path_length, temp_fpl); // need to do same 3 way min for end of path
          temp_endpt.y() = (mirrored) ? -1 * temp_endpt.y() : temp_endpt.y();
          path_i.closest_point = temp_endpt;
          path_i.obstruction = point;
        }
      }
      // calculate clearance
      else if((mag < r_1 || mag > r_2) && theta > 0) { 
        // handle negative clearance values? as is mentioned on the slide
        double temp_clear = (mag < r_1) ? abs(r_1 - mag)  : abs(mag - r_2);
        if (temp_clear > clearance_cap) 
          temp_clear = clearance_cap; // set clearance to c_max because we don't care at some point
        if (temp_clear < path_i.clearance)
          path_i.clearance = temp_clear;
      }
    }

    path_options.push_back(path_i);
  }

  // run thru all feasible paths, score them.
  if (debug_print) printf("\nIn Pick Arc, feasible arcs\n");
  
  for(unsigned i = 0; i < path_options.size(); i++) {
    visualization::DrawPathOption(path_options.at(i).curvature,
                                  path_options.at(i).free_path_length,
                                  path_options.at(i).clearance,
                                  0,
                                  false,
                                  local_viz_msg_);

    visualization::DrawCross(path_options.at(i).closest_point, .3, 0xab0065, local_viz_msg_);

    // use robot_rel_goal and optimal_endpt
    double dtgoal = magnitude(robot_rel_goal.x() - path_options.at(i).closest_point.x(), 
                              robot_rel_goal.y() - path_options.at(i).closest_point.y());

    path_score = (path_options.at(i).clearance * 15) + (path_options.at(i).free_path_length)  + (dtgoal * .01);

    if (debug_print) {
       printf("Arc %d, curvature = %f, fpl = %f, clearance = %f, dtg = %f, closest point = %f, %f\n", 
       i, path_options.at(i).curvature, path_options.at(i).free_path_length, path_options.at(i).clearance, dtgoal, path_options.at(i).closest_point.x(), path_options.at(i).closest_point.y());
    }

    if (path_score > best_path_score) {
      best_path_id = i;
      best_path_option = path_options.at(i);
      best_path_score = path_score;
    }
  }

  if (debug_print) printf("\n");
  if (debug_print) printf("\nAt the end of Pick Arc\n");
  if (debug_print) printf("Previous Score %f\n", prev_score);
  if (debug_print) printf("Current Score %f, Arc %d \n\n", best_path_score, best_path_id);

  if (prev_path_id == best_path_id) {
    return empty;
  } else {
    //already have removed all -ve fpl'd paths
    curr_score = best_path_score;
    curr_path = best_path_option;
    curr_path_id = best_path_id;
    return best_path_option;
  }
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
  if(debug_print) printf("d_travelled: %f, d_curr %f, d_max %f\n", d_travelled, d_curr, d_max);
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
    if(new_v_f > v_max) new_v_f = v_max;
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
  if (debug_print) printf("drive_msg_.velocity = %f\n", drive_msg_.velocity);

  // 5. save past state
  prev_loc = odom_loc_;
  // TODO: replace 0 with set curvature value
  toc_queue.Push(cycle_num, v_f - v_i, 0);
  if(debug_print) printf("pushed %f to queue with value %ld. now size %d.\n", v_f - v_i, cycle_num, toc_queue.Size());

  return;
}

}  // namespace navigation
