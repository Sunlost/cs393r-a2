#ifndef NAVIGATION_PARAMS_H
#define NAVIGATION_PARAMS_H

namespace navigation {
struct MotionLimits {
  float max_acceleration;
  float max_deceleration;
  float max_speed;

  MotionLimits(float max_acceleration, float max_deceleration, float max_speed)
      : max_acceleration(max_acceleration),
        max_deceleration(max_deceleration),
        max_speed(max_speed) {}
};

struct NavigationParams {
  NavigationParams()
      : dt(0.05),
        system_latency(0.1),
        linear_limits(1.0, 1.0, 1.0),
        angular_limits(0.5, 0.5, 0.5),
        max_curvature(1.0),
        max_path_length(0.5),
        max_clearance(0.5),
        goal_tolerance(1.0),
        robot_length(0.5),
        robot_width(0.5),
        robot_wheelbase(0.5),
        obstacle_margin(0.5),
        lidar_offset(0.5) {}

  double dt;
  double system_latency;
  
  // Parameters for the path sampler
  MotionLimits linear_limits;
  MotionLimits angular_limits;

  float max_curvature;
  float max_path_length;
  float max_clearance;

  // weights for the path evaluator
  float clearance_weight;
  float arc_length_weight;
  float distance_weight;

  double goal_tolerance;

  float robot_length;
  float robot_width;
  float robot_wheelbase;
  float obstacle_margin;
  float lidar_offset;

};
}  // namespace navigation
#endif  // NAVIGATION_PARAMS_H
