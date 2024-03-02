#ifndef ACKERMANN_MOTION_PRIMITIVES_H
#define ACKERMANN_MOTION_PRIMITIVES_H

#include <cmath>
#include <memory>
#include <vector>

#include "constant_curvature_arc.h"
#include "eigen3/Eigen/Dense"
#include "shared/math/math_util.h"
#include "navigation_params.h"

using std::fabs;
using std::max;
using std::min;
using math_util::Sq;
using math_util::Sign;
using math_util::Pow;

namespace motion_primitives {

class AckermannSampler {
 public:
  AckermannSampler() : linear_vel_(0), angular_vel_(0), local_target_(0, 0) {}
  AckermannSampler(const navigation::NavigationParams& nav_params)
      : nav_params_(nav_params),
        linear_vel_(0),
        angular_vel_(0),
        local_target_(0, 0) {}
  ~AckermannSampler() = default;

  void update(const Eigen::Vector2f& new_vel,
              const float new_ang_vel,
              const Eigen::Vector2f& new_local_target,
              const std::vector<Eigen::Vector2f>& new_point_cloud);

  void setPathLength(std::shared_ptr<ConstantCurvatureArc> path_ptr);
  void checkObstacles(std::shared_ptr<ConstantCurvatureArc> path_ptr);

  std::vector<std::shared_ptr<ConstantCurvatureArc>> getSamples(int n);

 private:
  navigation::NavigationParams nav_params_;

  float linear_vel_;
  float angular_vel_;
  Eigen::Vector2f local_target_;
  std::vector<Eigen::Vector2f> point_cloud_; // base_link frame
};

class AckermannEvaluator {
 public:
  AckermannEvaluator(const navigation::NavigationParams& nav_params)
      : nav_params_(nav_params) {}
  ~AckermannEvaluator() = default;

  void update(const Eigen::Vector2f& new_local_target);

  std::shared_ptr<ConstantCurvatureArc> findBestPath(
      std::vector<std::shared_ptr<ConstantCurvatureArc>>& samples);

  float evaluatePath(std::shared_ptr<ConstantCurvatureArc> path_ptr);

 private:
  navigation::NavigationParams nav_params_;

  Eigen::Vector2f local_target_;
};
}  // namespace motion_primitives

#endif  // ACKERMANN_MOTION_PRIMITIVES_H