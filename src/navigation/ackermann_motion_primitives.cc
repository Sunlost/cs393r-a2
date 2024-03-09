#include "ackermann_motion_primitives.h"

#include <iostream>

using std::cout;
using std::endl;

namespace motion_primitives {

void AckermannSampler::update(const Eigen::Vector2f& new_vel,
                              const float new_ang_vel,
                              const Eigen::Vector2f& new_local_target,
                              const std::vector<Eigen::Vector2f>& new_point_cloud) {
  linear_vel_ = new_vel.norm();
  angular_vel_ = new_ang_vel;
  local_target_ = new_local_target;
  point_cloud_ = new_point_cloud;

  if (FLAGS_v > 1) {
    cout << "================= [Navigation Sampler] Update ================" << endl;
    cout << "Linear velocity: " << linear_vel_ << endl;
    cout << "Angular velocity: " << angular_vel_ << endl;
    cout << "Local target: " << local_target_.transpose() << endl;
    cout << "==============================================================\n" << endl;
  }
}

std::vector<std::shared_ptr<ConstantCurvatureArc>> AckermannSampler::getSamples(int n) {
  std::vector<std::shared_ptr<ConstantCurvatureArc>> samples;
  // Whatever curvature we can achieve while accelerating we can achieve while
  // decelerating, so only consider acceleration
  const float max_dtheta_dot =  // Find angular velocity limit
      nav_params_.linear_limits.max_acceleration * nav_params_.max_curvature *
      nav_params_.dt;
  const float max_ds_dot = nav_params_.linear_limits.max_acceleration *
                           nav_params_.dt;  // Find maximum velocity

  float cmax = nav_params_.max_curvature;
  float cmin = -nav_params_.max_curvature;
  // Curvature limits in params only apply when we have 0 velocity
  if (linear_vel_ > max_ds_dot) {
    cmin = max(cmin, (angular_vel_ - max_dtheta_dot) / (linear_vel_ + max_ds_dot));
    cmax = min(cmax, (angular_vel_ + max_dtheta_dot) / (linear_vel_ + max_ds_dot));  //
  }
  const float dc = (cmax - cmin) / (n - 1);

  // 1 Take dynamics into for min and max curvature if necessary
  for (float c = cmin; c <= cmax; c += dc) {
    // 2 Compute the arc length
    auto sample = std::make_shared<ConstantCurvatureArc>(
        c, nav_params_.max_path_length, nav_params_.max_clearance);
    setPathLength(sample);
    checkObstacles(sample);

    samples.push_back(sample);
  }

  if (FLAGS_v > 1) {
    cout << "==================== [Navigation] Sampler ====================" << endl;
    cout << "Number of samples: " << samples.size() << endl;
    cout << "cmax: " << cmax << ", cmin: " << cmin << ", dc: " << dc << endl;
    cout << "Linear velocity: " << linear_vel_ << endl;
    cout << "Angular velocity: " << angular_vel_ << endl;
    cout << "range Angular velocity: " << angular_vel_ - max_dtheta_dot << " to "
         << angular_vel_ + max_dtheta_dot << endl;
    cout << "Local target: " << local_target_.transpose() << endl;
    cout << "==============================================================\n" << endl;
  }

  return samples;
}

void AckermannSampler::setPathLength(std::shared_ptr<ConstantCurvatureArc> path_ptr) {
  // Linear motion
  if (fabs(path_ptr->curvature()) < 1e-5) {
    path_ptr->set_arc_length(fmin(local_target_.x(), nav_params_.max_path_length));
    return;
  }

  const float radius = 1 / path_ptr->curvature();
  Eigen::Vector2f instant_center(0, radius);
  Eigen::Vector2f instant_center_to_goal =
      fabs(radius) * (local_target_ - instant_center).normalized();
  const float theta =
      atan2(fabs(instant_center_to_goal.x()), fabs(instant_center_to_goal.y()));
  const float arc_length = fabs(radius * theta);

  path_ptr->set_arc_length(fmin(arc_length, nav_params_.max_path_length));
}

void AckermannSampler::checkObstacles(std::shared_ptr<ConstantCurvatureArc> path_ptr) {
  static const bool kDebug = false;

  const float l = nav_params_.robot_length + 2 * nav_params_.obstacle_margin;
  const float w = nav_params_.robot_width + 2 * nav_params_.obstacle_margin;
  const float l_f = l - (l - nav_params_.robot_wheelbase) / 2;  // base to front
  const float l_r = l - l_f;                                    // base to rear

  // Add special case to handle when car is driving nearly straight
  if (fabs(path_ptr->curvature()) < 1e-5) {
    // Only check if points are within the width of the car
    const float min_y = -w / 2.0f;
    const float max_y = w / 2.0f;

    for (const auto& point : point_cloud_) {
      // Point is outside the lateral swept volume of the car
      if (point.y() < min_y || point.y() > max_y) {
        const float clearance = fmin(fabs(point.y() - min_y), fabs(point.y() - max_y));
        path_ptr->set_clearance(fmin(clearance, path_ptr->clearance()));
        continue;
      }

      // Point is outside the forward swept volume of the car
      if (point.x() > path_ptr->arc_length() + l_f || point.x() < -l_r) {
        continue;
      }

      // Point in inside the swept volume of the car
      float arc_length = fmax(point.x() - l_f, 0);
      path_ptr->set_arc_length(fmin(arc_length, path_ptr->arc_length()));
      path_ptr->set_clearance(0);
      break;
    }
    return;
  }

  // compute volume swept by the car during a single time (depends on curvature)
  const float r = 1 / path_ptr->curvature();
  const float r_base_min = fabs(r) - (w / 2.0f);
  const float r_base_max = fabs(r) + (w / 2.0f);
  const float r_front_min = sqrt(Sq(r_base_min) + Sq(l_f));
  const float r_front_max = sqrt(Sq(r_base_max) + Sq(l_f));
  const float r_rear_max = sqrt(Sq(r_base_max) + Sq(l_r));

  if (kDebug) {
    cout << "Checking obstacles for curvature " << path_ptr->curvature() << endl;
  }

  const Eigen::Vector2f instant_center(0, r);
  for (const auto& point : point_cloud_) {
    if (point.x() < -l_r) {
      // Case 1: Point is behind the car
      continue;
    }

    const float r_p = (point - instant_center).norm();
    if (r_p < r_base_min || r_p > r_front_max) {
      // Case 2: Point is out of sweep volume
      const float clearance = fmin(fabs(r_p - r_base_min), fabs(r_p - r_front_max));
      path_ptr->set_clearance(fmin(clearance, path_ptr->clearance()));
      continue;
    }

    int phase = 0;
    Eigen::Vector2f point_to_collision(0, 0);
    if (r_p < r_front_min) {
      // Case 3: Point is in swept volume in inside of car
      point_to_collision.y() = Sign(r) * (w / 2.0f);
      point_to_collision.x() = sqrt(Sq(r_p) - Sq(point_to_collision.y() - r));
      if (kDebug) {
        cout << point.transpose() << " is inner area of weep volume" << endl;
        phase = 3;
      }
    } else if (r_p > r_base_max && r_p < r_rear_max && fabs(point.x()) < l_r) {
      // Case 4: Point is in sweapt volum e between base_link and rear bumper of car
      point_to_collision.y() = -Sign(r) * (w / 2.0f);
      point_to_collision.x() = sqrt(Sq(r_p) - Sq(point_to_collision.y() - r));
      if (kDebug) {
        cout << point.transpose() << " is in rear area of sweep volume" << endl;
        phase = 4;
      }
    } else if (r_p < r_front_max && point.x() > l_f) {
      // Case 5: Point is in swept volume in front of car
      point_to_collision.x() = l_f;
      point_to_collision.y() = r - Sign(r) * sqrt(Sq(r_p) - Sq(point_to_collision.x()));
      if (kDebug) {
        cout << point.transpose() << " is in front area of sweep volume" << endl;
        phase = 5;
      }
    } else {
      // Case 6: Point is out of swept volume
      const float clearance = fabs(r_rear_max - r_p);
      path_ptr->set_clearance(fmin(clearance, path_ptr->clearance()));
      if (kDebug) {
        cout << point.transpose() << " is out of sweep volume" << endl;
        phase = 6;
      }
      continue;
    }

    // Compute the arc length until the collision
    const Eigen::Vector2f collision_radial = point_to_collision - instant_center;
    const Eigen::Vector2f point_radial = point - instant_center;
    const float swap_theta = acos(point_radial.dot(collision_radial) /
                                  (point_radial.norm() * collision_radial.norm()));
    const float arc_length = fabs(r) * swap_theta;

    if (kDebug) {
      cout << "Curvature: " << path_ptr->curvature()
           << ", Arc length to collision: " << arc_length << endl;
      cout << "Point to collision: " << point_to_collision.transpose()
           << ", r_p: " << r_p << ", r: " << r << ", swap_theta: " << swap_theta
           << ", phase: " << phase << endl;
      cout << phase << endl;
    }

    // Penalize clearance as well if we have a collision before arc length is traversed
    if (arc_length < path_ptr->arc_length()) {
      path_ptr->set_arc_length(arc_length);
      path_ptr->set_clearance(0);
    }
  }
}

}  // namespace motion_primitives

namespace motion_primitives {
void AckermannEvaluator::update(const Eigen::Vector2f& new_local_target) {
  local_target_ = new_local_target;
}

std::shared_ptr<ConstantCurvatureArc> AckermannEvaluator::findBestPath(
    std::vector<std::shared_ptr<ConstantCurvatureArc>>& samples) {
  if (samples.empty()) return nullptr;

  // return path that have the highest score
  std::shared_ptr<ConstantCurvatureArc> best_path = nullptr;

  if (FLAGS_v > 1) {
    cout << "==================== [Navigation Evaluator] All ===================="
         << endl;
    cout << "Local target: " << local_target_.transpose() << endl;
  }
  float best_score = 0;
  for (const auto& path : samples) {
    const float score = evaluatePath(path);

    if (best_path == nullptr || score > best_score) {
      best_path = path;
      best_score = score;
    }
  }

  if (FLAGS_v > 1) {
    cout << "================= [Navigation Evaluator] Best =================" << endl;
    cout << "Best path" << endl;
    cout << "Curvature: " << best_path->curvature()
         << ", Arc length: " << best_path->arc_length()
         << ", Clearance: " << best_path->clearance() << ", Score: " << best_score
         << endl;
    cout << "==============================================================\n" << endl;
  }

  return best_path;
}  // namespace motion_primitives

float AckermannEvaluator::evaluatePath(std::shared_ptr<ConstantCurvatureArc> path) {
  const float clearance = path->clearance();                                  // good
  const float arc_length = path->arc_length();                                // good
  const float dist_to_target = (local_target_ - path->getEndPoint()).norm();  // bad

  const float score = nav_params_.clearance_weight * clearance +
                      nav_params_.arc_length_weight * arc_length -
                      nav_params_.distance_weight * dist_to_target;
  // const float score = nav_params_.arc_length_weight * arc_length;

  if (FLAGS_v > 1) {
    cout << "Curvature: " << path->curvature() << ", Arc length: " << path->arc_length()
         << ", End point: (" << path->getEndPoint().transpose() << ")"
         << ", Clearance: " << path->clearance()
         << ", Distance to target: " << dist_to_target << ", Score: " << score << endl;
  }
  return score;
}

}  // namespace motion_primitives
