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
\file    particle-filter.cc
\brief   Particle Filter Starter Code
\author  Joydeep Biswas, (C) 2019
*/
//========================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "shared/math/geometry.h"
#include "shared/math/line2d.h"
#include "shared/math/math_util.h"
#include "shared/util/timer.h"

#include "config_reader/config_reader.h"
#include "particle_filter.h"

#include "vector_map/vector_map.h"

using geometry::line2f;
using std::cout;
using std::endl;
using std::string;
using std::swap;
using std::vector;
using Eigen::Vector2f;
using Eigen::Vector2i;
using math_util::DegToRad;
using math_util::RadToDeg;

using vector_map::VectorMap;

DEFINE_double(num_particles, 40, "Number of particles");

namespace particle_filter {

config_reader::ConfigReader config_reader_({"config/particle_filter.lua"});

ParticleFilter::ParticleFilter() :
    prev_odom_loc_(0, 0),
    prev_odom_angle_(0),
    odom_initialized_(false),
    sum_weight(FLAGS_num_particles),
    num_valid_particles(FLAGS_num_particles),
    num_updates_done(0),
    num_updates_reqd_for_resample(3),
    ith_ray(10),
    deg_offset(ith_ray * laser_ang_deg_res),
    dist_travelled(0),
    debug_print(false) {}

void ParticleFilter::GetParticles(vector<Particle>* particles) const {
  *particles = particles_;
}







void ParticleFilter::GetPredictedPointCloud(const Vector2f& loc,
                                            const float angle,
                                            int num_ranges,
                                            float range_min,
                                            float range_max,
                                            float angle_min,
                                            float angle_max,
                                            std::vector<Eigen::Vector2f>* scan_ptr) {

  std::vector<Eigen::Vector2f>& scan = *scan_ptr;
  // Eigen::Vector2f laser_loc(loc.x() +  0.2, loc.y());
  // printf("[GPPC]: laser_loc.x() = %f, laser_loc.y() = %f\n",laser_loc.x(), laser_loc.y());
  Eigen::Vector2f laser_loc(loc.x(), loc.y());
  // Compute what the predicted point cloud would be, if the car was at the pose
  // loc, angle, with the sensor characteristics defined by the provided
  // parameters.
  // This is NOT the motion model predict step: it is the prediction of the
  // expected observations, to be used for the update step.

  //printf("Scan ptr size before %ld\n", scan.size());
  scan.resize(num_ranges / ith_ray);
  //printf("Scan ptr size before %ld\n\n", scan.size());

  vector<float> scan_min_dists;
  for(int i = 0; i < num_ranges / ith_ray; i++) scan_min_dists.push_back(std::numeric_limits<float>::max());
 
  float degtoradresult = math_util::DegToRad(deg_offset);
  // iterate through scan to set points for each ray
  for (size_t j = 0; j < map_.lines.size(); ++j) {
    // iterate through map to check for collisions
    float alpha = angle + angle_min - degtoradresult;
    Eigen::Vector2f rm_pt(0, 0);
    for (size_t i = 0; i < scan.size(); ++i) {
      // to check for collisions, construct a line2f from range_min to range_max, in the direction of the ray, centered around laser pose
      const line2f map_line = map_.lines[j];
      // need to use math to calculate endpoint of the line based on range_max. treat laser location as point 0
      // 10 is a magic number rn describing the angle increment. we should tune that (and by extension num_ranges)
      alpha += degtoradresult;
      // the range max point
      // Eigen::Vector2f rm_pt(range_max * sin(alpha) + laser_loc.x(), );
      rm_pt.x() = range_max * cos(alpha) + laser_loc.x();
      rm_pt.y() = range_max * sin(alpha) + laser_loc.y();
      line2f my_line(laser_loc.x(), laser_loc.y(), rm_pt.x(), rm_pt.y());

      // check for intersection with this line and the map line
      Vector2f intersection_point; // Return variable
      bool intersects = map_line.Intersection(my_line, &intersection_point);
      
      if (intersects) {
        // if intersection exists, "first" collision wins
        float curr_dist = pow(intersection_point.x() - laser_loc.x(), 2) + pow(intersection_point.y() - laser_loc.y(), 2);
        if (curr_dist < scan_min_dists[i]) {
          scan_min_dists[i] = curr_dist;
          scan[i] = intersection_point;
        }
      } else if(scan_min_dists[i] == std::numeric_limits<float>::max()) {
        scan[i] = rm_pt;
        scan_min_dists[i] = std::numeric_limits<float>::max() - 1;
      }
    }

  }
}








float calc_variance(const vector<float>& ranges) {
  float avg = 0;
  for(size_t i = 0; i < ranges.size(); ++i) {
    avg += ranges[i];
  }
  avg /= ranges.size();
  float variance = 0;
  for(size_t i = 0; i < ranges.size(); ++i) {
    variance += (ranges[i] - avg) * (ranges[i] - avg);
  }
  
  return variance /= ranges.size();
}







void ParticleFilter::Update(const vector<float>& ranges,
                            float range_min,
                            float range_max,
                            float angle_min,
                            float angle_max,
                            Particle* p_ptr) {

  // Implement the update step of the particle filter here.
  // You will have to use the `GetPredictedPointCloud` to predict the expected
  // observations for each particle, and assign weights to the particles based
  // on the observation likelihood computed by relating the observation to the
  // predicted point cloud.

  // tunable param: num_ranges->aka the 10 thing.
  int num_ranges = ranges.size() / ith_ray;
  std::vector<Eigen::Vector2f> scan_ptr(num_ranges);
  for(int i = 0; i < num_ranges; i++) {
    scan_ptr[i] = Vector2f(0, 0);
  }
  // pass in robot's location and angle for loc and angle

  // use map-relative location actually
  // num_ranges should equal something like (angle_max - angle_min) / 10, so every 10 degrees we use the lidar range
  // use GetPredictedPointCloud to predict expected observations for particle conditioned on the map
  // GetPredictedPointCloud((Eigen::Vector2f const) loc, angle, num_ranges, range_min, range_max, angle_min, angle_max, &scan_ptr);
  // printf("[UPDATE]: p_loc x: %f, p_loc y: %f, p_loc angle: %f\n", 
  //     p_ptr->loc.x(), p_ptr->loc.y(), p_ptr->angle);
  GetPredictedPointCloud(p_ptr->loc, p_ptr->angle, ranges.size(), range_min, range_max, angle_min, angle_max, &scan_ptr);
  // compare particle observation to prediction
  double log_lik = 0.0;
  // tunable param: sd_squared
  // divisor == (laser std dev^2 == 0.01)
  float divisor = 0.01;
  // robustification will be on 14 - Expecting The Unexpected slide 28

  float d_short = 1;
  float d_short_squared = pow(d_short, 2);
  float d_long = 8;
  float d_long_squared = pow(d_long, 2);
  float delta = 0;
  for (size_t i = 0; i < scan_ptr.size(); ++i) {
    // s_hat is (dist btwn laser and scan[i] points) - range_min
    // TUNABLE: check if this should be sqnorm instead of norm if particle filter is slow
    double s_hat_dist = sqrt(pow(scan_ptr[i].x() - p_ptr->loc.x(), 2) + pow(scan_ptr[i].y() - p_ptr->loc.y(), 2));
    float s_hat = s_hat_dist - range_min;
    // s is the range, aka dist from laser to endpoint of observed
    float s = ranges[i * ith_ray]; // TUNABLE: every 10th laser?
    if (s < range_min || s > range_max) {
      continue;
    } else if (s < s_hat - d_short) {
      delta = d_short_squared / divisor;
    } else if (s > s_hat + d_long) {
      delta = d_long_squared / divisor;
    } else {
      delta = pow((s - s_hat), 2) / divisor;
    }

    // printf("[UPDATE]: s_hat: %f,s: %f\n", 
    //   s_hat, s);
    if(delta > 1000000) {
      printf("SCARILY LARGE S DELTA DETECTED?\n");
      printf("scan_ptr x: %f, scan_ptr y: %f, pptr x: %f, pptr y: %f, s_hat: %f, s: %f, delta: %f\n", scan_ptr[i].x(), scan_ptr[i].y(), p_ptr->loc.x(), p_ptr->loc.y(), s_hat, s, delta);
    }
    log_lik += delta;
    // if(debug_print) printf("[UPDATE]: %ld s_hat: %f s: %f log_lik: %f\n",
    //     i, s_hat, s, log_lik);
  }
  // assign weight to particle, times gamma
  p_ptr->weight = log_lik * 2 / num_ranges * -1;
  if(p_ptr->weight > 1000000) {
      printf("RESULTING IN CONCERNINGLY LARGE PARTICLE WEIGHT?\n");
      printf("weight: %f, particle x: %f, particle y: %f, particle angle: %f\n", p_ptr->weight, p_ptr->loc.x(), p_ptr->loc.y(), p_ptr->angle);
      printf("log_lik: %f", log_lik);
  }
}







void ParticleFilter::Resample() {
  // printf("[RESAMPLE]\n");
  vector<Particle> new_particles;
  sum_weight = 0;
  for (size_t i = 0; i < FLAGS_num_particles; ++i) {
    if (particles_[i].weight == 0) continue;
    particles_[i].weight = pow(particles_[i].weight, 2.718);
    sum_weight += particles_[i].weight;
  }

  double step_size = sum_weight / FLAGS_num_particles;
  double new_weight = 1 / FLAGS_num_particles;
  double last = rng_.UniformRandom(0, step_size) - step_size;
  int index = 0;
  double sum = 0;


  while(new_particles.size() < FLAGS_num_particles) {
    if(particles_[index].weight == 0) {
      index++;
      continue;
    }
    sum += particles_[index].weight;
    while(last + step_size < sum) {
      Particle p;
      p.loc.x() = particles_[index].loc.x();
      p.loc.y() = particles_[index].loc.y();
      p.angle = particles_[index].angle;
      p.weight = new_weight;
      new_particles.push_back(p);
      last += step_size;
    }
    index++;
  }

  // After resampling:
  particles_ = new_particles;
  num_valid_particles = FLAGS_num_particles;
  num_updates_done = 0;
}








void ParticleFilter::ObserveLaser(const vector<float>& ranges,
                                  float range_min,
                                  float range_max,
                                  float angle_min,
                                  float angle_max) {
  // A new laser scan observation is available (in the laser frame)
  // Call the Update and Resample steps as necessary.
  // printf("[OBSERVELASER INVOKED]\n");
  // TODO STEP 1: figure out how to call update
  if (dist_travelled < .15) return;
  dist_travelled = 0;
  // tunable param: d
  // if we have traveled at least distance d and sensor data is available
      // use for normalization later
  
  sum_weight = 0;

  double max_likelihood = 0;
  // loop over particle vector
  for (size_t i = 0; i < FLAGS_num_particles; ++i) {
    if (particles_[i].weight == 0) continue;
    ParticleFilter::Update(ranges, range_min, range_max, angle_min, angle_max, &(particles_[i]));
    double likelihood = particles_[i].weight;
    if(max_likelihood > likelihood) max_likelihood = likelihood;
  }
  // normalize all weights (see 16 - Problems with Particle Filters slide 10)
  for (size_t i = 0; i < FLAGS_num_particles; ++i) {
    // printf("[OBSERVELASER]: max_like: %f weight: %f reset weight: %f\n", max_likelihood, particles_[i].weight, (particles_[i].weight - max_likelihood));
    if(particles_[i].weight == 0) continue;
    particles_[i].weight = abs((particles_[i].weight - max_likelihood));
    if(particles_[i].weight == 0) num_valid_particles--;
    sum_weight += particles_[i].weight;
  }
  num_updates_done++;

  // if it has been n updates since our last resample
  if (num_updates_done >= num_updates_reqd_for_resample) ParticleFilter::Resample();
}








double magnitude(float x, float y) {
  return sqrt(pow(x, 2) + pow(y, 2));
}







void ParticleFilter::Predict(const Vector2f& odom_loc,
                             const float odom_angle) {
  // Implement the predict step of the particle filter here.
  // A new odometry value is available (in the odom frame)
  // Implement the motion model predict step here, to propagate the particles
  // forward based on odometry.

  if(debug_print) printf("[PREDICT] GIVEN LOC // x: %f, y: %f, angle: %f\n", odom_loc.x(), odom_loc.y(), odom_angle);
  if (!odom_initialized_) {
    if(debug_print) printf("[PREDICT] FIRST PREDICT -- CANCELLED EARLY\n");
    prev_odom_loc_.x() = odom_loc.x();
    prev_odom_loc_.y() = odom_loc.y();
    prev_odom_angle_ = odom_angle;
    odom_initialized_ = true;
    return;
  }

  // You will need to use the Gaussian random number generator provided. For
  // example, to generate a random number from a Gaussian with mean 0, and
  // standard deviation 2:
  // float x = rng_.Gaussian(0.0, 2.0);
  // if(debug_print) printf("Random number drawn from Gaussian distribution with 0 mean and "
  //        "standard deviation of 2 : %f\n", x);


  // calc x_hat, y_hat, and theta_hat from odom_loc and odom_angle
  double x_hat = odom_loc.x() - prev_odom_loc_.x();
  double y_hat = odom_loc.y() - prev_odom_loc_.y();
  Eigen::Rotation2Df r_prev_odom(-prev_odom_angle_);

  Vector2f T_odom(x_hat, y_hat);
  dist_travelled += pow(x_hat, 2) + pow(y_hat, 2);
  Vector2f T_delta_bl = r_prev_odom * T_odom;

  double theta_hat = odom_angle - prev_odom_angle_;
  // if(debug_print) printf("[PREDICT] x_hat: %f = %f - %f\n y_hat: %f = %f - %f\n theta_hat: %f = %f - %f\n",
      // x_hat, odom_loc.x(), prev_odom_loc_.x(), y_hat, odom_loc.y(), prev_odom_loc_.y(), 
      // theta_hat, odom_angle, prev_odom_angle_);
  // so for some reason, x_hat, y_hat, and theta_hat are kinda crazy. 
  // if(debug_print) printf("[PREDICT] x: %f y: %f theta: %f\n", 
      // odom_loc.x(), odom_loc.y(), odom_angle);
  // if(debug_print) printf("[PREDICT] prevx: %f prevy: %f prevtheta: %f\n", 
      // prev_odom_loc_.x(), prev_odom_loc_.y(), prev_odom_angle_);
  
  // tunable parameters
  float k_1 = 0.0001; //   x,y stddev's   mag(x,y) weight
  float k_2 = 0.0001; //   x,y stddev's   mag(theta) weight
  float k_3 = 0.001; // theta stddev's   mag(x,y) weight
  float k_4 = 0.00001; // theta stddev's   mag(theta) weight

  // e_x, e_y drawn from N(0, k1*sqrt(d_x^2 + d_y^2) + k2*||d_theta||)
  // e_theta drawn from N(0, k3*sqrt(d_x^2 + d_y^2) + k4*||d_theta||)
  // same distrib except for TUNABLE params k1+k2/k3+k4

  // update particles
  for (size_t i = 0; i < FLAGS_num_particles; ++i){
    if(particles_[i].weight == 0) continue;
    Eigen::Rotation2Df r_map(particles_[i].angle);
    Vector2f T_map_one(particles_[i].loc.x(), particles_[i].loc.y());
    Vector2f T_map = T_map_one + r_map * T_delta_bl;

    double xy_stddev = k_1*magnitude(T_map.x(), T_map.y()) + k_2*abs(theta_hat);
    double theta_stddev = k_3*magnitude(T_map.x(), T_map.y()) + k_4*abs(theta_hat);

    float e_x = rng_.Gaussian(0.0, xy_stddev);
    float e_y = rng_.Gaussian(0.0, xy_stddev);
    float e_theta = rng_.Gaussian(0.0, theta_stddev);

    // if(i == 0) {
    //   // have one noiseless particle to preserve the dataset in the event we pass through a narrow gap.
    //   // this one noiseless particle will successfully navigate the gap and will repopulate 
    //   // the particles via resample
    //   e_x = 0;
    //   e_y = 0;
    //   e_theta = 0;
    // }

    particles_[i].loc.x() = T_map.x() + e_x;
    particles_[i].loc.y() = T_map.y() + e_y;
    particles_[i].angle += theta_hat + e_theta;
    if(debug_print) printf("[PREDICT] e x: %f e y: %f e angle: %f\n", 
        e_x, e_y, e_theta);
    if(debug_print) printf("[PREDICT] x hat: %f y hat: %f angle hat: %f\n", 
        x_hat, y_hat, theta_hat);
    if(debug_print) printf("[PREDICT] particle x: %f particle y: %f particle angle: %f\n", 
        particles_[i].loc.x(), particles_[i].loc.y(), particles_[i].angle);
    
    line2f particle_line(T_map_one.x(), T_map_one.y(), particles_[i].loc.x(), particles_[i].loc.y());

    for (size_t j = 0; j < map_.lines.size(); ++j) {
      const line2f map_line = map_.lines[j];
      Vector2f intersection_point; // Return variable
      bool intersects = map_line.Intersection(particle_line, &intersection_point);
      if (intersects) {
        particles_[i].weight = 0;
        num_valid_particles--;
      } 
    }
  }
  prev_odom_loc_ = odom_loc;
  prev_odom_angle_ = odom_angle;
}







void ParticleFilter::Initialize(const string& map_file,
                                const Vector2f& loc,
                                const float angle) {
  // The "set_pose" button on the GUI was clicked, or an initialization message
  // was received from the log. Initialize the particles accordingly, e.g. with
  // some distribution around the provided location and angle.

  if(debug_print) printf("INITIALIZED!\n");

  map_.Load(map_file);
  
  vector<Particle> new_particles;

  if(debug_print) printf("INITIALIZE DATA: x: %f, y: %f, angle: %f\n", loc.x(), loc.y(), angle);
  num_valid_particles = FLAGS_num_particles;
  odom_initialized_ = false;
  num_updates_done = 0;
  dist_travelled = 0;
  double new_weights = 1 / FLAGS_num_particles;
  // initialize vector of particles with GetParticles
  for (size_t i = 0; i < FLAGS_num_particles; ++i){
    Particle p = Particle();
    // init in Gaussian distribution around loc and angle
    p.loc.x() = loc.x() + rng_.Gaussian(0.0, 0.001);
    p.loc.y() = loc.y() + rng_.Gaussian(0.0, 0.001);
    p.angle = angle + rng_.Gaussian(0.0, 0.02);
    p.weight = new_weights;
    sum_weight += p.weight;
    new_particles.push_back(p);
  }

  particles_ = new_particles;
}







void ParticleFilter::GetLocation(Eigen::Vector2f* loc_ptr, 
                                 float* angle_ptr) const {
  Vector2f& loc = *loc_ptr;
  float& angle = *angle_ptr;
  // Compute the best estimate of the robot's location based on the current set
  // of particles. The computed values must be set to the `loc` and `angle`
  // variables to return them. Modify the following assignments:

  double x_locs = 0.0;
  double y_locs = 0.0;
  double sines = 0.0;
  double cosines = 0.0;

  for (size_t i = 0; i < FLAGS_num_particles; ++i){
    if(particles_[i].weight == 0) {
      // printf("getloc skip %ld\n", i);
      continue;
    }
    // printf("[GETLOCATION] xwhy: %f  ywhy: %f  anglewhy: %f  weight: %f\n", 
    //     particles_[i].loc.x(), particles_[i].loc.y(), particles_[i].angle, particles_[i].weight);
    // temporarily have removed the weighted mean calculation
    x_locs += particles_[i].loc.x();
    y_locs += particles_[i].loc.y();
    sines += sin(particles_[i].angle) / num_valid_particles;
    cosines += cos(particles_[i].angle) / num_valid_particles;
  }
  // printf("[GETLOCATION] x sum: %f  y sum: %f  sines: %f  cosines: %f  sum_weight: %f\n", 
  //     x_locs, y_locs, sines, cosines, sum_weight);
  loc.x() = x_locs / num_valid_particles;
  loc.y() = y_locs / num_valid_particles;
  angle = atan2(sines, cosines);
  // printf("[GETLOCATION] predicted x: %f   predicted y: %f   predicted angle: %f\n\n", 
  //     loc.x(), loc.y(), angle);
}


}  // namespace particle_filter
