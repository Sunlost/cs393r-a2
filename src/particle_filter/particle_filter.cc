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
using vector_map::VectorMap;

DEFINE_double(num_particles, 8, "Number of particles");

namespace particle_filter {

config_reader::ConfigReader config_reader_({"config/particle_filter.lua"});

ParticleFilter::ParticleFilter() :
    prev_odom_loc_(0, 0),
    prev_odom_angle_(0),
    odom_initialized_(false) {}

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

  // TODO: figure out what frame of reference to use lawl
  // right now I'm using car's
  std::vector<Eigen::Vector2f>& scan = *scan_ptr;
  Eigen::Vector2f laser_loc(loc.x() + 0, loc.y() + 0.2);
  // Compute what the predicted point cloud would be, if the car was at the pose
  // loc, angle, with the sensor characteristics defined by the provided
  // parameters.
  // This is NOT the motion model predict step: it is the prediction of the
  // expected observations, to be used for the update step.

  scan.resize(num_ranges);
  // iterate through scan to set points for each ray
  for (size_t i = 0; i < scan.size(); ++i) {
    // iterate through map to check for collisions
    for (size_t j = 0; j < map_.lines.size(); ++j) {
      // to check for collisions, construct a line2f from range_min to range_max, in the direction of the ray, centered around laser pose
      const line2f map_line = map_.lines[j];
      // need to use math to calculate endpoint of the line based on range_max. treat laser location as point 0
      // 10 is a magic number rn describing the angle increment. we should tune that (and by extension num_ranges)
      float alpha = angle_min + 10*i;
      // the range max point
      Eigen::Vector2f rm_pt(range_max * sin(alpha), range_max * cos(alpha));
      line2f my_line(laser_loc.x(), laser_loc.y(), rm_pt.x(), rm_pt.y());
      printf("\n[GETPPC]: rmminusparticle x: %f rmminusparticle y: %f\n",
        rm_pt.x() - loc.x(), rm_pt.y() - loc.y());
        printf("[GETPPC]: particleloc x: %f particleloc y: %f\n",
    loc.x(), loc.y());
      // check for intersection with this line and the map line
      
      Vector2f intersection_point; // Return variable
      bool intersects = map_line.Intersection(my_line, &intersection_point);
      
      if (intersects) {
        // if intersection exists, "first" collision wins
        scan[i] = intersection_point;
        printf("[GETPPC]: intersection_point x: %f intersection_point y: %f\n",
        intersection_point.x(), intersection_point.y());
        break;
      } else {
        // else if no collision, set scan[i] to the point at range_max
        scan[i] = rm_pt;
        printf("%ld [GETPPC]: rm_pt x: %f rm_pt y: %f\n",
    i, rm_pt.x(), rm_pt.y());
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
    int num_ranges = ranges.size() / 10;
    std::vector<Eigen::Vector2f> scan_ptr(num_ranges);
    // pass in robot's location and angle for loc and angle
    // Eigen::Vector2f loc;
    // float angle;

    // Eigen::Vector2f laser_loc(loc.x() + 0, loc.y() + 0.2);

    // GetLocation(&loc, &angle); // TODO: doubt this is correct -sun
    // use map-relative location actually
    // num_ranges should equal something like (angle_max - angle_min) / 10, so every 10 degrees we use the lidar range
    // use GetPredictedPointCloud to predict expected observations for particle conditioned on the map
    // GetPredictedPointCloud((Eigen::Vector2f const) loc, angle, num_ranges, range_min, range_max, angle_min, angle_max, &scan_ptr);
    printf("[UPDATE]: p_loc x: %f, p_loc y: %f, p_loc angle: %f\n", p_ptr->loc.x(), p_ptr->loc.y(), p_ptr->angle);
    GetPredictedPointCloud(p_ptr->loc, p_ptr->angle, num_ranges, range_min, range_max, angle_min, angle_max, &scan_ptr);
  // compare particle observation to prediction
  double log_lik = 0.0;
  // tunable param: sd_squared
  float variance = calc_variance(ranges);
  // robustification will be on 14 - Expecting The Unexpected slide 28
  for (size_t i = 0; i < scan_ptr.size(); ++i) {
    // s_hat is (dist btwn laser and scan[i] points) - range_min
    // TUNABLE: check if this should be sqnorm instead of norm if particle filter is slow
    // TODO: check if I need to subtract laser_loc.x and .y from here (and also figure out how to calc that)
    double s_hat_dist = sqrt(pow(scan_ptr[i].x(), 2) + pow(scan_ptr[i].y(), 2));
    float s_hat = s_hat_dist - range_min;
    // s is the range, aka dist from laser to endpoint of observed
    float s = ranges[i * 10]; // TUNABLE: every 10th laser?
    log_lik -= pow((s - s_hat), 2) / variance;
     printf("%ld [UPDATE]: s_hat: %f s: %f log_lik: %f\n",
    i, s_hat, s, log_lik);

  }

  // assign weight to particle
  p_ptr->weight = log_lik;
}

void ParticleFilter::Resample() {
  vector<Particle> new_particles;

  double step_size = 1.0 / FLAGS_num_particles;
  double last = rng_.UniformRandom(0, step_size) - step_size;
  int index = 0;
  double sum = 0;

  while(new_particles.size() < FLAGS_num_particles) {
    sum += particles_[index].weight;
    while(last + step_size < sum) {
      Particle p;
      p.loc.x() = particles_[index].loc.x();
      p.loc.y() = particles_[index].loc.y();
      p.angle = particles_[index].angle;
      // assuming we do not need to copy weight through...
      p.weight = particles_[index].weight;
      new_particles.push_back(p);
      last += step_size;
    }
    index++;
  }

  // After resampling:
  particles_ = new_particles;
}

void ParticleFilter::ObserveLaser(const vector<float>& ranges,
                                  float range_min,
                                  float range_max,
                                  float angle_min,
                                  float angle_max) {
  // A new laser scan observation is available (in the laser frame)
  // Call the Update and Resample steps as necessary.
  cout << "me when I observe laser" << endl;
  // TODO STEP 1: figure out how to call update
  // init a num_updates variable
  int num_updates = 0;
  // tunable param: d
  // if we have traveled at least distance d and sensor data is available
      // use for normalization later
      double max_likelihood = 0;
      // loop over particle vector
        for (size_t i = 0; i < FLAGS_num_particles; ++i) {
          ParticleFilter::Update(ranges, range_min, range_max, angle_min, angle_max, &(particles_[i]));
          double likelihood = particles_[i].weight;
          if(max_likelihood < likelihood) max_likelihood = likelihood;
        }
        // normalize all weights (see 16 - Problems with Particle Filters slide 10)
        for (size_t i = 0; i < FLAGS_num_particles; ++i) {
          particles_[i].weight = abs((particles_[i].weight - max_likelihood));
        }
        num_updates++;

  // TODO STEP 2: figure out how to call resample
  // tunable param: n
  int n = 10; // FIX/Move somewhere else later
  // if it has been n updates since our last resample
  if (num_updates >= n) ParticleFilter::Resample();
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


  // You will need to use the Gaussian random number generator provided. For
  // example, to generate a random number from a Gaussian with mean 0, and
  // standard deviation 2:
  // float x = rng_.Gaussian(0.0, 2.0);
  // printf("Random number drawn from Gaussian distribution with 0 mean and "
  //        "standard deviation of 2 : %f\n", x);

  // calc x_hat, y_hat, and theta_hat from odom_loc and odom_angle
  double x_hat = odom_loc.x() - prev_odom_loc_.x();
  double y_hat = odom_loc.y() - prev_odom_loc_.y();
  double theta_hat = odom_angle - prev_odom_angle_;
  // printf("[PREDICT]: x_hat: %f = %f - %f\n y_hat: %f = %f - %f\n theta_hat: %f = %f - %f\n",
  //   x_hat, odom_loc.x(), prev_odom_loc_.x(), y_hat, odom_loc.y(), prev_odom_loc_.y(), theta_hat, odom_angle, prev_odom_angle_);
  // so for some reason, x_hat, y_hat, and theta_hat are kinda crazy. 
  // printf("x: %f y: %f theta: %f\n", 
  //       odom_loc.x(),
  //       odom_loc.y(),
  //       odom_angle);
  // printf("prevx: %f prevy: %f prevtheta: %f\n", 
  // prev_odom_loc_.x(),
  // prev_odom_loc_.y(),
  // prev_odom_angle_);
  
  // tunable parameters
  float k_1 = 0.01; //   x,y stddev's   mag(x,y) weight
  float k_2 = 0.01; //   x,y stddev's   mag(theta) weight
  float k_3 = 0.01; // theta stddev's   mag(x,y) weight
  float k_4 = 0.01; // theta stddev's   mag(theta) weight

  // e_x, e_y drawn from N(0, k1*sqrt(d_x^2 + d_y^2) + k2*||d_theta||)
  // e_theta drawn from N(0, k3*sqrt(d_x^2 + d_y^2) + k4*||d_theta||)
  // same distrib except for TUNABLE params k1+k2/k3+k4
  double xy_stddev = k_1*magnitude(x_hat, y_hat) + k_2*abs(theta_hat);
  double theta_stddev = k_3*magnitude(x_hat, y_hat) + k_4*abs(theta_hat);
  if (!odom_initialized_) {
      xy_stddev = 0;
      theta_stddev = 0;
  }
  // update particles
  for (size_t i = 0; i < FLAGS_num_particles; ++i){
    float e_x = rng_.Gaussian(0.0, xy_stddev);
    float e_y = rng_.Gaussian(0.0, xy_stddev);
    float e_theta = rng_.Gaussian(0.0, theta_stddev);
    
    particles_[i].loc.x() += x_hat + e_x;
    particles_[i].loc.y() += y_hat + e_y;
    particles_[i].angle += theta_hat + e_theta;
    printf("particle x: %f particle y: %f particle angle: %f\n", 
        particles_[i].loc.x(),
        particles_[i].loc.y(),
        particles_[i].angle);
    printf("e x: %f e y: %f e angle: %f\n", 
        e_x,
        e_y,
        e_theta);
      printf("x hat: %f y hat: %f angle hat: %f\n", 
      x_hat,
      y_hat,
      theta_hat);
  }
  odom_initialized_ = true;
  prev_odom_loc_ = odom_loc;
  prev_odom_angle_ = odom_angle;
}

void ParticleFilter::Initialize(const string& map_file,
                                const Vector2f& loc,
                                const float angle) {
  // The "set_pose" button on the GUI was clicked, or an initialization message
  // was received from the log. Initialize the particles accordingly, e.g. with
  // some distribution around the provided location and angle.
  map_.Load(map_file);
  
  // ? differences are still big so let's presume that ... well hold on let's see xhat first
  prev_odom_loc_.x() = loc.x();
  prev_odom_loc_.y() = loc.y();
  prev_odom_angle_ = angle;

  // problem is not in initialize bc these seem fine.
  // printf("xxxxx: %f yyyyyy: %f aaaaangle: %f\n", 
  //       loc.x(),
  //       loc.y(),
  //       angle);

// printf("prevx: %f prevy: %f prevangle: %f\n", 
//         prev_odom_loc_.x(),
//         prev_odom_loc_.y(),
//         prev_odom_angle_);

  odom_initialized_ = false;
  // initialize vector of particles with GetParticles
  for (size_t i = 0; i < FLAGS_num_particles; ++i){
    Particle p = Particle();
    // init in Gaussian distribution around loc and angle
    p.loc.x() = loc.x() + rng_.Gaussian(0.0, 1.0);
    p.loc.y() = loc.y() + rng_.Gaussian(0.0, 1.0);
    p.angle = angle + rng_.Gaussian(0.0, 1.0);
    p.weight = abs(rng_.Gaussian(0.0, 1.0));
    particles_.push_back(p);
  }
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
    // printf("xwhy: %f ywhy: %f anglewhy: %f\n", 
    //     particles_[i].loc.x(),
    //     particles_[i].loc.y(),
    //     particles_[i].angle);
    
    x_locs += particles_[i].loc.x() * particles_[i].weight;
    y_locs += particles_[i].loc.y() * particles_[i].weight;
    sines += sin(particles_[i].angle);
    cosines += cos(particles_[i].angle);
  }
  // printf("xlocs: %f ylocs: %f sines: %f\n", 
  //       x_locs,
  //       y_locs,
  //       sines);
  // I have temporarily removed the *100 bc the numbers look normal but what do I know
  loc.x() = x_locs;
  loc.y() = y_locs;
  angle = atan2(sines / FLAGS_num_particles, cosines / FLAGS_num_particles);
  printf("getloc x: %f getlocy: %f getlocangle: %f\n", 
        loc.x(),
        loc.y(),
        angle);
}


}  // namespace particle_filter
