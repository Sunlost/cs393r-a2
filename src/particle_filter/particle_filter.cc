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

DEFINE_double(num_particles, 50, "Number of particles");

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
      // Access the end points using `.p0` and `.p1` members:
      printf("P0: %f, %f P1: %f,%f\n", 
            my_line.p0.x(),
            my_line.p0.y(),
            my_line.p1.x(),
            my_line.p1.y());
      // check for intersection with this line and the map line
      Vector2f intersection_point; // Return variable
      bool intersects = map_line.Intersection(my_line, &intersection_point);
      if (intersects) {
        // if intersection exists, "first" collision wins
        scan[i] = intersection_point;
        break;
      } else {
        // else if no collision, set scan[i] to the point at range_max
        scan[i] = rm_pt;
      }
    }
  }
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

  // TODO: STEP 1
  // use GetPredictedPointCloud to predict expected observations for particle conditioned on the map
    // init a scan_ptr vector of points... size should be num_ranges
    // option to make scan_ptr the length of ranges or to make them like every 10th ray
    // tunable param: num_ranges->aka the 10 thing.
    int num_ranges = (angle_max - angle_min) / 10;
    std::vector<Eigen::Vector2f> scan_ptr(num_ranges);
    // pass in robot's location and angle for loc and angle
    Eigen::Vector2f loc;
    float angle;

    Eigen::Vector2f laser_loc(loc.x() + 0, loc.y() + 0.2);

    ParticleFilter::GetLocation(&loc, &angle); // TODO: doubt this is correct -sun
    // use map-relative location actually
    // num_ranges should equal something like (angle_max - angle_min) / 10, so every 10 degrees we use the lidar range
    ParticleFilter::GetPredictedPointCloud((Eigen::Vector2f const) loc, angle, num_ranges, range_min, range_max, angle_min, angle_max, &scan_ptr);

  // TODO: STEP 2
  // compare particle observation to prediction
  // ranges holds our observation s, scan_ptr (kinda) holds our prediction s_hat
  double log_lik = 0.0;
  // tunable param: sd_squared
  float sd_squared = 0.0025; // TODO: did the slides say to start at 0.05 std dev? unsure.
  // robustification will be on 14 - Expecting The Unexpected slide 28
  for (size_t i = 0; i < scan_ptr.size(); ++i) {
    // s_hat is (dist btwn laser and scan[i] points) - range_min
    // TUNABLE: check if this should be sqnorm instead of norm if particle filter is slow
    double s_hat_dist = sqrt(pow(scan_ptr[i].x() - laser_loc.x(), 2) + pow(scan_ptr[i].y() - laser_loc.y(), 2));
    float s_hat = s_hat_dist - range_min;
    // s is the range, aka dist from laser to endpoint of observed
    float s = ranges[i * 10]; // TUNABLE: every 10th laser?
    log_lik -= pow((s_hat - s) / sd_squared, 2);
  }

  // TODO: STEP 3
  // assign weight to particle
  // particle struct has field weight, set that to log_lik?
  /* Dr. Biswas's answer from 2/14
     These [point-cloud point location] estimations are based on the particles that we are modeling, 
     trying to guess what the point cloud would look like from that particle's perspective 
     
     If the particle's predicted point_cloud matches with the "ground truth" (what is this, ask Amanda)
     that is considered a "good" particle, and "bad" vice-versa. This will inform our weightings on the 
     particles. 
     
     need more clarification on what this "ground truth" is and how we use it...
  */
  p_ptr->weight = log_lik;
}

void ParticleFilter::Resample() {
  // Resample the particles, proportional to their weights.
  // The current particles are in the `particles_` variable. 
  // Create a variable to store the new particles, and when done, replace the
  // old set of particles:
  // vector<Particle> new_particles';
  // During resampling: 
  //    new_particles.push_back(...)
  // After resampling:
  // particles_ = new_particles;

  /*
    ED: Either we keep some tunable absolute threshold for the weights for discarding
    old unlikely particles 
    OR
    could figure out the range of of weights we have just calculated, figure out threshold and remove particles. 
    Don't filter out too many particles, maybe start by dropping 30% of the range of weights? 
    
    The latter sounds more correct.

    also need to clarify what they mean about "duplication", whether that is copying over an old point 
    from the old particles to the new ones 
    or 
    actually duplicating the particles s.t where there was one particle present in the old set, now there would be two. 
    
    The former sounds like the correct way.
  */ 

  // You will need to use the uniform random number generator provided. For
  // example, to generate a random number between 0 and 1:
  float x = rng_.UniformRandom(0, 1);
  printf("Random number drawn from uniform distribution between 0 and 1: %f\n",
         x);
}

void ParticleFilter::ObserveLaser(const vector<float>& ranges,
                                  float range_min,
                                  float range_max,
                                  float angle_min,
                                  float angle_max) {
  // A new laser scan observation is available (in the laser frame)
  // Call the Update and Resample steps as necessary.

  // TODO STEP 1: figure out how to call update
  // init a num_updates variable
  int num_updates = 0;
  // tunable param: d
  // if we have traveled at least distance d and sensor data is available
      // use for normalization later
      double max_likelihood = 0;
      // loop over particle vector
        for (size_t i = 0; i < FLAGS_num_particles; ++i) {
          ParticleFilter::Update(ranges, range_min, range_max, angle_min, angle_max, &particles_[i]);
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
  // if it has been n updates since our last resample
    // ParticleFilter::Resample()
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
  
  // tunable parameters
  float k_1 = 1.0; //   x,y stddev's   mag(x,y) weight
  float k_2 = 1.0; //   x,y stddev's   mag(theta) weight
  float k_3 = 1.0; // theta stddev's   mag(x,y) weight
  float k_4 = 1.0; // theta stddev's   mag(theta) weight

  // e_x, e_y drawn from N(0, k1*sqrt(d_x^2 + d_y^2) + k2*||d_theta||)
  // e_theta drawn from N(0, k3*sqrt(d_x^2 + d_y^2) + k4*||d_theta||)
  // same distrib except for TUNABLE params k1+k2/k3+k4
  double xy_stddev = k_1*magnitude(x_hat, y_hat) + k_2*abs(theta_hat);
  double theta_stddev = k_3*magnitude(x_hat, y_hat) + k_4*abs(theta_hat);
  // update particles
  for (size_t i = 0; i < FLAGS_num_particles; ++i){
    float e_x = rng_.Gaussian(0.0, xy_stddev);
    float e_y = rng_.Gaussian(0.0, xy_stddev);
    float e_theta = rng_.Gaussian(0.0, theta_stddev);
    particles_[i].loc.x() += x_hat + e_x;
    particles_[i].loc.y() += y_hat + e_y;
    particles_[i].angle += theta_hat + e_theta;
  }

}

void ParticleFilter::Initialize(const string& map_file,
                                const Vector2f& loc,
                                const float angle) {
  // The "set_pose" button on the GUI was clicked, or an initialization message
  // was received from the log. Initialize the particles accordingly, e.g. with
  // some distribution around the provided location and angle.
  map_.Load(map_file);

  // initialize vector of particles with GetParticles
  for (size_t i = 0; i < FLAGS_num_particles; ++i){
    Particle *p = new Particle();
    // init in Gaussian distribution around loc and angle
    p.loc.x() = loc.x() + rng_.Gaussian(0.0, );
    p.loc.y() = loc.y() + rng_.Gaussian(0.0, );
    p.angle = angle + rng_.Gaussian(0.0, );
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
  loc = Vector2f(0, 0);
  angle = 0;

  double x_locs = 0.0;
  double y_locs = 0.0;
  double sines = 0.0;
  double cosines = 0.0;

  for (size_t i = 0; i < FLAGS_num_particles; ++i){
    x_locs += particles_[i].loc.x() * particles_[i].weight;
    y_locs += particles_[i].loc.y() * particles_[i].weight;
    sines += sin(angle);
    cosines += cos(angle);
  }
  loc.x() = x_locs * 100 / FLAGS_num_particles;
  loc.y() = y_locs * 100 / FLAGS_num_particles;
  angle = atan2(sines / FLAGS_num_particles, cosines / FLAGS_num_particles);
}


}  // namespace particle_filter
