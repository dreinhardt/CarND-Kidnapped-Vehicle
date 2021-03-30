/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using namespace std;

using std::string;
using std::vector;
using std::normal_distribution;


void ParticleFilter::init(double gps_x, double gps_y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 30;  // TODO: Set the number of particles
  // with around 15 particles, the algo will not complete the round!
  
  // Random Value Generator
  std::default_random_engine gen;
  
  normal_distribution<double> dist_x(gps_x, std[0]);  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_y(gps_y, std[1]);  // Create normal distributions for y
  normal_distribution<double> dist_theta(theta, std[2]);  // Create normal distributions for theta

  for (int i = 0; i < num_particles; ++i) {

    // TODO: Sample from these normal distributions like this: 
    //   sample_x = dist_x(gen);
    //   where "gen" is the random engine initialized earlier.
    Particle particle;
    
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.id = i;
    particle.weight = 1.0;

    weights.push_back(particle.weight);
    particles.push_back(particle);
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  // See Lesson 5 - 9. Calculate Prediction...

  std::default_random_engine gen;

  for (int i = 0; i < num_particles; ++i) {
    if(yaw_rate == 0){
       throw std::exception(); }
    
    // We add the measurements to each particle
    Particle p_f; // particle "future"
    p_f.x = particles[i].x + (velocity / yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
    p_f.y = particles[i].y + (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
    p_f.theta = particles[i].theta + (yaw_rate * delta_t);
  
    // We add random Gaussian noise to each particle
    normal_distribution<double> dist_x(p_f.x, std_pos[0]);
    normal_distribution<double> dist_y(p_f.y, std_pos[1]);
    normal_distribution<double> dist_theta(p_f.theta, std_pos[2]);
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }
}



void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  // See overview in Lesson 5, chapter 15
  
  weights.clear();	// flush weights at the beginning


  for(int i=0; i<num_particles; i++){

    // Transform coordinates from car coords to map coords - see lesson 5, chapter 16
    vector<LandmarkObs> predictions; // my mapmark predictions (already transformed because of Map)
    for(unsigned int mm=0; mm < map_landmarks.landmark_list.size(); mm++){    
      if ( dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[mm].x_f, map_landmarks.landmark_list[mm].y_f ) < sensor_range){
        LandmarkObs mapmarks;  
        mapmarks.x = map_landmarks.landmark_list[mm].x_f;
        mapmarks.y = map_landmarks.landmark_list[mm].y_f;
        mapmarks.id = map_landmarks.landmark_list[mm].id_i;
        predictions.push_back(mapmarks);
      }
    }

    vector<LandmarkObs> t_observations; // my "transformed" landmark observations
    for(unsigned int lm=0; lm < observations.size(); lm++){
      LandmarkObs landmarks;
      landmarks.x = particles[i].x + (cos(particles[i].theta) * observations[lm].x) - (sin(particles[i].theta) * observations[lm].y);
      landmarks.y = particles[i].y + (sin(particles[i].theta) * observations[lm].x) + (cos(particles[i].theta) * observations[lm].y);
      landmarks.id = lm;
      t_observations.push_back(landmarks);
    }


    // Update weights
    particles[i].weight = 1.0;		// reset weight to have a starting neutral number 1 to be multiplied

    for (unsigned int o = 0; o < predictions.size(); o++) {
      int m_ID = -1;	// NaN
      double m_distance = 4*sensor_range;	// large number, 4 times the sensor_range

      for (unsigned int p = 0; p < t_observations.size(); p++) {
        double curr_dist = dist(t_observations[p].x, t_observations[p].y, predictions[o].x, predictions[o].y);

        if (curr_dist< m_distance){		// dataAssociation, search the shortest distance
          m_distance = curr_dist;		// always fit to the closest
          m_ID = p;						// the save best observation
        } 
      }

      if (m_ID != -1){
        double w = ( 1/(2*M_PI*std_landmark[0]*std_landmark[1])) * 
          exp( -( pow(t_observations[m_ID].x-predictions[o].x,2) / (2*pow(std_landmark[0], 2))
              + ( pow(t_observations[m_ID].y-predictions[o].y,2) / (2*pow(std_landmark[1], 2))) ) );		// Formla from the lessons
        particles[i].weight = particles[i].weight * w;		// multiply all weights with each other
      }
    }

    weights.push_back(particles[i].weight);		// Save the last weight!
    
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  // see lesson 4, chapter 20
  vector<Particle> random_particles;  // buffer weights
  std::default_random_engine gen;

  double mw = *max_element(weights.begin(), weights.end());  // Get largest weight
  double beta = 0.0; // my starting beta

  std::discrete_distribution<int> ddist(0, num_particles-1);
  int index = ddist(gen);   // my random starting index

  for (int i=0; i<num_particles; i++){
    std::uniform_real_distribution<double> udist(0.0, mw);
    beta = beta + udist(gen)* 2.0 * mw;
    
    while (weights[index] < beta) {    // the wheel
        beta = beta - weights[index];
        index = (index + 1) % num_particles;
      }
    random_particles.push_back(particles[index]);
  }
  
  particles = random_particles;   // set new particles
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}