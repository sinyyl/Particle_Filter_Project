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
#include <limits>
    
#include "helper_functions.h"

using std::string;
using std::vector;
using std::sin;
using std::cos;
using std::pow;

using std::endl;
using std::cout;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  // checl if it is initialized
  if(is_initialized) {
	return;
  } 
  // create random generator for gaussian distributation
  std::normal_distribution<double> distribution_x(x, std[0]);
  std::normal_distribution<double> distribution_y(y, std[1]);
  std::normal_distribution<double> distribution_theta(theta, std[2]);
  //std::default_random_engine generator;
  
  num_particles = 150;  // TODO: Set the number of particles
  // assign default values to the particles
  for(int i=0; i<num_particles; i++){
    // temperory variable of particle
    Particle thisParticle;
    thisParticle.id = i;
    // add gaussian noise
    thisParticle.x = distribution_x(generator);
    thisParticle.y = distribution_y(generator);
    thisParticle.theta = distribution_theta(generator);
    thisParticle.weight = 1.0;  // default weight is 1.0
    particles.push_back(thisParticle); // add to the particles list
  }
  //cout << "Initialized" << endl;
  // set the initialie flag to true
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
    //std::default_random_engine generator;
    // iterate through all particles
    //cout << "Before Predicted" << endl;
    for(int i=0; i<num_particles; i++){
      //Particle thisParticle = particles.at(i);

      // get temp variables
      double this_x = particles[i].x;
      double this_y = particles[i].y;
      double this_theta = particles[i].theta;
      // calculate the predicted x, y and theta
      double new_x;
      double new_y;
      // for low yaw rate
      if(abs(yaw_rate) < 0.0001){
        new_x = this_x +  velocity * delta_t * cos(this_theta);
        new_y = this_y +  velocity * delta_t * sin(this_theta);
      }else{
        // normal situation
        new_x = this_x +  velocity / yaw_rate * (sin(this_theta + yaw_rate * delta_t) - sin(this_theta));
        new_y = this_y +  velocity / yaw_rate * (cos(this_theta) - cos(this_theta + yaw_rate * delta_t));
      }
      double new_theta = this_theta + yaw_rate * delta_t;
      // add gaussian noise to the predicted data
      std::normal_distribution<double> distribution_x(new_x, std_pos[0]);
      std::normal_distribution<double> distribution_y(new_y, std_pos[1]);
      std::normal_distribution<double> distribution_theta(new_theta, std_pos[2]);
      // assign the data back to particle
      particles[i].x = distribution_x(generator);
      particles[i].y = distribution_y(generator);
      particles[i].theta = distribution_theta(generator);
    }
    //cout << "Predicted" << endl;
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
  // iterate through two vetors to find the nearest data points and associate them together
  // observations is the landmark truth position
  for(int i=0; i<predicted.size(); i++){
    double x_obs = predicted[i].x;  //read in the predicted x, y coordinates
    double y_obs = predicted[i].y;
    float minDistance = std::numeric_limits<float>::infinity();
    int associatedID;
    //vector<LandmarkObs> Newbservations;
    for(int j=0; j<observations.size(); j++){
      double x_pre = observations[j].x;  //read in the observations x, y coordinates
      double y_pre = observations[j].y;
      double thisDistance = dist(x_obs, y_obs, x_pre, y_pre);

      if(thisDistance <= minDistance){
        minDistance = thisDistance;    // update the new min
        associatedID = j;
      }
    }
    // use the index of landmark for wasier access later
    //Newbservations.push_back();
    observations[associatedID].id = predicted[i].id;
  }
  //observations = Newbservations;
  //cout << "Associated" << endl;
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
  // maximum weight
  double maxWeight = -1.0;
  vector<double> thisWeights;

  // iterate through all particles
  for(int i=0; i<num_particles; i++){
    // transfer the coordinates system from vehicle to map
    double thisX = particles[i].x;
    double thisY = particles[i].y;
    double thisTheta = particles[i].theta;


    // find the predicted landmarks in the measurement range
    vector<LandmarkObs> predictedLandmarks;
    for(int k=0; k<map_landmarks.landmark_list.size(); k++){
      double landmark_x = map_landmarks.landmark_list[k].x_f;
      double landmark_y = map_landmarks.landmark_list[k].y_f;
      // calculated the distance between landmark and predicted location
      double distance = dist(landmark_x, landmark_y, thisX, thisY);
      // if within range
      if(distance <= sensor_range){
        LandmarkObs thisPredictedLandmark;
        thisPredictedLandmark.x = landmark_x;
        thisPredictedLandmark.y = landmark_y;
        thisPredictedLandmark.id = map_landmarks.landmark_list[k].id_i;
        predictedLandmarks.push_back(thisPredictedLandmark);
      }
    }


    // iterate each observation data, transform and store in observedLandmarks
    vector<LandmarkObs> observedLandmarks;
    for(int j=0; j<observations.size(); j++){
      LandmarkObs thisLandmark;
      double thisObs_x = observations[j].x;
      double thisObs_y = observations[j].y;

      // calculate the transformation of coordinates
      double mapX = thisX + cos(thisTheta) * thisObs_x - sin(thisTheta) * thisObs_y;
      double mapY = thisY + sin(thisTheta) * thisObs_x + cos(thisTheta) * thisObs_y;
      thisLandmark.x = mapX;
      thisLandmark.y = mapY;
      // store the transformed landmark into a new list
      observedLandmarks.push_back(thisLandmark);
    }


    // using helper function to associate ovservation to landmatk truth
    dataAssociation(predictedLandmarks, observedLandmarks);


    // loop through the observations to add the association into particle struct
    vector<double> thisSense_x;
    vector<double> thisSense_y;
    vector<int> thisAssociation;
    for(int l=0; l<observedLandmarks.size(); l++){
      if(observedLandmarks[l].id != NULL){
        thisAssociation.push_back(observedLandmarks[l].id);
        thisSense_x.push_back(observedLandmarks[l].x); 
        thisSense_y.push_back(observedLandmarks[l].y);
      }
    }
    particles[i].associations = thisAssociation;
    particles[i].sense_x = thisSense_x; 
    particles[i].sense_y = thisSense_y;
    

    // calculate the guassian distribution
    // find the weight for each predicted landmark
    for(int j=0; j<predictedLandmarks.size(); j++){
      double mu_x = predictedLandmarks[j].x;
      double mu_y = predictedLandmarks[j].y;
        for(int k=0; k<particles[i].associations.size(); k++){
          if(particles[i].associations[k] == predictedLandmarks[j].id){
            double predX = particles[i].sense_x[k];
            double predY = particles[i].sense_y[k];
            particles[i].weight *= multiv_prob(std_landmark[0], std_landmark[1], predX, predY, mu_x, mu_y);
          }
        }
      }
    // add the weight into instance variable weights
    thisWeights.push_back(particles[i].weight);

    if(particles[i].weight > maxWeight){
      maxWeight = particles[i].weight;
    }
  }

    // change the weights instance variable to the current weight
    weights = thisWeights;
  
  
  // noamalize weight according to the maximum weight
  for(int i=0; i<num_particles; i++){
    // change both the weight in each particle and weights instance variable
    weights[i] = weights[i] / maxWeight;
    particles[i].weight = particles[i].weight / maxWeight;
  }
  maxWeight = 1.0;  // assign the max weight to 1

  //cout << "Updated" << endl;

}


void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine generator;
  vector<Particle> newParticles;
  std::discrete_distribution<> resampleDist(weights.begin(), weights.end());
  for(int i=0; i<num_particles; i++){
    int index = resampleDist(generator);
    newParticles.push_back(particles[index]);
  }
  // assign the new particles to the particles instance variable
  particles = newParticles;

  //cout << "Resampled" << endl;
  

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



// adopted from Quiz earlier to calculate multivariable gaussian distribution
double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
  //cout << "Multivariable_Gaussianed" << endl;
  return weight;
}
