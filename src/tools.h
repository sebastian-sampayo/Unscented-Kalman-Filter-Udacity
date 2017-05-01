/****************************************************************************\
 * Udacity Nanodegree: Self-Driving Car Engineering - December cohort
 * Project 6: Extended Kalman Filter
 * Date: 16th April 2017
 * 
 * Author: Sebastián Lucas Sampayo
 * e-mail: sebisampayo@gmail.com
 * file: tools.h
 * Description: Some tools for calculation of Jacobian matrix, RMSE and 
 * radar to model state conversion.
\****************************************************************************/

#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>

#include "Eigen/Dense"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate Jacobians.
  * @param x_state Model state vector in the format: {p_x, p_y, v_x, v_y}
  * @return Jacobian matrix
  */
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

  /**
  * A helper method to calculate RMSE.
  * @param estimations A vector of estimations
  * @param ground_truth A vector of the ground truth values for each estimation
  * @return The Root Mean Squared Error = \sqrt{1/N \sum{(x_i^e - x_i^t)^2}}
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /**
  * A helper method to convert from radar measurement to model state x 
  * (Polar to Cartesian coordinates).
  * @param radar_measurement Polar coordinates in the format: {ro, phi, ro_dot}
  * @return Cartesian coordinates in the format: {p_x, p_y, v_x, v_y}
  */
  Eigen::VectorXd ConvertRadar2Cartesian(const Eigen::VectorXd& radar_measurement);
  
private:
  //A parameter to tune the initial velocity when we have a radar measurement
  float init_velocity_scale = 0.2;
};

#endif /* TOOLS_H_ */
