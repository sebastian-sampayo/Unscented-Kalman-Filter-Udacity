#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  
  // State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  
  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd::Zero(n_aug_, 2*n_aug_ + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 4;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  // time when the state is true, in us
  time_us_ = 0;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // the current NIS for radar
  NIS_radar_ = 0;

  // the current NIS for laser
  NIS_laser_ = 0;
  
  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_ + 1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  double weight_i = 0.5/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    weights_(i) = weight_i;
  }
}

UKF::~UKF() {}

// ----------------------------------------------------------------------------
/**
 * AugmentedSigmaPoints
 * @param[out] Xsig_out The sigma points as columns in a matrix
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  
  //create augmented mean state
  x_aug.head(n_x_) = x_;
  for(int i = 0; i < n_aug_ - n_x_; ++i)
    x_aug(n_x_+i) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  (*Xsig_out).col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    (*Xsig_out).col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    (*Xsig_out).col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
}

// ----------------------------------------------------------------------------
/**
 * GenerateSigmaPoints
 * @param[out] Xsig_out The sigma points as columns in a matrix
 */
void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  
  //set first column of sigma point matrix
  (*Xsig_out).col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    (*Xsig_out).col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
    (*Xsig_out).col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
  }
  
}

// ----------------------------------------------------------------------------
/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
}

// ----------------------------------------------------------------------------
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

// ----------------------------------------------------------------------------
/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

// ----------------------------------------------------------------------------
/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
