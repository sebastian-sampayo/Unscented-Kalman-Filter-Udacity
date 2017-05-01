#include "ukf.h"

#include <iostream>

#include "Eigen/Dense"
#include "tools.h"

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
  x_ = VectorXd::Zero(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);
  
  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2*n_aug_ + 1);

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
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

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
    (*Xsig_out).col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
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
 * @param {MeasurementPackage} measurement_pack The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
#ifdef DEBUG
    cout << "---------- New measurement ---------" << endl;
#endif

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  //Initialize the state x_ with the first measurement.
  //TODO: Try different P matrix initializations. However, the final result 
  // (after several iterations) shouldn't be very dependent on this.
  if (!is_initialized_) {
    // first measurement
#ifdef DEBUG
    cout << "x State init" << endl;
#endif

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      Tools tools;
      VectorXd x_aux = tools.ConvertRadar2Cartesian(measurement_pack.raw_measurements_);
      x_(0) = x_aux(0); // px
      x_(1) = x_aux(1); // py
      x_(2) = sqrt(x_aux(2)*x_aux(2) + x_aux(3)*x_aux(3)); // v = sqrt(vx^2 + vy^2)
      x_(3) = 0; // yaw
      x_(4) = 0; // yaw dot
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x_ << measurement_pack.raw_measurements_[0], 
            measurement_pack.raw_measurements_[1],
            0, 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    time_us_ = measurement_pack.timestamp_;  
    return;
  }
  
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  //Compute the time elapsed between the current and previous measurements
  //Time is measured in seconds and timestamps are in microseconds.
  const float dt = (measurement_pack.timestamp_ - time_us_) / 1000000.0;
  assert(dt > 0); // if dt <= 0, the input data is corrupted
  time_us_ = measurement_pack.timestamp_;
  Prediction(dt);

#ifdef DEBUG
  // print the output
  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl << P_ << endl;
#endif
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
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Generate augmented sigma points, predict new sigma points using the model
  // and predict mean con covariance of the state.
  AugmentedSigmaPoints(&Xsig_aug);
  SigmaPointPrediction(delta_t, Xsig_aug, &Xsig_pred_);
  PredictMeanAndCovariance(Xsig_pred_);
}

// ----------------------------------------------------------------------------
/**
 * PredictMeanAndCovariance Predicts mean and covariance of the state
 * @param[in] Xsig_in Predicted sigma points of size [n_x_, (2*n_aug_+1)]
 */
void UKF::PredictMeanAndCovariance(const MatrixXd &Xsig_in) {
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_in.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_in.col(i) - x_;
    //angle normalization
    while (x_diff(3) >  M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

// ----------------------------------------------------------------------------
/**
 * SigmaPointPrediction Predicts sigma points from augmented sigma points
 * @param delta_t Time between k and k+1 in s
 * @param[in] Xsig_in Augmented sigma points generated with the  last 
 *                    a posteriori estimation.
 * @param[out] Xsig_out Predicted sigma points using the non-linear functions 
 *                      of the CVTR model.
 */
void UKF::SigmaPointPrediction(const double delta_t
                             , const MatrixXd &Xsig_in
                             , MatrixXd* Xsig_out) {
  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_in(0,i);
    double p_y = Xsig_in(1,i);
    double v = Xsig_in(2,i);
    double yaw = Xsig_in(3,i);
    double yawd = Xsig_in(4,i);
    double nu_a = Xsig_in(5,i);
    double nu_yawdd = Xsig_in(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    (*Xsig_out)(0,i) = px_p;
    (*Xsig_out)(1,i) = py_p;
    (*Xsig_out)(2,i) = v_p;
    (*Xsig_out)(3,i) = yaw_p;
    (*Xsig_out)(4,i) = yawd_p;
  }
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
