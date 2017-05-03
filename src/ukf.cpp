#include "ukf.h"

#include <iostream>

#include "Eigen/Dense"
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

// ----------------------------------------------------------------------------
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
  std_a_ = 0.5; // Analysing the data 2.5 m/s seems to be the max acceleration. However, lower values than 1.25 work better regarding the NIS values.

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.25; // This value is more difficult to get but 0.25 seems to work, as the yaw rate doesn't seem to be larger than 0.5.

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

// ----------------------------------------------------------------------------
UKF::~UKF() {}

// ----------------------------------------------------------------------------
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
  
  /*****************************************************************************
  *  Update
  ****************************************************************************/
  //Update measurement matrix and measurement covariance matrix in each case, 
  //then update the state, depending on the sensor type.
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    // Radar updates
    UpdateRadar(measurement_pack);
  } else if (use_laser_) {
    // Laser updates
    assert(measurement_pack.sensor_type_ == MeasurementPackage::LASER);
    UpdateLidar(measurement_pack);
  }

#ifdef DEBUG
  // print the output
  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl << P_ << endl;
#endif
}

// ----------------------------------------------------------------------------
void UKF::Prediction(double delta_t) {
  /**
  TODO: DONE

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Generate augmented sigma points, predict new sigma points using the model
  // and predict mean con covariance of the state.
  AugmentedSigmaPoints(&Xsig_aug);
  SigmaPointPrediction(delta_t, Xsig_aug);
  PredictMeanAndCovariance();
}

// ----------------------------------------------------------------------------
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
void UKF::SigmaPointPrediction(const double delta_t, const MatrixXd &Xsig_in) {
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
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

// ----------------------------------------------------------------------------
void UKF::PredictMeanAndCovariance() {
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) >  M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

// ----------------------------------------------------------------------------
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
  /**
  TODO: DONE?

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  // Here I use the standard Kalman Filter update step, as the Lidar equations
  // are linear.
  MatrixXd H = MatrixXd::Zero(2, n_x_);
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;
  const VectorXd z_pred = H * x_;
  const VectorXd z = measurement_pack.raw_measurements_;
  MatrixXd R = MatrixXd(2,2);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  //Calculate the kalman gain matrix
  VectorXd y = z - z_pred;
  const MatrixXd Ht = H.transpose();
  const MatrixXd S = H * P_ * Ht + R;
  const MatrixXd Si = S.inverse();
  const MatrixXd PHt = P_ * Ht;
  const MatrixXd K = PHt * Si;

  //angle normalization
  while(y(1) >  M_PI) y(1) -= 2.*M_PI;
  while(y(1) < -M_PI) y(1) += 2.*M_PI;

  //new estimates
  x_ = x_ + (K * y);
  const long x_size = x_.size();
  const MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;
  
  // Calculate NIS: Normalized Innovation Squared for consistency check
  NIS_laser_ = y.transpose() * Si * y;
}

// ----------------------------------------------------------------------------
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
  /**
  TODO: DONE

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // Number of components of the measurement vector: rho, phi, rhodot
  const int n_z = 3;
  //create a vector for the measurement prediction, covariance prediction and 
  // sigma points in the measurement space
  VectorXd z_pred = VectorXd::Zero(n_z);
  MatrixXd S =  MatrixXd::Zero(n_z, n_z);
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  PredictRadarMeasurement(z_pred, S, Zsig);
  
  // radar measurement
  VectorXd z = measurement_pack.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1) >  M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) >  M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1) >  M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
  // Calculate NIS: Normalized Innovation Squared for consistency check
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

// ----------------------------------------------------------------------------
void UKF::PredictRadarMeasurement(VectorXd &z_pred, MatrixXd &S, MatrixXd &Zsig) {
  //number of measurement components
  const int n_z = z_pred.size();
  
    //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  // VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  // MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1) >  M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
}

// ----------------------------------------------------------------------------
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
