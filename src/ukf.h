#ifndef UKF_H
#define UKF_H

#include <vector>
#include <string>
#include <fstream>

#include "Eigen/Dense"
#include "measurement_package.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param measurement_pack The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage measurement_pack);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * AugmentedSigmaPoints
   * @param[out] Xsig_out The sigma points as columns in a matrix
   */
  void AugmentedSigmaPoints(MatrixXd* Xsig_out);

  /**
   * SigmaPointPrediction Predicts sigma points from augmented sigma points 
   *                      using the non-linear functions of the CVTR model.
   * @param[in] delta_t Time between k and k+1 in s
   * @param[in] Xsig_in Augmented sigma points generated with the  last 
   *                    a posteriori estimation.
   */
  void SigmaPointPrediction(const double delta_t, const MatrixXd &Xsig_in);

  /**
   * PredictMeanAndCovariance Predicts mean and covariance of the state
   */
  void PredictMeanAndCovariance();

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param measurement_pack The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage measurement_pack);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param measurement_pack The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage measurement_pack);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param[out] z_out The measurement prediction at k+1
   * @param[out] S_out The measurement covariance matrix prediction at k+1
   * @param[out] Zsig The measurement sigma points as colums in a matrix
   */
  void PredictRadarMeasurement(VectorXd& z_out, MatrixXd& S_out, MatrixXd& Zsig);

  /**
   * GenerateSigmaPoints
   * @param[out] Xsig_out The sigma points as columns in a matrix
   */
  void GenerateSigmaPoints(MatrixXd* Xsig_out);
};

#endif /* UKF_H */
