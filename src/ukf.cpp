#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  


  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  // Initialize Process covariance matrix
  P_ << 0.15, 0, 0, 0, 0,
        0, 0.15, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  previous_timestamp_ = 0.0;
   //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    // initial weights vector
  weights_ = VectorXd(2 * n_aug_ + 1);

  // Initialize Lidar covariance matrix
  R_lidar_ = MatrixXd::Zero(2, 2);
  R_lidar_(0,0) = std_laspx_ * std_laspx_;
  R_lidar_(1,1) = std_laspy_ * std_laspy_;

  // Initialize Radar covariance matrix 
  R_radar_ = MatrixXd::Zero(3, 3);
  R_radar_(0,0) = std_radr_ * std_radr_;
  R_radar_(1,1) = std_radphi_ * std_radphi_;
  R_radar_(2,2) = std_radrd_ * std_radrd_;
}

UKF::~UKF() {}

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
    /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_) ||
     (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_))
       return;
  if (!is_initialized_){
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "UKF: " << endl;

    x_ << 1, 1, 1, 1, 0.1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        double rho = meas_package.raw_measurements_[0];
        double phi = meas_package.raw_measurements_[1];
        x_[0] = rho*cos(phi);
        x_[1] = rho*sin(phi);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      /**
      Initialize state.
      */
      x_[0] = meas_package.raw_measurements_[0];
      x_[1] = meas_package.raw_measurements_[1];
    }
    cout << x_ << endl;
    previous_timestamp_ = meas_package.timestamp_;
    weights_.fill(0.0);
    cout << weights_ << endl;
    for(unsigned int i=0; i<(2*n_aug_+1); i++){
      double wi;
      if(i==0){
          wi = lambda_/(lambda_ + n_aug_);
      }else{
          wi = 1/(2*(lambda_ + n_aug_));
      }
      weights_[i] = wi;
    }
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;
  Prediction(dt);
  cout << "The predicted state = " << endl;
  // print the output
  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl << P_ << endl;
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    // Radar updates
    cout << "Calling UpdateRadar with observation = " << endl << meas_package.raw_measurements_ << endl;
    UpdateRadar(meas_package);
    cout << "Lidar NIS = " << NIS_lidar_ << endl;
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    cout << "Calling UpdateLidar with observation = " << endl << meas_package.raw_measurements_ << endl;
    UpdateLidar(meas_package);
    cout << "Radar NIS = " << NIS_radar_ << endl;
  }
  cout << "Updating the state to:" << endl;
  // print the output
  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl << P_ << endl;
}

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
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
 
  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  { 
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
 
 
 cout << Xsig_aug << endl;
 for(unsigned int i=0; i < (2 * n_aug_ + 1); i++){
      double px = Xsig_aug(0, i);
      double py = Xsig_aug(1, i);
      double v = Xsig_aug(2, i);
      double yaw = Xsig_aug(3, i);
      double yaw_d = Xsig_aug(4, i);
      double v_d = Xsig_aug(5, i);
      double yaw_dd = Xsig_aug(6, i);
      
      VectorXd m = VectorXd(n_x_);
      
      if(yaw_d >= 0.001){
          m(0) = (v/yaw_d)*(sin(yaw+yaw_d*delta_t) - sin(yaw));
          m(1) = (v/yaw_d)*(-cos(yaw+yaw_d*delta_t) + cos(yaw));
      }else{
          m(0) = v*cos(yaw)*delta_t;
          m(1) = v*sin(yaw)*delta_t;
      }
      m(2) = 0.0;
      m(3) = yaw_d*delta_t;
      m(4) = 0.0;
      
      VectorXd nu = VectorXd(n_x_);
      nu(0) = 0.5*delta_t*delta_t*cos(yaw)*v_d;
      nu(1) = 0.5*delta_t*delta_t*sin(yaw)*v_d;
      nu(2) = v_d*delta_t;
      nu(3) = 0.5*delta_t*delta_t*yaw_dd;
      nu(4) = delta_t*yaw_dd;
      VectorXd xk = VectorXd(n_x_);
      xk << px, py, v, yaw, yaw_d;
      Xsig_pred_.col(i) = (xk + m + nu);
  }
   x_.fill(0.0);
   P_.fill(0.0);
  for(unsigned int i=0; i<(2*n_aug_+1); i++){
      x_ += (weights_[i] * Xsig_pred_.col(i));
  }
  //set weights
  //predict state mean
  //predict state covariance matrix
  for(unsigned int i=0; i<(2*n_aug_+1); i++){
      VectorXd Xdiff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (Xdiff(3)> M_PI) Xdiff(3)-=2.*M_PI;
    while (Xdiff(3)<-M_PI) Xdiff(3)+=2.*M_PI;
    P_ += (weights_[i] * Xdiff * Xdiff.transpose());
  }

}

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
  int n_z = 2;

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for(unsigned int i=0;i<(2*n_aug_+1);i++){
      VectorXd xdiff = Xsig_pred_.col(i) - x_;
      VectorXd zdiff = xdiff.head(2);
      S += weights_(i) * zdiff * zdiff.transpose();
      Tc = Tc + weights_(i) * xdiff * zdiff.transpose();
  }

  S += R_lidar_;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // Actual Radar measurement
  VectorXd z = meas_package.raw_measurements_;
  //residual
  VectorXd z_resid = z - x_.head(2);


  //calculate Lidar NIS
  NIS_lidar_ = z_resid.transpose() * S.inverse() * z_resid;

  //update state mean and covariance matrix
  x_ = x_ + K * z_resid;
  P_ = P_ - K*S*K.transpose();
}
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
  int n_z = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);


   for(unsigned int i=0;i < (2*n_aug_+1);i++){
      double px = Xsig_pred_(0, i);
      double py = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double yaw = Xsig_pred_(3,i);
      VectorXd zk_i = VectorXd(n_z);
      double rho = sqrt(px*px + py*py);
      Zsig(0,i) = rho;
      Zsig(1,i) = atan2(py, px);
      Zsig(2,i) = ((fabs(rho) > 0.0001) ? (px*cos(yaw)*v + py*sin(yaw)*v)/rho : 0.0 );
      z_pred += weights_(i) * Zsig.col(i);
  }
    //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for(unsigned int i=0;i<(2*n_aug_+1);i++){
      VectorXd zdiff = Zsig.col(i) - z_pred;

      //angle normalization
      while (zdiff(1)> M_PI) zdiff(1)-=2.*M_PI;
      while (zdiff(1)<-M_PI) zdiff(1)+=2.*M_PI;

      S += weights_(i) * zdiff * zdiff.transpose();
  }
  S += R_radar_;


  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);


  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd zdiff = Zsig.col(i) - z_pred;
    //angle normalization
    while (zdiff(1)> M_PI) zdiff(1)-=2.*M_PI;
    while (zdiff(1)<-M_PI) zdiff(1)+=2.*M_PI;

    // state difference
    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (xdiff(3)> M_PI) xdiff(3)-=2.*M_PI;
    while (xdiff(3)<-M_PI) xdiff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * xdiff * zdiff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // Actual Radar measurement
  VectorXd z = meas_package.raw_measurements_;
  //residual
  VectorXd z_resid = z - z_pred;

  //angle normalization
  while (z_resid(1)> M_PI) z_resid(1)-=2.*M_PI;
  while (z_resid(1)<-M_PI) z_resid(1)+=2.*M_PI;

  //calculate Radar NIS
  NIS_radar_ = z_resid.transpose() * S.inverse() * z_resid;

  //update state mean and covariance matrix
  x_ = x_ + K * z_resid;
  P_ = P_ - K*S*K.transpose();
}
