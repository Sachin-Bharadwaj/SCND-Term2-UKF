#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define FAILURE -1

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 1;//30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;//30;

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
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;
  P_ << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0,
       0, 0, 1, 0, 0,
       0, 0, 0, 1, 0,
       0, 0, 0, 0, 1;
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0.5/(lambda_+n_aug_));
  weights_(0) = lambda_/(lambda_+n_aug_);
  H_laser_ = MatrixXd(2, n_x_);
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;
  R_laser_ = MatrixXd(2,2);
  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0, std_radrd_*std_radrd_;
  // open outfiles for dumping NIS data
  fh_NIS_r.open(out_fname_NIS_radar);
  fh_NIS_l.open(out_fname_NIS_laser);
  // check if file can be opened
  if(!fh_NIS_r.is_open())
  {
    cerr << "Cannot open output file: " << out_fname_NIS_radar << endl;
    exit(FAILURE);
  }
  if(!fh_NIS_l.is_open())
  {
    cerr << "Cannot open output file: " << out_fname_NIS_laser << endl;
    exit(FAILURE);
  }
}

UKF::~UKF()
{
  // close the file handlers
  if(fh_NIS_r.is_open())
  {
    fh_NIS_r.close();
  }
  if(fh_NIS_l.is_open())
  {
    fh_NIS_l.close();
  }
}

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
  if(!is_initialized_){
    cout << "Initializing with First measurements" << endl;
    x_ << 1, 1, 0, 0, 0;
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      cout<< "Initializing with Radar measurement" << endl;
      x_(0) = meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]);
      x_(1) = meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]);       	
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      cout << " Initializing with LASER measurement" << endl;
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1]; 
    }
  is_initialized_ = true;
  time_us_ = meas_package.timestamp_;
  cout << "UKF initialized with state vec :"<<endl;
  cout << x_ << endl;
  return;
  } // end of is_initialized
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; // in sec
  time_us_ = meas_package.timestamp_;
  // call the prediction step
  Prediction(delta_t);
  //cout <<"Prediction clean..." << endl;
  // call measurement update
  if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    UpdateLidar(meas_package);
    //cout << "LASER Update clean.." << endl;
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
    //cout << "Radar Update clean.." << endl;
  }   
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
  // Step1: Generate Sigma Points
  VectorXd x_aug_ = VectorXd(n_aug_); // augmented state vec
  MatrixXd x_sig_aug_ = MatrixXd(n_aug_, 2*n_aug_ + 1); // sigma point matrix
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_); // augmented state cov matrix
  // create augmented state vector
  x_aug_.head(n_x_) = x_;
  x_aug_(5) = 0.0; // zero mean noise
  x_aug_(6) = 0.0; // zero mean noise
  // create augmented state cov Mat
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(5, 5) = std_a_ * std_a_;
  P_aug_(6, 6) = std_yawdd_ * std_yawdd_;
  // compute square root using cholesky
  MatrixXd P_aug_sqrt_ = P_aug_.llt().matrixL();
  // compute sigma point
  x_sig_aug_.col(0) = x_aug_; // mean vector   
  for(unsigned int i=0; i<n_aug_; i++)
  {
    x_sig_aug_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_)*P_aug_sqrt_.col(i);
    x_sig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_)*P_aug_sqrt_.col(i); 
    
  }
  // cout << "Inside Prediction... " << endl;
  // End of sigma point generation

  // Step2: Predict Sigma Point by running it through process model
  for(unsigned int i=0; i< 2*n_aug_+1; i++)
  {
    // extract values for better readability
    double p_x = x_sig_aug_(0,i);
    double p_y = x_sig_aug_(1,i);
    double v   = x_sig_aug_(2,i);
    double yaw = x_sig_aug_(3,i);
    double yawd = x_sig_aug_(4,i);
    double nu_a = x_sig_aug_(5,i);
    double nu_yawdd = x_sig_aug_(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if(fabs(yawd) > 0.001)
    {
      px_p = p_x + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
    }
    else
    {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    } 
    
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise contribution
    px_p = px_p + 0.5*nu_a*delta_t*delta_t*cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t*sin(yaw);
    v_p  = v_p + nu_a*delta_t;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma points
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  // End of Sigma Point Prediction
  //cout << "Inside Predction.. end of step2" << endl;
  // Step3: Predict Mean and Cov from Precited Sigma Points
  x_ = Xsig_pred_ * weights_;
  //cout << "Here.." << endl;
  //cout << x_ << endl;
  P_.fill(0.0); // Imp: Reset the P_ otherwise it will get accumulated inside the loop
  for(unsigned int i=0; i<2*n_aug_; i++)
  {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    if(x_diff(3) >= M_PI) x_diff(3) -= 2.*M_PI;
    if(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
    // set the state cov Mat (sum of outer products)
    P_ = P_ + weights_(i)* x_diff * x_diff.transpose();
  }

  //cout << "Inside Prediction.. end of step3" << endl; 	
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
  // State vector and LIDAR measurements have linear relationship
  // compute the residual
  VectorXd res = meas_package.raw_measurements_ - H_laser_ * x_;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si =  S.inverse();
  MatrixXd K = P_ * Ht * Si;
  // update state and cov Mat
  x_ = x_ + K* res;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K*H_laser_) * P_;

  // update NIS for laser
  NIS_laser_ = res.transpose()*Si*res;
  fh_NIS_l << NIS_laser_ << endl;
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
  // Radar Measurements are non-linear function of state variable, so need sigma points, take short cut, already computed sigma points while using update step
  // transform sigma points into measurement space
  int n_z = 3;
  MatrixXd Zsig_ = MatrixXd(n_z, 2*n_aug_ + 1);
  VectorXd z_pred_ = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z, n_z); 
  MatrixXd Tc_ = MatrixXd(n_x_, n_z);
  for(unsigned int i=0; i<2*n_aug_+1; i++)
  {
    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = v*cos(yaw);
    double v2 = v*sin(yaw);

    // measurement model
    Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);
    Zsig_(1,i) = atan2(p_y, p_x);
    Zsig_(2,i) = (p_x*v1 + p_y*v2) / Zsig_(0,i);
   
  }
  // calculate mean predicted measurement
  z_pred_ = Zsig_ * weights_;
  //cout << "z_pred_" << endl;
  //cout << z_pred_<< endl;
  // calculate measurement cov mat S
  S.fill(0.0);
  for(unsigned int i=0; i<2*n_aug_+1; i++)
  {
    // compute diff
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    // angle normalization
    if(z_diff(1) >= M_PI) z_diff(1) -=2.*M_PI;
    if(z_diff(1) < -M_PI) z_diff(1) +=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  // add the measurement noise uncertainty
  S = S + R_radar_;
  //cout << "Final S..." << endl;
  //cout << S << endl;
  // compute cross correlation matrix Tc
  Tc_.fill(0.0);
  for(unsigned int i=0; i<2*n_aug_+1; i++)
  {
    // centered predicted measurement
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    // angle normalization
    if(z_diff(1) >= M_PI) z_diff(1) -= 2.*M_PI;
    if(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    if(x_diff(3) >= M_PI) x_diff(3) -= 2.*M_PI;
    if(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc_ = Tc_ + weights_(i) * x_diff * z_diff.transpose();   
  }  
  // calculate Kalman Gain
  MatrixXd Si = S.inverse();
  MatrixXd K = Tc_ * Si;
  // update state mean and cov Mat
  VectorXd res = meas_package.raw_measurements_ - z_pred_;
  // angle normalization
  if(res(1) >= M_PI) res(1) -= 2.*M_PI;
  if(res(1) < -M_PI) res(1) += 2.*M_PI;
  // state meas and cov update
  x_ = x_ + K * res;
  P_ = P_ - K * S * K.transpose(); 
  // update NIS for radar
  NIS_radar_ = res.transpose() * Si * res;
  fh_NIS_r << NIS_radar_ << endl; 
}
