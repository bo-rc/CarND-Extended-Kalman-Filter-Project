#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  VectorXd x_in(4);
  x_in << 0, 0, 0, 0;

  MatrixXd P_in(4,4);
  P_in << 1, 0, 1, 0,
	  0, 1, 0, 1,
	  0, 0, 1000, 0,
	  0, 0, 0, 1000;

  MatrixXd F_in(4,4);
  F_in << 1, 0, 1, 0,
	  0, 1, 0, 1,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  MatrixXd H_in(2, 4);
  H_in << 1, 0, 0, 0,
	  0, 1, 0, 0;

  MatrixXd Hj_in(3,4);
  Hj_in << 1, 0, 0, 0,
	   0, 1, 0, 0,
	   0, 0, 1, 0;

  MatrixXd Q_in(4,4);
  Q_in << 1, 0, 1, 0,
	  0, 1, 0, 1,
	  1, 0, 1, 0,
	  0, 1, 0, 1;

  ekf_.Init(x_in, P_in, F_in, H_in, Hj_in, R_laser_, R_radar_, Q_in);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

// convert phi to [-PI, PI]
double FusionEKF::_convert_phi(const double & phi) {

    if (phi >= -PI && phi <= PI) {
	return phi;
    }

    if (phi > PI) {
	return _convert_phi(phi - 2.*PI);
    }

    if (phi < -PI) {
	return _convert_phi(phi + 2.*PI);
    }
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	auto rho = measurement_pack.raw_measurements_[0];
	auto phi = measurement_pack.raw_measurements_[1];
	auto rho_dot = measurement_pack.raw_measurements_[2];

	// converting phi to [-pi, pi)
	phi = _convert_phi(phi);

	auto px = rho * sin(phi);
	auto py = rho * cos(phi);
	auto vx = rho_dot * sin(phi);
	auto vy = rho_dot * cos(phi);

	ekf_.x_ << px, py, vx, vy;

	previous_timestamp_ = measurement_pack.timestamp_;
	is_initialized_ = true;
	return;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	//set the state with the initial location and zero velocity
	ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

	previous_timestamp_ = measurement_pack.timestamp_;
	is_initialized_ = true;
	return;
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
  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  float noise_ax = 9;
  float noise_ay = 9;

  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
	      0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
	      dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
	      0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      ekf_.Hj_ = tools.CalculateJacobian(ekf_.x_);
      VectorXd z_norm(3);
      z_norm << measurement_pack.raw_measurements_[0], _convert_phi(measurement_pack.raw_measurements_[1]), measurement_pack.raw_measurements_[2];
      ekf_.UpdateEKF(z_norm);
  } else {
    // Laser updates
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
