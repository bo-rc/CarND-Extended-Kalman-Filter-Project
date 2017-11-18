#include "kalman_filter.h"
#include <iostream>
using std::cout;
using std::endl;

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
			MatrixXd &H_in, MatrixXd &Hj_in, MatrixXd &R_laser_in,
			MatrixXd &R_radar_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  Hj_ = Hj_in;
  R_laser_ = R_laser_in;
  R_radar_ = R_radar_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
    /**
     * predict the state
     */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    /**
     * update the state by using Kalman Filter equations
     */
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_laser_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
   /**
    * update the state by using Extended Kalman Filter equations
    */
    double px = x_[0];
    double py = x_[1];
    double vx = x_[2];
    double vy = x_[3];

    double rho = sqrt(px*px + py*py);

    double phi = 0;
    if (fabs(rho) > 0.00001) {
	phi = asin(py/rho);

	if (px >= 0) {
	    if (py >= 0) {
		phi = phi;
	    } else {
		phi = phi;
	    }
	} else {
	    if (py >= 0) {
		phi = PI - phi;
	    } else {
		phi = - (PI + phi);
	    }
	}
    }

    double rho_dot = 0;
    if (fabs(rho) > 0.00001) {
	rho_dot = (px*vx + py*vy) / rho;
    }

    VectorXd z_pred(3);
    z_pred << rho, phi, rho_dot;

    VectorXd y = z - z_pred;

    // normalize phi difference to [-PI, PI]
    while (y[1] > PI) y[1] -= 2.*PI;
    while (y[1] < -PI) y[1] += 2.*PI;

    // same as above, which is more elegant.
//    // this handles discontinuity of phi at -PI/PI transition.
//    if (fabs(y[1]) > PI) {
//	if (z[1] > 0) {
//	    y[1] = -(2*PI - fabs(y[1])); // clock-wise difference
//	} else {
//	    y[1] = (2*PI - fabs(y[1])); // counter clock-wise difference
//	}
//    }

    MatrixXd Hjt = Hj_.transpose();
    MatrixXd S = Hj_ * P_ * Hjt + R_radar_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Hjt;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * Hj_) * P_;
}
