#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <iostream>
#include <math.h>
const float pi = 3.14159265;

using namespace std;



KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

	MatrixXd I = MatrixXd::Identity(4, 4);

	VectorXd y = z - H_ * x_;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K_ =  P_ * Ht * Si;

	//new state
	x_ = x_ + (K_ * y);
	P_ = (I - K_ * H_) * P_;


}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	MatrixXd I = MatrixXd::Identity(4, 4);

	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);

	// prevent a divide by zero error
	float rho = sqrt(px*px + py*py);
	if (rho < 0.00001){
		rho = 0.00001;
	}
	float theta = atan2(py, px);
	float ro_dot = (px*vx + py*vy)/rho;

	VectorXd z_pred = VectorXd(3);
	z_pred << rho, theta, ro_dot;

	VectorXd zm = z;
	// make sure the radius, rho, is positive
	if (zm(0) < 0){
		zm(0) = -1 * zm(0);
	}

	VectorXd y =  z - z_pred;
	// make sure theta is between -pi and pi
	while (y(1) < -pi){
		y(1) += 2 * pi;
	}

	while (y(1) > pi){
		y(1) -= 2 * pi;
	}

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K_ =  P_ * Ht * Si;
	//new state
	x_ = x_ + (K_ * y);
	P_ = (I - K_ * H_) * P_;

}
