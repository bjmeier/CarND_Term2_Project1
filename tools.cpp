#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size()==0){
	    cout << "no estimations" << endl;
		    return rmse;
  }
  if( estimations.size() != ground_truth.size()){
	  cout << "inconsistant estimation and ground truth points" << endl;
      return rmse;
  }
  for(int i=0; i < estimations.size(); ++i){
    if (estimations[i].rows() != ground_truth[i].rows()){
	    cout << "inconsistant estimation and ground truth rows" << endl;
	    return rmse;
	}
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
      // ... your code here
	  VectorXd res = (estimations[i]-ground_truth[i]).array().square();
	  rmse = rmse + res;
  }

  rmse = rmse.array()/estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //check division by zer0
  float px2py2, H00, H01, H10, H11, H20, H21, H22, H23;
  px2py2 = px*px + py*py;
  if (px2py2 < 0.0001){
      px2py2 = 0.0001;
  }
  H00 = px * pow(px2py2,-0.5);
  H01 = py * pow(px2py2,-0.5);
  H10 = -py / px2py2;
  H11 =  px / px2py2;
  H20 = py * (vx*py - vy*px) * pow(px2py2, -1.5);
  H21 = px * (vy*px - vx*py) * pow(px2py2, -1.5);
  H22 = H00;
  H23 = H01;
  Hj <<    H00, H01, 0, 0,
           H10, H11, 0, 0,
           H20, H21, H22,H23;

  return Hj;

}
