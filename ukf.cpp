#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


//DISCLAIMER: code heavily relies on codes shown in lectures and questions asked and answered on the Udacity forums. 

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
  x_ << 0.1, 0.1, 0.1, 0.1, 0.1;

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.1; // Most cars reach 60 mph in roughly 10 sec, this corresponds roughtly to 2.5 m/s*s. However, most of the time, drives tend to accelarate much more slowly.

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.1; // This means roughly expecting that acceleration is most of the time less than 60 degrees in a second. That sounds conservative enough.
  // (1, 0.4), (2,2), (0.5, 0.5), (5, 1.5), (5, 2)=(11, 33), (3,3)=(11, 32), (5,5)=(12, 39)
  // (3, 0.1) = (16, 45), (2,1)=(11.67, 29.12), (2,0.5)=(12, 29), (2,0.1)=(16, 44)
  //  (2,4)=(11.62, 32.55), (2, 0.7)
  // (2, 1) = (11, 11, 29, 22), (1,2)=(11,11,29,24)
  // (1,1)=(11,11, 28, 20), 
  // (1,0.5) = (0.1170, 0.1107, 0.2813, 0.2036)
  // (0.5, 0.5) = (0/1158, 0.1117, 0.2784, 0.1915)
  // (0.5, 0.3) = (0.1184, 0.1118, 0.2891, 0.2037)
  // (0.3, 0.5) = (0.1156, 0.1120, 0.2809, 0.1899)
  // (0.1, 0.5) = (0.1156, 0.1120, 0.2809, 0.1899)
  // (0.1, 0.1) = ()
  // (0.2, 0.2) = ()
  // (1.5, 1) = (11,11,28, 21)
  // Laser m easurement noise standard deviation position1 in m
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

  P_ << 1.0, 0, 0, 0, 0,
	  0, 1, 0, 0, 0,
	  0, 0, 1, 0, 0,
	  0, 0, 0,1, 0,
	  0, 0, 0, 0, 1;

  /*P_ << 2, 0, 0, 0, 0,
	  0, 4, 0, 0, 0,
	  0, 0, 1, 0, 0,
	  0, 0, 0, 0.5, 0,
	  0, 0, 0, 0, 0.5;*/




  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  // dim measurement vector of radar
  n_z = 3;

  // dim measurement vector of lidar
  n_l_ = 2;


  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;


  for (int i = 1; i< 2 * n_aug_ + 1; i++) { 
	  double weight = 0.5 / (n_aug_ + lambda_);
	  weights_(i) = weight;
  }

  ///* initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;


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
	if (!is_initialized_) {
		//long long previous_timestamp_;
		previous_timestamp_ = meas_package.timestamp_;
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

			float rho = meas_package.raw_measurements_(0); // this is the distance to the pedestrian
			float phi = meas_package.raw_measurements_(1); // angle between the direction of vehicle motion and the object we tracking
			float rhodot = meas_package.raw_measurements_(2); // change rate of rho
			float px = rho * cos(phi);
			float py = rho * sin(phi);
			//      float vx = rhodot * cos(phi);
			//      float vy = rhodot * sin(phi);
			//      float v = sqrt(vx*vx + vy*vy);
			//if (px < 0.001 && px>-0.001) px = 0.001;
			//if (py < 0.001 && py>-0.001) py = 0.001;
			float v = 0;
			//      float yaw = cos(vx/vy);
			float yaw = 0;
			float yawdot = 0;
			x_ << px, py, v, yaw, yawdot;

		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			float px = meas_package.raw_measurements_(0);
			float py = meas_package.raw_measurements_(1);
			//if (px < 0.001 && px>-0.001) px = 0.001;
			//if (py < 0.001 && py>-0.001) py = 0.001;
			x_ << px, py, 0, 0, 0;
		}
		is_initialized_ = true;
		//cout << "initialization OK" << endl;
		return;
	}

	float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
	//cout << dt << endl;
	previous_timestamp_ = meas_package.timestamp_;
	//cout << "Run prediction " << endl;
	Prediction(dt);
	//cout << "Prediciton done " << endl;
	//cout << "Xsig_pred_ in Proces measurement" << Xsig_pred_ << endl;
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		//cout << "Update radar" << endl;
		UpdateRadar(meas_package);

	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		//cout << "Update lidar" << endl;
		UpdateLidar(meas_package);
	}

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	//cout << " PREDICTION" << endl;
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */


//Sigma points 

//create augmented mean vector

	//// AUGMENTATION ASSIGNMENT
	VectorXd x_aug = VectorXd(7);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);



	//create sigma point matrix
	//MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

	//calculate square root of P
	//MatrixXd A = P.llt().matrixL();

	//create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_*std_a_;
	P_aug(6, 6) = std_yawdd_*std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
	}



	//// SIGMA POINT PREDICTION ASSIGMENT
	//create matrix with predicted sigma points as columns
	//MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	//double delta_t = 0.1; //time diff in sec


						  //predict sigma points
	for (int i = 0; i< 2 * n_aug_ + 1; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			//px_p = p_x + v / 0.001 * (sin(yaw + 0.001*delta_t) - sin(yaw));
			py_p = p_y + v*delta_t*sin(yaw);
			//py_p = p_y + v / 0.001 * (cos(yaw) - cos(yaw + 0.001*delta_t));
		}

		//if (px < 0.001 && px>-0.001) px = 0.001;
		//if (py < 0.001 && py>-0.001) py = 0.001;

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
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}

 //// Predicted Mean And Covariance Assigments
 //create vector for predicted state
	VectorXd x = VectorXd(n_x_);

	//create covariance matrix for prediction
	MatrixXd P = MatrixXd(n_x_, n_x_);

	// set weights
	double weight_0 = lambda_ / (lambda_ + n_aug_);
	weights_(0) = weight_0;
	for (int i = 1; i<2 * n_aug_ + 1; i++) {  //2n+1 weights
		double weight = 0.5 / (n_aug_ + lambda_);
		weights_(i) = weight;
	}

	//predicted state mean
	x.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x = x + weights_(i) * Xsig_pred_.col(i);
	}

	//predicted state covariance matrix
	P.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

											   // state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		//If x_diff(3) > PI then
			//x_diff(3) = ((x_diff(3) - M_PI) % 2 * M_PI) - M_PI

			//If x_diff(3) < -PI then
			//x_diff(3) = ((x_diff(3) - M_PI) % 2 * M_PI) + M_PI

		P = P + weights_(i) * x_diff * x_diff.transpose();
	}

	x_ = x;
	P_ = P;

	//cout << "Xsig_pred_ in PREDICTION" << Xsig_pred_<< endl;

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
  //create matrix for sigma points in measurement space
	MatrixXd Lsig = MatrixXd(n_l_, 2 * n_aug_ + 1);


	//mean predicted measurement
	VectorXd l_pred = VectorXd(n_l_);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_l_, n_l_);

	Lsig.fill(0.0);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

											   // extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		//double v = Xsig_pred_(2, i);
		//double yaw = Xsig_pred_(3, i);

		//double v1 = cos(yaw)*v;
		//double v2 = sin(yaw)*v;

		// measurement model
		Lsig(0, i) = p_x;
		Lsig(1, i) = p_y;
	}


	l_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		l_pred = l_pred + weights_(i) * Lsig.col(i);
	}


	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
											   //residual
		VectorXd l_diff = Lsig.col(i) - l_pred;

		S = S + weights_(i) * l_diff * l_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_l_, n_l_);
	R << std_laspx_ * std_laspx_, 0,
		0, std_laspy_ * std_laspy_;
	S = S + R;

	VectorXd l = meas_package.raw_measurements_;

	//// UKF Update assignment

	MatrixXd Tc = MatrixXd(n_x_, n_l_);

	/*******************************************************************************
	* Student part begin
	******************************************************************************/

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

											   //residual
		VectorXd l_diff = Lsig.col(i) - l_pred;


		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;


		Tc = Tc + weights_(i) * x_diff * l_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd l_diff = l - l_pred;


	//update state mean and covariance matrix
	x_ = x_ + K * l_diff;
	P_ = P_ - K*S*K.transpose();

	// Calculate NIS
	//NIS_laser_ = l_diff.transpose() * S.inverse() * l_diff;



}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	//cout << "Updating radar" << endl;
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */



	//// PREDICT RADAR MEASUREMENT ASSIGNMENT
  //create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	//cout << "radar 1" << endl;

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	//cout << "radar 2" << endl;

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	//cout << "radar 3" << endl;

	Zsig.fill(0.0);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		//cout << "radar 4.0.0" << endl;
		//cout << "Xsig_pred_ "  << Xsig_pred_ << endl;								   // extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		//cout << "radar 4.0.0" << endl;
		double p_y = Xsig_pred_(1, i);
		//cout << "radar 4.0.1" << endl;
		double v = Xsig_pred_(2, i);
		//cout << "radar 4.0.2" << endl;
		double yaw = Xsig_pred_(3, i);
		//cout << "radar 4.0.3" << endl;

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                        //r
		Zsig(1, i) = atan2(p_y, p_x);                                 //phi
		Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	//cout << "radar 4" << endl;

	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//cout << "radar 5" << endl;
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
											   //residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//cout << "radar 6" << endl;

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_radr_*std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;
	S = S + R;

	//cout << "radar 7" << endl;

	//// UKF Update assignment

	MatrixXd Tc = MatrixXd(n_x_, n_z);

	VectorXd z = meas_package.raw_measurements_;

	/*******************************************************************************
	* Student part begin
	******************************************************************************/

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

											   //residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}


	//cout << "radar 8" << endl;
	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();
	//cout << "radar end" << endl;
}
