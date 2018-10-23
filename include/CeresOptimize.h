//
// Created by mitc on 18-10-22.
//

#ifndef CAMERACALIBRATION_CERESOPTIMIZE_H
#define CAMERACALIBRATION_CERESOPTIMIZE_H
#include <iostream>
#include <fstream>
#include <ceres/ceres.h>
#include <opencv2/core/core.hpp>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include <ceres/rotation.h>
using namespace std;
using namespace cv;
void rotationMatrixToEulerAngles(double* R, double* r);
void OptimizeARt(const vector<Point3d>& pts_3d, const vector<Point2d>& pts_2d, double* H);
void Optimize(const vector<Point3d>& pts_3d, const vector<Point2d>& pts_2d, double* rot,
              double* cere_tran, double* cere_intrinsic, double* cere_distCoeff);
#endif //CAMERACALIBRATION_CERESOPTIMIZE_H
