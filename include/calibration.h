//
// Created by mitc on 18-10-19.
//

#ifndef CAMERACALIBRATION_CALIBRATION_H
#define CAMERACALIBRATION_CALIBRATION_H

#include "opencv2/core/core.hpp"

using namespace std;

void calibrate(std::string InputFileName);
void CalHomograph(int nPoints, double* objectPoints, double* imagePoints, double* H);
void solveAx0(int m, int n, double* A, double* x);
void CalB(int nH, double* H, double* B);
int  CalA(double* B, double* A);
int CalInternal(int nHomographies, double* H, double N[9], double A[9]);
int CalExternal(double H[9], double A[9], double K[12], double N[9]);
double CalDistortion(int nPts, double * imgPts, double* imgPtsProjected, double A[9], double* k);
void Reprojection(int nPoints, double* H, double* objectPoints, double* imgPtRepro);
#endif //CAMERACALIBRATION_CALIBRATION_H
