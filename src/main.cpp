//
// Created by mitc on 18-10-18.
//
#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <fstream>

#include "matrixSolve.h"
#include "calibration.h"

using namespace std;
using namespace cv;

int mycalibrate()
{
    string InputFileName="../calibrationdata.txt";
    calibrate(InputFileName);
    return 0;
}

int main(){
    mycalibrate();
    return 0;
}