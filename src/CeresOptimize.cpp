//
// Created by mitc on 18-10-21.
//
#include "CeresOptimize.h"
using namespace std;
using namespace ceres;
using namespace cv;

struct ProjectionErrorARt{
    ProjectionErrorARt(Point3d objectPoint,Point2d imagePoint)
            : _objectPoint(objectPoint),_imagePoint(imagePoint){}

    template <typename T>
    bool operator()(const T* H, T* residuals) const {
        T M[3];
        M[0]=T(_objectPoint.x);
        M[1]=T(_objectPoint.y);
        M[2]=T(_objectPoint.z);
        T imageM[3];
        imageM[2]=H[6]*M[0]+H[7]*M[1]+H[8];
        imageM[0]=(H[0]*M[0]+H[1]*M[1]+H[2]*M[2])/imageM[2];
        imageM[1]=(H[3]*M[0]+H[4]*M[2]+H[5])/imageM[2];


        residuals[0]=_imagePoint.x - imageM[0];
        residuals[1]=_imagePoint.y - imageM[1];

        return true;
    }

private:
    Point3d _objectPoint;
    Point2d _imagePoint;
};

struct ProjectionError{
    ProjectionError(Point3d objectPoint,Point2d imagePoint)
            : _objectPoint(objectPoint),_imagePoint(imagePoint){}

    template <typename T>
    bool operator()(const T* RotationEuler, const T* TranslationMatrix,
                    const T* IntrinsicMatrix, const
                    T* DistCoeff,
                    T* residuals) const {
        T M[3];
        M[0]=T(_objectPoint.x);
        M[1]=T(_objectPoint.y);
        M[2]=T(_objectPoint.z);
        T MC[3];
        AngleAxisRotatePoint(RotationEuler,M,MC);
        MC[0]+=TranslationMatrix[0];
        MC[1]+=TranslationMatrix[1];
        MC[2]+=TranslationMatrix[2];

        T imageM[2];
        imageM[0]=(IntrinsicMatrix[0]*MC[0]+IntrinsicMatrix[1]*MC[1]+IntrinsicMatrix[3]*MC[2])/MC[2];
        imageM[1]=(IntrinsicMatrix[3]*MC[0]+IntrinsicMatrix[4]*MC[2])/MC[2];

        T predictPoint[2];
        T r2=imageM[0]*imageM[0]+imageM[1]*imageM[1];
        predictPoint[0]=imageM[0]*(T(1)+DistCoeff[0]*r2+DistCoeff[1]*r2*r2+DistCoeff[2]*r2*r2*r2)
                        +T(2)*DistCoeff[3]*imageM[0]*imageM[1]+DistCoeff[4]*(r2+T(2)*imageM[0]*imageM[0]);
        predictPoint[1]=imageM[1]*(T(1)+DistCoeff[0]*r2+DistCoeff[1]*r2*r2+DistCoeff[2]*r2*r2*r2)
                        +DistCoeff[3]*(r2+T(2)*imageM[0]*imageM[0])+T(2)*DistCoeff[4]*imageM[0]*imageM[1];


        residuals[0]=_imagePoint.x - predictPoint[0];
        residuals[1]=_imagePoint.y - predictPoint[1];

        return true;
    }

private:
    Point3d _objectPoint;
    Point2d _imagePoint;
};

void rotationMatrixToEulerAngles(double* R, double* r)
{

    float sy = sqrt(R[0] * R[0] +  R[3] * R[3] );

    bool singular = sy < 1e-6;

    float x, y, z;
    if (!singular)
    {
        x = atan2(R[6], R[7]);
        y = atan2(-R[6], sy);
        z = atan2(R[3], R[0]);
    }
    else
    {
        x = atan2(-R[5], R[4]);
        y = atan2(-R[6], sy);
        z = 0;
    }
    r[0]=x;
    r[1]=y;
    r[2]=z;
}
void OptimizeARt(const vector<Point3d>& pts_3d, const vector<Point2d>& pts_2d,double* H)
{
    int i;
    Problem problem;
    for(i=0;i<pts_3d.size();i++)
    {
        CostFunction* costfunction=new ceres::AutoDiffCostFunction<ProjectionErrorARt,2,9>(new ProjectionErrorARt(pts_3d[i],pts_2d[i]));
        problem.AddResidualBlock(costfunction,NULL,H);
    }
    cout << "generate ART options" << endl;
    Solver::Options options;
    options.linear_solver_type=DENSE_QR;
    options.minimizer_progress_to_stdout=false;
    Solver::Summary summary;
    cout << "ART start solve ..." << endl;
    Solve(options,&problem,&summary);

    cout << "ART solve is done !" << endl;
}

void Optimize(const vector<Point3d>& pts_3d, const vector<Point2d>& pts_2d,double* rot,
                 double* cere_tran, double* cere_intrinsic, double* cere_distCoeff)
{
    double* cere_rot;
    rotationMatrixToEulerAngles(rot, cere_rot);
    int i;
    Problem problem;
    for(i=0;i<pts_3d.size();i++)
    {
        CostFunction* costfunction=new ceres::AutoDiffCostFunction<ProjectionError,2,3,3,5,5>(new ProjectionError(pts_3d[i],pts_2d[i]));
        problem.AddResidualBlock(costfunction,NULL,cere_rot,cere_tran,cere_intrinsic,cere_distCoeff);
    }
    cout << "generate PARAMETER options" << endl;
    Solver::Options options;
    options.linear_solver_type=DENSE_QR;
    options.minimizer_progress_to_stdout=false;
    Solver::Summary summary;
    cout << "start PARAMETER solve ..." << endl;
    Solve(options,&problem,&summary);

    cout << "PARAMETER solve is done !" << endl;
}