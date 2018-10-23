//
// Created by mitc on 18-10-19.
//
#include "calibration.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include <iostream>
#include <fstream>

#include "matrixSolve.h"
#include "CeresOptimize.h"
using namespace cv;


void calibrate(string InputFileName)
{
    ifstream fin(InputFileName);

    Mat img;
    cout<<"start extracting corner point ..." << endl;

    int image_count = 0;
    Size image_size;
    Size board_size = Size(7, 7);
    vector<Point2f> image_point_buf;
    vector<vector<Point2f>> image_points_seq;
    string filename;
    while(getline(fin,filename))
    {
        image_count++;
        Mat imageInput = imread(filename);
        if (image_count == 1)
        {
            image_size.width=imageInput.cols;
            image_size.height=imageInput.rows;

        }
        if (0 == findChessboardCorners(imageInput,board_size,image_point_buf))
        {
            cout << filename << "can not find chessboard corners!" << endl;
            image_count--;
        }
        else
        {
            Mat view_gray;
            cvtColor(imageInput,view_gray,CV_RGB2GRAY);
            find4QuadCornerSubpix(view_gray,image_point_buf,Size(11,11));
            image_points_seq.push_back(image_point_buf);
            drawChessboardCorners(view_gray,board_size,image_point_buf,true);

        }
    }
    int total = image_points_seq.size();
    int CornerNum = board_size.width * board_size.height;
    for (int ii = 0; ii < total; ii++)
    {
        if (0 == ii % CornerNum)
        {
            int i = -1;
            i = ii / CornerNum;
            int j = i + 1;
        }
    }
    int nPoints=image_count*image_points_seq[0].size();
    cout << "corner extraction is done !"<<endl;

    cout << "Begin calibration ... " << endl;

    Size square_size = Size(21, 21);// board 3D information
    vector<vector<Point3f>> object_points;
    int i, j, t, ii, jj;
    for (t = 0; t < image_count; t++)
    {
        vector<Point3f> tempPointSet;
        for (i = 0; i < board_size.height ; i++)
        {
            for (j = 0; j< board_size.width ; j++)
            {
                Point3f realPoint;
                realPoint.x = i * square_size.width;
                realPoint.y = j * square_size.height;
                realPoint.z = 0;
                tempPointSet.push_back(realPoint);

            }
        }
        object_points.push_back(tempPointSet);
    }
    double imagePointsCopy[2*nPoints];
    double imagePoints[2*nPoints];
    double imagePointsProjected[2*nPoints];
    double K_distortion[5] = {0.0,0.0,0.0,0.0,0.0};
    double K_distortion_res[5] = {0.0,0.0,0.0,0.0,0.0};

    for( i=0; i < image_count; i++)
    {
        for (j=0; j < image_points_seq[i].size(); j++)
        {
            imagePoints[2*(i*image_points_seq[i].size())+j]=image_points_seq[i][j].x;
            imagePoints[2*(i*image_points_seq[i].size())+j+1]=image_points_seq[i][j].y;
            imagePointsCopy[2*(i*image_points_seq[i].size())+j]=image_points_seq[i][j].x;
            imagePointsCopy[2*(i*image_points_seq[i].size())+j+1]=image_points_seq[i][j].y;
        }
    }
    double wholeH[image_count*3*3], A[3*3], A_res[3*3], N[3*3], K[image_count*3*4],err,tmperr;
    double *_Hpt, *_Kpt, *imgPtsProj;
    double Ncopy[3*3];
    cout << "begin calhomograph" << endl;
    vector<vector<Point3d>> object_points3d;
    vector<vector<Point2d>> image_points_seq2d;
    for(ii=0; ii < image_count; ii++)
    {
        vector<Point3d> tmp_object3d;
        vector<Point2d> tmp_image2d;
        for(jj=0; jj < image_points_seq[ii].size(); jj++)
        {
            Point3d tmp3d;
            Point2d tmp2d;

            tmp3d.x=(double)object_points[ii][jj].x;
            tmp3d.y=(double)object_points[ii][jj].y;
            tmp3d.z=(double)object_points[ii][jj].z;
            tmp_object3d.push_back(tmp3d);

            tmp2d.x=(double)image_points_seq[ii][jj].x;
            tmp2d.y=(double)image_points_seq[ii][jj].y;
            tmp_image2d.push_back(tmp2d);

        }
        object_points3d.push_back(tmp_object3d);
        image_points_seq2d.push_back(tmp_image2d);
    }
    for(i=0;i<image_count;i++) {
        double image_points_seq_double[image_points_seq[i].size() * 2];
        double object_points_double[object_points[i].size() * 3];
        for (int k = 0; k < image_points_seq[i].size(); k++) {
            image_points_seq_double[2 * k] = image_points_seq[i][k].x;
            image_points_seq_double[2 * k + 1] = image_points_seq[i][k].y;

            object_points_double[3 * k] = object_points[i][k].x;
            object_points_double[3 * k + 1] = object_points[i][k].y;
            object_points_double[3 * k + 2] = object_points[i][k].z;

        }
        double *H = &wholeH[9 * i];
        CalHomograph(image_points_seq[0].size(), object_points_double, image_points_seq_double, H);
        OptimizeARt(object_points3d[i],image_points_seq2d[i],H);
        cout << endl << "------the--" << i << "th--image--Homograph------" << endl;
        for (int j = 0; j < 9; j++) {
            if (j % 3 == 0)
                cout << endl;
            cout << H[j] << " ";
        }
        imgPtsProj = &imagePointsProjected[2*image_points_seq[i].size()*i];
        Reprojection(image_points_seq[i].size(),H,object_points_double,imgPtsProj);
    }
    if(CalInternal(image_count, wholeH, N, A)) {
        for (i = 0; i < image_count; i++) {
            _Hpt = &wholeH[9 * i];
            _Kpt = &K[12 * i];

            CalExternal(_Hpt, A, _Kpt, Ncopy);
            tmperr = CalDistortion(image_points_seq[0].size() * 5, imagePointsCopy,
                          imagePointsProjected, A, K_distortion);
            if(i == 0)
            {
                err = tmperr;
                memcpy(A_res,A,sizeof(double)*9);
                memcpy(K_distortion_res,K_distortion,sizeof(double)*5);
            } else{
                if(tmperr < err)
                {
                    err = tmperr;
                    memcpy(A_res,A,sizeof(double)*9);
                    memcpy(K_distortion_res,K_distortion,sizeof(double)*5);
                }
            }
        }


    }
    cout << endl << "-------ins-------" << endl;
    for (int j = 0; j < 9; j++) {
        if (j % 3 == 0 && j != 0)
            cout << endl;
        cout << "A " << A[j] << " ";
    }

    cout << endl << "-------disto------" << endl;
    for (int j = 0; j < 5; j++) {
        if (j % 5 == 0 && j != 0)
            cout << endl;
        cout << "k " << K_distortion[j] << " ";
    }
}

void solveAx0(int m, int n, double* A, double* x)
{
    double* U   =(double*) malloc(sizeof(double)*m*n);
    double* D   =(double*) malloc(sizeof(double)*n);
    double* V   =(double*) malloc(sizeof(double)*n*n);
    double* tmp =(double*) malloc(sizeof(double)*m);

    matrixSVDAx0(A,m,n,x,U,D,V,tmp);
    if (1) {
        matrixScaVec(D,n,1/D[0],D);
    }
    free(tmp);
    free(V);
    free(D);
    free(U);
}

void CalHomograph(int nPoints, double* objectPoints, double* imagePoints, double* H)
{
    int k;

    double* L=(double*) malloc(2*nPoints*9*sizeof(double));

    for(k=0; k<nPoints; k++) {
        double X = objectPoints[3*k+0];
        double Y = objectPoints[3*k+1];
        double W = objectPoints[3*k+2];
        double u = imagePoints[2*k+0];
        double v = imagePoints[2*k+1];
        int    i=2*k;
        L[9*i+0] =    X; L[9*i+1] =    Y; L[9*i+2] =    W;
        L[9*i+3] =    0; L[9*i+4] =    0; L[9*i+5] =    0;
        L[9*i+6] = -u*X; L[9*i+7] = -u*Y; L[9*i+8] = -u*W;
        i=2*k+1;
        L[9*i+0] =    0; L[9*i+1] =    0; L[9*i+2] =    0;
        L[9*i+3] =    X; L[9*i+4] =    Y; L[9*i+5] =    W;
        L[9*i+6] = -v*X; L[9*i+7] = -v*Y; L[9*i+8] = -v*W;
    }

    solveAx0(2*nPoints,9,L,H);
    matrixScaMat(H, 1/1e-9,3,3,H );

    free(L);

}

void CalB(int nH, double* H, double* B)
{

    int m = 2*nH;
    int n = 4;
    double* V =(double*) calloc(sizeof(double),m*n);
    double* x =(double*) malloc(sizeof(double)*n);
    int i;

    double Uo = 0.0;
    double Vo = 0.0;
    double a,b,c,d,e,f;

    for(i=0;i<nH;i++)
    {
        double* h = &H[9*i];
        int line1 = 2*n*i;
        int line2 = line1+4;

        a =  h[0]*h[1];
        b =  (h[0]*h[4] + h[3]*h[1]);
        c =  h[3]*h[4];
        d =  (h[0]*h[7]+h[6]*h[1]);
        e =  (h[3]*h[7]+h[6]*h[4]);
        f =  h[6]*h[7];

        V[line1+0]= a - d*Uo + f*Uo*Uo;
        V[line1+1]= b - d*Vo - e*Uo + 2.0*f*Uo*Vo;
        V[line1+2]= c - e*Vo + f*Vo*Vo;
        V[line1+3]= f;

        a = h[0]*h[0] - h[1]*h[1];
        b = 2*(h[0]*h[3] - h[1]*h[4]);
        c = h[3]*h[3] - h[4]*h[4];
        d = 2*(h[0]*h[6] - h[1]*h[7]);
        e = 2*(h[3]*h[6] - h[4]*h[7]);
        f = h[6]*h[6] - h[7]*h[7];

        V[line2+0]= a - d*Uo + f*Uo*Uo;
        V[line2+1]= b - d*Vo - e*Uo + 2.0*f*Uo*Vo;
        V[line2+2]= c - e*Vo + f*Vo*Vo;
        V[line2+3]= f;

        matrixNormalizeVector(4,&V[line1]);
        matrixNormalizeVector(4,&V[line2]);

    }

    solveAx0(m,n,V,x);

    x[0]/=x[3];
    x[1]/=x[3];
    x[2]/=x[3];
    x[3]/=x[3];

    i=0;
    B[i++]=x[0];               B[i++]=x[1];               B[i++]= -Vo*x[1] - Uo*x[0];
    B[i++]=x[1];               B[i++]=x[2];               B[i++]= -Vo*x[2] - Uo*x[1];
    B[i++]=-Vo*x[1] - Uo*x[0]; B[i++]=-Vo*x[2] - Uo*x[1]; B[i++]= Vo*Vo*x[2] + 2.0*Vo*Uo*x[1] + Uo*Uo*x[0] + x[3];

    free(x);
    free(V);
}

int  CalA(double* B, double* A)
{
    double alpha,betha,gamma,u0,v0,lambda;
    double den = B[0]*B[4]-B[1]*B[1];

    if  (fabs(den) < 1e-10 )
    {
        return 0;
    }
    v0 = (B[1]*B[2]-B[0]*B[5]) / den;
    if (fabs(B[0]) < 1e-10)
    {
        return 0;
    }
    lambda = B[8]-(B[2]*B[2]+v0*(B[1]*B[2]-B[0]*B[5]))/B[0];
    if (lambda/B[0] < 0)
    {
        return 0;
    }
    alpha=sqrt(lambda / B[0]);
    if ((lambda*B[0] / den) < 0)
    {
        return 0;
    }
    betha = sqrt(lambda*B[0]/den);
    gamma = - B[1]*alpha*alpha*betha/lambda;
    u0  = gamma*v0/betha-B[2]*alpha*alpha/lambda;

    A[0]=alpha; A[1]=gamma; A[2]=u0;
    A[3]=0;     A[4]=betha; A[5]=v0;
    A[6]=0;     A[7]=0;     A[8]=1;

    return 1;
}

void correctA(double* AL, double* N, double* A)
{
    double Ninv[3*3];
    double det=matrixInv(N,Ninv);
    memcpy(N,Ninv,9*sizeof(double));
    matrixMultAB(Ninv,AL,A);
}

int CalInternal(int nHomographies, double* H, double N[9], double A[9])
{
    double B[3*3], AL[3*3];
    CalB(nHomographies,H,B);
    if(CalA(B,AL))
    {
        correctA(AL,N,A);
        return 1;
    }
    return 0;
}

int CalExternal(double H[9], double A[9], double K[12], double N[9])
{
    double Q[9],R[9],H_[9];
    double lambda1,lambda2,lambda3;
    double Ainv[3*3],hi[3], Ninv[3*3];
    double r1[3],r2[3],r3[3],t[3];
    int i;

    double lbda1,lbda2,lbda3;

    matrixInv(N,Ninv);
    matrixInv(A,Ainv);

    matrixAB(Ninv,H,3,3,3,H_);


    matrixCol(H_,0,3,3,hi);//h1
    matrixMultAb(Ainv,hi,r1);

    lbda1 = matrixNormalizeVector(3,r1);
    lambda1 = 1/lbda1;

    matrixCol(H_,1,3,3,hi);//h2
    matrixMultAb(Ainv,hi,r2);

    lbda2 = matrixNormalizeVector(3,r2);
    lambda2 = 1/lbda2;

    matrixCross(r1,r2,r3);

    matrixCol(H_,2,3,3,hi);//h3
    matrixMultAb(Ainv,hi,t);

    lbda3 = (lbda1 + lbda2)/2;
    lbda3 = 1./lbda3;

    lambda3 = (lambda1+lambda2)/2;
    matrixScaVec(t,3,lambda3,t);

    i=0;
    Q[i++]=r1[0]; Q[i++]=r2[0]; Q[i++]=r3[0];
    Q[i++]=r1[1]; Q[i++]=r2[1]; Q[i++]=r3[1];
    Q[i++]=r1[2]; Q[i++]=r2[2]; Q[i++]=r3[2];


    for(i=0;i<9;i++){
        R[i]= Q[i];
    };

    i=0;
    K[i++]=R[0]; K[i++]=R[1]; K[i++]=R[2]; K[i++]=t[0];
    K[i++]=R[3]; K[i++]=R[4]; K[i++]=R[5]; K[i++]=t[1];
    K[i++]=R[6]; K[i++]=R[7]; K[i++]=R[8]; K[i++]=t[2];

    return 1;
}

double CalDistortion(int nPts, double * imgPts,
                    double* imgPtsProjected, double A[9], double* k)
{
    int i;
    double u,v,r2,err;
    double *_imgPts_uv,*_imgPts_uv_proj;
    double* D =(double*) malloc(sizeof(double)*((nPts)*2 *5));
    double* d =(double*) malloc(sizeof(double)*((nPts)*2));

    for (i=0;i<nPts;i++) {

        u=imgPts[2*i];
        v=imgPts[2*i+1];
        r2=u*u+v*v;
        int k=2*i;
        D[5*k+0] = u*r2; D[5*k+1] = u*r2*r2; D[5*k+2] = u*r2*r2*r2; D[5*k+3] =    2*u*v; D[5*k+4] = r2+2*u*u;
        k=2*i+1;
        D[5*k+0] = v*r2; D[5*k+1] = v*r2*r2; D[5*k+2] = v*r2*r2*r2; D[5*k+3] = r2+2*v*v; D[5*k+4] = 2*u*v;

        _imgPts_uv = &imgPts[i*2];
        _imgPts_uv_proj = &imgPtsProjected[i*2];

        d[2*i] = _imgPts_uv[0] - _imgPts_uv_proj[0];
        d[2*i+1] = _imgPts_uv[1] - _imgPts_uv_proj[1];

        err += d[2*i]*d[2*i];
        err += d[2*i+1]*d[2*i+1];
    };

    matrixAxb(D,((nPts)*2), 2, d, k);
//    free(D);
//    free(d);
    return err;

}

void Reprojection(int nPoints, double* H, double* objectPoints, double* imgPtRepro)
{
    double res[3],object[4];
    int  k;

    for(k=0; k< nPoints; k++)
    {
        object[0] = objectPoints[3*k];
        object[1] = objectPoints[(3*k)+1];
        object[2] = 1.0;

        matrixAb(H,object,3,3,res);

        imgPtRepro[2*k] = res[0]/(res[2]+1e-10);
        imgPtRepro[(2*k)+1] = res[1]/(res[2]+1e-10);
    }
}

