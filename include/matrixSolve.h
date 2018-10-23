//
// Created by mitc on 18-10-20.
//

#ifndef CAMERACALIBRATION_MATRIXSOLVE_H
#define CAMERACALIBRATION_MATRIXSOLVE_H

using namespace std;
void matrixSVD(double* a, int m, int n, double* u, double* d, double* v, double* tmp);
double matrixSVDAx0(double* a, int m, int n, double* x, double* u, double* d, double* v, double* tmp);
void matrixSVDAxb(double* a, int m, int n, double* x, double* b, double* u, double* d, double* v, double* tmp);
void matrixSVDbacksub(double* u, double* d, double* v, int m, int n, double* b, double* x, double* tmp);
double matrixNormalizeVector(int n, double* v);
double matrixDotProd(double* u, double* v, int m);//dot product of vectors -> s={v}.{u}  Dimensions: {v}=m and  {u}=m
double matrixDet(double* A);
double matrixInv( double* A, double* Ainv );//computes the inverse,[Ainv], of the matrix [A]
void matrixAb(double* a, double* b, int m, int n, double* x);//computes the matrix product {x}=[A]{b}
void matrixAB(double* a, double* b, int m, int p, int n, double* x);//computes the matrix product [X]=[A][B]
void matrixCol(double* A,int col,  int m, int n, double* x);//gets a column -> {x}=col(A)
void matrixMultAb(double* A, double* b, double* x);//computes the transformation {x}=[A]{b}
void matrixCross(double *a,double *b,double *c);// computes the cross product {c} = {a}x{b}
void matrixScaVec(double* u, int m,  double s, double* x);//multiplies a vector by a scalar -> {x}=s{u}  Dimensions: {u}=m and {x}=m
void matrixMultAB(double* A, double* B, double* AB);//computes the matrix product [AB]=[A][B]
void matrixAt(double* a, int m, int n, double* x);//computes the transpose [X]=[A]T
void matrixAxb(double *A, int m, int n, double *b, double *x);//find solution  x to the system  Ax=b
void matrixScaMat(double* a, double s,int m, int n, double* x);//multiplies a matrix by a scalar
#endif //CAMERACALIBRATION_MATRIXSOLVE_H
