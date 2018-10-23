//
// Created by mitc on 18-10-20.
//
#define SIGN(a,b) ((b) >= 0. ? fabs(a) : -fabs(a))
#include "matrixSolve.h"
#include <math.h>
#include <iostream>

void matrixSVD(double*a, int m, int n, double* u, double* d, double* v, double* tmp)
{
    int flag,i,its,j,jj,k,l,nm;
    double anorm,c,f,g,h,s,scale,x,y,z;

    for(i=0;i<m;i++)
        for(j=0;j<n;j++)
            u[i*n+j]=a[i*n+j];

    g=scale=anorm=0.;
    for (i=0;i<n;i++) {
        l=i+2;
        tmp[i]=scale*g;
        g=s=scale=0.;
        if (i < m) {
            for (k=i;k<m;k++) scale +=  fabs(u[k*n+i]);
            if (scale != 0.) {
                for (k=i;k<m;k++) {
                    u[k*n+i] /= (scale+1e-10);
                    s += u[k*n+i]*u[k*n+i];
                }
                f=u[i*n+i];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                u[i*n+i]=f-g;
                for (j=l-1;j<n;j++) {
                    for (s=0.,k=i;k<m;k++) s += u[k*n+i]*u[k*n+j];
                    f=s/(h+1e-10);
                    for (k=i;k<m;k++) u[k*n+j] += f*u[k*n+i];
                }
                for (k=i;k<m;k++) u[k*n+i] *= scale;
            }
        }
        d[i]=scale *g;
        g=s=scale=0.;
        if (i+1 <= m && i+1 != n) {
            for (k=l-1;k<n;k++) scale += fabs(u[i*n+k]);
            if (scale!=0.) {
                for (k=l-1;k<n;k++) {
                    u[i*n+k] /= (scale+1e-10);
                    s += u[i*n+k]*u[i*n+k];
                }
                f=u[i*n+l-1];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                u[i*n+l-1]=f-g;
                for (k=l-1;k<n;k++) tmp[k]=u[i*n+k]/(h+1e-10);
                for (j=l-1;j<m;j++) {
                    for (s=0.,k=l-1;k<n;k++) s += u[j*n+k]*u[i*n+k];
                    for (k=l-1;k<n;k++) u[j*n+k] += s*tmp[k];
                }
                for (k=l-1;k<n;k++) u[i*n+k] *= scale;
            }
        }
        anorm=(anorm>(fabs(d[i])+fabs(tmp[i]))?anorm:(fabs(d[i])+fabs(tmp[i])));
    }
    for (i=n-1;i>=0;i--) {
        if (i < n-1) {
            if (g!=0.) {
                for (j=l;j<n;j++)
                    v[j*n+i]=(u[i*n+j]/(u[i*n+l]+1e-10))/(g+1e-10);
                for (j=l;j<n;j++) {
                    for (s=0.,k=l;k<n;k++) s += u[i*n+k]*v[k*n+j];
                    for (k=l;k<n;k++) v[k*n+j] += s*v[k*n+i];
                }
            }
            for (j=l;j<n;j++) v[i*n+j]=v[j*n+i]=0.;
        }
        v[i*n+i]=1.;
        g=tmp[i];
        l=i;
    }
    for (i=(m<n?m:n)-1;i>=0;i--) {
        l=i+1;
        g=d[i];
        for (j=l;j<n;j++) u[i*n+j]=0.;
        if (g != 0.) {
            g=1./g;
            for (j=l;j<n;j++) {
                for (s=0.,k=l;k<m;k++) s += u[k*n+i]*u[k*n+j];
                f=(s/u[i*n+i])*g;
                for (k=i;k<m;k++) u[k*n+j] += f*u[k*n+i];
            }
            for (j=i;j<m;j++) u[j*n+i] *= g;
        } else for (j=i;j<m;j++) u[j*n+i]=0.;
        ++u[i*n+i];
    }
    for (k=n-1;k>=0;k--) {
        for (its=0;its<30;its++) {
            flag=1;
            for (l=k;l>=0;l--) {
                nm=l-1;
                if ((fabs(tmp[l])+anorm) == anorm) {
                    flag=0;
                    break;
                }
                if ((fabs(d[nm])+anorm) == anorm) break;
            }
            if (flag) {
                c=0.;
                s=1.;
                for (i=l;i<k+1;i++) {
                    f=s*tmp[i];
                    tmp[i]=c*tmp[i];
                    if ((fabs(f)+anorm) == anorm) break;
                    g=d[i];
                    h=sqrt(f*f+g*g);
                    d[i]=h;
                    h=1./(h+1e-10);
                    c=g*h;
                    s = -f*h;
                    for (j=0;j<m;j++) {
                        y=u[j*n+nm];
                        z=u[j*n+i];
                        u[j*n+nm]=y*c+z*s;
                        u[j*n+i]=z*c-y*s;
                    }
                }
            }
            z=d[k];
            if (l == k) {
                if (z < 0.) {
                    d[k] = -z;
                    for (j=0;j<n;j++) v[j*n+k] = -v[j*n+k];
                }
                break;
            }
            if (its == 49) return ;
            x=d[l];
            nm=k-1;
            y=d[nm];
            g=tmp[nm];
            h=tmp[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0f*h*y);
            g=sqrt(f*f+1.);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/(x+1e-10);
            c=s=1.;
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=tmp[i];
                y=d[i];
                h=s*g;
                g=c*g;
                z=sqrt(f*f+h*h);
                tmp[j]=z;
                c=f/(z+1e-10);
                s=h/(z+1e-10);
                f=x*c+g*s;
                g = g*c-x*s;
                h=y*s;
                y *= c;
                for (jj=0;jj<n;jj++) {
                    x=v[jj*n+j];
                    z=v[jj*n+i];
                    v[jj*n+j]=x*c+z*s;
                    v[jj*n+i]=z*c-x*s;
                }
                z=sqrt(f*f+h*h);
                d[j]=z;
                if (z) {
                    z=1.0f/(z+1e-10);
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=0;jj<m;jj++) {
                    y=u[jj*n+j];
                    z=u[jj*n+i];
                    u[jj*n+j]=y*c+z*s;
                    u[jj*n+i]=z*c-y*s;
                }
            }
            tmp[l]=0.;
            tmp[k]=f;
            d[k]=x;
        }
    }
}

double matrixSVDAx0(double* a, int m, int n, double* x, double* u, double* d, double* v, double* tmp)
{
    double wmax,wmin,wmin2;
    int i,j,jmin;

    matrixSVD(a,m,n,u,d,v,tmp);

    wmax=d[0]; wmin = d[0];
    for (j=1;j<n;j++){
        if (d[j] < wmin && d[j] >1e-3) { wmin=d[j]; jmin =j; }
        if (d[j] > wmax) wmax=d[j];
    }

    wmin2=wmax;
    for (j=0;j<n;j++){
        if (j==jmin) continue;
        if (d[j] < wmin2 && d[j]>1e-3) wmin2=d[j];
    }

    for (i=0;i<n;i++)
        x[i] = v[i*n+jmin];

    return (wmin/wmax);
}

void matrixSVDAxb(double* a, int m, int n, double* x, double* b, double* u, double* d, double* v, double* tmp)
{
    double wmax,wmin;
    int j;

    matrixSVD(a,m,n,u,d,v,tmp);

    wmax=d[0];
    for (j=1;j<n;j++)
        if (d[j] > wmax) wmax = d[j];

    wmin=1.e-6*wmax;
    for (j=0;j<n;j++)
        if (d[j]<wmin) d[j]=0;

    matrixSVDbacksub(u,d,v,m,n,b,x,tmp);
}

#define TINY 1.0e-20

void matrixSVDbacksub(double* u, double* d, double* v, int m, int n, double* b, double* x, double* tmp)
{
    int j,i;
    double s;

    for (j=0;j<n;j++) {
        s=0.0;
        if (d[j]<TINY) {
            for (i=0;i<m;i++) s += u[i*n+j]*b[i]; /* computes [U]T{b} */
            s /= d[j];   /* multiply by [d]-1 */
        }
        tmp[j]=s;
    }
    for (i=0;i<n;i++) {
        s=0.0;
        for (j=0;j<n;j++) s += v[i*n+j]*tmp[j];  /* computes [V]{tmp} */
        x[i]=s;
    }
}

double matrixNormalizeVector(int n, double* v)
{
    double norm=matrixDotProd(v,v,n);
    int i;
    if (norm>1e-10) {
        norm = sqrt(norm);
        for (i=0;i<n;i++) v[i]/=norm;
    }
    return norm;
}

double matrixDotProd(double* u, double* v, int m)
{
    double tmp=0;
    int i;
    for (i=0;i<m;i++)
        tmp+=u[i]*v[i];
    return tmp;
}
double matrixDet(double* A)
{
    return A[0]*A[4]*A[8]+A[1]*A[5]*A[6]+A[2]*A[3]*A[7]-A[6]*A[4]*A[2]-A[7]*A[5]*A[0]-A[8]*A[3]*A[1];
}

double matrixInv( double* A, double* Ainv )
{
    double det = matrixDet ( A );

    if ( fabs(det)<1e-10 ) return det;

    Ainv[0] =  (A[4]*A[8]-A[7]*A[5])/det;
    Ainv[1] = -(A[1]*A[8]-A[7]*A[2])/det;
    Ainv[2] =  (A[1]*A[5]-A[4]*A[2])/det;

    Ainv[3] = -(A[3]*A[8]-A[6]*A[5])/det;
    Ainv[4] =  (A[0]*A[8]-A[6]*A[2])/det;
    Ainv[5] = -(A[0]*A[5]-A[3]*A[2])/det;

    Ainv[6] =  (A[3]*A[7]-A[6]*A[4])/det;
    Ainv[7] = -(A[0]*A[7]-A[6]*A[1])/det;
    Ainv[8]=   (A[0]*A[4]-A[3]*A[1])/det;

    return det;
}

void matrixAB(double* a, double* b, int m, int p, int n, double* x)
{
    int i,j,k;
    for (i=0;i<m;i++) {
        for (j=0;j<n;j++) {
            x[i*n+j]=0.;
            for (k=0;k<p;k++)
                x[i*n+j]+=a[i*p+k]*b[k*n+j];
        }
    }
}

void matrixCol(double* A,int col,  int m, int n, double* x)
{
    int i;
    for (i=0;i<m;i++){
        x[i]=A[i*n+col];
    }
}

void matrixAb(double* a, double* b, int m, int n, double* x)
{
    int i,j;
    for (i=0;i<m;i++) {
        x[i]=0.;
        for (j=0;j<n;j++)
            x[i]+=a[i*n+j]*b[j];
    }
}

void matrixMultAb(double* A, double* b, double* x)
{
    x[0]=A[0]*b[0]+A[1]*b[1]+A[2]*b[2];
    x[1]=A[3]*b[0]+A[4]*b[1]+A[5]*b[2];
    x[2]=A[6]*b[0]+A[7]*b[1]+A[8]*b[2];
}

void matrixCross(double *a,double *b,double *c)
{
    c[0]= a[1]*b[2] - a[2]*b[1];
    c[1]= a[2]*b[0] - a[0]*b[2];
    c[2]= a[0]*b[1] - a[1]*b[0];
}

void matrixScaVec(double* u, int m,  double s, double* x)
{
    int i;
    for (i=0;i<m;i++) x[i]=s*u[i];
}

void matrixMultAB(double* A, double* B, double* AB)
{
    AB[0]=A[0]*B[0]+A[1]*B[3]+A[2]*B[6];
    AB[1]=A[0]*B[1]+A[1]*B[4]+A[2]*B[7];
    AB[2]=A[0]*B[2]+A[1]*B[5]+A[2]*B[8];

    AB[3]=A[3]*B[0]+A[4]*B[3]+A[5]*B[6];
    AB[4]=A[3]*B[1]+A[4]*B[4]+A[5]*B[7];
    AB[5]=A[3]*B[2]+A[4]*B[5]+A[5]*B[8];

    AB[6]=A[6]*B[0]+A[7]*B[3]+A[8]*B[6];
    AB[7]=A[6]*B[1]+A[7]*B[4]+A[8]*B[7];
    AB[8]=A[6]*B[2]+A[7]*B[5]+A[8]*B[8];
}

void matrixAt(double* a, int m, int n, double* x)
{
    int i,j;
    for (i=0;i<m;i++)
        for (j=0;j<n;j++)
            x[j*m+i]=a[i*n+j];
}

/*
*						[A]x = b
*				    [At][A]x = [At]b
*	 inv([At][A])	[At][A]x =	inv([At][A])[At]b
*						[I]x =  inv([At][A])[At]b
*/
void matrixAxb(double *A, int m, int n, double *b, double *x)
{
    double *At, *AtA, *invAtA, *invAtA_At;

    At = (double*) malloc(sizeof(double)*n*m);
    AtA = (double*) malloc(sizeof(double)*m*m);
    invAtA = (double*) malloc(sizeof(double)*n*n);
    invAtA_At = (double*) malloc(sizeof(double)*n*m);

    matrixAt(A, m, n, At);
    matrixAB(At,A,n,m,n,AtA);
//    matrixInv(AtA,invAtA);
    matrixAB(invAtA,At,n,n,m,invAtA_At);
    matrixAb(invAtA_At,b,n,m,x);

    return;
}

void matrixScaMat(double* a, double s,int m, int n, double* x)
{
    int i,j;
    for (i=0;i<m;i++)
        for (j=0;j<n;j++)
            x[i*n+j]=s*a[i*n+j];
}