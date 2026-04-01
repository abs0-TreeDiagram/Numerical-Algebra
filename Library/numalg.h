#ifndef NUMALG_H
#define NUMALG_H

//basic matrix operation
void rowTrans_swap(int rowCount,int colCount,double ** pmat,int i,int j);
void rowTrans_scalarMult(int rowCount,int colCount,double ** pmat,int i,double k);
void rowTrans_add(int rowCount,int colCount,double ** pmat,int i1,double k,int i2);
void matScalarMult(int i,int j,double ** p,double k);
double ** matCreate(int i,int j);
void matFree(double ** p,int i,int j);
double ** matPlus(int i,int j,double ** p,int k,int l,double ** q);
double ** matMult(int i,int j,double ** p,int k,int l,double ** q);
void matDisp(int i,int j,double ** mat);
double ** matCopy(int i,int j,double ** mat);
double ** matTransfer(int i,int j,double * p);
double norm2(int n,double ** p);
double vectError2(int n,double ** p,double ** q);
double matmaxAbs(int i,int j,double ** p);
double matmaxAbsLoc(int i,int j,double ** p);

//solving linear equation set
double ** U_solve(int n,double ** U,double ** b);
double ** LESsolve_Gauss(int n,double ** A,double ** b);

double specRadius_powerIter(int n,double ** mat,double ** x0,int N,double epsilon);

#endif