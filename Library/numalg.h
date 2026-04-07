#ifndef NUMALG_H
#define NUMALG_H

//basic matrix operation
void rowTrans_swap(int rowCount, int colCount, double **pmat, int i, int j);
void rowTrans_scalarMult(int rowCount, int colCount, double **pmat, int i, double k);
void rowTrans_add(int rowCount, int colCount, double **pmat, int i1, double k, int i2);
void colTrans_add(int rowCount, int colCount, double **pmat, int j1, double k, int j2);
void colTrans_scalarMult(int rowCount, int colCount, double **pmat, int j, double k);
void matScalarMult(int i, int j, double **p, double k);
double ** matCreate(int i, int j);
void matFree(double **p, int i, int j);
void matPlus(int i, int j, double **p, int k, int l, double **q);
double ** matMult(int i, int j, double **p, int k, int l, double **q);
void matDisp(int i, int j, double **mat);
double ** matCopy(int i, int j, double **mat);
double ** matTransfer(int i, int j, double *p);
double norm2(int n, double **p);
double vectError2(int n, double **p, double **q1);
double matmaxAbs(int i, int j, double **p);
double matmaxAbsLoc(int i, int j, double **p);
void matTranspose(int i, int j, double ***p);
double ** matCol(int i, int j, double **p, int k);
void colJoint(int i, int j, double ***pA, int k, double ***b);
void rowJoint(int i, int j, double ***pA, int k, double ***B);
double vectInProd(int n, double **p1, double **q);
void matFill(int i, int j, double **p, double k);
double ** eye(int n);
void errorElim(int i, int j, double **A);

//solving linear equation set
double ** U_solve(int n, double **U1, double **b);
double ** L_solve(int n, double **L1, double **b);
double *** decomp_simpleLU(int n, double **A);
double ** LESsolve_Gauss(int n, double **A, double **b);
double ** LESsolve_simpleLU(int n, double **A, double **b);

//eigvalue and eigvector
double maxEigvalue_powerIter(int n, double **mat, double **x0, int N, double epsilon);
double maxEigvalue_powerIter_AtikenAcc(int n, double **mat, double **x0, int N, double epsilon);
double *** decomp_QR(int n, double **p1);
double ** eigValue_basicQR(int n, double **A, int N, double epsilon);
double ** normalHouseholderTrans(int n, double **vect, int i);
double ** semiUpTriangle_Householder(int n, double **A);
double minEigvalue_invPowerIter(int n, double **mat, double apprx, double **x0, int N, double epsilon);

#endif
