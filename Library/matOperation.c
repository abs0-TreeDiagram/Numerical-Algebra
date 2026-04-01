#include "numalg.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

void rowTrans_swap(int rowCount,int colCount,double ** pmat,int i,int j)
{
    double * temp;
    temp=pmat[i];
    pmat[i]=pmat[j];
    pmat[j]=temp;
}

void rowTrans_scalarMult(int rowCount,int colCount,double ** pmat,int i,double k)

{
    double * p=&pmat[i][0];
    for (int j=0;j<colCount;j++)
    {
        *p=(*p)*k;
        p++;
    }
}

void rowTrans_add(int rowCount,int colCount,double ** pmat,int i1,double k,int i2)
{
    double * p,* q;
    p=&pmat[i1][0];
    q=&pmat[i2][0];
    for (int m=0;m<colCount;m++)
    {
        *q+=k*(*p);
        q++;
        p++;
    }
}

void matScalarMult(int i,int j,double ** p,double k)
{
    for (int m=0;m<i;m++)
    {
        rowTrans_scalarMult(i,j,p,m,k);
    }
}

double ** matCreate(int i,int j)
{
    double ** p=(double **)malloc(i*sizeof(double *));
    for (int m=0;m<i;m++)
    {
        p[m]=(double *)malloc(j*sizeof(double));
    }
    return p;
}

void matFree(double ** p,int i,int j)
{
    for (int m=0;m<i;m++)
    {
        free(p[m]);
    }
    free(p);
}

double ** matPlus(int i,int j,double ** p,int k,int l,double ** q)
{
    if ((i!=k)||(j!=l))
    {
        printf("ERROR while performing matrix addtion: unmatched matrix size.\n");
        exit(0);
    }
    else
    {
        double ** r=matCreate(i,j);
        for (int m=0;m<i;m++)
        {
            for (int n=0;n<j;n++)
            {
                r[m][n]=p[m][n]+q[m][n];
            }
        }
        return r;
    }
}

double ** matMult(int i,int j,double ** p,int k,int l,double ** q)
{
    if (j!=k)
    {
        printf("ERROR while performing matrix multiplication: unmatched matrix size.\n");
        exit(0);
    }
    else
    {
        double ** r=matCreate(i,l);
        for (int m=0;m<i;m++)
        {
            for (int n=0;n<l;n++)
            {
                r[m][n]=0;
                for (int o=0;o<j;o++)
                {
                    r[m][n]+=p[m][o]*q[o][n];
                }
            }
        }
        return r;
    }
}

void matDisp(int i,int j,double ** mat)
{
	for (int m=0;m<i;m++)
	{
		for (int n=0;n<j;n++)
		{
			printf("%lf",mat[m][n]);
		}
		printf("\n");
	}
}

double ** matCopy(int i,int j,double ** mat)
{
	double ** p=matCreate(i,j);
    for (int m=0;m<i;m++)
    {
        for (int n=0;n<j;n++)
        {
            p[m][n]=mat[m][n];
        }
    }
    return p;
}

double ** matTransfer(int i,int j,double * p)
{
    double ** ans=matCreate(i,j);
    for (int m=0;m<i;m++)
    {
        for (int n=0;n<j;n++)
        {
            ans[m][n]=p[m*j+n];
        }
    }
    return ans;
}

double norm2(int n,double ** p)
{
    double ans=0;
    for (int i=0;i<n;i++)
    {
        ans+=(p[i][0])*(p[i][0]);
    }
    ans=sqrt(ans);
    return ans;
}

double vectError2(int n,double ** p,double ** q)
{
    matScalarMult(n,1,q,-1);
    double ** r=matPlus(n,1,p,n,1,q);
    return norm2(n,r);
}

double matmaxAbs(int i,int j,double ** p)
{
    double ans=0;
    for (int m=0;m<i;m++)
    {
        for (int n=0;n<j;n++)
        {
            if (fabs(p[m][n])> 0)
            {
                ans=fabs(p[m][n]);
            }
        }
    }
    return ans;
}

double matmaxAbsLoc(int i,int j,double ** p)
{
    double ans=0,ans1[2]={0,0};
    for (int m=0;m<i;m++)
    {
        for (int n=0;n<j;n++)
        {
            if (fabs(p[m][n])> 0)
            {
                ans1[0]=m;
                ans1[1]=n;
            }
        }
    }
    return ans1[0]*j+ans1[1];
}