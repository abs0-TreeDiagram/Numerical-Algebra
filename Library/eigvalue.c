#include"numalg.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double specRadius_powerIter(int n,double ** mat,double ** x0,int N,double epsilon)
{
    double ** v1=matCopy(n,1,x0);
    double ** v2;
    double l1,l2;
    for (int count=0;;count++)
    {
        
        v2=matMult(n,n,mat,n,1,v1);
        if (count>N)
        {
            printf("Iteration FAILED in specRadius_powerIter\n");
            return -1;
        }
        
        int maxLoc=matmaxAbsLoc(n,1,v2);
        l2=(fabs(v2[maxLoc][0]/v1[maxLoc][0]));
        if (count!=0)
        {
            if (fabs(l1-l2)<epsilon)
            {
                return l2;
            }
        }
        
        matScalarMult(n,1,v2,1/matmaxAbs(n,1,v2));
    
        matFree(v1,n,1);
        v1=matCopy(n,1,v2);
        l1=l2;
    }
}

double specRadius_powerIter_AtikenAcc(int n,double ** mat,double ** x0,int N,double epsilon)
{
	double ** v1, ** v2;//iter vect
	int maxLoc;
	double a[3];//direct quotient
	double l1,l2;//iter eigvalue
	
	//initialise: fill a[3]
	v1=matCopy(n,1,x0);
	v2=matCopy(n,1,x0);
	for (int i=0;i<3;i++)
	{
		maxLoc=matmaxAbsLoc(n,1,v1);
		matScalarMult(n,1,v1,1/v1[maxLoc][0]);
		matFree(v2,n,1);
		v2=matMult(n,n,mat,n,1,v1);
		a[i]=v2[maxLoc][0];
		matFree(v1,n,1);
		v1=matCopy(n,1,v2);
	}
	//iter
	l1=a[0]-((a[1]-a[0])*(a[1]-a[0]))/(a[2]-2*a[1]+a[0]);
	a[0]=a[1];
	a[1]=a[2];
	maxLoc=matmaxAbsLoc(n,1,v1);
	matScalarMult(n,1,v1,1/v1[maxLoc][0]);
	matFree(v2,n,1);
	v2=matMult(n,n,mat,n,1,v1);
	a[2]=v2[maxLoc][0];
	matFree(v1,n,1);
	v1=matCopy(n,1,v2);
	l2=a[0]-((a[1]-a[0])*(a[1]-a[0]))/(a[2]-2*a[1]+a[0]);
	for (int count;;count++)
	{
		if (count>N)
		{
			printf("Iteration FAILED in specRadius_powerIter_AtikenAcc\n");
			return -1;
		}
		else if(fabs(l2-l1)<epsilon)
		{
			return l2;
		}
		l1=l2;
		a[0]=a[1];
		a[1]=a[2];
		maxLoc=matmaxAbsLoc(n,1,v1);
		matScalarMult(n,1,v1,1/v1[maxLoc][0]);
		matFree(v2,n,1);
		v2=matMult(n,n,mat,n,1,v1);
		a[2]=v2[maxLoc][0];
		matFree(v1,n,1);
		v1=matCopy(n,1,v2);
		l2=a[0]-((a[1]-a[0])*(a[1]-a[0]))/(a[2]-2*a[1]+a[0]);
	}
}

double *** QR_decomp(int n,double ** p1)
{
    double ** p=matCopy(n,n,p1);
    double *** b=(double ***)malloc(n*sizeof(double**));

    for(int m=0;m<n;m++)
    {
        b[m]=matCol(n,n,p,m);
    }
    double ** R=matCreate(n,n);
    matFill(n,n,R,0);
    for (int m=0;m<n;m++)
    {
        R[m][m]=1;
    }
    
    double ** temp1,** temp2,temp3;
    for (int m=1;m<n;m++)
    {
        for (int i=0;i<m;i++)
        {
            temp1=matCopy(n,1,b[i]);
            temp3=vectInProd(n,temp1,b[m]);
            matScalarMult(n,1,temp1,-temp3);
            matPlus(n,1,b[m],n,1,temp1);
            R[i][m]=temp3;
        }
    }

    for (int m=0;m<n;m++)
    {
        temp3=norm2(n,b[m]);
        matScalarMult(n,1,b[m],1/temp3);
        rowTrans_scalarMult(n,n,R,m,temp3);
    }

    for (int m=1;m<n;m++)
    {
        colJoint(n,m,&b[0],b[m]);
    }

    double *** ans=(double ***)malloc(2*sizeof(double **));
    ans[0]=matCopy(n,n,b[0]);
    ans[1]=matCopy(n,n,R);
    for (int m=0;m<n;m++)
    {
        matFree(b[m],n,1);
    }
    free(b);
    matFree(R,n,n);
    return ans;
}
