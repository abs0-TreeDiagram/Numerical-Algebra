#include"numalg.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double maxEigvalue_powerIter(int n,double ** mat,double ** x0,int N,double epsilon)
{
    double ** v1=matCopy(n,1,x0);
    double ** v2;
    double l1,l2;
    for (int count=0;;count++)
    {
        
        v2=matMult(n,n,mat,n,1,v1);
        if (count>N)
        {
            printf("Iteration FAILED in maxEigvalue_powerIter\n");
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

double maxEigvalue_powerIter_AtikenAcc(int n,double ** mat,double ** x0,int N,double epsilon)
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
			printf("Iteration FAILED in maxEigvalue_powerIter_AtikenAcc\n");
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

double *** decomp_QR(int n,double ** p1)
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
    double k,nm;
    for (int m=1;m<n;m++)
    {
        for (int i=0;i<m;i++)
        {
            //calc k
            nm=norm2(n,b[i]);
            k=vectInProd(n,b[m],b[i])/(nm*nm);
            //perform col trans
            colTrans_add(n,n,p,i,-k,m);
            rowTrans_add(n,n,R,m,k,i);
        }
        matFree(b[m],n,1);
        b[m]=matCol(n,n,p,m);
    }
    
    for (int m=0;m<n;m++)
    {
        nm=norm2(n,b[m]);
        colTrans_scalarMult(n,n,p,m,1.0/nm);
        rowTrans_scalarMult(n,n,R,m,nm);
    }
    double *** ans=(double ***)malloc(2*sizeof(double **));
    ans[0]=p;
    ans[1]=R;
    return ans;

}

double ** eigValue_basicQR(int n,double ** A,int N,double epsilon)
{
    double ** A1, ** A2, ***temp, err,**ans;
    A1=matCopy(n,n,A);
    temp=decomp_QR(n,A1);
    A2=matMult(n,n,temp[1],n,n,temp[0]);
    
    for (int count=0;;count++)
    {//matDisp(n,n,A2);
        if (count>N)
        {
            printf("Iteration FAILED in eigValue_basicQR.\n");
            return NULL;
        }
        else
        {
            err=0;
            for (int i=0;i<n;i++)
            {
                err+=fabs(A1[i][i]-A2[i][i]);
            }
            if (err<epsilon)
            {
                ans=matCreate(1,n);
                matFill(1,n,ans,0);
                for (int i=0;i<n;i++)
                {
                    ans[0][i]=A2[i][i];
                }
                return ans;
            }
        }
        matFree(A1,n,n);
        A1=matCopy(n,n,A2);
        matFree(temp[0],n,n);
        matFree(temp[1],n,n);
        temp=decomp_QR(n,A1);
        matFree(A2,n,n);
        A2=matMult(n,n,temp[1],n,n,temp[0]);
    }
}

double ** normalHouseholderTrans(int n,double ** vect,int i)
{
    double ** w=matCopy(n,1,vect);
    w[i][0]=w[i][0]+norm2(n,vect);
    matScalarMult(n,1,w,1.0/norm2(n,w));
    double ** H=eye(n);
    double ** wT=matCopy(n,1,w);
    matTranspose(n,1,&wT);
    double ** wwT2=matMult(n,1,w,1,n,wT);
    matScalarMult(n,n,wwT2,-2);
    matPlus(n,n,H,n,n,wwT2);
    matScalarMult(n,n,H,-1);
    return H;
}

double ** semiUpTriangle_Householder(int n,double ** A)
{
    if (n<=2)
    {
        return A;
    }
    double ** A1=matCopy(n,n,A);
    double ** elim=matCreate(n-1,1);
    for (int i=0;i<n-1;i++)
    {
        elim[i][0]=A1[i+1][0];
    }
    double ** G=normalHouseholderTrans(n-1,elim,0);
    colJoint(n-1,1,NULL,n-1,&G);
    rowJoint(1,n,NULL,n-1,&G);
    G[0][0]=1;

    double ** H=matCopy(n,n,G);
    double ** temp=matMult(n,n,H,n,n,A1);

    double ** A2=matMult(n,n,temp,n,n,H);
    if (n==3)
    {
        return A2;
    }

    for (int i=1;i<=n-3;i++)
    {
        matFree(A1,n,n);
        A1=matCopy(n,n,A2);
        matFree(A2,n,n);

        matFree(elim,n-i,1);
        elim=matCreate(n-i-1,1);
        for (int j=0;j<n-i-1;j++)
        {
            elim[j][0]=A1[i+1+j][i];
        }
        
        matFree(G,n-i,n-i);
        G=normalHouseholderTrans(n-i-1,elim,0);

        H=eye(i+1);
        colJoint(i+1,i+1,&H,n-i-1,NULL);
        colJoint(i+1,n-i-1,NULL,n-i-1,&G);
        rowJoint(i+1,n,&H,n-i-1,&G);

        matFree(temp,n,n);
        temp=matMult(n,n,H,n,n,A1);
        A2=matMult(n,n,temp,n,n,H);
    }
    return A2;
}



/*
double ** eigValue_HouseholderAcc(int n,double ** A,int N,double epsilon)
{
    double ** U=semiUpTriangle_Householder(n,A);
    


    return ans;
}*/

double minEigvalue_invPowerIter(int n,double ** mat,double apprx,double ** x0,int N,double epsilon)
{
    double ** x1,** x2;
    x1=matCopy(n,1,x0);
    int r=matmaxAbsLoc(n,1,x1);
    matScalarMult(n,1,x1,1.0/matmaxAbs(n,1,x1));
    double ** mat1=matCopy(n,n,mat);
    for (int i=0;i<n;i++)
    {
        mat1[i][i]=mat1[i][i]-apprx;
    }
    double *** LU=decomp_simpleLU(n,mat1);
    double ** y=L_solve(n,LU[0],x1);
    x2=U_solve(n,LU[1],y);
    double l1,l2;
    l1=x1[r][0]/x2[r][0];

    matFree(x1,n,1);
    x1=matCopy(n,1,x2);
    r=matmaxAbsLoc(n,1,x1);
    matScalarMult(n,1,x1,1.0/matmaxAbs(n,1,x1));
    y=L_solve(n,LU[0],x1);
    x2=U_solve(n,LU[1],y);
    l2=x1[r][0]/x2[r][0];

    for (int count=1;;count++)
    {
        if (count>N)
        {
            printf("Iteration FAILED in minEigvalue_invPowerIter.\n");
            return -1;
        }
        else if(fabs(l2-l1)<epsilon)
        {
            return l2+apprx;
        }
        l1=l2;
        matFree(x1,n,1);
        x1=matCopy(n,1,x2);
        r=matmaxAbsLoc(n,1,x1);
        matScalarMult(n,1,x1,1.0/matmaxAbs(n,1,x1));
        y=L_solve(n,LU[0],x1);
        x2=U_solve(n,LU[1],y);
        l2=x1[r][0]/x2[r][0];
    }
}
