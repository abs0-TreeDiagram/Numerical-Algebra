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

//---------------------------以上为matTrans以下为直接法LES----------------------

double ** U_solve(int n,double ** U,double ** b)
{
    //check
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<=i;j++)
        {
            if (((j<i)&&(U[i][j]!=0)))
            {
                printf("ERROR while performing U-solving: U should be upper-triangle.\n");
                exit(0);
            }
        }
    }
    //perform solving
    double ** ans=matCreate(n,1);
    for (int m=n-1;m>=0;m--)
    {
        ans[m][0]=b[m][0];
        for (int p=m+1;p<n;p++)
        {
            ans[m][0]-=U[m][p]*ans[p][0];
        }
        ans[m][0]/=U[m][m];
    }
    return ans;
}

double ** LESsolve_Gauss(int n,double ** A,double ** b)
{    
    double ** mat=matCopy(n,n,A);
    double ** B=matCopy(n,1,b);
    for (int t=0;t<=n-1;t++)//主元行
    {
        if (mat[t][t]==0)
        {
            printf("ERROR while performing simple Gaussian elimination:Zero pivot.\n");
            exit(0);
        }
        for (int u=t+1;u<=n-1;u++)//主元行下每一行
        {
            rowTrans_add(n,1,B,t,-mat[u][t]/mat[t][t],u);
            rowTrans_add(n,n,mat,t,-mat[u][t]/mat[t][t],u);
            
        }
    }
    double ** ans=U_solve(n,mat,B);

    matFree(mat,n,n);
    matFree(B,n,1);
    return ans;
}

//----------------------------以上为直接法LES以下为特征值特征根

double specRadius_powerIter(int n,double ** mat,double ** x0,int N,double epsilon)
{
    double ** v1=matCopy(n,1,x0);
    double ** v2;
    double l1,l2;double tt;
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
            tt=fabs(l1-l2);
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

int main(int argc, char *argv[]) 
{
	double A1[9]={
		2,-1,0,
		0,2,-1,
		0,-1,2
	};
	double b1[3]=
	{
		0,
		0,
		1
	};
	double ** A=matTransfer(3,3,&A1[0]);
	double ** b=matTransfer(3,1,&b1[0]);
	printf("%lf\n",specRadius_powerIter(3,A,b,1000,0.00001));
}