#include"numalg.h"
#include<stdlib.h>
#include<stdio.h>

double ** U_solve(int n,double ** U1,double ** b)
{
    double ** U=matCopy(n,n,U1);
    errorElim(n,n,U);
    //check
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<=i;j++)
        {
            if (((j<i)&&((U[i][j]!=0)&&(U[i][j]!=-0))))
            {
                printf("ERROR while performing U-solving: U should be upper-triangle.\n");
                exit(0);
            }
            if ((j==i)&&(U[i][j]==0))
            {
                printf("ERROR while performing U-solving: Zero on diagonal.\n");
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

double ** L_solve(int n,double ** L1,double ** b)
{
    double ** L=matCopy(n,n,L1);
    errorElim(n,n,L);
    for (int i=0;i<n;i++)
    {
        for (int j=i;j<n;j++)
        {
            if (((j>i)&&((L[i][j]!=0)&&(L[i][j]!=-0))))
            {
                printf("ERROR while performing L-solving: L should be lower-triangle.\n");
                exit(0);
            }
            if ((j==i)&&(L[i][j]==0))
            {
                printf("ERROR while performing L-solving: Zero on diagonal.\n");
                exit(0);
            }
        }
    }
    double ** ans=matCreate(n,1);
    ans[0][0]=b[0][0];
    for (int i=1;i<n;i++)
    {
        ans[i][0]=b[i][0];
        for (int j=0;j<i;j++)
        {
            ans[i][0]=ans[i][0]-L[i][j]*ans[j][0];
        }
        ans[i][0]=ans[i][0]/L[i][i];
    }
    return ans;
}

double *** decomp_simpleLU(int n,double ** A)
{
    double ** L=eye(n);
    double ** U=matCopy(n,n,A);
    double k;

    for (int i=0;i<=n-1;i++)
    {
        if (U[i][i]==0)
        {
            printf("ERROR while performing simple LU decomposition: Zero pivot.\n");
            exit(0);
        }
        for (int j=i+1;j<n;j++)
        {
            k=U[j][i]/U[i][i];
            rowTrans_add(n,n,U,i,-k,j);
            colTrans_add(n,n,L,j,k,i);
        }
    }
    double *** ans=(double ***)malloc(2*sizeof(double **));
    ans[0]=L;
    ans[1]=U;
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

double ** LESsolve_simpleLU(int n,double ** A,double ** b)
{
    double *** LU=decomp_simpleLU(n,A);
    double ** Y=L_solve(n,LU[0],b);
    double ** X=U_solve(n,LU[1],Y);
    return X;
}
