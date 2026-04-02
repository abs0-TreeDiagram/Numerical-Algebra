#include"numalg.h"
#include<stdlib.h>
#include<stdio.h>

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
