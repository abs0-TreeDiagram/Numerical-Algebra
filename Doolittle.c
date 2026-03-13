#include<stdio.h>
#include<math.h>

//---------------------------以下为编辑区---------------------------

//行列数和增广矩阵

#define N 3//增广矩阵的行数

double mat[N][N]=//矩阵
{
    {2,2,3},
    {4,7,7},
    {-2,4,5}
};

//---------------------------以上为编辑区---------------------------

void row_act1(int i,int j)
{
    double temp;
    for (int k=0;k<=N-1;k++)
    {
        temp=mat[i][k];
        mat[i][k]=mat[j][k];
        mat[j][k]=temp;
    }
    printf("行变换：交换%d、%d行\n",i,j);
}

void row_act2(int i,double k)
{
    for (int j=0;j<=N-1;j++)
    {
        mat[i][j]*=k;
    }
}

void row_act3(int i1,double k,int i2)
{
    for (int j=0;j<=N-1;j++)
    {
        mat[i2][j]+=mat[i1][j]*k;
    }
    printf("行变换：%d行的%lf倍加到%d行\n",i1,k,12);
}

int main()
{
    //初始化L矩阵
    double L[N][N];
    for (int i=0;i<N;i++)
    {
        for (int j=0;j<N;j++)
        {
            L[i][j]=0;
        }
    }
    for (int i=0;i<N;i++)
    {
        L[i][i]=1;
    }


    for (int t=0;t<=N-1;t++)//主元行
    {
        if (mat[t][t]==0)
        {
            printf("ERROR\n");
            return -1;
        }
        for (int u=t+1;u<=N-1;u++)//主元行下每一行
        {
            L[u][t]=mat[u][t]/mat[t][t];
            row_act3(t,-mat[u][t]/mat[t][t],u);
        }
    }

    //输出

    printf("A=LU\n");
    printf("L=\n");

    for(int s=0;s<=N-1;s++)
    {
        for (int r=0;r<=N-1;r++)
        {
            printf("%lf ",L[s][r]);
        }
        printf("\n");
    }

    printf("U=\n");

    for(int s=0;s<=N-1;s++)
    {
        for (int r=0;r<=N-1;r++)
        {
            printf("%lf ",mat[s][r]);
        }
        printf("\n");
    }   
}