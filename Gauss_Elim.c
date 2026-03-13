#include<stdio.h>



//---------------------------以下为编辑区---------------------------

//行列数和增广矩阵

#define N 3//增广矩阵的行数

double mat[N][N+1]=//矩阵
{
    {1,1,1,0},
    {1,2,1,0},
    {1,1,3,0}
};

//---------------------------以上为编辑区---------------------------

void row_act1(int i,int j)
{
    double temp;
    for (int k=0;k<=N;k++)
    {
        temp=mat[i][k];
        mat[i][k]=mat[j][k];
        mat[j][k]=temp;
    }
}

void row_act2(int i,double k)
{
    for (int j=0;j<=N;j++)
    {
        mat[i][j]*=k;
    }
}

void row_act3(int i1,double k,int i2)
{
    for (int j=0;j<=N;j++)
    {
        mat[i2][j]+=mat[i1][j]*k;
    }
}

int main()
{
    for (int t=0;t<=N-1;t++)//主元行
    {
        if (mat[t][t]==0)
        {
            printf("ERROR\n");
            return -1;
        }
        for (int u=t+1;u<=N-1;u++)//主元行下每一行
        {
            row_act3(t,-mat[u][t]/mat[t][t],u);
        }
    }

    //输出

    for(int s=0;s<=N-1;s++)
    {
        for (int r=0;r<=N;r++)
        {
            printf("%lf ",mat[s][r]);
        }
        printf("\n");
    }
}