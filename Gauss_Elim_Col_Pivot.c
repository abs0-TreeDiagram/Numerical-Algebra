#include<stdio.h>
#include<math.h>


//---------------------------以下为编辑区---------------------------

//行列数和增广矩阵

#define N 3//增广矩阵的行数

double mat[N][N+1]=//矩阵
{
    {1,1,1,6},
    {12,-3,3,15},
    {-18,3,-1,-15}
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
    printf("行变换：交换%d、%d行\n",i,j);
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
    printf("行变换：%d行的%lf倍加到%d行\n",i1,k,12);
}

int main()
{
    for (int t=0;t<=N-1;t++)//主元行
    {
        //寻找最大行
        int imax=t;
        for (int u=t+1;u<=N-1;u++)
        {
            if (fabs(mat[imax][t])<fabs(mat[u][t]))
            {
                imax=u;
            }
        }
        //交换
        row_act1(t,imax);
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