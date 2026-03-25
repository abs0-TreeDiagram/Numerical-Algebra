#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//--------------------以下为编辑区--------------------

#define N 3//增广矩阵的行数

#define MaxIterCount 10000//最大迭代次数

#define Epsilon 0.000001//迭代准出误差

double mat[N][N+1]={//必须对称正定
    {2,-1,0,1},
    {-1,2,-1,0},
    {0,-1,2,1.8}
};

double x0[N]={0,0,0};//迭代初始点

//--------------------以上为编辑区--------------------

double xy(double a[N],double b[N])//向量内积
{
    double ans=0;
    for (int i=0;i<=N-1;i++)
    {
        ans+=a[i]*b[i];
    }
    return ans;
}

double * Ay(double A[N][N],double y[N])//形如Ay的矩阵与向量乘法
{
    double * ans=(double *)malloc(N*sizeof(double));
    for (int i=0;i<=N-1;i++)
    {
        ans[i]=0;
        for (int j=0;j<=N-1;j++)
        {
            ans[i]+=A[i][j]*y[j];
        }
    }
    return ans;
}

double xAy(double x[N],double A[N][N],double y[N])//形如xAy的两步乘法
{
    double * temp=Ay(A,y);
    double ans=0;
    ans=xy(x,temp);
    free(temp);
    return ans;
}

double Vect_Error_Norm(double x1[N], double x2[N])
{
    int n=N;
    double ans=0;
    for (int i=0;i<=n-1;i++)
    {
        ans+=(x1[i]-x2[i])*(x1[i]-x2[i]);
    }
    ans=sqrt(ans);
    return ans;
}

int main()
{
    //截取系数矩阵与右端向量
    double A[N][N],A1[N][N],b[N];
    for (int i=0;i<=N-1;i++)
    {
        for (int j=0;j<=N-1;j++)
        {
            A[i][j]=mat[i][j];
            A1[i][j]=mat[i][j]/2;
        }
        b[i]=-mat[i][N];
    }
    //初始化迭代数组
    double x1[N],x2[N];
    for (int i=0;i<=N-1;i++)
    {
        x1[i]=x0[i];
    }
    //每一步迭代，x1为迭代起点，x2为迭代结果。
    //每步结束后，先判敛，若未达准出条件，使x1=x2，进行下一步迭代
    //开始迭代
    double D[N];//梯度
    double t,app[2];//步长系数和关于t的二次多项式的系数
    double * temp;
    for (int count=0;;count++)
    {
        //计算梯度
        temp=Ay(A,x1);
        for (int i=0;i<=N-1;i++)
        {
            D[i]=temp[i];
        }
        free(temp);
        for (int i=0;i<=N-1;i++)
        {
            D[i]+=b[i];
        }
        //计算步长系数
        app[0]=xAy(D,A1,D);
        app[1]=2*xAy(x1,A1,D)+xy(D,b);
        if (app[0]==0)
        {
            t=0;
        }
        else {
        t=-app[1]/(2*app[0]);}
        //执行迭代
        for (int i=0;i<=N-1;i++)
        {
            x2[i]=x1[i]+t*D[i];
        }
        //输出每一步
/*
        for (int i=0;i<=N-1;i++)
            {
                printf("%lf ",x2[i]);
            }
            printf("\n");
*/
        //检查退出条件并返回值
        if (Vect_Error_Norm(x1,x2)<=Epsilon)
        {
            printf("Iteration SUCCEEDED.\nIteration Count:%d\nResult:\n",count);
            for (int i=0;i<=N-1;i++)
            {
                printf("%lf ",x2[i]);
            }
            printf("\n");
            break;
        }
        if (count>MaxIterCount)
        {
            printf("Iteration FAILED with iteration number exceeding the limit.");
            break;
        }
        //未退出
        for (int i=0;i<=N-1;i++)
        {
            x1[i]=x2[i];
        }
    }
}