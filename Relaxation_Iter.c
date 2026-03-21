#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//--------------------以下为编辑区--------------------

#define N 3 //增广矩阵的行数

#define MaxIterCount 10000//最大迭代次数

#define Epsilon 0.00000001//迭代准出误差

#define Omega 0.7//松弛因子

double mat[N][N+1]=
{
    {10,-1,-2,72},
    {-1,10,-2,83},
    {-1,-1,5,42}
};

//--------------------以上为编辑区--------------------

void row_scalar_mult(double * A, int i, int j, int i1, double k)
{
    //矩阵A为i×j的，第i1行数乘k
    if(i1 < 0 || i1 >= i) {
        printf("ERROR: Row index out of bounds.\n");
        exit(-1);
    }
    for(int col = 0; col < j; col++) {
        A[i1 * j + col] *= k;
    }
}

double Vect_Error_Norm(double *x1, double *x2, int n)
{
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
    //迭代公式：x=Bx+g
    //提取对角阵D供后续使用
    double D[N];
    for (int i=0;i<N;i++)
    {
        D[i]=mat[i][i];
    }
    //求B和g
    double * pmat=&mat[0][0];
    for (int i=0;i<=N-1;i++)
    {
        row_scalar_mult(pmat,N,N+1,i,1.0/mat[i][i]);
    }
    //此时mat非增广部分为-B+I，增广部分为g
    double B[N*N];
    for (int i=0;i<N;i++)
    {
        for (int j=0;j<N;j++)
        {
            B[i*N+j]=-mat[i][j];
        }
    }
    for (int i=0;i<N;i++)
    {
        B[i*N+i]=0;
    }

    double g[N];
    for (int i=0;i<N;i++)
    {
        g[i]=mat[i][N];
    }

    double x1[N] = {0}, x2[N]= {0};
    for (int count=0;;count++)
    {
        for (int i=0;i<=N;i++)//计算第i个分量
        {
            x2[i]=0;
            for (int j=0;j<=i-1;j++)
            {
                x2[i]+=Omega*B[i * N+j]*x2[j];
            }
            for (int j=i;j<=N-1;j++)
            {
                x2[i]+=Omega*B[i*N+j]*x1[j];
            }
            x2[i]+=Omega*g[i];
            x2[i]+=(1-Omega)*x1[i];
        }
        //检查准出条件
        if (Vect_Error_Norm(x1, x2, N)<=Epsilon)
        {
            printf("Iteration succeed.\n");
            printf("Cycle count: %d\n",count);
            printf("Result: \n");
            for (int i=0;i<N;i++)
            {
                printf("%f\n",x2[i]);
            }
            break;
        }
        if (count>MaxIterCount)
        {
            printf("Iteration failed with cycle count beyond limit.\n");
            break;
        }
        //未退出，继续下一轮
        for (int i=0;i<N;i++)
        {
            x1[i]=x2[i];
        }
    }
}