#include"numalg.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

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