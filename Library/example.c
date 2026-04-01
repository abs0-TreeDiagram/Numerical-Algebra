#include<stdlib.h>
#include<stdio.h>
#include<math.h>


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
