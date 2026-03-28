#include<stdbool.h>
#include<stdlib.h>
#include<stdio.h>

struct element
{
    bool isZero;
    double * valuePointer;
};


struct spsMatrix
{
    int rowSize;
    int columnSize;
    struct element * elementList;
};

struct spsMatrix * createSpsMat(int rowSize,int columnSize)
{
    struct spsMatrix * mat = (struct spsMatrix *)malloc(sizeof(struct spsMatrix));
    mat->rowSize = rowSize;
    mat->columnSize = columnSize;
    mat->elementList = (struct element *)malloc(sizeof(struct element) * rowSize * columnSize);
    for(int i = 0;i < rowSize * columnSize;i++)
    {
        mat->elementList[i].isZero = true;
        mat->elementList[i].valuePointer = NULL;
    }
    return mat;

};

void destroySpsMat(struct spsMatrix * mat)
{
    free(mat->elementList);
    free(mat);
}

void spsMatInput(struct spsMatrix * mat,int i,int j,double value)
{
    mat->elementList[i * mat->columnSize + j].isZero = false;
    mat->elementList[i * mat->columnSize + j].valuePointer = (double *)malloc(sizeof(double));
    *(mat->elementList[i * mat->columnSize + j].valuePointer) = value;
}

double spsMatGet(struct spsMatrix * mat,int i,int j)
{
    if(mat->elementList[i * mat->columnSize + j].isZero)
    {
        return 0;
    }
    else
    {
        return *(mat->elementList[i * mat->columnSize + j].valuePointer);
    }
}

void main()
{
    struct spsMatrix * mat = createSpsMat(20000,20000);
    for (int i=0;i < 10000;i++)
    {
        spsMatInput(mat,i,i,2.0);
    }
    for (int i=0;i < 10000;i++)
    {
        spsMatInput(mat,i,i+1,-1.0);
    }
    for (int i=0;i < 10000;i++)
    {
        spsMatInput(mat,i+1,i,-1.0);
    }
    printf("-----\n");
    printf("%lf\n",spsMatGet(mat,0,0));
}