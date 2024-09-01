// digitRecognition_224101059.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<conio.h>
#include "string.h"
#include "stdio.h"
#include<stdlib.h>
#include<cstdlib>
#include<time.h>
#include<math.h>
#define framesize 320
#define sample_count 20
#define p 12
#define pii 22.0/7.0
#define size 30000
typedef long double lld;

int maxx(int a,int b)
{
	if(a>b)
	return a;
	return b;
}

void dc_shift(lld arr[])
{
	lld avg,count=0;
	int n=sizeof(arr);
	for(int i=0;i<n;i++)
		count+=arr[i];
	avg=count/n;
	for(int i=0;i<n;i++)
		arr[i]=arr[i]-avg;
}

void normalisation(lld arr[])
{
	lld max_val=0,min_val=0;
	int n=sizeof(arr);
	for(int i=0;i<n;i++)
	{
		if(arr[i]>max_val)
		max_val=arr[i];
		if(arr[i]<min_val)
		min_val=arr[i];
	}
	max_val=maxx(max_val,0-min_val);
	for(int i=0;i<n;i++)
	arr[i]=(arr[i]*5000)/max_val;
}

void calculateRi(lld arr[],lld R[])
{
	for(int i=0;i<=p;i++)
	{
		lld temp=0;
		for(int j=0;j<framesize-i;j++)
			temp+=arr[j]*arr[j+i];
		R[i]=temp;
	}
}

void calculateAi(lld R[],lld A[])
{
	lld E[p+1],k[p+1],alpha[p+1][p+1];
	E[0]=R[0];
	lld summation;
	for(int i=1;i<=p;i++)
	{
		summation=0;
		for(int j=1;j<=i-1;j++)
			summation+=alpha[j][i-1]*R[i-j];
		k[i]=(R[i]-summation)/E[i-1];
		alpha[i][i]=k[i];
		for(int j=1;j<=i-1;j++)
		alpha[j][i]=alpha[j][i-1]-k[i]*alpha[i-j][i-1];
		E[i]=(1-(k[i]*k[i]))*E[i-1];
	}
	for(int i=1;i<=12;i++)
		A[i]=alpha[i][p];
}

void calculateCi(lld R[],lld A[],lld C[])
{
	C[0]=log(R[0]*R[0]);
	lld summation;
	for(int i=1;i<=12;i++)
	{
		summation=0;
		for(int k=1;k<=i-1;k++)
			summation+=(k*C[k]*A[i-k])/i;
		C[i]=A[i]+summation;
	}
}

void hammingWindow(lld arr[])
{
	lld w;
	for(int i=0;i<framesize;i++)
	{
		w=0.54-(0.46*cos((2*pii*i)/(framesize-1)));
		arr[i]=arr[i]*w;
	}
}

void raisedsinewindow(lld arr[])
{
	lld w;
	for(int i=1;i<=p;i++)
	{
		w=1+((p/2)*sin(pii*i/p));
		arr[i]=w*arr[i];
	}
}

lld Energy(lld arr[])
{
	lld e=0;
	for(int i=0;i<sizeof(arr);i++)
		e+=arr[i]*arr[i];
	return e;
}


void create_universe()
{
	char filename[40];
	lld arr[40000],frame_arr[framesize+1];
	int index,z;
	lld count;
	lld R[p+1],A[p+1],C[p+1];
	FILE* fptr=fopen("universe.csv","w");
	if(fptr==NULL)
	{
		printf("Unable to open reference file");
		return;
	}
	printf("\nCreating Universe Vector");
	for(int i=0;i<=9;i++)
	{
		for(int j=1;j<=sample_count;j++)
		{
			sprintf(filename,"training_sample/234101066_E_%d_%d.txt",i,j);
			FILE *fp=fopen(filename,"r");
			if(fp==NULL)
			{
				printf("unable to access training sample\n");
				return;
			}
			int k=0;
			while(fscanf(fp,"%Lf",&arr[k])!=EOF)
				k++;
			fclose(fp);
			dc_shift(arr);
			normalisation(arr);
			for(int x=0;x+framesize<k;x+=80)
			{
				index=0;
				for(int y=x;y<x+framesize;y++)
				{
					frame_arr[index]=arr[y];
					index++;
				}
				hammingWindow(frame_arr);
				calculateRi(frame_arr,R);
				calculateAi(R,A);
				calculateCi(R,A,C);
				raisedsinewindow(C);
				//dumping Ci values
				for(z=1;z<=p;z++)
					fprintf(fptr,"%Lf,",C[z]);
				fprintf(fptr,"\n");
			}
			printf(".");
		}
	}
	fclose(fptr);
	printf("\nUniverse Creation Successfull!\n");
}