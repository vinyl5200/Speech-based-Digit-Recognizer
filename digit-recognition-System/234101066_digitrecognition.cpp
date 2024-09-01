// digitRecognition_234101066.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<conio.h>
#include "string.h"
#include "stdio.h"
#include<stdlib.h>
#include<cstdlib>
#include<time.h>
#include<math.h>
#include<windows.h>
#include<iostream>
#pragma comment(lib,"Winmm.lib")
#include "create_universe.h"
#include "create_codebook.h"
#include "training.h"
#include "testing.h"
#include "live_testing.h"
using namespace std;
int _tmain(int argc, _TCHAR* argv[])
{
	int x;
	while(1)
	{
		printf("\nEnter 1:Create Universe \nEnter 2:Create Codebook \nEnter 3:Model Training \nEnter 4:Prerecorded Model Testing  \n0:Exit\n");
		scanf("%d",&x);
		switch(x)
		{
			case 1:
				create_universe();
				break;
			case 2:
				create_codebook();
				break;
			case 3:
				training();
				break;
			case 4:
				testing();
				break;
			case 0:
				return 0;
			default:
				printf("Invalid input\n");
		}
	}
	getch();
	return 0;
}

