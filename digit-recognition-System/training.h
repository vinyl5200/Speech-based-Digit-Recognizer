#define random_size 1000
#define N 5
#define M 32
long double alpha[random_size][N+1];
long double beta[random_size][N+1];
long double gamma[random_size][N+1];
long double delta[random_size][N+1];
long double zai[random_size][N+1][N+1];
int shai[random_size][N+1];
int qstar[random_size];
long double A[N+1][N+1];
long double B[N+1][M+1];
long double pi[N+1];
long double A_comp[N+1][N+1];
long double B_comp[N+1][M+1];
long double pi_comp[N+1];
int O[random_size];
int T;
void load_codebook()
{
	FILE *fp;
	fp=fopen("codebook.txt","r");
	for(int i=1;i<=32;i++)
	{
		for(int j=1;j<=p;j++)
		{
			fscanf(fp,"%Lf ",&codebook[i][j]);
		}
	}
}
int tokhura_distance_index(lld C[])
{
	lld min_dis=INT_MAX,dis;
	int index=1;
	for(int i=1;i<=32;i++)
	{
		dis=0;
		for(int j=1;j<=p;j++)
		dis+=tokhura_weight[j]*pow((C[j]-codebook[i][j]),2);
		if(dis<min_dis)
		{
			min_dis=dis;
			index=i;
		}
	}
	return index;
}
lld forward_procedure()
{
	//Step 1: Initialisation
	for(int i=1;i<=N;i++)
		alpha[1][i]=pi[i]*B[i][O[1]];
	//Step 2: Induction
	for(int t=1;t<=T-1;t++)
	{
		for(int j=1;j<=N;j++)
		{
			long double temp=0;
			for(int i=1;i<=N;i++)
				temp+=alpha[t][i]*A[i][j];
			alpha[t+1][j]=temp*B[j][O[t+1]];
		}
	}
	//Step 3: Termination
	long double probability=0;
	for(int i=1;i<=N;i++)
		probability+=alpha[T][i];
	return probability;
}
void backward_procedure()
{
	//Initialisation
	for(int i=1;i<=N;i++)
		beta[T][i]=1;
	//Induction
	for(int t=T-1;t>=1;t--)
	{
		for(int i=1;i<=N;i++)
		{
			long double temp=0;
			for(int j=1;j<=N;j++)
			{
				temp+=A[i][j]*B[j][O[t+1]]*beta[t+1][j];
			}
			beta[t][i]=temp;
		}
	}
}
void soln_problem2()
{
	for(int i=1;i<=N;i++)
	{
		for(int t=1;t<=T;t++)
		{
			long double temp=0;
			for(int j=1;j<=N;j++)
			{
				temp+=alpha[t][j]*beta[t][j];
			}
			gamma[t][i]=(alpha[t][i]*beta[t][i])/temp;
		}
	}
}
lld viterbi_algo()
{
	//Step 1:Initialization
	for(int i=1;i<=N;i++)
	{
		delta[1][i]=pi[i]*B[i][O[1]];
		shai[1][i]=0;
	}
	//Step 2: Recursion
	long double max,arg_max,temp;
	int state=-1;
	for(int t=2;t<=T;t++)
	{
		for(int j=1;j<=N;j++)
		{
			max=0;
			arg_max=0;
			for(int i=1;i<=N;i++)
			{
				temp=delta[t-1][i]*A[i][j];
				if(temp>arg_max)
				{
					arg_max=temp;
					state=i;
				}
				temp*=B[j][O[t]];
				if(temp>max)
					max=temp;
			}
			delta[t][j]=max;
			shai[t][j]=state;
		}
	}
	//Step 3:Termination
	long double Pstar;
	max=0,arg_max=0;
	for(int i=1;i<=N;i++)
	{
		if(delta[T][i]>max)
			max=delta[T][i];
		if(delta[T][i]>arg_max)
		{
			arg_max=delta[T][i];
			state=i;
		}
	}
	Pstar=max;
	qstar[T]=state;
	return Pstar;
}
void re_estimation()
{
	long double numerator,denominator;
	//Re-estimation of pi array
	for(int i=1;i<=N;i++)
	{
		pi[i]=gamma[1][i];
	}
	//Re-estimation of A matrix
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			numerator=0,denominator=0;
			for(int t=1;t<=T-1;t++)
			{
				numerator+=zai[t][i][j];
				denominator+=gamma[t][i];
			}
			A[i][j]=numerator/denominator;
		}
	}
	//Re-estimation of B matrix
	for(int j=1;j<=N;j++)
	{
		denominator=0;
		for(int t=1;t<=T;t++)
		{
			denominator+=gamma[t][j];
		}
		for(int k=1;k<=M;k++)
		{
			numerator=0;
			for(int t=1;t<=T;t++)
			{
				if(O[t]==k)
					numerator+=gamma[t][j];
			}
			B[j][k]=numerator/denominator;
		}
	}
}
void maintain_stochastic()
{
	lld row_sum=0,max=0;
	int index;
	//pi values
	for(int i=1;i<=N;i++)
	{
		if(pi[i]>max)
		{
			max=pi[i];
			index=i;
		}
		row_sum+=pi[i];
	}
	if(row_sum<1)
		pi[index]+=1-row_sum;
	else if(row_sum>1)
		pi[index]-=row_sum-1;
	//A values
	for(int i=1;i<=N;i++)
	{
		row_sum=0,max=0;
		for(int j=1;j<=N;j++)
		{
			if(A[i][j]>max)
			{
				max=A[i][j];
				index=j;
			}
			row_sum+=A[i][j];		
		}
		if(row_sum<1)
			A[i][index]+=1-row_sum;
		else if(row_sum>1)
			A[i][index]-=row_sum-1;
	}
	//B values
	for(int i=1;i<=N;i++)
	{
		row_sum=0,max=0;
		for(int j=1;j<=M;j++)
		{
			if(B[i][j]<(1e-30))
				B[i][j]=1e-30;
			if(B[i][j]>max)
			{
				max=B[i][j];
				index=j;
			}
			row_sum+=B[i][j];		
		}
		if(row_sum<1)
			B[i][index]+=1-row_sum;
		else if(row_sum>1)
			B[i][index]-=row_sum-1;
	}
}
void soln_problem3()
{
	long double numerator,denominator;
	for(int t=1;t<=T-1;t++)
	{
		for(int i=1;i<=N;i++)
		{
			for(int j=1;j<=N;j++)
			{
				numerator=alpha[t][i]*A[i][j]*B[j][O[t+1]]*beta[t+1][j];
				denominator=0;
				for(int x=1;x<=N;x++)
				{
					for(int y=1;y<=N;y++)
					{
						denominator+=alpha[t][x]*A[x][y]*B[y][O[t+1]]*beta[t+1][y];
					}
				}
				zai[t][i][j]=numerator/denominator;
			}
		}
	}
	re_estimation(); 
	maintain_stochastic();
}
void initialize_lambda()
{
	for(int i=1;i<=N;i++)
	{
		if(i==1)
			pi[i]=1.0;
		else
			pi[i]=0.0;
	}
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			if(j==i && j+1<=N)
			{
				A[i][j]=0.8;
				A[i][j+1]=0.2;
			}
			else if(j==i && j+1>N)
				A[i][j]=1;
		}
	}
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=M;j++)
		{
			B[i][j]=1.0/32;
		}
	}
}
void initialize_lambda_comp()
{
	for(int i=1;i<=N;i++)
		pi_comp[i]=0;
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			A_comp[i][j]=0;
		}
	}
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=M;j++)
		{
			B_comp[i][j]=0;
		}
	}
}
void read_lambda_values(int digit)
{
	FILE *fpt;
	char filename[40];
	sprintf(filename,"lambda_values/234101066_E_%d_pi.txt",digit);
	fpt=fopen(filename,"r");
	for(int i=1;i<=N;i++)
		fscanf(fpt,"%lf",&pi[i]);
	fclose(fpt);
	sprintf(filename,"lambda_values/234101066_E_%d_A.txt",digit);
	fpt=fopen(filename,"r");
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
			fscanf(fpt,"%lf",&A[i][j]);
	}
	fclose(fpt);
	sprintf(filename,"lambda_values/234101066_E_%d_B.txt",digit);
	fpt=fopen(filename,"r");
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=M;j++)
			fscanf(fpt,"%lf",&B[i][j]);
	}
	fclose(fpt);
}
void converge_model(int train_count,int digit)
{
	if(train_count==1)
		initialize_lambda();
	else
		read_lambda_values(digit);
	lld old_Pstar,new_Pstar=10;
	int countt=0;
	do
	{
		countt++;
		old_Pstar=new_Pstar;
		forward_procedure();	
		backward_procedure();
		soln_problem2();
		new_Pstar=viterbi_algo();
		soln_problem3();
	}while(countt!=1000 && old_Pstar/new_Pstar!=1.0);
}

void add_lambda_values()
{
	for(int i=1;i<=N;i++)
		pi_comp[i]+=gamma[1][i];
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			A_comp[i][j]+=A[i][j];
		}
	}
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=M;j++)
		{
			B_comp[i][j]+=B[i][j];
		}
	}
}

void takeAvg_and_dumb_lambda_values(int digit)
{
	for(int i=1;i<=N;i++)
		pi_comp[i]/=sample_count;
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
			A_comp[i][j]/=sample_count;
	}
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=M;j++)
			B_comp[i][j]/=sample_count;
	}
	//print in file
	FILE *fpp;
	char filename[40];
	sprintf(filename,"lambda_values/234101066_E_%d_pi.txt",digit);
	fpp=fopen(filename,"w");
	for(int i=1;i<=N;i++)
		fprintf(fpp,"%e ",pi_comp[i]);
	fclose(fpp);
	sprintf(filename,"lambda_values/234101066_E_%d_A.txt",digit);
	fpp=fopen(filename,"w");
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=N;j++)
		{
			fprintf(fpp,"%e ",A_comp[i][j]);
		}
		fprintf(fpp,"\n");
	}
	fclose(fpp);
	sprintf(filename,"lambda_values/234101066_E_%d_B.txt",digit);
	fpp=fopen(filename,"w");
	for(int i=1;i<=N;i++)
	{
		for(int j=1;j<=M;j++)
		{
			fprintf(fpp,"%e ",B_comp[i][j]);
		}
		fprintf(fpp,"\n");
	}
	fclose(fpp);
}
void training()
{
	char filename[40];
	lld arr[40000],frame_arr[framesize+1];
	int index,framecount;
	lld count;
	lld R[p+1],A[p+1],C[p+1];
	load_codebook();
	int train_count=1;
	printf("\n--------------------------------------MODEL TRAINING STARTED---------------------------------------\n");
	while(train_count<=2)
	{
	for(int i=0;i<=9;i++)
	{
		initialize_lambda_comp();
		printf("training_sample/234101066_E_%d\n",i);
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
			framecount=0;
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
				framecount++;
				O[framecount]=tokhura_distance_index(C);
			}
			
			T=framecount;
			converge_model(train_count,i);
			add_lambda_values();
		}
		takeAvg_and_dumb_lambda_values(i);
	}
	printf("\nTraining phase %d completed!\n\n",train_count);
	train_count++;
	}
	printf("\n--------------------------------------MODEL TRAINED---------------------------------------\n");
}