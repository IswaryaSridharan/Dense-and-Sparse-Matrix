#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include <stdint.h>

#define MASTER 0
#define FROM_MASTER 1
#define FROM_SLAVE 9

#define ROWS 1375
#define COLUMNS 1375

double** matrixA = NULL;
double** matrixB = NULL;
double matrixR[8607][3];
double *XA = NULL;
double *XIA = NULL;
double *XJA = NULL;
double *YA = NULL;
double *YIA = NULL;
double *YJA = NULL;
double *ZA = NULL;
double *ZIA = NULL;
double *ZJA = NULL;
double *TA = NULL;
double *TIA = NULL;
double *TJA = NULL;

double *TB = NULL;
double *TIB = NULL;
double *TJB = NULL;



double *MA = NULL;
double *CA = NULL;
int Xlen=0,Ylen=0,Zlen = 0;
int rank;
int size;
double start_time;
double end_time;
FILE *file;
//int ROWS,COLUMNS;

int generateRandomNumber()
{
	return rand()%10 + 1;
	//return 1;
}

void initMatrixSparse(double** a, double** b)
{
	int i,j;
	for(i=0; i<ROWS; i++)
	{
	    	for(j=0; j<COLUMNS; j++)
		{
			if(i == j)
			{
		    		a[i][j] = generateRandomNumber();
				b[i][j] = generateRandomNumber();
			}
			else
			{
		    		a[i][j] = 0;
				b[i][j] = 0;
			}

	    	}
      	}
}

void printMatrix(double** mat, int rows, int cols)
{
	int i,j;
	for(i = 0; i < rows; ++i) 
	{
		for (int j = 0; j < cols; ++j) 
		{
			printf(" %lf ", mat[i][j]);
		}
		printf("\n");
	}
}


void readMatrix()
{

	int rows = 8607,columns = 3;

	if(file == NULL)
	{
		printf(" File doesn't exist ");
		exit(0);
	}
	else
	{

		int i,j;
		for(i=0;i<rows; i++)
		{
			for (j=0; j <columns; j++) {
				fscanf(file,"%lf",&matrixR[i][j]);
				//printf("%lf \t",matrixR[i][j]);
			}
		}
		fclose(file);
		//ROWS = matrixR[0][0];
		//COLUMNS = matrixR[0][1];
		int m,n; double v;
		for(i=1; i<rows; i++)
		{
			m = matrixR[i][0];
			n = matrixR[i][1];
			v = matrixR[i][2];
			matrixA[m][n] = v; matrixB[m][n] = v;
		}
	//	printMatrix(matrixA,ROWS,COLUMNS);
	//	printMatrix(matrixB,ROWS,COLUMNS);

	

		
	}
}


double **memory_allocate_2d_contiguous() {
    int i;
    double *value = (double *)malloc(ROWS * COLUMNS * sizeof(double));
    double **a = (double **)malloc(ROWS * sizeof(double*));
    for (i=0; i<ROWS; i++)
    {
        a[i] = &(value[COLUMNS*i]);
    }
    return a;
}

void memoryAllocate()
{
	matrixA = memory_allocate_2d_contiguous();
	matrixB = memory_allocate_2d_contiguous();

	XA = (double *)malloc(ROWS*COLUMNS*sizeof(double));
	XIA = (double *)malloc((ROWS*(COLUMNS + 1)) *sizeof(double));
	XJA = (double *)malloc(ROWS*COLUMNS*sizeof(double));
	YA = (double *)malloc(ROWS*COLUMNS*sizeof(double));
	YIA = (double *)malloc((ROWS*(COLUMNS +1)) * sizeof(double));
	YJA = (double *)malloc(ROWS*COLUMNS*sizeof(double));
	ZA = (double *)malloc(ROWS*COLUMNS*sizeof(double));
	ZIA = (double *)malloc((ROWS*(COLUMNS +1)) * sizeof(double));
	ZJA = (double *)malloc(ROWS*COLUMNS*sizeof(double));
	TA = (double *)malloc(ROWS*COLUMNS*sizeof(double));
	TIA = (double *)malloc((ROWS*(COLUMNS +1)) * sizeof(double));
	TJA = (double *)malloc(ROWS*COLUMNS*sizeof(double));
	TB = (double *)malloc(ROWS*COLUMNS*sizeof(double));
	TIB = (double *)malloc((ROWS*(COLUMNS +1)) * sizeof(double));
	TJB = (double *)malloc(ROWS*COLUMNS*sizeof(double));
	MA = (double *)malloc((ROWS*(COLUMNS +1))*sizeof(double));
	CA = (double *)malloc((ROWS*(COLUMNS + 1))*sizeof(double));
	file = fopen("./nnc1374.mtx","r");


}

void free_mem(double *a)
{
	if (a!= NULL)free(a);
}

void free_2d_mem(double **a)
{
	int i;
	for (i=0; i<ROWS; i++)
	{
        	free_mem(a[i]);
	}
	if (a != NULL)free(a);
}

void freeMemory()
{
	free_2d_mem(matrixA);
	free_2d_mem(matrixB);
	free_mem(XA);
	free_mem(XIA);
	free_mem(XJA);
	free_mem(YA);
	free_mem(YIA);
	free_mem(YJA);
	free_mem(ZA);
	free_mem(ZIA);
	free_mem(ZJA);
	free_mem(TA);
	free_mem(TIA);
	free_mem(TJA);
	free_mem(TB);
	free_mem(TIB);
	free_mem(TJB);
	free_mem(MA);
	free_mem(CA);
}

void SparseMultiply(double *XA, double *XIA, double *XJA, double *YA, double *YIA, double *YJA, int row_count)
{
	int i,j,k,Aval,Bval,z=0;
	int Aindex = 0,Bindex = 0,Cindex = 0;
	int Arow,Acol,Brow,Bcol,rescol;

	// Loop through each rows
	for(i = 1; i<row_count + 1; i++)
	{
		// Initialize result matrix row count array
		TIA[i] = 0;

		// Find number of elements in row i
		Aval = XIA[i] - XIA[i-1];
	
		// Loop through each non zero elements in row i
		for(j=0; j<Aval; j++)
		{
			MA[j] = XA[Aindex]; //Stores the value of XA
			CA[j] = XJA[Aindex]; //Stores the column number for each values in row i

			Aindex++;
		}

		// calculate each element in column rescol of the result matrix row i
		for (rescol = 1; rescol <= COLUMNS; rescol++)
		{		
			int sum = 0;
			for(j=0; j<Aval; j++)
			{
				int x=CA[j];
				int Yrowcnt = YIA[x] - YIA[x-1];
				if (Yrowcnt > 0)
				{
					for(k=0; k<Yrowcnt; k++)
					{
						Bindex = YIA[x-1];
						if (YJA[Bindex + k] ==  rescol) sum+=YA[Bindex + k] * MA[j];
					}
				}
			}
			if (sum > 0)
			{
				TA[Cindex] = sum;
				TJA[Cindex] = rescol;
				Cindex++;
			}
		}
		TIA[i-1] = Cindex;
	}
	Zlen = Cindex;
}

void MatrixMultiplySparse(double **matrixA, double **matrixB, int from, int to)
{
	int i, j, m=0, n=0,a=0,b=0,countA=0,countB=0,x=0,y=0;
	XIA[a] = 0;
	YIA[b] = 0;

	for(i=from; i<to; i++)
	{
		for(j=0;j<COLUMNS; j++)
		{
			if(matrixA[i][j] != 0)
			{
				XA[m] = matrixA[i][j];
				m++; countA++;
				XJA[x] = j+1; x++;
			}
		}
		a++;
		XIA[a] = countA;

	}
	for(i=0; i<ROWS; i++)
	{
		for(j=0; j<COLUMNS; j++)
		{
			if(matrixB[i][j] != 0)
			{
				YA[n] = matrixB[i][j];
				n++; countB++;
				YJA[y] = j+1; y++;
			}
		}
		b++;
		YIA[b] = countB;

	}
	Xlen = m;
	Ylen = n;

	SparseMultiply(XA,XIA,XJA,YA,YIA,YJA, to-from);
}	



void printEntries(int lastValue, int Xl)
{

	int i;
	printf("\n ************** XA ***************** \n");		
	for(i=0; i<Xl; i++)
	{
		printf("%lf \t",TB[i]);
	}

	printf(" \n ************** XIA ***************** \n");		
	for(i=0; i<ROWS+1; i++)
	{
		printf("%lf \t",TIB[i]);
	}

	printf(" \n************** XJA ***************** \n");		
	for(i=0; i<Xl; i++)
	{
		printf("%lf \t",TJB[i]);
	}

	printf(" \n************** YA ***************** \n");		
	for(i=0; i<Ylen; i++)
	{
		printf("%lf \t",YA[i]);
	}

	printf(" \n************** YIA ***************** \n");		
	for(i=0; i<ROWS+1; i++)
	{
		printf("%lf \t",YIA[i]);
	}

	printf(" \n************** YJA ***************** \n");		
	for(i=0; i<Ylen; i++)
	{
		printf("%lf \t",YJA[i]);
	}

	printf(" \n************** ZA ***************** \n");		
	for(i=0; i<lastValue; i++)
	{
		printf("%lf \t",ZA[i]);
	}

	printf(" \n************** ZIA ***************** \n");		
	for(i=0; i<ROWS+1; i++)
	{
		printf("%lf \t",ZIA[i]);
	}

	printf(" \n************** ZJA ***************** \n");		
	for(i=0; i<lastValue; i++)
	{
		printf("%lf \t",ZJA[i]);
	}
}

void printEntriesForOne()
{

	int i;
	printf("\n ************** XA ***************** \n");		
	for(i=0; i<Xlen; i++)
	{
		printf("%lf \t",XA[i]);
	}

	printf(" \n ************** XIA ***************** \n");		
	for(i=0; i<ROWS+1; i++)
	{
		printf("%lf \t",XIA[i]);
	}

	printf(" \n************** XJA ***************** \n");		
	for(i=0; i<Xlen; i++)
	{
		printf("%lf \t",XJA[i]);
	}

	printf(" \n************** YA ***************** \n");		
	for(i=0; i<Ylen; i++)
	{
		printf("%lf \t",YA[i]);
	}

	printf(" \n************** YIA ***************** \n");		
	for(i=0; i<ROWS+1; i++)
	{
		printf("%lf \t",YIA[i]);
	}

	printf(" \n************** YJA ***************** \n");		
	for(i=0; i<Ylen; i++)
	{
		printf("%lf \t",YJA[i]);
	}

	printf(" \n************** ZA ***************** \n");		
	for(i=0; i<Zlen; i++)
	{
		printf("%lf \t",TA[i]);
	}


	printf(" \n************** ZIA ***************** \n");		
	for(i=0; i<ROWS+1; i++)
	{
		printf("%lf \t",TIA[i-1]);
	}

	printf(" \n************** ZJA ***************** \n");		
	for(i=0; i<Zlen; i++)
	{
		printf("%lf \t",TJA[i]);
	}
}

int main(int argc, char *argv[])
{
	int from;
	int to;
	int count_RC;

	MPI_Init(&argc, &argv); //initialize MPI operations 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get the rank
        MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processes

	memoryAllocate();

	//printf("\n\n Rank is : %d \n\n ",rank);
	//printf("\n\n Size is : %d \n\n ",size);
	if(size == 1)
	{
		int i=0;
		start_time = MPI_Wtime();

		//initMatrixSparse(matrixA, matrixB);
		readMatrix();

		//printMatrix(matrixA, ROWS, COLUMNS);
		//printMatrix(matrixB, ROWS, COLUMNS);
		printf("\nMatrixSparse 1374 : The size of the process is : %d \n ",size);
			
		from = (i * ROWS)/size;
		to  = ((i+1) *ROWS)/size;
	        MatrixMultiplySparse(matrixA,matrixB,from,to);
		end_time = MPI_Wtime();
		//printEntriesForOne();
		printf("\n\n Running Time : %lf \n\n ", end_time - start_time);
		MPI_Finalize();
	}

	if(size > 1)
	{

		// Split the work to slaves
		if(rank == 0)
		{
			printf("\n MatrixSparse 1374 : The size of the process is : %d \n ",size);
			start_time = MPI_Wtime();


		//	initMatrixSparse(matrixA, matrixB);
			readMatrix();


			//printMatrix(matrixA, ROWS, COLUMNS);
			//printMatrix(matrixB, ROWS, COLUMNS);
			//printf("\n\n The size of the process is : %d \n\n ",size);
			
			for(int i=1; i<size; i++)
			{
				//from = (i *ROWS)/size;
				//to  = ((i+1) *ROWS)/size;
				//count_RC = to - from;
                        	//MPI_Isend(&from, 1, MPI_INT, i, FROM_MASTER, MPI_COMM_WORLD, &request);
                        	//MPI_Isend(&to, 1, MPI_INT, i, FROM_MASTER +1, MPI_COMM_WORLD, &request);

                        	MPI_Send(&(matrixA[0][0]), ROWS*COLUMNS, MPI_DOUBLE, i, FROM_MASTER, MPI_COMM_WORLD);
                        	MPI_Send(&(matrixB[0][0]), ROWS*COLUMNS, MPI_DOUBLE, i, FROM_MASTER + 1, MPI_COMM_WORLD);
			}
		}

		if(rank > 0)
		{

			from = ((rank-1) *ROWS)/(size-1);
			to  = ((rank) *ROWS)/(size-1);
			count_RC = to - from;
		
			//MPI_Recv(&from, 1, MPI_INT, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
			//MPI_Recv(&to, 1, MPI_INT, 0, FROM_MASTER, MPI_COMM_WORLD, &status);

                	MPI_Recv(&(matrixA[0][0]), ROWS * COLUMNS, MPI_DOUBLE, 0, FROM_MASTER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                	MPI_Recv(&(matrixB[0][0]), ROWS * COLUMNS, MPI_DOUBLE, 0, FROM_MASTER + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


	        	MatrixMultiplySparse(matrixA,matrixB,from,to);



                	//MPI_Isend(&from, 1, MPI_INT, 0, FROM_SLAVE, MPI_COMM_WORLD, &request);
                	//MPI_Isend(&to, 1, MPI_INT, 0, FROM_SLAVE+1, MPI_COMM_WORLD, &request);

                	MPI_Send(&Zlen, 1, MPI_INT, 0, FROM_SLAVE, MPI_COMM_WORLD);
                	MPI_Send(&count_RC, 1, MPI_INT, 0, FROM_SLAVE + 1, MPI_COMM_WORLD);

                	MPI_Send(&(TA[0]), Zlen, MPI_DOUBLE, 0, FROM_SLAVE + 2, MPI_COMM_WORLD);
                	MPI_Send(&(TIA[0]), count_RC, MPI_DOUBLE, 0, FROM_SLAVE + 3, MPI_COMM_WORLD);
                	MPI_Send(&(TJA[0]), Zlen, MPI_DOUBLE, 0, FROM_SLAVE + 4, MPI_COMM_WORLD);

                	MPI_Send(&Xlen, 1, MPI_INT, 0, FROM_SLAVE + 5, MPI_COMM_WORLD);
                	MPI_Send(&(XA[0]), Xlen, MPI_DOUBLE, 0, FROM_SLAVE + 6, MPI_COMM_WORLD);
                	MPI_Send(&(XIA[0]), count_RC + 1, MPI_DOUBLE, 0, FROM_SLAVE + 7, MPI_COMM_WORLD);
                	MPI_Send(&(XJA[0]), Xlen, MPI_DOUBLE, 0, FROM_SLAVE + 8, MPI_COMM_WORLD);


			if(rank == 1)
			{
                		MPI_Send(&Ylen, 1, MPI_INT, 0, FROM_SLAVE + 9, MPI_COMM_WORLD);
                		MPI_Send(&(YA[0]), Ylen, MPI_DOUBLE, 0, FROM_SLAVE + 10, MPI_COMM_WORLD);
                		MPI_Send(&(YIA[0]), ROWS+1, MPI_DOUBLE, 0, FROM_SLAVE + 11, MPI_COMM_WORLD);
                		MPI_Send(&(YJA[0]), Ylen, MPI_DOUBLE, 0, FROM_SLAVE + 12, MPI_COMM_WORLD);

			}



		}

		
		if(rank == 0)
		{
			ZIA[0] =0;TIB[0] = 0;
			int k,Xl=0,lastValue=0;
			int ZrowCount, addZrows=1;
			int XrowCount, addXrows=1;
			int x=0,y=0;
			for(int i=1; i<size; i++)
			{
				//MPI_Recv(&from, 1, MPI_INT, i, FROM_SLAVE, MPI_COMM_WORLD, &status);
        	    		//MPI_Recv(&to, 1, MPI_INT, i, FROM_SLAVE + 1, MPI_COMM_WORLD, &status);

        	    		MPI_Recv(&Zlen, 1, MPI_INT, i, FROM_SLAVE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	    		MPI_Recv(&count_RC, 1, MPI_INT, i, FROM_SLAVE + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		


        	    		MPI_Recv(&(TA[0]), Zlen, MPI_DOUBLE, i, FROM_SLAVE + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	    		MPI_Recv(&(TIA[0]), count_RC, MPI_DOUBLE, i, FROM_SLAVE + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	    		MPI_Recv(&(TJA[0]), Zlen, MPI_DOUBLE, i, FROM_SLAVE + 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        	    		MPI_Recv(&Xlen, 1, MPI_INT, i, FROM_SLAVE + 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	    		MPI_Recv(&(XA[0]), Xlen, MPI_DOUBLE, i, FROM_SLAVE + 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	    		MPI_Recv(&(XIA[0]), count_RC +1, MPI_DOUBLE, i, FROM_SLAVE + 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	    		MPI_Recv(&(XJA[0]), Xlen, MPI_DOUBLE, i, FROM_SLAVE + 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				if(i == 1)
				{
        	    			MPI_Recv(&Ylen, 1, MPI_INT, i, FROM_SLAVE + 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	    			MPI_Recv(&(YA[0]), Ylen, MPI_DOUBLE, i, FROM_SLAVE + 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	    			MPI_Recv(&(YIA[0]), ROWS+1, MPI_DOUBLE, i, FROM_SLAVE + 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	    			MPI_Recv(&(YJA[0]), Ylen, MPI_DOUBLE, i, FROM_SLAVE + 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				
				for(k=lastValue; k<(lastValue+Zlen); k++)
				{
					ZA[k] = TA[k-lastValue];
					ZJA[k] = TJA[k-lastValue];

				}

				for(k=Xl; k<(Xl+Xlen); k++)
				{
					TB[k] = XA[k-Xl];
					TJB[k] = XJA[k-Xl];

				}
				lastValue += Zlen;
				Xl += Xlen;

				for(ZrowCount=addZrows; ZrowCount <=(addZrows + count_RC); ZrowCount++)
				{
					ZIA[ZrowCount] = TIA[ZrowCount-addZrows] + x;
				}
				x += TIA[count_RC -1];
				addZrows += count_RC;

				for(XrowCount=addXrows; XrowCount <=(addXrows + count_RC); XrowCount++)
				{
					TIB[XrowCount] = XIA[XrowCount-addXrows +1] + y;
				}
				y += XIA[count_RC];
				addXrows += count_RC;

			}
			end_time = MPI_Wtime();
			//printEntries(lastValue,Xl);
			printf("\n\n Running Time : %lf \n\n ", end_time - start_time);
		}
//		freeMemory();
		int p = MPI_Finalize();

	}
		
	return 0;
}
