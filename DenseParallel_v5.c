#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include <stdint.h>

#define MASTER 0
#define FROM_MASTER 1
#define FROM_SLAVE 9

#define ROWS 10000
#define COLUMNS 10000

double** matrixA = NULL;
double** matrixB = NULL;
double** matrixC = NULL;
double** matrixR = NULL;

int rank;
int size;
double start_time;
double end_time;


int generateRandomNumber()
{
	return rand()%10 + 1;
	//return 1;
}

void initMatrixDense(double** a, double** b)
{
	int i,j;
	for(i=0; i<ROWS; i++)
	{
	    	for(j=0; j<COLUMNS; j++)
		{
		    	a[i][j] = generateRandomNumber();
			b[i][j] = generateRandomNumber();


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
	matrixC = memory_allocate_2d_contiguous();
	matrixR = memory_allocate_2d_contiguous();
}

void MatrixMultiplyDense(double **matrixA, double **matrixB, int from, int to)
{
	int i,j,k;
	double temp;

	for(i=from; i<to; i++)
	{
		for(j=0; j<COLUMNS; j++)
		{
			temp=0;
			for(k=0;k < ROWS-7; k+=8)
			{
				temp += matrixA[i][k] * matrixB[k][j];
				temp += matrixA[i][k+1] * matrixB[k+1][j];
				temp += matrixA[i][k+2] * matrixB[k+2][j];
				temp += matrixA[i][k+3] * matrixB[k+3][j];
				temp += matrixA[i][k+4] * matrixB[k+4][j];
				temp += matrixA[i][k+5] * matrixB[k+5][j];
				temp += matrixA[i][k+6] * matrixB[k+6][j];
				temp += matrixA[i][k+7] * matrixB[k+7][j];
			}
			for(;k<ROWS;k++)
			{
				temp += matrixA[i][k] * matrixB[k][j];
			}
			matrixC[i][j]=temp;
		}
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

		initMatrixDense(matrixA, matrixB);

		//printMatrix(matrixA, ROWS, COLUMNS);
		//printMatrix(matrixB, ROWS, COLUMNS);
		printf("\nDense Parallel : The size of the process is : %d \n ",size);
			
		from = (i * ROWS)/size;
		to  = ((i+1) *ROWS)/size;
	        MatrixMultiplyDense(matrixA,matrixB,from,to);
		end_time = MPI_Wtime();

		//printMatrix(matrixC,ROWS,COLUMNS);
		printf("\n\n Running Time : %lf \n\n ", end_time - start_time);
		MPI_Finalize();
	}

	if(size > 1)
	{
		// Split the work to slaves
		if(rank == 0)
		{
			printf("\nDense Parallel : The size of the process is : %d \n ",size);
			start_time = MPI_Wtime();

			initMatrixDense(matrixA, matrixB);

			//printMatrix(matrixA, ROWS, COLUMNS);
			//printMatrix(matrixB, ROWS, COLUMNS);
			//printf("\n\n The size of the process is : %d \n\n ",size);
			
			for(int i=1; i<size; i++)
			{

                        	MPI_Send(&(matrixA[0][0]), ROWS*COLUMNS, MPI_DOUBLE, i, FROM_MASTER, MPI_COMM_WORLD);
                        	MPI_Send(&(matrixB[0][0]), ROWS*COLUMNS, MPI_DOUBLE, i, FROM_MASTER + 1, MPI_COMM_WORLD);
			}
		}

		if(rank > 0)
		{

			from = ((rank-1) *ROWS)/(size-1);
			to  = (rank *ROWS)/(size-1);
			count_RC = to - from;
		
			//MPI_Recv(&from, 1, MPI_INT, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
			//MPI_Recv(&to, 1, MPI_INT, 0, FROM_MASTER, MPI_COMM_WORLD, &status);

                	MPI_Recv(&(matrixA[0][0]), ROWS * COLUMNS, MPI_DOUBLE, 0, FROM_MASTER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                	MPI_Recv(&(matrixB[0][0]), ROWS * COLUMNS, MPI_DOUBLE, 0, FROM_MASTER + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


	        	MatrixMultiplyDense(matrixA,matrixB,from,to);
	               	MPI_Send(&count_RC, 1, MPI_INT, 0, FROM_SLAVE, MPI_COMM_WORLD);
	               	MPI_Send(&from, 1, MPI_INT, 0, FROM_SLAVE + 1, MPI_COMM_WORLD);
                	MPI_Send(&(matrixC[from][0]), (count_RC * COLUMNS), MPI_DOUBLE, 0, FROM_SLAVE + 2, MPI_COMM_WORLD);

	               	MPI_Send(&to, 1, MPI_INT, 0, FROM_SLAVE + 3, MPI_COMM_WORLD);

		}

		if(rank == 0)
		{

			int r,count=0;

			for(int i=1; i<size; i++)
			{
			
                		MPI_Recv(&count_RC, 1, MPI_INT, i, FROM_SLAVE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                		MPI_Recv(&from, 1, MPI_INT, i, FROM_SLAVE + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                		MPI_Recv(&(matrixC[from][0]), (count_RC * COLUMNS), MPI_DOUBLE, i, FROM_SLAVE +2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                		MPI_Recv(&to, 1, MPI_INT, i, FROM_SLAVE + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			}


			end_time = MPI_Wtime();
			//printMatrix(matrixC,ROWS,COLUMNS);
			printf("\n\n Running Time : %lf \n\n ", end_time - start_time);

		}


		int p = MPI_Finalize();

	}
		
	return 0;
}





