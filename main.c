#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

double** make_matrix(int size, double init_value){
	double** matrix = malloc(sizeof(double*)*size);
	for(int y = 0; y < size; y++){
		matrix[y] = malloc(sizeof(double)*size);
		for (int x = 0; x < size; x++){
			if (y == 0 || x == 0){
				matrix[y][x] = init_value;
			}
		}
	}
	return matrix;
}

void print_matrix(double** matrix, int size){
	for(int y = 0; y < size; y++){
		for(int x = 0; x < size; x++){
			printf("%f ", matrix[y][x]);
		}
		printf("\n");
	}
}

void print_array(double* array, int size){
	for(int x = 0; x < size; x++){
		printf("%f, ", array[x]);
		if ((x+1) % 5 == 0){
			printf("\n");
		}
	}
	printf("\n");
}

double* convert_matrix(double** matrix, int size){
	double* p_array = malloc(sizeof(double)*((size-2)*(size-2))*5);
	int i = 0;
	for(int y = 1; y < size-1; y++){
		for (int x = 1; x < size-1; x++){
			p_array[i] = matrix[y][x-1]; // Left
			p_array[++i] = matrix[y-1][x]; // Top
			p_array[++i] = matrix[y][x+1]; // Right
			p_array[++i] = matrix[y+1][x]; // Bottom
			p_array[++i] = matrix[y][x]; // Centre
		}
	}

	return p_array;
}

int main(int argc, char** argv){
	int rc, myrank, nproc, namelen;

	int matrix_size = 5;
	int n_procs = 4;
	double** matrix = make_matrix(matrix_size, 1.0);
	print_matrix(matrix, matrix_size);
	printf("\n");
	double* p_array = convert_matrix(matrix, matrix_size);
	print_array(p_array, ((matrix_size-2)*(matrix_size-2))*5);
	char name[MPI_MAX_PROCESSOR_NAME];

	rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS){
		printf("Error starting MPI program\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	if (myrank == 0){
		printf("main reports %d procs\n", nproc);
		
	}

	int* counts = (int*)malloc(sizeof(int)*n_procs);
	int* strides = (int*)malloc(sizeof(int)*n_procs);
	int* displace = (int*)malloc(sizeof(int)*n_procs);
	int prop, rem, offset;
	prop = ((matrix_size-2)*(matrix_size-2))/n_procs;
	rem = ((matrix_size-2)*(matrix_size-2))%n_procs;
	double* recieve = malloc(sizeof(double)*((matrix_size-2)*(matrix_size-2))*5);
	for (int x = 0; x < n_procs; x++){
		counts[x] = prop*5;
		strides[x] = prop*5;
		displace[x] = offset;
		offset += strides[x];
	}

	counts[n_procs-1] += rem*5;
	strides[n_procs-1] += rem*5;

	

	MPI_Scatterv(p_array, counts, displace, MPI_DOUBLE, recieve, prop+rem, MPI_DOUBLE, 0, MPI_COMM_WORLD);




	// int p_size = ((matrix_size-1)*(matrix_size-1))/nproc;
	// int rem = ((matrix_size-1)*(matrix_size-1))%nproc;
	// double* p_matrix = malloc(sizeof(double)*p_size);
	// int start_point = (myrank*p_size)+matrix_size+1;
	// for(int x = start_point; x < p_size+start_point; x++){

	// 	double left_value = matrix[x-1];
	// 	double right_value = matrix[x-1];
	// 	double top_value = matrix[x-matrix_size];
	// 	double bottom_value = matrix[x+matrix_size];
	// 	double avg = (left_value+right_value+top_value+bottom_value)/4.0;
	// 	p_matrix[x-start_point] = avg;

	// }

	namelen = MPI_MAX_PROCESSOR_NAME;
	MPI_Get_processor_name(name, &namelen);
	printf("hello world %d from %s\n", myrank, name);

	MPI_Finalize();

	return 0;
}