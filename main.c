#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

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

void print_int_array(int* array, int size){
	for(int x = 0; x < size; x++){
		printf("%d\t", array[x]);
	}
	printf("\n");
}

double* convert_matrix(double** matrix, int size){
	double* p_array = malloc(sizeof(double)*((size-2)*(size-2))*5);
	int i = 0;
	for(int y = 1; y < size-1; y++){
		for (int x = 1; x < size-1; x++){
			p_array[i++] = matrix[y][x-1]; // Left
			p_array[i++] = matrix[y-1][x]; // Top
			p_array[i++] = matrix[y][x+1]; // Right
			p_array[i++] = matrix[y+1][x]; // Bottom
			p_array[i++] = matrix[y][x]; // Centre
		}
	}

	return p_array;
}

// Find a string in an array of strings.
// Returns the index of the string if it found, -1 if not found.
int str_array_find(char* arr[], int size, char* string){
	for(int x = 0; x < size; x++){
		if(strcmp(arr[x], string) == 0){
			return x;
		}
	}
	return -1;
}

int main(int argc, char** argv){
	

	int matrix_size = 5;
	double default_value = 1.0;
	double precision = 0.01;

	// Check to see if we have any arguments.
	int dim_found = str_array_find(argv, argc, "-d");
	int def_val_found = str_array_find(argv, argc, "-v");
	int prec_found = str_array_find(argv, argc, "-p");

	char *eptr;

	// If we have arguments, start setting up our variables.
	if (dim_found > -1){
		matrix_size = atoi(argv[dim_found+1]);
	}
	if (def_val_found > 1){
		default_value = strtod(argv[def_val_found+1], &eptr);
	}
	if (prec_found > -1){
		precision = strtod(argv[prec_found+1], &eptr);
	}

	int rc, myrank, nproc;

	double** matrix = make_matrix(matrix_size, default_value);
	double* p_array = convert_matrix(matrix, matrix_size);

	//char name[MPI_MAX_PROCESSOR_NAME];

	rc = MPI_Init(&argc, &argv);

	if (rc != MPI_SUCCESS){
		printf("Error starting MPI program\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	double t1, t2;
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	int* counts = (int*)malloc(sizeof(int)*nproc );
	int* strides = (int*)malloc(sizeof(int)*nproc );
	int* displace = (int*)malloc(sizeof(int)*nproc );
	int* recieve_counts = (int*)malloc(sizeof(int)*nproc );
	int* recieve_displace = (int*)malloc(sizeof(int)*nproc );
	int prop, rem;
	int offset = 0;
	int recieve_offset = 0;
	prop = ((matrix_size-2)*(matrix_size-2))/nproc;
	rem = ((matrix_size-2)*(matrix_size-2))%nproc;
	double* recieve;
	for (int x = 0; x < nproc ; x++){
		counts[x] = prop*5;
		strides[x] = prop*5;
		displace[x] = offset;
		offset += strides[x];
		recieve_counts[x] = prop;
		recieve_displace[x] = recieve_offset;
		recieve_offset += prop;
	}

	counts[nproc-1] += rem*5;
	strides[nproc-1] += rem*5;
	recieve_counts[nproc-1] += rem;
	
	int recieve_count = prop*5;
	if (myrank == nproc-1){
		recieve_count += rem*5;
	}

	// if (myrank == 0){
	// 	printf("Counts: ");
	// 	print_int_array(counts, nproc );
	// 	printf("Strides: ");
	// 	print_int_array(strides, nproc );
	// 	printf("Displace: ");
	// 	print_int_array(displace, nproc );
	// 	printf("Recieve Counts: ");
	// 	print_int_array(recieve_counts, nproc );
	// 	printf("Recieve Displace: ");
	// 	print_int_array(recieve_displace, nproc );

	// }

	recieve = malloc(sizeof(double)*recieve_count);
	int* proc_status = (int*)calloc(nproc , sizeof(int));
	int status;
	int processing = 1;

	while (processing){
		
		status = 1;
		MPI_Scatterv(p_array, counts, displace, MPI_DOUBLE, recieve, recieve_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);


		double* calculated = malloc(sizeof(double)*(recieve_count/5));
		for (int x = 0; x < recieve_count; x+=5){
			calculated[x/5] = (recieve[x] + recieve[x+1] + recieve[x+2] + recieve[x+3])/4.0;
			//printf("Comparing %f and %f\n", calculated[x/5], recieve[x+4]);
			if ((calculated[x/5] - recieve[x+4]) >= precision){
				
				status = 0;
			}
		}

		double* root_retrieve = malloc(sizeof(double)*((matrix_size-2)*(matrix_size-2)));

		MPI_Gatherv(calculated, recieve_count/5, MPI_DOUBLE, root_retrieve, recieve_counts, recieve_displace, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (myrank == 0){

			int retrieve_count = 0;
			for(int y = 1; y < matrix_size-1; y++){
				for (int x = 1; x < matrix_size-1; x++){
					matrix[y][x] = root_retrieve[retrieve_count++];
				}
			}

			//printf("Matrix: \n");
			//print_matrix(matrix, matrix_size);
		}
		free(calculated);
		free(root_retrieve);
		free(p_array);

		p_array = convert_matrix(matrix, matrix_size);

		MPI_Allgather(&status, 1, MPI_INT, proc_status, 1, MPI_INT, MPI_COMM_WORLD);
		//printf("Statuses: \n");
		//print_int_array(proc_status, nproc );

		processing = 0;
		for(int x = 0; x < nproc ; x++){
			if (!proc_status[x]){
				processing = 1;
				x = nproc ;
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	t2 = MPI_Wtime();
	if(myrank == 0){
		if (matrix_size < 11){
			print_matrix(matrix, matrix_size);
		}
		
		printf("\nMatrix size: %d\n", matrix_size);
		printf("Initial value: %f\n", default_value);
		printf("Target precision: %f\n", precision);
		printf("Number of processors: %d\n", nproc );
		printf("Processing time: %f seconds\n", (t2-t1));
	}

	MPI_Finalize();

	return 0;
}