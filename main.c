#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

// Make an initial square matrix that takes a size and an initial value for the 
// left and top values
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

// Print out the matrix for debugging purposes
void print_matrix(double** matrix, int size){
	for(int y = 0; y < size; y++){
		for(int x = 0; x < size; x++){
			printf("%f ", matrix[y][x]);
		}
		printf("\n");
	}
}

// Print out an array of doubles for debugging purposes
void print_array(double* array, int size){
	for(int x = 0; x < size; x++){
		printf("%f, ", array[x]);
		if ((x+1) % 5 == 0){
			printf("\n");
		}
	}
	printf("\n");
}

// Print out an array of integers for debugging purposes
void print_int_array(int* array, int size){
	for(int x = 0; x < size; x++){
		printf("%d\t", array[x]);
	}
	printf("\n");
}

// Convert a 2D matrix into a problem array that takes all the calculations
// needed into one long array
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
	
	// Default values
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

	//char name[MPI_MAX_PROCESSOR_NAME];

	rc = MPI_Init(&argc, &argv);

	//MPI_Bcast(p_array, ((matrix_size-2)*(matrix_size-2)*5), MPI_DOUBLE, 0, MPI_COMM_WORLD);



	// matrix = make_matrix(matrix_size, default_value);

	// 	// Create our problem array, this will be split via scatter to share the
	// 	// calculations.
	// 	p_array = convert_matrix(matrix, matrix_size);

	if (rc != MPI_SUCCESS){
		printf("Error starting MPI program\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	// Init the timer
	double t1, t2;
	// Sync all processors before starting the timer, because each processor has
	// their own version of the timer, so we decide the total time based on the
	// slowest process.
	MPI_Barrier(MPI_COMM_WORLD); 
	t1 = MPI_Wtime();

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	// Check that that correct number of processors has been initialised
	if (myrank == 0){
		printf("main reports %d procs \n", nproc);
	}

	double** matrix;
	double* p_array;

	// We only need to populate this for rank 0
	if (myrank == 0){
		
		matrix = make_matrix(matrix_size, default_value);

		// Create our problem array, this will be split via scatter to share the
		// calculations.
		p_array = convert_matrix(matrix, matrix_size);


	} else {
		p_array = NULL;
	}

	// Here we need arrays of ints that store specific information for scatterv
	// and gatherv

	// Scatterv
	// No. elements to send to each processor
	int* counts = (int*)malloc(sizeof(int)*nproc);

	// The starting indecies of each proportion
	int* displace = (int*)malloc(sizeof(int)*nproc);

	// The number of elements each process will recieve
	int* recieve_counts = (int*)malloc(sizeof(int)*nproc);
	// The starting point of each proportion in the recieve buffer
	int* recieve_displace = (int*)malloc(sizeof(int)*nproc);

	int prop, rem;
	int offset = 0;
	int recieve_offset = 0;

	// Number of elements to give to each process
	prop = ((matrix_size-2)*(matrix_size-2))/nproc;
	// The remainder to share equally among the processes
	rem = ((matrix_size-2)*(matrix_size-2))%nproc;
	// Populate the arrays to be passed into MPI functions
	for (int x = 0; x < nproc ; x++){
		// Proportion*5 because we have 5 elements for each point in the matrix 
		// that needs to be calculated
		counts[x] = prop*5;

		displace[x] = offset;
		// If the processor rank is less than the remainder, then add an
		// additional point to calculate
		if (x < rem){
			counts[x] += 5;
		}
		offset += counts[x];

		// Apply same logic for the recieve counts
		recieve_counts[x] = prop;
		recieve_displace[x] = recieve_offset;
		
		if (x < rem){
			recieve_counts[x] += 1;	
		}
		recieve_offset += recieve_counts[x];
	}

	// The recieve amount for a particular process, will vary for each process
	
	int recieve_count = recieve_counts[myrank]*5;

	// Uncomment for extra information about how data is split
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

	// Receive buffer
	double* recieve = malloc(sizeof(double)*recieve_count);
	// Status of all processes, we can use this to check if all the processes
	// are done calculating up to the precision.
	int* proc_status = (int*)calloc(nproc , sizeof(int));
	int status; // The status of a particular process
	int processing = 1; // The processing state of a process

	while (processing){
		// Set status to 1 indicates that we have reached the goal state, this
		// will be disproved when we check the calculations against the previous
		status = 1;

		// Scatter the problem array equally among the different processers.
		MPI_Scatterv(p_array, counts, displace, MPI_DOUBLE, recieve, 
			recieve_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// Init the array that will take all the results so that it can be
		// gathered
		double* calculated = malloc(sizeof(double)*(recieve_count/5));

		// Calculate each point using the elements of what we have recieved, the
		// first four elements will be avereged, we will then use the last elem
		// to compare to see if it is below the precision
		for (int x = 0; x < recieve_count; x+=5){
			calculated[x/5] = 
				(recieve[x] + recieve[x+1] + recieve[x+2] + recieve[x+3])/4.0;
			//printf("Comparing %f and %f\n", calculated[x/5], recieve[x+4]);

			// Fabs is used to account for the possibility of negative values
			if ((fabs(calculated[x/5]) - fabs(recieve[x+4])) >= precision){
				// We can see that we still require calculating if at least one
				// point doesn't meet the requirements.
				status = 0;
			}
		}

		// The array that will take back all the calculations into one array
		double* root_retrieve = 
			malloc(sizeof(double)*((matrix_size-2)*(matrix_size-2)));

		// Gather all the calculations and place them into root_retrieve
		MPI_Gatherv(calculated, recieve_count/5, MPI_DOUBLE, root_retrieve, 
			recieve_counts, recieve_displace, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// The root will then take those values and put it back into the matrix
		if (myrank == 0){

			int retrieve_count = 0;
			for(int y = 1; y < matrix_size-1; y++){
				for (int x = 1; x < matrix_size-1; x++){
					matrix[y][x] = root_retrieve[retrieve_count++];
				}
			}

			// Uncomment if you want to see each iteration of the matrix
			// printf("Matrix: \n");
			// print_matrix(matrix, matrix_size);
		}

		// Free calculated as that is no longer in use at this point, it will be
		// re-allocated when we continue to process
		free(calculated);
		// Free root_retrieve as we have already made use of the values and put
		// them back into the matrix
		free(root_retrieve);

		// Gather all the statuses for so we can stop processing, this needs to
		// be communicated to each process so that they can all stop
		// individually.
		MPI_Allgather(&status, 1, MPI_INT, proc_status, 1, MPI_INT, 
			MPI_COMM_WORLD);

		// Stop processing
		processing = 0;
		for(int x = 0; x < nproc ; x++){
			// If one processor is still not reached the goal state, then
			// continue processing
			if (!proc_status[x]){
				processing = 1;
				x = nproc ;
			}
		}

		// If we do continue processing, then reconvert the matrix back into a
		// problem array from root
		if (processing && myrank == 0){
			// Check timings of how long it takes to convert a matrix to a
			// problem array
			double c_time1, c_time2;
			c_time1 = MPI_Wtime();
			free(p_array);
			p_array = convert_matrix(matrix, matrix_size);
			c_time2 = MPI_Wtime();
			printf("Conversion time: %f secs\n", (c_time2-c_time1));
		}
	}

	// Sync the processors again to get the final time
	MPI_Barrier(MPI_COMM_WORLD);
	t2 = MPI_Wtime();
	// Output the stats from root
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