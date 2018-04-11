#include "header.h"

int main(int argc,char *argv[])
{
	int i = 0, j = 0, k = 2, iterations = 0, count = 0, chunk_count = 0, stop_inner_loop = FALSE, stop_external_loop = FALSE;
	double x, y, diameter = 0, calc_qm = 0, QM = 0, answer = 0, start, end;
	int N = 0, MAX = 0, LIMIT = 0; 
	point_t *points_arr, *chunk;
	cluster_t *clusters_arr; 

	int  namelen, numprocs, myid;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);	

	MPI_Get_processor_name(processor_name,&namelen);

	MPI_Status status;

	// create new MPI types:
	MPI_Datatype MPI_POINT = mpiPoint();
	MPI_Datatype MPI_CLUSTER = mpiCluster();


	if(myid == MASTER)
	{
		points_arr = (point_t*)calloc(MAX_NUM_POINTS, sizeof(point_t)); // before reading amount of points from file (N) - allocate array to maximum.

		readFromFile(&N, &MAX, &LIMIT, &QM, points_arr);
		printf("N = %d, MAX = %d, LIMIT = %d, QM = %.3f\n", N, MAX, LIMIT, QM);fflush(stdout);

		clusters_arr = (cluster_t*)calloc(MAX, sizeof(cluster_t));
		chunk = (point_t*)calloc(CHUNK_SIZE, sizeof(point_t));

		start = MPI_Wtime(); // get start time

		// k-means algorithem
		for(k = 2; k <= MAX; k++)  // stopping conditions (master's external loop) - k = max
		{
			calc_qm = 0;

			// determine clusters for the first time
			initializeClusters(k, N, clusters_arr, points_arr);

			for(iterations = 0; iterations < LIMIT; iterations++) //stopping conditions (master's inner loop) - max num of iterations.
			{
				count = groupPointsToClusters(N, k, clusters_arr, points_arr); 

				calculateNewClustersCenters(N, k, clusters_arr, points_arr);

				if(count == N) // stopping conditions (master's inner loop) - No point has changed it's cluster.
				{
					break;
				}
			}

			// find two points in cluster with the largest distance from each other:
			for(i = 0; i < k; i++) // clusters array
			{
				diameter = 0;
				chunk_count = 0;

				// stop values for slaves:
				stop_inner_loop = FALSE;
				stop_external_loop = FALSE;

				// current cluster coordinates:
				x = clusters_arr[i].x;
				y = clusters_arr[i].y;

				for(j = 1; j < numprocs; j++) // send to slaves:
				{
					// current cluster
					MPI_Send(&x, 1, MPI_DOUBLE, j , 0, MPI_COMM_WORLD);
					MPI_Send(&y, 1, MPI_DOUBLE, j , 0, MPI_COMM_WORLD);

					// first job
					insertPointsToChunk(chunk, points_arr, chunk_count);
					MPI_Send(chunk, CHUNK_SIZE, MPI_POINT, j , 0, MPI_COMM_WORLD);
					chunk_count++;
				}

				for(j = numprocs-1; j < N/CHUNK_SIZE; j++) // send the remaining jobs
				{
					MPI_Recv(&answer, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status); // get max_dist for this chunk
					if(answer > diameter) // if greater than diameter - change it
					{
						diameter = answer;
					}

					// send next chunk
					insertPointsToChunk(chunk, points_arr, chunk_count);
					MPI_Send(chunk, CHUNK_SIZE, MPI_POINT, status.MPI_SOURCE , 0, MPI_COMM_WORLD);
					chunk_count++;

					MPI_Send(&stop_inner_loop, 1, MPI_INT, status.MPI_SOURCE , 0, MPI_COMM_WORLD); // sends stop value of slave's inner loop 
				}

				stop_inner_loop = TRUE;
				for(j = 1; j < numprocs; j++) // recv two last answers
				{
					MPI_Recv(&answer, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
					if(answer > diameter)
					{
						diameter = answer;
					}

					chunk[MASTER].id = STOP; // sending a junk value in that chunk (we already sent all jobs)
					MPI_Send(chunk, CHUNK_SIZE, MPI_POINT, status.MPI_SOURCE , 0, MPI_COMM_WORLD);

					MPI_Send(&stop_inner_loop, 1, MPI_INT, status.MPI_SOURCE , 0, MPI_COMM_WORLD);
				}

				if(diameter != 0)
				{
					calc_qm += calculateQM(k, i, diameter, clusters_arr);
				}


				for(j = 1; j < numprocs; j++) // sends stop value of slave's external loop
				{
					MPI_Send(&stop_external_loop, 1, MPI_INT, j , 0, MPI_COMM_WORLD);
				}


			}	

			if(calc_qm < QM) // stopping conditions (External loop) - qm < QM
			{
				break;
			}

		}

		stop_external_loop = TRUE;
		for(j = 1; j < numprocs; j++) // stop slaves - outer loop:
		{
			// current cluster
			MPI_Send(&x, 1, MPI_DOUBLE, j , 0, MPI_COMM_WORLD);
			MPI_Send(&y, 1, MPI_DOUBLE, j , 0, MPI_COMM_WORLD);

			// first job
			chunk[TRUE].id = STOP; // sending a junk value in that chunk - id in TRUE = 1 is STOP = -1
			MPI_Send(chunk, CHUNK_SIZE, MPI_POINT, j , 0, MPI_COMM_WORLD);

			MPI_Send(&stop_external_loop, 1, MPI_INT, j , 0, MPI_COMM_WORLD); // sends stop value of slave's external loop
		}

		end = MPI_Wtime(); // get end time
		printf("Time = %.3f\n", end - start );fflush(stdout);

		writeToFile(k, calc_qm, clusters_arr); // write results to output file

		free(points_arr);
		free(clusters_arr);
		free(chunk);

	}

	else // SLAVES
	{
		double max_dist;
		chunk = (point_t*)calloc(CHUNK_SIZE, sizeof(point_t));
		stop_external_loop = FALSE;

		while(stop_external_loop != TRUE)
		{
			// get current cluster:
			MPI_Recv(&x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&y, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

			// get first chunk
			MPI_Recv(chunk, CHUNK_SIZE, MPI_POINT, 0, 0, MPI_COMM_WORLD, &status);

			if(chunk[TRUE].id != STOP) // if the external loop is not done - enter the inner loop
			{
				stop_inner_loop = FALSE;
			}

			while(stop_inner_loop != TRUE)
			{
				if(chunk[MASTER].id != STOP) // if it is not junk value - enter the function
				{
					max_dist = calculateDiameter(x, y, chunk);
				}

				MPI_Send(&max_dist, 1, MPI_DOUBLE, 0 , 0, MPI_COMM_WORLD); // send answer to master

				MPI_Recv(chunk, CHUNK_SIZE, MPI_POINT, 0, 0, MPI_COMM_WORLD, &status); // get new chunk

				MPI_Recv(&stop_inner_loop, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status); // get stop value 
			}

			MPI_Recv(&stop_external_loop, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status); // get stop value
		}

		free(chunk);
	}

	MPI_Finalize();
	return 0;
}
