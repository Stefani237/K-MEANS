# include "header.h"

MPI_Datatype mpiCluster()
{
	// create new MPI_Datatype for cluster struct
	Cluster cluster;
	int block_lengths[CLUSTER_FIELDS] = {1, 1, 1}; // cluster contains 3 fields - each field is sized 1
	MPI_Aint offsets[CLUSTER_FIELDS];
	MPI_Datatype types[CLUSTER_FIELDS] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT}; // field's types
	MPI_Datatype MPI_CLUSTER;

	offsets[0] = (char*)&(cluster.x) - (char*)&(cluster);
	offsets[1] = (char*)&(cluster.y) - (char*)&(cluster);
	offsets[2] = (char*)&(cluster.num_of_points) - (char*)&(cluster);

	MPI_Type_create_struct(CLUSTER_FIELDS, block_lengths, offsets, types, &MPI_CLUSTER);
	MPI_Type_commit(&MPI_CLUSTER);
	return MPI_CLUSTER;
}


MPI_Datatype mpiPoint()
{
	// create new MPI_Datatype for point struct
	Point point;
	int block_lengths[POINT_FIELDS] = {1, 1, 1, 1, 1};  // point contains 5 fields - each field is sized 1
	MPI_Aint offsets[POINT_FIELDS];
	MPI_Datatype types[POINT_FIELDS] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; // field's types
	MPI_Datatype MPI_POINT;

	offsets[0] = (char*)&(point.id) - (char*)&(point);
	offsets[1] = (char*)&(point.x) - (char*)&(point);
	offsets[2] = (char*)&(point.y) - (char*)&(point);
	offsets[3] = (char*)&(point.cluster_x) - (char*)&(point);
	offsets[4] = (char*)&(point.cluster_y) - (char*)&(point);


	MPI_Type_create_struct(POINT_FIELDS, block_lengths, offsets, types, &MPI_POINT);
	MPI_Type_commit(&MPI_POINT);
	return MPI_POINT;
}

void distance(double x1, double y1, double x2, double y2, double *dist)
{	

	double dx = 0, dy = 0;
	dx = x1 - x2;
	dy = y1 - y2;
	*dist = sqrt(dx*dx + dy*dy); 
}


void readFromFile(int *N, int *MAX, int *LIMIT, double *QM, point_t *points_arr)
{

	int i;
	FILE *file;

	file = fopen(INPUT, "r");

	if (file == NULL) // opening file failed 
	{
		printf("Failed opening the file. Exiting!\n");
		return;
	}

	// reads the first rows from file:
	fscanf(file, "%d   %d   %d   %lf", N, MAX, LIMIT, QM);

	// read all points from file and insert them into an array:
	for(i = 0; i < *N; i++)
	{
		fscanf(file, "%d %lf %lf", &points_arr[i].id, &points_arr[i].x, &points_arr[i].y);	 
	}

	fclose(file); // close file after we finished reading all of it
	
}

void writeToFile(int k, double qm, cluster_t *clusters_arr)
{

	int i;
	FILE *file;

	file = fopen(OUTPUT, "w"); 
	if (file == NULL) // opening file failed 
	{
		printf("Failed opening the file. Exiting!\n");
		return;
	}
	fseek(file, 0, SEEK_SET);

	// write results to file:
	fputs("Number of best measure clusters \n", file);
	fprintf(file, "K = %d  qm = %lf\n", k, qm);
	fputs("Centers of the clusters: \n", file);

	for (i = 0; i < k; i++)
	{
		fprintf(file, "x = %lf, y = %lf\n", clusters_arr[i].x, clusters_arr[i].y);
	}

	fclose(file); //close file
	
}

// This function initializes k clusters to be the first k points read from the file.
void initializeClusters(int k, int N, cluster_t *clusters_arr, point_t *points_arr)
{
	int i;

#pragma omp parallel for
	for(i = 0; i < k; i++) 
	{
		clusters_arr[i].x = points_arr[i].x;
		clusters_arr[i].y = points_arr[i].y;
	}

}

// This function inserts each time 1000 diffrent points into chunk. 
void insertPointsToChunk(point_t *chunk, point_t *points_arr, const int counter)
{
	int i = 0,j;

	for(j = CHUNK_SIZE*counter; j < CHUNK_SIZE*(counter+1); j++)
	{
		chunk[i].id = points_arr[j].id;
		chunk[i].x = points_arr[j].x;
		chunk[i].y = points_arr[j].y;
		chunk[i].cluster_x = points_arr[j].cluster_x;
		chunk[i].cluster_y = points_arr[j].cluster_y;	

		i++;
	}
}

// Every time the points are re-assigned to the clusters, the amount of points of all clusters must be reset.
void resetPointsAmount(int k, cluster_t *clusters_arr)
{
	int i;

#pragma omp parallel for
	for(i = 0; i < k; i++) // clusters array
	{
		clusters_arr[i].num_of_points = 0;
	}
}

// for each point - check it's distance from every cluster center to find the closest:
int groupPointsToClusters(int N, int k, cluster_t *clusters_arr, point_t *points_arr)
{
	int i, j, cluster, count = 0;
	double min_dist = MAX_NUM_POINTS, dist, old_x, old_y;

	resetPointsAmount(k, clusters_arr);

	for(j = 0; j < N; j++) // points array
	{
		for(i = 0; i < k; i++) // clusters array
		{
			// calculate the distance from current point to each cluster center
			distance(clusters_arr[i].x, clusters_arr[i].y, points_arr[j].x, points_arr[j].y, &dist);
			
			 
			if(dist < min_dist) // if current distance is smaller then min_dist - update min_dist and save clusters's index.
			{
				min_dist = dist;
				cluster = i;
			}
		}

		// save old cluster center
		old_x = points_arr[j].cluster_x;
		old_y = points_arr[j].cluster_y;

		if(min_dist != MAX_NUM_POINTS) // only if it's not equal to initial value
		{
			points_arr[j].cluster_x = clusters_arr[cluster].x;
			points_arr[j].cluster_y = clusters_arr[cluster].y;
			clusters_arr[cluster].num_of_points++; // count how many points are in this cluster.

			// if this point didn't change her cluster 
			if(old_x == points_arr[j].cluster_x && old_y == points_arr[j].cluster_y)
			{
				count++;
			}

		}

		min_dist = MAX_NUM_POINTS;
	}
	return count;
}

void calculateNewClustersCenters(int N, int k, cluster_t *clusters_arr, point_t *points_arr)
{
	int i, j;
	double x, y, sum_x = 0, sum_y = 0;


	for(i = 0; i < k; i++) // clusters array
	{
		if(clusters_arr[i].num_of_points != FALSE) // empty cluster won't change it's center
		{
			x = clusters_arr[i].x;
			y = clusters_arr[i].y;
#pragma omp parallel for private (j), reduction (+: sum_x), reduction (+: sum_y)
			for(j = 0; j < N; j++) // points array
			{
				
				if(points_arr[j].cluster_x == x && points_arr[j].cluster_y == y) // if point is in this cluster
				{

					sum_x += points_arr[j].x;
					sum_y += points_arr[j].y;
				}

			}
			// the average will give us the new cluster centers
			clusters_arr[i].x = sum_x/clusters_arr[i].num_of_points;
			clusters_arr[i].y = sum_y/clusters_arr[i].num_of_points;

		}
		sum_x = 0;
		sum_y = 0;
	}
}

double calculateDiameter(double x, double y, point_t *chunk)
{
	int j, i;
	double max_dist = 0, dist, x1, x2, y1, y2;

#pragma omp parallel for private (i, x1, x2, y1, y2, dist)
	for(j = 0; j < CHUNK_SIZE; j++) // points array
	{
		for(i = 0 ; i < CHUNK_SIZE; i++) // points array
		{
			
			if(i != j) // if this is not the same point
			{
				x1 = chunk[j].cluster_x;
				y1 = chunk[j].cluster_y;
				x2 = chunk[i].cluster_x;
				y2 = chunk[i].cluster_y;

				if(x == x1  && x == x2 && y == y1 && y == y2) // if both points are from current cluster
				{
					// calculate distance:
					distance(chunk[j].x, chunk[j].y, chunk[i].x, chunk[i].y, &dist);

					if(dist > max_dist) // if current calculated distance is greather then max_dist, replace it.
					{
						max_dist = dist;
					}
				}

			}
		}

	}
	
	
	return max_dist; // diameter for this clustr

}

double calculateQM(int k, int i, double diameter, cluster_t *clusters_arr)
{
	int idx;
	double dist, qm = 0;
#pragma omp parallel for reduction (+: qm)
	for(idx = 0; idx < k; idx++) // clusters array
	{
		if(i != idx) // no need to check distance of center with itself
		{
			distance(clusters_arr[idx].x, clusters_arr[idx].y, clusters_arr[i].x, clusters_arr[i].y, &dist);
			qm += (diameter / dist);
		}

	}
	return qm;
}