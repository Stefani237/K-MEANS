#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define MAX_NUM_POINTS 500000
#define CHUNK_SIZE 1000
#define POINT_FIELDS 5
#define CLUSTER_FIELDS 3
#define STOP -1
#define TRUE 1
#define FALSE 0
#define MASTER 0
#define INPUT "C:\\InitialMPIproject\\Initial\\data.txt"
#define OUTPUT "C:\\InitialMPIproject\\Initial\\outputData.txt"

struct Point
{
	int id;
	double x, y, cluster_x, cluster_y;
} typedef point_t;

struct Cluster
{
	double x, y;
	int num_of_points;
} typedef cluster_t;


MPI_Datatype mpiCluster();
MPI_Datatype mpiPoint();
void distance(double x1, double y1, double x2, double y2, double *dist);
void printCenters(int k, double qm, cluster_t *clusters_arr);
void readFromFile(int *N, int *MAX, int *LIMIT, double *QM, point_t *points_arr);
void writeToFile(int k, double qm, cluster_t *clusters_arr);
void initializeClusters(int k, int N, cluster_t *clusters_arr, point_t *points_arr);
void insertPointsToChunk(point_t *chunk, point_t *points_arr, const int counter);
void resetPointsAmount(int k, cluster_t *clusters_arr);
int groupPointsToClusters(int N, int k, cluster_t *clusters_arr, point_t *points_arr);
void calculateNewClustersCenters(int N, int k, cluster_t *clusters_arr, point_t *points_arr);
double calculateDiameter(double x, double y, point_t *chunk);
double calculateQM(int k, int i, double diameter, cluster_t *clusters_arr);