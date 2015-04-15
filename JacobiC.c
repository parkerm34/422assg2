/* Jacobi iteration using pthreads and semaphores

     gcc JacobiC.c -o JacobiC -lpthread
     JacobiC N numProcs [L] [T] [R] [B] [E]
  Author: Parker Mathewson
  Date: 3/10/15
  Program: JacobiC.c

  Details:
	This file uses many of the alorithms found in the book, as well as functions names
	closely related to the onces in the book, because all of them are based off of the book.
	The barrier is bassed on a dissemination barrier but is a sort of own barrier I have
	made out of sempahores(...I think), one that counts to the number of processes, and the other to lock
	at the barrier until that number is reached.
*/

#define _REENTRANT
#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/times.h>
#include <time.h>
#include <limits.h>
#define SHARED 1
#define MAXGRID 258   /* maximum grid size, including boundaries */
#define MAXWORKERS 16  /* maximum number of worker threads */

void *Worker(void *);
void InitializeGrids();
void Barrier();

struct timeval startTime, endTime;

clock_t start, stop;
double time_spent;

sem_t bar;
sem_t e;
pthread_mutex_t barrier;  /* mutex semaphore for the barrier */
pthread_cond_t go;        /* condition variable for leaving */
int numArrived = 0;       /* count of the number who have arrived */

float EPSILON = 0.1;
float TOP = 1.0;
float RIGHT = 80.0;
float BOTTOM = 80.0;
float LEFT = 1.0;

int N, numProcs, MAXITERS, stripSize;
double maxDiff[MAXWORKERS];
double grid[MAXGRID][MAXGRID], newGrid[MAXGRID][MAXGRID];


void usage()
{
    printf("JacobiC Usage\n\n");
    printf("JacobiC N numProcs [L] [T] [R] [B] [E]\n");
    printf("  N is the size of the grid that will be NxN (required)\n");
    printf("  numProcs is the number of processes (threads) to create to solve the problem ");
    printf("(required, bounded by minimum of 1, maximum of 16)\n");
    printf("  L is the fixed value for the left edge of the grid [optional, default = 1.0]\n");
    printf("  T is the fixed value for the top edge of the grid [optional, default = 1.0]\n");
    printf("  R is the fixed value for the right edge of the grid [optional, default = 80.0]\n");
    printf("  B is the fixed value for the bottom edge of the grid [optional, default = 80.0]\n");
    printf("  E is the epsilon value used to determine when to stop the computation [optional, default = 0.1]\n\n");
    printf("This program will print out the settings to a file called JacobiC.log\n");
    printf("as well as the final grid contents.\n\n");
    printf("The number of processes, iterations, and computation time will be printed out\n");
    printf("to standard output.\n\n");

    exit(1);
}


int main(int argc, char *argv[])
{
  pthread_t workerid[MAXWORKERS];
  pthread_attr_t attr;
  int i, j;
  double maxdiff = 0.0;
  FILE *results;
  float le = 0.0, to = 0.0, ri = 0.0, bo = 0.0, ep = 0.0;
  int seconds;
  double microseconds;
  int ugh;
  int ugh2;

  ugh = sem_init(&e, 1, 16);
  ugh2 = sem_init(&bar, 1, 1);

  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  pthread_mutex_init(&barrier, NULL);
  pthread_cond_init(&go, NULL);

  if(argc < 3)
    usage();
  N = atoi(argv[1]);
  numProcs = atoi(argv[2]);
  if(argc > 3)
  {
    if(argc != 8)
    {
      usage();
    }
    if(argv[3] != NULL)
        LEFT = atof(argv[3]);
    if(argv[4] != NULL)
        TOP = atof(argv[4]);
    if(argv[5] != NULL)
        RIGHT = atof(argv[5]);
    if(argv[6] != NULL)
        BOTTOM = atof(argv[6]);
    if(argv[7] != NULL)
        EPSILON = atof(argv[7]);
  }
  MAXITERS = 20000;
  stripSize = N / numProcs;
  InitializeGrids();

  start = clock();
  gettimeofday(&startTime, NULL);

  for (i = 0; i < numProcs; i++)
    pthread_create(&workerid[i], &attr, Worker, (void *) i);
  for (i = 0; i < numProcs; i++)
    pthread_join(workerid[i], NULL);

  gettimeofday(&endTime, NULL);
  stop = clock();
  time_spent = (double)(stop - start) / CLOCKS_PER_SEC;
  seconds = endTime.tv_sec - startTime.tv_sec;
  if(seconds > 0)
    microseconds = startTime.tv_usec - endTime.tv_usec;
  else
    microseconds = endTime.tv_usec - startTime.tv_usec;

  for (i = 0; i < numProcs; i++)
  {
    if (maxdiff < maxDiff[i])
      maxdiff = maxDiff[i];
  }
  printf("main: numProcs = %d, N = %d\n", numProcs, N);
  printf("execution time:  %d seconds, %.0f microseconds\n", seconds, microseconds);
  results = fopen("JacobiC.log", "w");

  fprintf(results, "Grid\t = %dx%d\nnumProcs = %d\nleft\t = %f\ntop\t = %f\n", N, N, numProcs, LEFT, TOP);
  fprintf(results, "right\t = %f\nbottom\t = %f\nepsilon\t = %f\n", RIGHT, BOTTOM, EPSILON);
  fprintf(results, "execution time = %d seconds, %.0f microseconds\n\n", seconds, microseconds);
  for (i = 0; i <= N+1; i++)
  {
    fprintf(results, "\t");
    for (j = 0; j <= N + 1; j ++)
    {
      fprintf(results, "%07.4f ", newGrid[i][j]);
    }
    fprintf(results, "\n");
  }
}


void *Worker(void *arg) {
  int myid = (int) arg;
  double maxdiff, temp;
  int i, j, iters;
  int first, last;

  first = myid * stripSize + 1;
  last = first + stripSize - 1;

  for (iters = 1; iters <= MAXITERS; iters ++)
  {
    for (i = first; i <= last; i ++)
    {
      for (j = 1; j <= N; j ++)
        newGrid[i][j] = (grid[i - 1][j] + grid[i + 1][j] + grid[i][j - 1] + grid[i][j + 1]) * 0.25;
    }
    Barrier();

    for (i = first; i <= last; i ++)
    {
      for (j = 1; j <= N; j ++)
      {
        grid[i][j] = (newGrid[i - 1][j] + newGrid[i + 1][j] + newGrid[i][j - 1] + newGrid[i][j + 1]) * 0.25;

        if(abs(grid[i][j] - newGrid[i][j]) > maxdiff)
	  maxdiff = abs(grid[i][j] - newGrid[i][j]);
      }
    }
    maxDiff[myid] = maxdiff;
    if(maxdiff < EPSILON)
      break;
    maxdiff = 0.0;
    Barrier();
  }
}

void InitializeGrids()
{
  int i, j;
  for (i = 0; i <= N + 1; i ++)
  {
    for (j = 0; j <= N + 1; j ++)
    {
      grid[i][j] = 0.0;
      newGrid[i][j] = 0.0;
    }
  }
  for (i = 0; i <= N + 1; i ++)
  {
    grid[i][0] = LEFT;
    newGrid[i][0] = LEFT;
    grid[i][N + 1] = RIGHT;
    newGrid[i][N + 1] = RIGHT;
  }
  for (j = 0; j <= N + 1; j ++) 
  {
    grid[0][j] = TOP;
    newGrid[0][j] = TOP;
    grid[N + 1][j] = BOTTOM;
    newGrid[N + 1][j] = BOTTOM;
  }
}

void Barrier()
{
  int result, result2;
  int i = 0;

  result2 = sem_wait(&e);
  result = sem_wait(&bar);

  numArrived++;
  if(numArrived == numProcs)
  {
    numArrived = 0;
    
    for(i = 0; i < numProcs; i++ )
    {
      result2 = sem_post(&e);
    }
  }

  result = sem_post(&bar);

}
