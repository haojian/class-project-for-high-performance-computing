/* Compute/draw mandelbrot set using MPI/MPE commands
 * Written Winter, 1998, W. David Laverell.
 * Simplified Winter 2002, Joel Adams. 
 */

//#define USE_MPE 1
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#include "display.h"
#endif


/* compute the mandelbrot-set function for a given
 *  point in the complex plane.
 * Receive: doubles x and y,
 *          complex c.
 * Modify: doubles ans_x and ans_y.
 * POST: ans_x and ans_y contain the results of the mandelbrot-set
 *        function for x, y, and c.
 */
void compute(double x, double y, double c_real, double c_imag,
        double *ans_x, double *ans_y)
{
        *ans_x = x*x - y*y + c_real;
        *ans_y = 2*x*y + c_imag;
}

/* compute the 'distance' between x and y.
 * Receive: doubles x and y.
 * Return: x^2 + y^2.
 */
double distance(double x, double y)
{
        return(x*x + y*y);
}


void taskSchedule(int process_number, int window_size, int **task_pool)
{
	int maxtask_process= (int)(window_size/process_number) + 1;
	int i;
	for(i = 0; i<process_number; i++)
	{
		task_pool[i] = (int*)malloc(sizeof(int)*maxtask_process);
	}
	int *taskIndex_process = (int*)malloc(sizeof(int)*process_number);
	for(i=1; i<=window_size; i++)
	{
		int target_process= i % process_number;
		task_pool[target_process][taskIndex_process[target_process]] = i;
		taskIndex_process[target_process]++;
	}
}

void free_taskpool(int process_number, int **task_pool)
{
	int i;
	for(i=0; i< process_number; i++)
		free(task_pool[i]);
	free(task_pool);
}

int main(int argc, char* argv[])
{
    	const int  WINDOW_SIZE = 1024;

    	int        n,ix, iy, i;
    	double     spacing=.005, x, y, c_real, c_imag, x_center = 0.0, y_center = 0.0;

	int cur_rank, total_process, message;
	int tag = 199;
	MPI_Status status;

	double starttime, endtime;

	
#ifdef USE_MPE   
    MPE_XGraph graph;
#endif
    MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &cur_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &total_process);

#ifdef USE_MPE
    MPE_Open_graphics( &graph, MPI_COMM_WORLD, 
                         getDisplay(),
                         -1, -1,
                         WINDOW_SIZE, WINDOW_SIZE, 0 );
#endif						
		
	i=0;
	MPI_Send(&cur_rank,1, MPI_INT, cur_rank, tag, MPI_COMM_WORLD);

	while(1)
	{
		MPI_Recv(&message, 1, MPI_INT, cur_rank, tag, MPI_COMM_WORLD, &status);

		iy = message;
		for (ix = 0; ix < WINDOW_SIZE; ix++)
        	{
         		c_real=(ix - 400) * spacing - x_center;
         		c_imag=(iy - 400) * spacing - y_center;
         		x = y = 0.0;
         		n = 0;

        		 while (n < 50 && distance(x,y) < 4.0)
         		{
           		 compute(x,y,c_real,c_imag,&x,&y);
          		  n++;
        		 }

	 		if (n < 50) {
#ifdef USE_MPE
          		MPE_Draw_point(graph,ix,iy,MPE_RED);
#else
	 		//fprintf (stdout, "%d,%d,RED\n", ix, iy);
#endif
	 		} else {
#ifdef USE_MPE
	   		MPE_Draw_point(graph,ix,iy,MPE_BLACK);
#else
	   		//fprintf (stdout, "%d,%d,BLACK\n", ix, iy);
#endif
	 		}
       		}
		int nextRowIndex = message + total_process;
		if( nextRowIndex > WINDOW_SIZE)
			break;
		else
			MPI_Send(&nextRowIndex,1, MPI_INT, cur_rank, tag, MPI_COMM_WORLD);
	}
	getchar();

#ifdef USE_MPE
    MPE_Close_graphics( &graph );
#endif

    MPI_Finalize();
    return 0;
}

