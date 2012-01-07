/* Compute/draw mandelbrot set using MPI/MPE commands
 * Written Winter, 1998, W. David Laverell.
 * Simplified Winter 2002, Joel Adams. 
 */

//#define USE_MPE 1
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#include "display.h"
#endif


#ifdef USE_MPE   
    MPE_XGraph graph;
#endif

int total_process;

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

void drawRow(int rowIndex)
{
	    double     spacing=.005, x, y, c_real, c_imag, x_center = 0.0, y_center = 0.0;
		int ix, n;
	    const int  WINDOW_SIZE = 1024;
for (ix = 0; ix < WINDOW_SIZE; ix++)
        	{
         		c_real=(ix - 400) * spacing - x_center;
         		c_imag=(rowIndex - 400) * spacing - y_center;
         		x = y = 0.0;
         		n = 0;

        		 while (n < 50 && distance(x,y) < 4.0)
         		{
           		 compute(x,y,c_real,c_imag,&x,&y);
          		  n++;
        		 }

	 		if (n < 50) {
#ifdef USE_MPE
          		MPE_Draw_point(graph,ix,rowIndex,MPE_RED);
#else
	 		//fprintf (stdout, "%d,%d,RED\n", ix, rowIndex);
#endif
	 		} else {
#ifdef USE_MPE
	   		MPE_Draw_point(graph,ix,rowIndex,MPE_BLACK);
#else
	   		//fprintf (stdout, "%d,%d,BLACK\n", ix, rowIndex);
#endif
	 		}
       		}

}

void masterProcess()
{
	int i; 		//Loop index
	int b_is_done;
	int flag, source;
	MPI_Status status;
	int message;
	int cur_row_index = 0;
	for(i=1; i < total_process; i++)
	{
		MPI_Send(&cur_row_index, 1, MPI_INT, i, 199, MPI_COMM_WORLD);
		cur_row_index++;
	}

	b_is_done = 0;

	while(!b_is_done)
	{
		MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
		if(flag)
		{
			source = status.MPI_SOURCE;
	  		MPI_Recv(&message, 1, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if(cur_row_index == 1024)
			{
				b_is_done = 1;
				continue;
			}

			MPI_Send(&cur_row_index, 1, MPI_INT, source, 199, MPI_COMM_WORLD);
			cur_row_index++;
		}
	}
	for(i = 1; i<total_process; i++)
		  MPI_Send(NULL, 0, MPI_CHAR, i, 200, MPI_COMM_WORLD);
}


void slaveProcess()
{
	int b_is_done = 0;
	int source=0;        // 0 means source is master process
	MPI_Status status;
	int message;
	int  WINDOW_SIZE = 1024;

  while (!b_is_done)
    {
      MPI_Recv(&message, 1, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      switch (status.MPI_TAG)
	{
	case 199:
		drawRow(message);
	  	MPI_Send(&WINDOW_SIZE, 1, MPI_INT, 0, 199, MPI_COMM_WORLD);
	  break;

	case 200:    // tells slave to stop waiting for messages and terminate
	  b_is_done=1;
	  break;
	}
    }
}

int main(int argc, char* argv[])
{
	const int  WINDOW_SIZE = 1024;
	int cur_rank;
	
	double starttime, endtime;
    MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &cur_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &total_process);

#ifdef USE_MPE
    MPE_Open_graphics( &graph, MPI_COMM_WORLD, 
                         getDisplay(),
                         -1, -1,
                         WINDOW_SIZE, WINDOW_SIZE, 0 );
#endif	

	starttime = MPI_Wtime();

  	switch (cur_rank)
    	{
    		case 0:
      			masterProcess();
      			break;
    		default:
      			slaveProcess();
      			break;
    	}
	endtime = MPI_Wtime();
	printf("Time cost for process %d = %f\n", cur_rank, endtime - starttime );

	getchar();

#ifdef USE_MPE
    MPE_Close_graphics( &graph );
#endif

    	MPI_Finalize();
    	return 0;
}

