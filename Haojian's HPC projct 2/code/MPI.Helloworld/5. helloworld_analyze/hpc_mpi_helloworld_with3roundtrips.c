#include <mpi.h>
#include <stdio.h>
#include <time.h>

int  main(int argc, char **argv)
{
	MPI_Status status;
	int cur_rank, total_process, repeated_times, tag, next, from, message;
	double starttime, endtime;

	MPI_Init(&argc, &argv);
       MPI_Comm_rank (MPI_COMM_WORLD, &cur_rank);
       MPI_Comm_size (MPI_COMM_WORLD, &total_process);
        
	
	tag = 199;
	next = (cur_rank + 1) % total_process;
	from = (cur_rank + total_process -1) % total_process;
	int i;
	for(i=0; i< 4; i++)
	{
		tag += i;
		repeated_times =1;
		if(cur_rank == 0 )
		{
			//printf("Proces %d sending %d to %d\n", cur_rank, cur_rank, next);
			starttime = MPI_Wtime();
			MPI_Send(&cur_rank, 1 , MPI_INT, next, tag, MPI_COMM_WORLD);
		}

	while(repeated_times !=0)
	{
		MPI_Recv(&message, 1, MPI_INT, from, tag, MPI_COMM_WORLD, &status);
		
		repeated_times--;
		//printf("Process %d received %d \n", cur_rank, message);
		if(cur_rank == 0)
		{
			endtime = MPI_Wtime();
			printf("Time cost for message transmission = %f\n", endtime - starttime );
			break;
		}
		else
		{
			//printf("Proces %d sending %d to %d\n", cur_rank, cur_rank, next);
			MPI_Send(&cur_rank, 1 , MPI_INT, next, tag, MPI_COMM_WORLD);
		}
	}
	}
        MPI_Finalize();
	return 0;
}
