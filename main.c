#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
//#include <unistd.h>
#include <mpi.h>
#include "utility.h"

int main(){
    time_t t;
    int num_density=28;
    srand((unsigned) time(&t));

    int steps = 50000;
    double density[num_density];
    double temp[num_density];

    int i;
    double step = (1.0-0.01)/num_density;
    for (i = 0; i < num_density; i++) {
        density[i] = step*i + 0.01;
    }

    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


    int nstart = world_rank*num_density/world_size;
    int nend = nstart + num_density/world_size;

    //for(i=nstart;i<nend;i++){
     //   temp[i] = simulation(density[i],steps,world_rank);
        //printf("Processore: %d, densità: %lf, temperatura: %lf\n",world_rank,density[i],temp[i]);
    //}


    MPI_Finalize();

    return 0;
}