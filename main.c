#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>

#include <mpi.h>
#include "utility.h"

#define num_temp 100


double *simulation(double density, double temperature, int steps, double potentials[num_temp][3],int rank);

double deltatable(double x);

int main(){

	double density = 0.02;
	double temperature = 2.0;

	double *temps;

	int total_steps = 5000000;

	int steps,i;

	double potentials[num_temp][3];
	double final_potentials[num_temp][3];

	double cv[num_temp];

    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank==0){
    	generate_initial_configurations(world_size,density,temperature);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    srand(rank);
	steps = total_steps/world_size;
	if(total_steps%world_size!=0){
		int rest = total_steps%world_size;
		if(rank<rest){
			steps++;
		}
	}

	temps = simulation(density,temperature,steps,potentials,rank);

	
	double dpotentials[num_temp][2];
	for(i=0;i<num_temp;i++){
		dpotentials[i][0] = sqrt((1/steps)*(potentials[i][1]-potentials[i][0]*potentials[i][0]));
		dpotentials[i][1] = sqrt((1/steps)*(potentials[i][2]-potentials[i][1]*potentials[i][1]));
	}

	MPI_Allreduce(potentials,final_potentials,num_temp*3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	if(rank==0){
		char filename[80];
		sprintf(filename,"%lf.csv",density);
    	FILE *f = fopen(filename, "w+");
		for(i=0;i<num_temp;i++){
			final_potentials[i][0] /= world_size; 
			final_potentials[i][1] /= world_size; 
			final_potentials[i][2] /= world_size; 
			cv[i] = (final_potentials[i][1]-final_potentials[i][0]*final_potentials[i][0])/(temps[i]*temps[i]);
			//printf("Cv: %lf",cv[i]);
      		fprintf(f, "%lf,%lf\n",cv[i],temps[i]);
    		
    		
		}
		fclose(f);
	}

    MPI_Finalize();

    return 0;
}



double *simulation(double density, double temp_simul, int steps, double potentials[num_temp][3],int rank){ 

 	static double temp[num_temp];
	int i,t,k,j;
	double first_temp = 1.0;
	double last_temp = 4.0;
	double step = (last_temp-first_temp)/num_temp;
	for (i = 0; i < num_temp; i++) {
		temp[i] = step*i + first_temp;
	}

    double betas[num_temp];
    for (i = 0; i < num_temp; i++) {
		betas[i] = 1.0/temp[i];
	}
    double normalization[num_temp];
    
    double V = N/density;
    double L =pow(V,1.0/3.0);

    char filename[80];
    sprintf(filename,"./%lf/%d.csv",density,rank);
    double positions[N][3];
    readmatrix(positions,filename);
    printf("Attendi, trovo il migliore delta.... Processore %d\n",rank);
   
    double delta = deltatable(density);//find_delta2(density);//find_delta(temp_simul,L,0.01/(pow(density,1.0/3.0)),positions);//find_delta(temp_simul,L,0.02/(pow(density,1.0/3.0)));//0.05/(pow(density,1.0/3.0));//find_delta(temp_simul,L,0.01);
  
    printf("Processore %d Delta migliore trovato: %lf, adesso inizio la simulazione\n",rank,delta);

    int acceptance = 0;

    double beta = 1.0/temp_simul;
    

    //printmatrix(positions);
    double potential = calculate_potential(positions,L);
    //printf("Potenziale %lf\n",potential);
    double Vt_simul = potential;
    double Vt2_simul =  potential*potential;

    double new_positions[N][3];
    double new_potential;
    double p;
    double xi;
    
    srand((unsigned)rank);

    for(t=0;t<num_temp;t++)
    {
    	potentials[t][0] =  potential*exp((beta - betas[t])*Vt_simul);
        potentials[t][1] =  potential*potential*exp((beta - betas[t])*Vt_simul);
        potentials[t][2] =  pow(potential,4.0)*exp((beta - betas[t])*Vt_simul);
        normalization[t] = exp((beta - betas[t])*Vt_simul);

    }
    //loop principale della catena
    for(k=0;k<steps;k++)
    {
        for(j=0 ; j<N;j++){
        	for(i=0 ; i<3;i++){
        		new_positions[j][i] = positions[j][i]+myrandom(-delta/2,delta/2);
                new_positions[j][i] -= L*rint(new_positions[j][i]/L);
        	}
        }
        //printmatrix(new_positions);
        new_potential = calculate_potential(new_positions,L);
        //printf("Potenziale %lf\n",new_potential);
        //sleep(10);
        p = min(1,exp(-beta*(new_potential - potential)));
        xi = ((double)rand()/(double)(RAND_MAX));
        if (xi < p){
        	copy_matrix(new_positions,positions);
            potential = new_potential;
            acceptance+=1;
        }
        Vt_simul += potential;
        Vt2_simul+= potential*potential;
    	for(t=0;t<num_temp;t++)
    	{
    		potentials[t][0] +=  potential*exp((beta - betas[t])*potential);
       		potentials[t][1] +=  potential*potential*exp((beta - betas[t])*potential);
        	potentials[t][2] +=  pow(potential,4.0)*exp((beta - betas[t])*potential);
        	normalization[t] += exp((beta - betas[t])*potential);
    	}
    	
    	//sleep(1);
    }

    
    for(t=0;t<num_temp;t++)
    {
        potentials[t][0] /=  normalization[t];
        potentials[t][1] /=  normalization[t];  
        potentials[t][2] /=  normalization[t];     
    }
    printf("Processore: %d, Acceptance ratio: %lf step fatti: %d\n",rank,acceptance*100/(double)steps,k);

    return temp;
}
