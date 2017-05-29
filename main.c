#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <mpi.h>
#include "utility.h"

#define num_temp 50


double *simulation(double density, double temperature, int steps, double potentials[num_temp][3],int rank);

int main(){

	double density = 0.01;
	double temperature = 2.0;

	double *temps;

	int total_steps = 100000;

	int steps;

	double potentials[num_temp][3];
	double final_potentials[num_temp][3];

    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank==0){
    	generate_initial_configurations(world_size,density,temperature);
    }
    MPI_Barrier(MPI_COMM_WORLD);

	steps = total_steps/world_size;
	if(total_steps%world_size!=0){
		int rest = total_steps%world_size;
		if(rank<rest){
			steps++;
		}
	}

	temps = simulation(density,temperature,steps,potentials,rank);

	//double dpotentials[num_temp][2];
	//for(i=0;i<num_temp;i++){
	//	dpotentials[i][0] = sqrt((1/steps)*())
	//}

	MPI_Reduce(potentials,final_potentials,num_temp*3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	if(rank==0){
		printf("Potenziale medio: %lf, Potenziale medio quadro: %lf, Potenziale medio alla quarta: %lf\n",final_potentials[0][0]/world_size,final_potentials[0][1]/world_size,final_potentials[0][2]/world_size);
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
    sprintf(filename,"%d.csv",rank);
    double positions[N][3];
    readmatrix(positions,filename);
    printf("Attendi, trovo il migliore delta.... Processore %d\n",rank);
   
    double delta = find_delta(temp_simul,L,0.01/(pow(density,1.0/3.0)),positions);//find_delta(temp_simul,L,0.02/(pow(density,1.0/3.0)));//0.05/(pow(density,1.0/3.0));//find_delta(temp_simul,L,0.01);
  
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
        		new_positions[j][i] = positions[j][i]+myrandom(-delta/2,delta/2);;
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
            //printf("accettata!!\n");
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
    }

    
    for(t=0;t<num_temp;t++)
    {
        potentials[t][0] /=  normalization[t];
        potentials[t][1] /=  normalization[t];  
        potentials[t][2] /=  normalization[t];     
    }

    printf("Processore: %d, Acceptance ratio: %lf\n",rank,acceptance*100/(double)steps);

    return temp;
}