#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "utility.h"


void initialize(double L,double positions[N][3]){
    float i = 5;
    int index = 0;
    float j,k,m;
    for(j =0;j<i;j++){
    	for(k =0;k<i;k++){
    		for(m =0;m<i;m++){
    			positions[index][0]= j*(L/i) - L/2.0;
                positions[index][1] = k*(L/i) - L/2.0;
                positions[index][2] = m*(L/i) - L/2.0;
                index++;
			}
    	}
    }
}


double tmpfind_delta(double temperature,double L, double delta){
    int steps = 0;
    int acceptance = 0;
    float acceptance_ratio=0;
    double beta = 1.0/temperature;
    int i,j;
    double positions[N][3];
    double potential;
    double new_potential;
    double p;
    double xi;
    double new_positions[N][3];

    initialize(L,positions);
    potential = calculate_potential(positions,L);
    while(steps<1000){
    	//printmatrix(positions);
        for(j=0; j<N;j++){
        	for(i=0; i<3;i++){
        		new_positions[j][i] = positions[j][i]+myrandom(-delta/2,delta/2);
                new_positions[j][i] -= L*rint(new_positions[j][i]/L);
        	}
        }  
        //printmatrix(new_positions);
        new_potential = calculate_potential(new_positions,L);
        p = min(1,exp(-beta*(new_potential - potential)));
        xi = ((double)rand()/(double)(RAND_MAX));
        if (xi < p){
        	copy_matrix(new_positions,positions);
            potential = new_potential;
            acceptance+=1;
        }
        steps+=1;

        
    }
    acceptance_ratio = acceptance*100/(steps-1);
    //printf("Delta: %lf, Acceptance ratio: %lf\n",delta,acceptance_ratio);
    if(acceptance_ratio<60 && acceptance_ratio>40){
        return delta;
    }
    else if(acceptance_ratio == 0.0){
        return tmpfind_delta(temperature,L, delta/1.5);
    }
    else{
        return tmpfind_delta(temperature,L, delta*(acceptance_ratio/50)); //???? più graduale delta 
    }
}


void tmpsimulation(int n, double density,double temp){
 
    double V = N/density;
    double L = pow(V,1.0/3.0);

    //printf("Attendi, trovo il migliore delta.... Processore %d\n",rank);
   
    double delta = tmpfind_delta(temp,L,0.01/(pow(density,1.0/3.0)));//find_delta(temp_simul,L,0.02/(pow(density,1.0/3.0)));//0.05/(pow(density,1.0/3.0));//find_delta(temp_simul,L,0.01);
  
    //printf("Processore %d Delta migliore trovato: %lf, adesso inizio la simulazione\n",rank,delta);

    int acceptance = 0;

    double beta = 1.0/temp;
    //inizializzazione
    double positions[N][3];
    initialize(L,positions); 
    //printmatrix(positions);
    double potential = calculate_potential(positions,L);
    //printf("Potenziale %lf\n",potential);
   
    double new_positions[N][3];
    double new_potential;
    double p;
    double xi;
    
    int dadecidere = 10000;
    char filename[80];
    int i,j;
    int remaining = n;
    int k=0;
    //loop principale della catena
    while(remaining>0)
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
        p = min(1,exp(-beta*(new_potential - potential)));
        xi = ((double)rand()/(double)(RAND_MAX));
        if (xi < p){
        	copy_matrix(new_positions,positions);
            potential = new_potential;
            acceptance+=1;
            //printf("accettata!!\n");
        }
        
        if(k%dadecidere){
            sprintf(filename,"%d.csv",remaining-1);
            printmatrix_onfile(positions,filename);
            remaining--;
        }

        k++;
    }
}

void generate_initial_configurations(int n, double density, double temperature){
    time_t t;

    srand(23);

    tmpsimulation(n,density, temperature);
}
