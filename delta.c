#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utility.h"

#define N 125


double find_delta3(double temperature,double L, double delta,double initial_position[N][3]){
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

    copy_matrix(initial_position,positions);
    potential = calculate_potential(positions,L);
    while(steps<50000){
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
    if(acceptance_ratio<55 && acceptance_ratio>45){
        return delta;
    }
    else if(acceptance_ratio == 0.0){
        return find_delta(temperature,L, delta/1.5,initial_position);
    }
    else{
        return find_delta(temperature,L, delta*(acceptance_ratio/50),initial_position); //???? più graduale delta + (acceptance_ratio-50)%(L/2)
        //delta*(acceptance_ratio/50)
    }
} 


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


int main(){

	double temp = 2.0;
	int i;
	int num_dens = 50;
	double first = 0.01;
	double last = 1.0;
	double step = (last-first)/num_dens;
	double density[num_dens];
	double delta[num_dens];
	double positions[N][3];
	double V,L;
	for (i = 0; i < num_dens; i++) {
		density[i] = step*i + first;
	}

	for (i = 0; i < num_dens; i++) {
		V = N/density[i];
    	L =pow(V,1.0/3.0);
		initialize(L,positions);
		printf("Densità: %lf delta: %lf\n",density[i],find_delta3(temp,L,0.01/(pow(density[i],1.0/3.0)),positions));
	}

	return 0;
}