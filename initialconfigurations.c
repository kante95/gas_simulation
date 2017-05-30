#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>


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



void tmpsimulation(int n, double density,double temp){

    char fi[80];
    sprintf(fi,"%lf",density);
    mkdir(fi, 0700);
 
    double V = N/density;
    double L = pow(V,1.0/3.0);

    //inizializzazione
    double positions[N][3];
    initialize(L,positions); 

    //printf("Attendi, trovo il migliore delta.... Processore %d\n",rank);
   
    double delta = find_delta(temp,L,0.01/(pow(density,1.0/3.0)),positions);//find_delta(temp_simul,L,0.02/(pow(density,1.0/3.0)));//0.05/(pow(density,1.0/3.0));//find_delta(temp_simul,L,0.01);
  
    //printf("Processore %d Delta migliore trovato: %lf, adesso inizio la simulazione\n",rank,delta);

    int acceptance = 0;

    double beta = 1.0/temp;

    //printmatrix(positions);
    double potential = calculate_potential(positions,L);
    //printf("Potenziale %lf\n",potential);
   
    double new_positions[N][3];
    double new_potential;
    double p;
    double xi;
    
    int dadecidere = 5000;
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
            sprintf(filename,"./%lf/%d.csv",density,remaining-1);
            printmatrix_onfile(positions,filename);
            remaining--;
        }

        k++;
    }
}

void generate_initial_configurations(int n, double density, double temperature){
    char filename[80];
    int i;
    int ok=1;
    printf("Controllo se ci sono i file delle configurazioni iniziali...\n");
    for(i=0;i<n;i++){
        sprintf(filename,"./%lf/%d.csv",density,i);
        if(access(filename,F_OK) == -1){
            ok = 0;
            break;
        }
    }
    if(ok==0){
        printf("Configurazioni iniziali non trovate o incomplete, inizio a generarle...\n");
        srand(23);
        tmpsimulation(n,density, temperature);
        printf("Inizio simulazione\n");
    }
    else{
        printf("Ci sono tutti, inizio la simulazione...\n");
    }
}
