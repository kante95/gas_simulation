#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utility.h"

double min(double a,double b){
  double result = ((a > b)?b:a);
  return result;
}


void printmatrix(double matrix[N][3]){
	int i;
	for(i=0;i<N;i++){
		printf("Particella %d: [%lf] [%lf] [%lf]\n",i+1,matrix[i][0],matrix[i][1],matrix[i][2]);
	}
}

int find_maximum(double a[], int n) {
  int c, max, index;
 
  max = a[0];
  index = 0;
 
  for (c = 1; c < n; c++) {
    if (a[c] > max) {
       index = c;
       max = a[c];
    }
  } 
  return index;
}


void copy_matrix(double source[N][3],double destination[N][3]){
	int i,j;
	for(i=0;i<N;i++)
		for(j=0;j<3;j++)
			destination[i][j] = source[i][j];
}


double lennard_jones(double r){
  double r2 = r*r;
  double r6 = r2*r2*r2;
  double r12 = r6*r6;
  return 4*(1.0/r12-1.0/r6);
}

double calculate_potential(double positions[N][3],double L){

    double potential = 0;
    double r;
    int i,j,x;
    double d;
    for(i=0;i<N;i++){
      for(j=i+1;j<N;j++){
        r=0;
        for(x = 0; x<3;x++){
          d = positions[i][x] - positions[j][x]-L*rint((positions[i][x] - positions[j][x])/L);
          r+=d*d;
        }
        r = sqrt(r);
        if(r<L/2.0){
          potential+=2*lennard_jones(r);
        }
      }
    }

    return potential;
}


double myrandom(double from,double to){
  return ((double)rand()/(double)(RAND_MAX))*(to-from) + from;
}


void printmatrix_onfile(double matrix[N][3],char *filename){
    int i;
    FILE *f = fopen(filename, "w+");
    for(i=0;i<N;i++){
      fprintf(f, "%lf,%lf,%lf\n",matrix[i][0],matrix[i][1],matrix[i][2]);
    }
    fclose(f);
}


void readmatrix(double matrix[N][3],char *filename){
    int i=0,j=0;
    FILE *f = fopen(filename, "r");
    char *line = NULL;
    char *record;
    size_t len = 0;
    ssize_t read;

    while ((read = getline(&line, &len, f)) != -1) {
      j=0;
      record = strtok(line,",");
      while(record != NULL)
      {
        matrix[i][j] = atof(record);
        record = strtok(NULL,",");
        j++;
      }
      i++;
    }
    fclose(f);
}

double find_delta2(double density){
  return 0.03037587 + (0.4316096 - 0.03037587)/(1 + pow(density/0.008681893,2.252729));
}

double find_delta(double temperature,double L, double delta,double initial_position[N][3]){
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
    if(acceptance_ratio<70 && acceptance_ratio>45){
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
    float i = 8;
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
