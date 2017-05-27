#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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
        if(r<L/2.0 && r!=0.0){
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