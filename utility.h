#define N 512

double min(double a,double b);

void printmatrix(double matrix[N][3]);

int find_maximum(double a[], int n);

void copy_matrix(double source[N][3],double destination[N][3]);

double lennard_jones(double r);

double calculate_potential(double positions[N][3],double L);

double myrandom(double from,double to);

void printmatrix_onfile(double matrix[N][3],char *filename);

void generate_initial_configurations(int n, double density, double temperature);

void readmatrix(double matrix[N][3],char *filename);

double find_delta(double temperature,double L, double delta,double initial_position[N][3]);

double find_delta2(double density);

void initialize(double L,double positions[N][3]);

