#include <stdio.h>
#include <math.h>

struct ModelParams{
  double k;
  double c;
  double alpha;
  double d;
  double R;
  double L;
  double T;
  double q;
  double hT;
  double hZ;
};

double k = 0.59;
double c = 1.65;
double alpha = 0.003;
double d = 0.5;
double R = 10;
double L = 50;
double T = 200;
double q = 2 * alpha / pow(R, 2) * (R + d);

int calcIndex(int i, int j, int I){
    return j*(I + 1) + i;
}

double initCondition(int i, int I) {
    return 1 + cos(M_PI*i/I);
}

double leftCondition(double bottom, double bottom_right, const ModelParams params){
    return bottom - (params.q*params.hT)/params.c*bottom + (params.k*params.hT)/(params.c*pow(params.hZ, 2))*(bottom_right - bottom);
}

double innerNode(double bottom, double bottom_left, double bottom_right, const ModelParams params){
    return bottom - (params.q*params.hT)/params.c*bottom + (params.k*params.hT)/(params.c*pow(params.hZ, 2))*(bottom_right - 2 * bottom + bottom_left);
}


ModelParams prepareParams(int I, int J){
	ModelParams params;
    params.k = k;
    params.c = c;
    params.alpha = alpha;
    params.d = d;
    params.R = R;
    params.L = L;
    params.T = T;
    params.q = q;
    params.hT = T/J;
    params.hZ = L/I;

    return params;

}

void saveMatrix(double* matrix, int I, int J, char* filepath){
    FILE* fptr = fopen(filepath, "w");
    for(int j=0; j <= J; j++){
        for(int i=0; i<= I; i++){
            int idx = j*(I + 1) + i;
            fprintf(fptr, "%f ", matrix[idx]);
        }
        fprintf(fptr, "\n");
    } 
}