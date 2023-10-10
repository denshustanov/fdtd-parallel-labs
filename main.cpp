#include "heat_scheme.cpp"

using namespace std;


void calcMSecuential(double* m, int I, int J, const ModelParams params){
    double hT = params.T/J;
    double hZ = params.L/I;
    double bottom;
    double bottomLeft;
    double bottomRight;
    for(int i=0; i <= I; i++){
        m[calcIndex(i, 0, I)] = initCondition(i, I);
    }

    for(int j=0; j < J; j++){
        m[calcIndex(I, j, I)] = 0;
        bottom = m[calcIndex(0, j, I)];
        bottomRight = m[calcIndex(1, j, I)];
        m[calcIndex(0, j+1, I)] = leftCondition(bottom, bottomRight, params);
        for(int i=1; i<I; i++){
            bottom = m[calcIndex(i, j, I)];
            bottomRight = m[calcIndex(i + 1, j, I)];
            bottomLeft = m[calcIndex(i - 1, j, I)];
            m[calcIndex(i, j+1, I)] = innerNode(bottom, bottomLeft, bottomRight, params);
        }
    }
}


int main(){
	int I = 10000;
    int J = 1000;

   	ModelParams params = prepareParams(I, J);
    int vecSize = (I+1)*(J+1);

    double *m = (double *) calloc(vecSize, sizeof(double));
    calcMSecuential(m, I, J, params);

    saveMatrix(m, I, J, "output.txt");

	return 0;
}