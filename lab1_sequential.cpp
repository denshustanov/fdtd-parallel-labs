#include "heat_scheme.cpp"
#include <chrono>
#include <iostream>

using namespace std;


void calcM(double* m, int I, int J, const ModelParams params){
    double hT = params.T/J;
    double hZ = params.L/I;
    double bottom;
    double bottomLeft;
    double bottomRight;
    for(int i=0; i <= I; i++){
        m[calcIndex(i, 0, I)] = initCondition(i, I);
    }

    for(int j=1; j <= J; j++){
        bottom = m[calcIndex(0, j-1, I)];
        bottomRight = m[calcIndex(1, j-1, I)];
        m[calcIndex(0, j, I)] = leftCondition(bottom, bottomRight, params);
        for(int i=1; i<I; i++){
            bottom = m[calcIndex(i, j-1, I)];
            bottomRight = m[calcIndex(i + 1, j-1, I)];
            bottomLeft = m[calcIndex(i - 1, j-1, I)];
            m[calcIndex(i, j, I)] = innerNode(bottom, bottomLeft, bottomRight, params);
            cout << 1<<", "<< j << ", "<<i<<endl;
        }
        m[calcIndex(I, j, I)] = 0;
    }
}


int main(){
	int I = 100;
    int J = 1000;
    const int iterations = 1;
   	ModelParams params = prepareParams(I, J);
    int vecSize = (I+1)*(J+1);

    double *m = (double *) calloc(vecSize, sizeof(double));
    int64_t avg = 0;
    int64_t duration;
    for (int i = 0; i < iterations; i++)
    {
        auto start = chrono::high_resolution_clock::now();
        calcM(m, I, J, params);
        auto stop = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(stop-start).count();
		cout << "Iteration: " << i << "; Time: " << duration << ";\n";
        avg += duration/iterations;
    }
    
    cout << "Average sequential time: " << avg << " ms" << endl;
    
    const char* path = "output.txt";

    saveMatrix(m, I, J, path);

    free(m);

	return 0;
}