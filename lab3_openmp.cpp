#include "heat_scheme.cpp"
#include <omp.h>
#include<iostream>
#include <chrono>
using namespace std;

const int STREAMS = 10;

void barrier(){
	#pragma omp barrier
}

void calcM(double* m, int I, int J, const ModelParams params){
	double hT = params.T/J;
    double hZ = params.L/I;
    double bottom;
    double bottomLeft;
    double bottomRight;
    int i, j, ID, threads;
	threads = 1;


    #pragma omp parallel num_threads(1) private(i, j, ID) shared(m)
    {
    	// threads = omp_get_num_threads();
		ID = omp_get_thread_num();
	
    	for(i=ID; i <= I; i+=threads)
    	{
        	m[calcIndex(i, 0, I)] = initCondition(i, I);
			
    	}
    	#pragma omp barrier
		

		for (j = 1; j <= J; j++)
		{	
			for(i = ID; i<=I; i+=threads)
			{
				// #pragma omp critical
				// cout << ID << ", " << i << endl;
				// if(i==0)
    			// {
					
	        		bottom = m[calcIndex(i, j-1, I)];
        		// 	bottomRight = m[calcIndex(1, j-1, I)];
        		// 	m[calcIndex(0, j, I)] = leftCondition(bottom, bottomRight, params);
    			// } else if(i == I){
					m[calcIndex(i, j, I)] = bottom;
				// } else
    			// {
    			// 	bottom = m[calcIndex(i, j-1, I)];
            	// 	bottomRight = m[calcIndex(i + 1, j-1, I)];
        	    // 	bottomLeft = m[calcIndex(i - 1, j-1, I)];
		        //     m[calcIndex(i, j, I)] = innerNode(bottom, bottomLeft, bottomRight, params);
    			// }
			}
			#pragma omp barrier
			// if(ID == 0){
			// 	cout << j << endl;
			// }

		}
    }
}

int main(){
	int I = 100;
    int J = 1000;

   	ModelParams params = prepareParams(I, J);
    int vecSize = (I+1)*(J+1);

    double *m = (double *) calloc(vecSize, sizeof(double));

	const int iterations = 1;
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
	cout << "Average time: " << avg << ";\n";
	saveMatrix(m, I, J, "output_omp.txt");
	free(m);
	return 0;
}