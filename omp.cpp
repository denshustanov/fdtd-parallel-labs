#include "heat_scheme.cpp"
#include <omp.h>
#include<iostream>

const int STREAMS = 10;

void calcM(double* m, int I, int J, const ModelParams params){
	double hT = params.T/J;
    double hZ = params.L/I;
    double bottom;
    double bottomLeft;
    double bottomRight;
    int i, j;

    // omp_set_dynamic(0);
	// omp_set_num_threads(4);

    #pragma omp parallel num_threads(8)
    {
    	#pragma omp for schedule(static)
    	for(int i=0; i <= I; i++)
    	{
        	m[calcIndex(i, 0, I)] = initCondition(i, I);
    	}
    	
    	// std::cout<<"a\n";

    	#pragma omp for schedule(static)
    	for (int i = 0; i < I; ++i)
    	{
    		for(j=1; j < J; ++j)
    		{

    			if(i==0)
    			{
	        		bottom = m[calcIndex(0, j, I)];
	        		bottomRight = m[calcIndex(1, j, I)];
	        		m[calcIndex(0, j+1, I)] = leftCondition(bottom, bottomRight, params);
    			} else
    			{
    				bottom = m[calcIndex(i, j, I)];
	            	bottomRight = m[calcIndex(i + 1, j, I)];
	            	bottomLeft = m[calcIndex(i - 1, j, I)];
	            	m[calcIndex(i, j+1, I)] = innerNode(bottom, bottomLeft, bottomRight, params);
    			}
    		}
    	}
    }
}

int main(){
	int I = 10000;
    int J = 1000;

   	ModelParams params = prepareParams(I, J);
    int vecSize = (I+1)*(J+1);

    double *m = (double *) calloc(vecSize, sizeof(double));

	calcM(m, I, J, params);
	saveMatrix(m, I, J, "output_omp.txt");
	return 0;
}